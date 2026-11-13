/*
 * FlexibleContractorDense.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 */
#include "A2ADenseContract.hpp"
#include <Hadrons/Module.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace A2ADense;

namespace FlexibleContractorDense
{
    class GlobalPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(GlobalPar,
                                        TrajRange,    trajCounter,
                                        unsigned int, nt,
                                        std::string,  output);
    };

    class TermPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(TermPar,
                                        std::string,  name,
                                        std::string,  role,
                                        unsigned int, offset);
    };

    class ContractionPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(ContractionPar,
                                        std::string,          name,
                                        std::vector<TermPar>, terms,
                                        std::string,          translations,
                                        unsigned int,         min_t,
                                        unsigned int,         max_t,
                                        bool,                 translationAverage);
    };

    class CorrelatorResult: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(CorrelatorResult,
                                        std::vector<MatrixPar>, a2aMatrix,
                                        ContractionPar,         contraction,
                                        std::vector<ComplexD>,  correlator);
    };
}

using namespace FlexibleContractorDense;

struct FlexibleContractorDensePar
{
    GlobalPar                   global;
    std::vector<MatrixPar>      a2aMatrix;
    std::vector<ContractionPar> contraction;
};

struct Term
{
    const Field  *field;
    bool         sink;
    unsigned int offset;
};

struct ContractionPlan
{
    ContractionPar          par;
    std::vector<Term>       terms;
    std::vector<Connection> edges;
    std::set<unsigned int>  translations;
    std::vector<MatrixPar>  inputs;
    std::string             filename;
};

std::map<std::string, Field> loadInputs(const FlexibleContractorDensePar &par,
                                        const unsigned int traj)
{
    std::map<std::string, MatrixPar> byName;
    std::map<std::string, Field>     fields;
    std::set<std::string>            used;

    for (auto &m: par.a2aMatrix)
    {
        if (!byName.emplace(m.name, m).second)
        {
            HADRONS_ERROR(Definition, "duplicate a2aMatrix name '" + m.name + "'");
        }
    }
    for (auto &c: par.contraction)
    {
        for (auto &t: c.terms)
        {
            used.insert(t.name);
        }
    }
    for (auto &name: used)
    {
        auto it = byName.find(name);

        if (it == byName.end())
        {
            HADRONS_ERROR(Definition, "no a2aMatrix named '" + name + "'");
        }

        MatrixPar m = it->second;

        tokenReplace(m.file, "traj", traj);
        fields.emplace(name, loadField(m, par.global.nt));
    }

    return fields;
}

ContractionPlan makePlan(const ContractionPar &p, const std::map<std::string, Field> &fields,
                         const GlobalPar &global, const unsigned int traj)
{
    const unsigned int nt   = global.nt;
    const std::size_t  n    = p.terms.size();
    const std::string  what = "contraction '" + p.name + "'";
    ContractionPlan    plan;

    if (n < 2)
    {
        HADRONS_ERROR(Size, what + ": needs at least two terms");
    }
    if ((p.min_t > p.max_t) or (p.max_t >= nt))
    {
        HADRONS_ERROR(Range, what + ": need min_t <= max_t < nt");
    }
    plan.par          = p;
    plan.translations = parseTimeRange(p.translations, nt);
    if (plan.translations.empty())
    {
        HADRONS_ERROR(Argument, what + ": empty translations");
    }
    for (auto &t: p.terms)
    {
        if ((t.role != "source") and (t.role != "sink"))
        {
            HADRONS_ERROR(Argument, what + ": role of '" + t.name + "' is '" + t.role
                          + "', expected 'source' or 'sink'");
        }
        plan.terms.push_back({&fields.at(t.name), t.role == "sink", t.offset % nt});
    }
    if (plan.terms[0].sink)
    {
        HADRONS_ERROR(Argument, what + ": the first term must be a source");
    }
    for (std::size_t j = 0; j < n; ++j)
    {
        const std::size_t k = (j + 1) % n;

        if (plan.terms[j].sink == plan.terms[k].sink)
        {
            HADRONS_ERROR(Argument, what + ": neighbouring terms '" + p.terms[j].name
                          + "' and '" + p.terms[k].name + "' are both " + p.terms[j].role
                          + "s; a fixed-offset pair is a product, build it with A2APixPi");
        }
        plan.edges.push_back(makeConnection(plan.terms[j].field->col,
                                            plan.terms[k].field->row, nt,
                                            what + ", edge " + p.terms[j].name
                                            + " -> " + p.terms[k].name));
    }
    for (auto &t: plan.terms)
    {
        if (std::find(plan.inputs.begin(), plan.inputs.end(), t.field->par) == plan.inputs.end())
        {
            plan.inputs.push_back(t.field->par);
        }
    }
    plan.filename = global.output + "/" + ModuleBase::resultFilename(p.name, traj);
    if (fileExists(plan.filename))
    {
        HADRONS_ERROR(Io, what + ": output '" + plan.filename + "' already exists");
    }

    return plan;
}

void writeResult(ResultWriter &writer, const std::string &name, CorrelatorResult &result)
{
    for (auto &z: result.correlator)
    {
        if (!std::isfinite(z.real()) or !std::isfinite(z.imag()))
        {
            HADRONS_ERROR(Range, "non-finite correlator in '" + name + "'");
        }
    }
    write(writer, name, result);
}

void runContraction(const ContractionPlan &plan, const unsigned int nt)
{
    const std::size_t        n = plan.terms.size();
    Matrix                   prod, tmp, buf;
    TraceBuffers             tbuf;
    std::vector<std::size_t> ts(n);
    CorrelatorResult         result;
    double                   tMul = 0., tTrace = 0.;

    result.a2aMatrix = plan.inputs;
    result.contraction = plan.par;
    result.correlator.assign(nt, 0.);
    std::cout << "======== tr(";
    for (std::size_t j = 0; j < n; ++j)
    {
        const auto &t = plan.par.terms[j];

        std::cout << t.name << "[" << (plan.terms[j].sink ? "tLast" : "dt") << "+"
                  << plan.terms[j].offset << "]" << ((j + 1 < n) ? " " : ")");
    }
    std::cout << " -> " << plan.filename << std::endl;
    makeFileDir(plan.filename);

    ResultWriter writer(plan.filename);

    for (auto dt: plan.translations)
    {
        for (unsigned int s = plan.par.min_t; s <= plan.par.max_t; ++s)
        {
            const unsigned int tLast = (dt + s) % nt;
            const Matrix       *p;

            for (std::size_t j = 0; j < n; ++j)
            {
                ts[j] = ((plan.terms[j].sink ? tLast : dt) + plan.terms[j].offset) % nt;
            }
            p = &plan.terms[0].field->slice[ts[0]];

            auto start = std::chrono::steady_clock::now();

            for (std::size_t j = 1; j + 1 < n; ++j)
            {
                multiply(tmp, *p, ts[j - 1], plan.terms[j].field->slice[ts[j]], ts[j],
                         plan.edges[j - 1], buf);
                prod.swap(tmp);
                p = &prod;
            }
            tMul  += seconds(start);
            start  = std::chrono::steady_clock::now();
            accTrace(result.correlator[s], *p, ts[0], ts[n - 2],
                     plan.terms[n - 1].field->slice[ts[n - 1]], ts[n - 1],
                     plan.edges[n - 2], plan.edges[n - 1], tbuf);
            tTrace += seconds(start);
        }
        if (!plan.par.translationAverage)
        {
            writeResult(writer, plan.par.name + "_dt_" + std::to_string(dt), result);
            std::fill(result.correlator.begin(), result.correlator.end(), ComplexD(0.));
        }
    }
    if (plan.par.translationAverage)
    {
        for (auto &z: result.correlator)
        {
            z /= static_cast<double>(plan.translations.size());
        }
        writeResult(writer, plan.par.name, result);
    }
    std::cout << "  products " << tMul << " s, traces " << tTrace << " s" << std::endl;
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cerr << "usage: " << argv[0] << " <parameter file>" << std::endl;

        return EXIT_FAILURE;
    }
    H5NS::Exception::dontPrint();
    try
    {
        FlexibleContractorDensePar par;
        XmlReader                  reader(argv[1]);
        std::set<std::string>      names;

        read(reader, "global",      par.global);
        read(reader, "a2aMatrix",   par.a2aMatrix);
        read(reader, "contraction", par.contraction);
        for (auto &c: par.contraction)
        {
            if (!names.insert(c.name).second)
            {
                HADRONS_ERROR(Definition, "duplicate contraction name '" + c.name + "'");
            }
        }
        std::cout << "FlexibleContractorDense: " << Eigen::nbThreads() << " thread(s)"
                  << std::endl;

        const auto &range = par.global.trajCounter;

        for (unsigned int traj = range.start; traj < range.end; traj += range.step)
        {
            std::cout << ":::::::: Trajectory " << traj << std::endl;

            auto start  = std::chrono::steady_clock::now();
            auto fields = loadInputs(par, traj);

            std::cout << "Loaded " << fields.size() << " field(s) in " << seconds(start)
                      << " s" << std::endl;

            std::vector<ContractionPlan> plans;

            for (auto &c: par.contraction)
            {
                plans.push_back(makePlan(c, fields, par.global, traj));
            }
            for (auto &plan: plans)
            {
                runContraction(plan, par.global.nt);
            }
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "FlexibleContractorDense: " << e.what() << std::endl;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
