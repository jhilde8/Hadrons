/*
 * BubbleContractorDense.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 */
#include "A2ADenseContract.hpp"
#include <Hadrons/Module.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace A2ADense;

namespace BubbleContractorDense
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
                                        unsigned int, offset);
    };

    class BubblePar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(BubblePar,
                                        std::string,          name,
                                        std::vector<TermPar>, terms);
    };

    class BubbleResult: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(BubbleResult,
                                        std::vector<MatrixPar>, a2aMatrix,
                                        BubblePar,              bubble,
                                        std::vector<ComplexD>,  correlator);
    };
}

using namespace BubbleContractorDense;

struct BubbleContractorDensePar
{
    GlobalPar              global;
    std::vector<MatrixPar> a2aMatrix;
    std::vector<BubblePar> bubble;
};

struct Term
{
    const Field  *field;
    unsigned int offset;
};

struct BubblePlan
{
    BubblePar               par;
    std::vector<Term>       terms;
    std::vector<Connection> edges;
    std::vector<MatrixPar>  inputs;
    std::string             filename;
};

std::map<std::string, Field> loadInputs(const BubbleContractorDensePar &par,
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
    for (auto &b: par.bubble)
    {
        for (auto &t: b.terms)
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

BubblePlan makePlan(const BubblePar &p, const std::map<std::string, Field> &fields,
                    const GlobalPar &global, const unsigned int traj)
{
    const unsigned int nt   = global.nt;
    const std::size_t  n    = p.terms.size();
    const std::string  what = "bubble '" + p.name + "'";
    BubblePlan         plan;

    if ((n < 1) or (n > 2))
    {
        HADRONS_ERROR(Size, what + ": needs one or two terms; a longer fixed-offset "
                      "chain is a product, build it with A2APixPi");
    }
    plan.par = p;
    for (auto &t: p.terms)
    {
        plan.terms.push_back({&fields.at(t.name), t.offset % nt});
    }
    for (std::size_t j = 0; j < n; ++j)
    {
        const std::size_t k = (j + 1) % n;

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

void runBubble(const BubblePlan &plan, const unsigned int nt)
{
    const std::size_t n = plan.terms.size();
    TraceBuffers      tbuf;
    BubbleResult      result;

    result.a2aMatrix = plan.inputs;
    result.bubble    = plan.par;
    result.correlator.assign(nt, 0.);
    std::cout << "======== tr(";
    for (std::size_t j = 0; j < n; ++j)
    {
        std::cout << plan.par.terms[j].name << "[t+" << plan.terms[j].offset << "]"
                  << ((j + 1 < n) ? " " : ")");
    }
    std::cout << " -> " << plan.filename << std::endl;

    auto start = std::chrono::steady_clock::now();

    for (unsigned int t = 0; t < nt; ++t)
    {
        const unsigned int t0 = (t + plan.terms[0].offset) % nt;

        if (n == 1)
        {
            accSelfTrace(result.correlator[t], plan.terms[0].field->slice[t0], t0,
                         plan.edges[0]);
        }
        else
        {
            const unsigned int t1 = (t + plan.terms[1].offset) % nt;

            accTrace(result.correlator[t], plan.terms[0].field->slice[t0], t0, t0,
                     plan.terms[1].field->slice[t1], t1, plan.edges[0], plan.edges[1],
                     tbuf);
        }
    }
    std::cout << "  traces " << seconds(start) << " s" << std::endl;
    for (auto &z: result.correlator)
    {
        if (!std::isfinite(z.real()) or !std::isfinite(z.imag()))
        {
            HADRONS_ERROR(Range, "non-finite bubble in '" + plan.par.name + "'");
        }
    }
    makeFileDir(plan.filename);

    ResultWriter writer(plan.filename);

    write(writer, plan.par.name, result);
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
        BubbleContractorDensePar par;
        XmlReader                reader(argv[1]);
        std::set<std::string>    names;

        read(reader, "global",    par.global);
        read(reader, "a2aMatrix", par.a2aMatrix);
        read(reader, "bubble",    par.bubble);
        for (auto &b: par.bubble)
        {
            if (!names.insert(b.name).second)
            {
                HADRONS_ERROR(Definition, "duplicate bubble name '" + b.name + "'");
            }
        }
        std::cout << "BubbleContractorDense: " << Eigen::nbThreads() << " thread(s)"
                  << std::endl;

        const auto &range = par.global.trajCounter;

        for (unsigned int traj = range.start; traj < range.end; traj += range.step)
        {
            std::cout << ":::::::: Trajectory " << traj << std::endl;

            auto start  = std::chrono::steady_clock::now();
            auto fields = loadInputs(par, traj);

            std::cout << "Loaded " << fields.size() << " field(s) in " << seconds(start)
                      << " s" << std::endl;

            std::vector<BubblePlan> plans;

            for (auto &b: par.bubble)
            {
                plans.push_back(makePlan(b, fields, par.global, traj));
            }
            for (auto &plan: plans)
            {
                runBubble(plan, par.global.nt);
            }
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "BubbleContractorDense: " << e.what() << std::endl;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
