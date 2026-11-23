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

struct FixedProduct
{
    std::vector<Term>       factors;
    std::vector<Connection> edges;
    Field                   field;
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

const Field *fixedProduct(std::map<std::string, FixedProduct> &products,
                          const std::vector<Term> &factors,
                          const std::vector<std::string> &factorNames,
                          const unsigned int nt, const std::string &what)
{
    const std::size_t m = factors.size();
    std::string       key = "(";

    for (std::size_t j = 0; j < m; ++j)
    {
        key += factorNames[j] + "[+" + std::to_string(factors[j].offset) + "]"
               + ((j + 1 < m) ? " " : ")");
    }

    auto it = products.find(key);

    if (it != products.end())
    {
        return &it->second.field;
    }

    FixedProduct f;

    f.factors = factors;
    for (std::size_t j = 0; j + 1 < m; ++j)
    {
        f.edges.push_back(makeConnection(factors[j].field->col, factors[j + 1].field->row, nt,
                                         what + ", edge " + factorNames[j] + " -> "
                                         + factorNames[j + 1]));
    }
    f.field.par.name       = key;
    f.field.row            = factors.front().field->row;
    f.field.row.timeOffset = (f.field.row.timeOffset + factors.front().offset) % nt;
    f.field.col            = factors.back().field->col;
    f.field.col.timeOffset = (f.field.col.timeOffset + factors.back().offset) % nt;
    f.field.isProduct      = true;

    return &products.emplace(key, std::move(f)).first->second.field;
}

void computeFixed(FixedProduct &f, const unsigned int nt)
{
    const std::size_t          m = f.factors.size();
    std::vector<Matrix>        g(m);
    std::vector<const Matrix*> gp(m);
    std::vector<std::size_t>   ts(m);
    Matrix                     tmp, gtmp;

    f.field.slice.resize(nt);
    for (unsigned int t = 0; t < nt; ++t)
    {
        Matrix &out = f.field.slice[t];

        for (std::size_t j = 0; j < m; ++j)
        {
            ts[j] = (t + f.factors[j].offset) % nt;
        }
        for (std::size_t j = 0; j < m; ++j)
        {
            const bool inner = (j > 0), outer = (j + 1 < m);

            gp[j] = &gathered(g[j], f.factors[j].field->slice[ts[j]],
                              inner ? f.edges[j - 1] : openEdge, inner ? ts[j - 1] : 0,
                              outer ? f.edges[j] : openEdge, outer ? ts[j + 1] : 0, gtmp);
        }
        out.noalias() = (*gp[0])*(*gp[1]);
        for (std::size_t j = 2; j < m; ++j)
        {
            tmp.noalias() = out*(*gp[j]);
            out.swap(tmp);
        }
    }
    f.factors.clear();
    f.edges.clear();
}

ContractionPlan makePlan(const ContractionPar &p, const std::map<std::string, Field> &fields,
                         std::map<std::string, FixedProduct> &products,
                         const GlobalPar &global, const unsigned int traj)
{
    const unsigned int nt   = global.nt;
    const std::size_t  n    = p.terms.size();
    const std::string  what = "contraction '" + p.name + "'";
    ContractionPlan    plan;
    std::vector<Term>  parTerms;

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
        parTerms.push_back({&fields.at(t.name), t.role == "sink", t.offset % nt});
    }
    if (parTerms[0].sink)
    {
        HADRONS_ERROR(Argument, what + ": the first term must be a source");
    }

    std::size_t k = 0;

    while ((k < n) and (parTerms[k].sink or !parTerms[(k + n - 1) % n].sink))
    {
        ++k;
    }
    if (k == n)
    {
        HADRONS_ERROR(Argument, what + ": needs at least one sink term");
    }
    for (std::size_t j = 0; j < n;)
    {
        std::vector<Term>        factors;
        std::vector<std::string> factorNames;
        const bool               sink = parTerms[(k + j) % n].sink;

        while ((j < n) and (parTerms[(k + j) % n].sink == sink))
        {
            factors.push_back(parTerms[(k + j) % n]);
            factorNames.push_back(p.terms[(k + j) % n].name);
            ++j;
        }
        if (factors.size() == 1)
        {
            plan.terms.push_back(factors[0]);
        }
        else
        {
            plan.terms.push_back({fixedProduct(products, factors, factorNames, nt, what),
                                  sink, 0});
        }
    }

    const std::size_t ne = plan.terms.size();

    for (std::size_t j = 0; j < ne; ++j)
    {
        const Term &a = plan.terms[j], &b = plan.terms[(j + 1) % ne];

        plan.edges.push_back(makeConnection(a.field->col, b.field->row, nt,
                                            what + ", edge " + a.field->par.name
                                            + " -> " + b.field->par.name));
    }
    for (auto &t: parTerms)
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

// At each (dt, s) every term's expanded axes meet one W time, so all terms are
// gathered against their neighbours first and the chain multiplies dense-sized
// blocks only.
void runContraction(const ContractionPlan &plan, const unsigned int nt)
{
    const std::size_t          n = plan.terms.size();
    std::vector<Matrix>        g(n - 1);
    std::vector<const Matrix*> gp(n - 1);
    MatrixTr                   last;
    Matrix                     prod, tmp, gtmp;
    std::vector<std::size_t>   ts(n);
    CorrelatorResult           result;
    double                     tGather = 0., tMul = 0., tTrace = 0.;

    result.a2aMatrix = plan.inputs;
    result.contraction = plan.par;
    result.correlator.assign(nt, 0.);
    std::cout << "======== tr(";
    for (std::size_t j = 0; j < n; ++j)
    {
        std::cout << plan.terms[j].field->par.name << "["
                  << (plan.terms[j].sink ? "tLast" : "dt") << "+"
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

            auto start = std::chrono::steady_clock::now();

            for (std::size_t j = 0; j + 1 < n; ++j)
            {
                const std::size_t l = (j + n - 1) % n;

                gp[j] = &gathered(g[j], plan.terms[j].field->slice[ts[j]],
                                  plan.edges[l], ts[l], plan.edges[j], ts[j + 1], gtmp);
            }
            select(last, plan.terms[n - 1].field->slice[ts[n - 1]], plan.edges[n - 2],
                   ts[n - 2], plan.edges[n - 1], ts[0], gtmp);
            tGather += seconds(start);
            start    = std::chrono::steady_clock::now();
            p        = gp[0];
            for (std::size_t j = 1; j + 1 < n; ++j)
            {
                tmp.noalias() = (*p)*(*gp[j]);
                prod.swap(tmp);
                p = &prod;
            }
            tMul   += seconds(start);
            start   = std::chrono::steady_clock::now();
            A2AContraction::accTrMul(result.correlator[s], *p, last);
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
    std::cout << "  gathers " << tGather << " s, products " << tMul << " s, traces "
              << tTrace << " s" << std::endl;
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

            std::vector<ContractionPlan>        plans;
            std::map<std::string, FixedProduct> products;
            std::set<const Field *>             fieldChain;

            for (auto &c: par.contraction)
            {
                plans.push_back(makePlan(c, fields, products, par.global, traj));
            }
            start = std::chrono::steady_clock::now();
            for (auto &f: products)
            {
                computeFixed(f.second, par.global.nt);
            }
            if (!products.empty())
            {
                std::cout << "Fixed products: " << products.size() << " in "
                          << seconds(start) << " s" << std::endl;
            }
            for (auto &plan: plans)
            {
                for (auto &t: plan.terms)
                {
                    fieldChain.insert(t.field);
                }
            }
            for (auto it = fields.begin(); it != fields.end();)
            {
                it = fieldChain.count(&it->second) ? std::next(it) : fields.erase(it);
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
