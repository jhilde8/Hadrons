/*
 * K2pipiContractorDense.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 */
#include "A2ADenseContract.hpp"
#include <Hadrons/Module.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace A2ADense;

namespace K2pipiContractorDense
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

    class ContractionPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(ContractionPar,
                                        std::string,              name,
                                        std::vector<TermPar>,     terms,
                                        std::vector<std::string>, swept,
                                        std::string,              translations,
                                        unsigned int,             min_t,
                                        unsigned int,             max_t,
                                        bool,                     translationAverage);
    };

    class CorrelatorResult: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(CorrelatorResult,
                                        std::vector<MatrixPar>, a2aMatrix,
                                        ContractionPar,         contraction,
                                        std::string,            swept,
                                        std::vector<ComplexD>,  correlator);
    };
}

using namespace K2pipiContractorDense;

struct K2pipiContractorDensePar
{
    GlobalPar                   global;
    std::vector<MatrixPar>      a2aMatrix;
    std::vector<ContractionPar> contraction;
};

struct Term
{
    const Field  *field;
    unsigned int offset;
};

struct SweptEdges
{
    Connection inner, closing;
};

struct ContractionPlan
{
    ContractionPar                                  par;
    std::vector<Term>                               terms;
    std::vector<Connection>                         chain;
    std::vector<SweptEdges>                         swept;
    std::set<unsigned int>                          translations;
    std::vector<bool>                               inTranslations;
    std::vector<MatrixPar>                          inputs;
    std::string                                     filename;
    std::vector<Matrix>                             fixed;
    std::vector<std::size_t>                        rowTime, colTime;
    std::vector<std::vector<std::vector<ComplexD>>> corr;
};

std::map<std::string, MatrixPar> matricesByName(const K2pipiContractorDensePar &par,
                                                const unsigned int traj)
{
    std::map<std::string, MatrixPar> byName;

    for (auto m: par.a2aMatrix)
    {
        tokenReplace(m.file, "traj", traj);
        if (!byName.emplace(m.name, m).second)
        {
            HADRONS_ERROR(Definition, "duplicate a2aMatrix name '" + m.name + "'");
        }
    }

    return byName;
}

const MatrixPar &matrix(const std::map<std::string, MatrixPar> &byName, const std::string &name)
{
    auto it = byName.find(name);

    if (it == byName.end())
    {
        HADRONS_ERROR(Definition, "no a2aMatrix named '" + name + "'");
    }

    return it->second;
}

ContractionPlan makePlan(const ContractionPar &p, const std::map<std::string, Field> &fixed,
                         const std::map<std::string, SliceReader> &swept,
                         const GlobalPar &global, const unsigned int traj)
{
    const unsigned int    nt   = global.nt;
    const std::size_t     n    = p.terms.size();
    const std::string     what = "contraction '" + p.name + "'";
    ContractionPlan       plan;
    std::set<std::string> sweptNames;

    if (n < 1)
    {
        HADRONS_ERROR(Size, what + ": needs at least one fixed term");
    }
    if (p.swept.empty())
    {
        HADRONS_ERROR(Size, what + ": needs at least one swept field");
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
    plan.inTranslations.assign(nt, false);
    for (auto dt: plan.translations)
    {
        plan.inTranslations[dt] = true;
    }
    for (auto &t: p.terms)
    {
        const Field &f = fixed.at(t.name);

        plan.terms.push_back({&f, t.offset % nt});
        if (std::find(plan.inputs.begin(), plan.inputs.end(), f.par) == plan.inputs.end())
        {
            plan.inputs.push_back(f.par);
        }
    }
    for (std::size_t j = 0; j + 1 < n; ++j)
    {
        plan.chain.push_back(makeConnection(plan.terms[j].field->col,
                                            plan.terms[j + 1].field->row, nt,
                                            what + ", edge " + p.terms[j].name
                                            + " -> " + p.terms[j + 1].name));
    }
    for (auto &name: p.swept)
    {
        const SliceReader &e = swept.at(name);

        if (!sweptNames.insert(name).second)
        {
            HADRONS_ERROR(Definition, what + ": swept field '" + name + "' listed twice");
        }
        plan.swept.push_back({makeConnection(plan.terms[n - 1].field->col, e.row, nt,
                                             what + ", edge " + p.terms[n - 1].name
                                             + " -> " + name),
                              makeConnection(e.col, plan.terms[0].field->row, nt,
                                             what + ", edge " + name + " -> "
                                             + p.terms[0].name)});
    }
    plan.filename = global.output + "/" + ModuleBase::resultFilename(p.name, traj);

    return plan;
}

// Internal edges of the fixed chain meet one W time, so those axes are
// gathered before any product; the two outer axes stay open, since they
// close against the swept field at every tau.
void computeFixed(ContractionPlan &plan, const unsigned int nt)
{
    const std::size_t          n = plan.terms.size();
    std::vector<Matrix>        g(n);
    std::vector<const Matrix*> gp(n);
    Matrix                     prod, tmp, gtmp;
    std::vector<std::size_t>   ts(n);

    std::cout << "======== tr(";
    for (std::size_t j = 0; j < n; ++j)
    {
        std::cout << plan.par.terms[j].name << "[dt+" << plan.terms[j].offset << "] ";
    }
    std::cout << "swept[dt+tau]) for " << plan.par.swept.size() << " swept field(s) -> "
              << plan.filename << std::endl;
    plan.fixed.resize(nt);
    plan.rowTime.assign(nt, 0);
    plan.colTime.assign(nt, 0);
    plan.corr.assign(plan.swept.size(), std::vector<std::vector<ComplexD>>(
                         nt, std::vector<ComplexD>(nt, 0.)));
    for (auto dt: plan.translations)
    {
        const Matrix *p;

        for (std::size_t j = 0; j < n; ++j)
        {
            ts[j] = (dt + plan.terms[j].offset) % nt;
        }
        for (std::size_t j = 0; j < n; ++j)
        {
            const bool inner = (j > 0), outer = (j + 1 < n);

            gp[j] = &gathered(g[j], plan.terms[j].field->slice[ts[j]],
                              inner ? plan.chain[j - 1] : openEdge, inner ? ts[j - 1] : 0,
                              outer ? plan.chain[j] : openEdge, outer ? ts[j + 1] : 0, gtmp);
        }
        p = gp[0];
        for (std::size_t j = 1; j < n; ++j)
        {
            tmp.noalias() = (*p)*(*gp[j]);
            prod.swap(tmp);
            p = &prod;
        }
        plan.fixed[dt]   = *p;
        plan.rowTime[dt] = ts[0];
        plan.colTime[dt] = ts[n - 1];
    }
    plan.terms.clear();
}

void sweep(const SliceReader &reader, std::vector<ContractionPlan> &plans,
           const unsigned int nt)
{
    const std::string &name = reader.par.name;
    Matrix             e;
    TraceBuffers       tbuf;
    double             tRead = 0., tTrace = 0.;

    for (unsigned int top = 0; top < nt; ++top)
    {
        auto start = std::chrono::steady_clock::now();

        reader.read(e, top);
        tRead += seconds(start);
        start  = std::chrono::steady_clock::now();
        for (auto &plan: plans)
        {
            const auto &list = plan.par.swept;
            const auto it    = std::find(list.begin(), list.end(), name);

            if (it == list.end())
            {
                continue;
            }

            const std::size_t i      = it - list.begin();
            const SweptEdges  &edges = plan.swept[i];

            for (unsigned int tau = plan.par.min_t; tau <= plan.par.max_t; ++tau)
            {
                const unsigned int dt = (top + nt - tau) % nt;

                if (plan.inTranslations[dt])
                {
                    accTrace(plan.corr[i][dt][tau], plan.fixed[dt], plan.rowTime[dt],
                             plan.colTime[dt], e, top, edges.inner, edges.closing, tbuf);
                }
            }
        }
        tTrace += seconds(start);
    }
    std::cout << "  '" << name << "': read " << tRead << " s, traces " << tTrace << " s"
              << std::endl;
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

void writeResults(const ContractionPlan &plan, const std::map<std::string, SliceReader> &swept,
                  const unsigned int nt)
{
    makeFileDir(plan.filename);

    ResultWriter writer(plan.filename);

    for (std::size_t i = 0; i < plan.par.swept.size(); ++i)
    {
        const std::string &name = plan.par.swept[i];
        const std::string stem  = plan.par.name + "_" + name;
        CorrelatorResult  result;

        result.a2aMatrix = plan.inputs;
        result.a2aMatrix.push_back(swept.at(name).par);
        result.contraction = plan.par;
        result.swept       = name;
        if (plan.par.translationAverage)
        {
            result.correlator.assign(nt, 0.);
            for (auto dt: plan.translations)
            {
                for (unsigned int tau = 0; tau < nt; ++tau)
                {
                    result.correlator[tau] += plan.corr[i][dt][tau];
                }
            }
            for (auto &z: result.correlator)
            {
                z /= static_cast<double>(plan.translations.size());
            }
            writeResult(writer, stem, result);
        }
        else
        {
            for (auto dt: plan.translations)
            {
                result.correlator = plan.corr[i][dt];
                writeResult(writer, stem + "_dt_" + std::to_string(dt), result);
            }
        }
    }
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
        K2pipiContractorDensePar par;
        XmlReader                reader(argv[1]);
        std::set<std::string>    names;

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
        std::cout << "K2pipiContractorDense: " << Eigen::nbThreads() << " thread(s)"
                  << std::endl;

        const auto &range = par.global.trajCounter;
        const auto nt     = par.global.nt;

        for (unsigned int traj = range.start; traj < range.end; traj += range.step)
        {
            std::cout << ":::::::: Trajectory " << traj << std::endl;

            const auto                         byName = matricesByName(par, traj);
            const auto                         &cs    = par.contraction;
            std::map<std::string, Field>       fixed;
            std::map<std::string, SliceReader> swept;
            std::map<std::string, std::size_t> lastUse;
            std::vector<ContractionPlan>       plans;
            std::size_t                        nLoad = 0, peak = 0;
            double                             tLoad = 0., tFixed = 0.;

            for (std::size_t i = 0; i < cs.size(); ++i)
            {
                const auto filename = par.global.output + "/"
                                      + ModuleBase::resultFilename(cs[i].name, traj);

                if (fileExists(filename))
                {
                    HADRONS_ERROR(Io, "contraction '" + cs[i].name + "': output '"
                                      + filename + "' already exists");
                }
                for (auto &t: cs[i].terms)
                {
                    lastUse[t.name] = i;
                }
                for (auto &name: cs[i].swept)
                {
                    if (!swept.count(name))
                    {
                        swept.emplace(name, openSliceReader(matrix(byName, name), nt));
                    }
                }
            }
            for (std::size_t i = 0; i < cs.size(); ++i)
            {
                auto start = std::chrono::steady_clock::now();

                for (auto &t: cs[i].terms)
                {
                    if (!fixed.count(t.name))
                    {
                        fixed.emplace(t.name, loadField(matrix(byName, t.name), nt));
                        ++nLoad;
                    }
                }
                peak   = std::max(peak, fixed.size());
                tLoad += seconds(start);
                start  = std::chrono::steady_clock::now();
                plans.push_back(makePlan(cs[i], fixed, swept, par.global, traj));
                computeFixed(plans.back(), nt);
                tFixed += seconds(start);
                for (auto &t: cs[i].terms)
                {
                    if (lastUse.at(t.name) == i)
                    {
                        fixed.erase(t.name);
                    }
                }
            }
            std::cout << "Loaded " << nLoad << " fixed field(s), at most " << peak
                      << " resident, in " << tLoad << " s; opened " << swept.size()
                      << " swept field(s); fixed chains in " << tFixed << " s" << std::endl;
            for (auto &s: swept)
            {
                sweep(s.second, plans, nt);
            }
            for (auto &plan: plans)
            {
                writeResults(plan, swept, nt);
            }
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "K2pipiContractorDense: " << e.what() << std::endl;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
