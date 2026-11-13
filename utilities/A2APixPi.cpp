/*
 * A2APixPi.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 */
#include "A2ADenseContract.hpp"

using namespace Grid;
using namespace Hadrons;
using namespace A2ADense;

namespace A2APixPi
{
    class GlobalPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(GlobalPar,
                                        TrajRange,    trajCounter,
                                        unsigned int, nt);
    };

    class ProductPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(ProductPar,
                                        std::string,  left,
                                        std::string,  right,
                                        unsigned int, deltaL,
                                        unsigned int, deltaR,
                                        std::string,  output,
                                        std::string,  dataset);
    };
}

struct A2APixPiPar
{
    A2APixPi::GlobalPar               global;
    std::vector<MatrixPar>            a2aMatrix;
    std::vector<A2APixPi::ProductPar> product;
};

struct ProductPlan
{
    A2APixPi::ProductPar par;
    const Field          *left, *right;
    Connection           connection;
    ProductMetadata      metadata;
};

std::map<std::string, Field> loadInputs(const A2APixPiPar &par, const unsigned int traj)
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
    for (auto &p: par.product)
    {
        used.insert(p.left);
        used.insert(p.right);
    }
    for (auto &name: used)
    {
        auto it = byName.find(name);

        if (it == byName.end())
        {
            HADRONS_ERROR(Definition, "no a2aMatrix named '" + name + "'");
        }

        MatrixPar m = it->second;

        if ((m.row.timeOffset != 0) or (m.col.timeOffset != 0))
        {
            HADRONS_ERROR(Argument, "'" + name + "': inputs must have zero time offsets");
        }
        tokenReplace(m.file, "traj", traj);

        Field f = loadField(m, par.global.nt);

        if (f.isProduct)
        {
            HADRONS_ERROR(Argument, "'" + name + "' is a product; products are not "
                          "accepted as inputs");
        }
        fields.emplace(name, std::move(f));
    }

    return fields;
}

ProductPlan makePlan(const A2APixPi::ProductPar &p, const std::map<std::string, Field> &fields,
                     const unsigned int nt, const unsigned int traj)
{
    ProductPlan plan;
    const std::string what = "product '" + p.output + "'";

    if ((p.deltaL >= nt) or (p.deltaR >= nt) or ((p.deltaL != 0) and (p.deltaR != 0)))
    {
        HADRONS_ERROR(Range, what + ": need deltaL, deltaR < nt with at least one zero");
    }
    plan.par = p;
    tokenReplace(plan.par.output, "traj", traj);
    plan.left       = &fields.at(p.left);
    plan.right      = &fields.at(p.right);
    plan.connection = makeConnection(plan.left->col, plan.right->row, nt, what);
    plan.metadata.row            = plan.left->par.row;
    plan.metadata.row.timeOffset = (plan.left->par.row.timeOffset + p.deltaL) % nt;
    plan.metadata.col            = plan.right->par.col;
    plan.metadata.col.timeOffset = (plan.right->par.col.timeOffset + p.deltaR) % nt;
    plan.metadata.left           = plan.left->par;
    plan.metadata.right          = plan.right->par;
    plan.metadata.deltaL         = p.deltaL;
    plan.metadata.deltaR         = p.deltaR;
    plan.metadata.nt             = nt;
    for (unsigned int t = 0; t < nt; ++t)
    {
        if (fileExists(sliceFilename(plan.par.output, t)))
        {
            HADRONS_ERROR(Io, what + ": '" + sliceFilename(plan.par.output, t)
                          + "' already exists");
        }
    }

    return plan;
}

void runProduct(const ProductPlan &plan, const unsigned int nt)
{
    Matrix result, buf;
    double tMul = 0., tWrite = 0.;

    std::cout << "======== " << plan.par.left << "[t+" << plan.par.deltaL << "] * "
              << plan.par.right << "[t+" << plan.par.deltaR << "] -> "
              << plan.par.output << " (" << plan.left->row.size << "x"
              << plan.right->col.size << ", summed " << plan.connection.active
              << ")" << std::endl;
    makeFileDir(plan.par.output);
    for (unsigned int t = 0; t < nt; ++t)
    {
        const unsigned int tL = (t + plan.par.deltaL) % nt;
        const unsigned int tR = (t + plan.par.deltaR) % nt;
        auto               start = std::chrono::steady_clock::now();

        multiply(result, plan.left->slice[tL], tL, plan.right->slice[tR], tR,
                 plan.connection, buf);
        tMul += seconds(start);
        if (!result.allFinite())
        {
            HADRONS_ERROR(Range, "non-finite product at t = " + std::to_string(t));
        }
        start = std::chrono::steady_clock::now();
        writeSlice(sliceFilename(plan.par.output, t), plan.par.dataset, result,
                   plan.metadata);
        tWrite += seconds(start);
    }
    std::cout << "  gather+GEMM " << tMul << " s, write " << tWrite << " s" << std::endl;
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
        A2APixPiPar par;
        XmlReader   reader(argv[1]);

        read(reader, "global",    par.global);
        read(reader, "a2aMatrix", par.a2aMatrix);
        read(reader, "product",   par.product);
        std::cout << "A2APixPi: " << Eigen::nbThreads() << " thread(s)" << std::endl;

        const auto &range = par.global.trajCounter;

        for (unsigned int traj = range.start; traj < range.end; traj += range.step)
        {
            std::cout << ":::::::: Trajectory " << traj << std::endl;

            auto start  = std::chrono::steady_clock::now();
            auto fields = loadInputs(par, traj);

            std::cout << "Loaded " << fields.size() << " field(s) in " << seconds(start)
                      << " s" << std::endl;

            std::vector<ProductPlan> plans;

            for (auto &p: par.product)
            {
                plans.push_back(makePlan(p, fields, par.global.nt, traj));
            }
            for (auto &plan: plans)
            {
                runProduct(plan, par.global.nt);
            }
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "A2APixPi: " << e.what() << std::endl;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
