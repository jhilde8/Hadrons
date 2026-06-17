#ifndef TestMFShells_hpp
#define TestMFShells_hpp

#include <string>
#include <vector>

// Returns all momentum strings in shells 0..maxShell (cumulative).
// Each string is "px py pz" with components in {-1,0,1}.
//
//   0  -> {0,0,0}               (1 vector,  total 1)
//   1  -> + permutations of {+-1,0,0}    (6 vectors, total  7)
//   2  -> + permutations of {+-1,+-1,0}  (12 vectors, total 19)
//   3  -> + all signs of {+-1,+-1,+-1}   (8 vectors,  total 27)
inline std::vector<std::string> momentumShells(int maxShell)
{
    std::vector<std::string> mom;

    if (maxShell >= 0)
        mom.push_back("0 0 0");

    if (maxShell >= 1)
        for (auto s : {"1 0 0", "-1 0 0", "0 1 0", "0 -1 0", "0 0 1", "0 0 -1"})
            mom.push_back(s);

    if (maxShell >= 2)
        for (auto s : {"1 1 0", "-1 1 0", "1 -1 0", "-1 -1 0",
                       "1 0 1", "-1 0 1", "1 0 -1", "-1 0 -1",
                       "0 1 1", "0 -1 1", "0 1 -1", "0 -1 -1"})
            mom.push_back(s);

    if (maxShell >= 3)
        for (auto s : {"1 1 1", "-1 1 1", "1 -1 1", "-1 -1 1",
                       "1 1 -1", "-1 1 -1", "1 -1 -1", "-1 -1 -1"})
            mom.push_back(s);

    return mom;
}

#endif // TestMFShells_hpp
