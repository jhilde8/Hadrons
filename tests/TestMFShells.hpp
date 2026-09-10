#ifndef TestMFShells_hpp
#define TestMFShells_hpp

#include <string>
#include <vector>

// Turn a comma-separated command-line list into the whitespace-separated form
// the module Par structs expect.
//
// No Grid command-line option may contain whitespace: GridCmdOptionPayload
// returns a single argv element, so "--gammas Gamma5 Identity" yields only
// "Gamma5". A quoted payload survives direct invocation but not an unquoted
// variable expansion in a batch script, which is how it usually arrives. Hence
// the Grid convention of dots for coordinates (8.8.8.16) and commas for lists.
// The module side parses with strToVec and wants spaces, so the two meet here.
//
// A payload that already contains spaces passes through unchanged, so both
// forms work once the shell has been persuaded to deliver one argv element.
inline std::string cliListToPar(const std::string &csl)
{
    std::string par = csl;

    for (auto &c : par)
        if (c == ',') c = ' ';

    return par;
}

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
