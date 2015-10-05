#include "shell_vec_functions.h"

#include "../include/libint.h"

#include <algorithm>
#include <numeric>

namespace mpqc {
namespace basis {

int64_t max_am(ShellVec const &shell_vec) {
    int64_t am = 0;

    for (auto const &sh : shell_vec) {
        for (auto const &c : sh.contr) {
            am = std::max(int64_t(c.l), am);
        }
    }

    return am;
}

int64_t max_nprim(ShellVec const &shell_vec) {
    return int64_t(std::max_element(shell_vec.begin(), shell_vec.end(),
                                    [](Shell const &a, Shell const &b) {
                                        return std::max(int64_t(a.nprim()),
                                                        int64_t(b.nprim()));
                                    })->nprim());
}

int64_t nfunctions(ShellVec const &shell_vec) {
    return std::accumulate(shell_vec.begin(), shell_vec.end(), int64_t(0),
                           [](int64_t x, Shell const &s) {
        return x + s.size();
    });
}

} // namespace basis
} // namespace mpqc
