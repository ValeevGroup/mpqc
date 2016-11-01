#include <mpqc/chemistry/qc/basis/shell_vec_functions.h>

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
    std::size_t n = 0;
    for(auto const& shell : shell_vec){
        n = std::max(shell.nprim(), n);
    }
    return n;
}

int64_t nfunctions(ShellVec const &shell_vec) {
    return std::accumulate(shell_vec.begin(), shell_vec.end(), int64_t(0),
                           [](int64_t x, Shell const &s) {
        return x + s.size();
    });
}

std::vector<std::vector<libint2::Shell>>
reblock_basis(std::vector<libint2::Shell> shells, std::size_t blocksize){

    std::vector<std::vector<libint2::Shell>> result;

    std::vector<libint2::Shell> tmp;
    std::size_t tmp_size = 0;
    for (auto& shell : shells){

        std::size_t shell_size = shell.size();
        tmp_size += shell_size;

        if(4*tmp_size < 5*blocksize){
            tmp.push_back(shell);
        }
        else{
            result.push_back(tmp);
            tmp = std::vector<libint2::Shell>();
            tmp.push_back(shell);
            tmp_size = shell_size;
        }
    }

    // handle the boundary condition
    if(4*tmp_size < 5*blocksize){

        // if boundary is less than 2/3 of the block size
        // include it to previous block
        if(3*tmp_size < 2*blocksize){
            result.back().insert(result.back().end(),tmp.begin(),tmp.end());
        }
        else{
            result.push_back(tmp);
        }
    }
    return result;
}

} // namespace basis
} // namespace mpqc
