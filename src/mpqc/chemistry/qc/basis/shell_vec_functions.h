

#ifndef MPQC_BASIS_SHELLVECFUNCTIONS_H
#define MPQC_BASIS_SHELLVECFUNCTIONS_H


#include <mpqc/chemistry/qc/basis/basis.h>

#include <libint2/shell.h>

namespace mpqc {
namespace basis {

/*! \brief Returns the maximum angular momement of any shell in the vector. */
int64_t max_am(ShellVec const &);

/*! \brief Returns the maximum number of primatives of any shell in the vector.
 */
int64_t max_nprim(ShellVec const &);

/*! \brief Returns the maximum number of primatives of any shell in the vector.
 */
int64_t nfunctions(ShellVec const &);

// reblock based on blocksize
std::vector<std::vector<libint2::Shell>>
        reblock_basis(std::vector<libint2::Shell> shells, std::size_t blocksize);

} // namespace basis
} // namespace mpqc



#endif // MPQC_BASIS_SHELLVECFUNCTIONS_H
