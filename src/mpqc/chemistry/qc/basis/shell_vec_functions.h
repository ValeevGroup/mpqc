

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_BASIS_SHELL_VEC_FUNCTIONS_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_BASIS_SHELL_VEC_FUNCTIONS_H_

#include "mpqc/chemistry/qc/basis/basis.h"

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
std::vector<std::vector<libint2::Shell>> reblock_basis(
    std::vector<libint2::Shell> shells, std::size_t blocksize);

}  // namespace basis
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_BASIS_SHELL_VEC_FUNCTIONS_H_
