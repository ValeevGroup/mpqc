#pragma once

#ifndef MPQC_BASIS_SHELLVECFUNCTIONS_H
#define MPQC_BASIS_SHELLVECFUNCTIONS_H

#include "../common/typedefs.h"
#include "basis_fwd.h"

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

} // namespace basis
} // namespace mpqc



#endif // MPQC_BASIS_SHELLVECFUNCTIONS_H
