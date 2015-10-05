#pragma once
#ifndef MPQC_BASIS_BASISFWD_H
#define MPQC_BASIS_BASISFWD_H

#include "../common/typedefs.h"
#include <vector>

namespace mpqc {
namespace basis {

using ShellVec = std::vector<Shell>;

class BasisSet;
class AtomBasisSet;
class ClusterShells;

} // namespace basis
} // namespace tcc

#endif // MPQC_BASIS_BASISFWD_H
