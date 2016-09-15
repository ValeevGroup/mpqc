#pragma once
#ifndef TCC_COMMON_FWDDECL_H
#define TCC_COMMON_FWDDECL_H

#pragma GCC diagnostic push
#pragma GCC system_header
#include <tiledarray_fwd.h>
#pragma GCC diagnostic pop

/////////////////////////////////////////////////////////
// Madness FWD
/////////////////////////////////////////////////////////
namespace Madness {
    class World;
}

/////////////////////////////////////////////////////////
// TiledArray FWD
/////////////////////////////////////////////////////////
namespace TiledArray {
class Range;
class TiledRange;
class TiledRange1;

template<typename T>
class SparseShape;
} // namespace TiledArray

/////////////////////////////////////////////////////////
// libint FWD
/////////////////////////////////////////////////////////
namespace libint2 {
    struct Atom;
    struct Shell;
}

/////////////////////////////////////////////////////////
// mpqc FWD
/////////////////////////////////////////////////////////
namespace mpqc {
namespace tensor {

template <typename T>
class Tile;

template <typename T>
class DecomposedTensor;

} // namespace tensor

namespace integrals {

template <typename E>
class EnginePool;

} // namespace integrals
} // namespace mpqc

#endif // TCC_COMMON_FWDDECL_H
