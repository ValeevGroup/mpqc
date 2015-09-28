#pragma once
#ifndef TCC_COMMON_FWDDECL_H
#define TCC_COMMON_FWDDECL_H

#pragma GCC diagnostic push
#pragma GCC system_header
#include <tiledarray_fwd.h>
#pragma GCC diagnostic pop

namespace TiledArray {
class Range;
class TiledRange;
class TiledRange1;

template<typename T>
class SparseShape;
} // namespace TiledArray

namespace tcc {

namespace tensor {

template <typename T>
class TilePimpl;

template <typename T>
class Tile;

} // namespace tensor

namespace integrals {

template <typename E>
class EnginePool;

} // namespace integrals

} // namespace tcc

#endif // TCC_COMMON_FWDDECL_H
