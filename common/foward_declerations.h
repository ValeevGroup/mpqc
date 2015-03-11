#pragma once 
#ifndef TCC_COMMON_FWDDECL_H
#define TCC_COMMON_FWDDECL_H

#pragma GCC diagnostic push
#pragma GCC system_header
#include <tiledarray_fwd.h>
#pragma GCC diagnostic pop

namespace tcc {
namespace tensor {
    template <typename T> 
    class TilePimpl;
} // namespace tensor
} // namespace tcc

#endif // TCC_COMMON_FWDDECL_H
