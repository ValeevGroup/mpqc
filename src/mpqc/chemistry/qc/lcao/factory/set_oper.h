//
// Created by Chong Peng on 2/13/17.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_SET_OPER_H
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_SET_OPER_H

#include <tiledarray.h>

namespace mpqc {
namespace lcao {
namespace gaussian {
namespace detail {

/// interface template function
template <typename Tile>
inline void set_oper(std::function<Tile(TA::TensorD &&)> &op) {}

/// specilazation when Tile=TA::TensorD
template <>
inline void set_oper(std::function<TA::TensorD(TA::TensorD &&)> &op) {
  op = TA::detail::Noop<TA::TensorD, TA::TensorD, true>();
}

}  // namespace detail
}  // namespace gaussian
}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_SET_OPER_H
