//
// task_integrals_helper.h
//
// Copyright (C) 2015 Drew Lewis
// Maintainer Drew Lewis
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_TASK_INTEGRAL_KERNELS_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_TASK_INTEGRAL_KERNELS_H_

#include <libint2/engine.h>
#include <tiledarray.h>

#include "mpqc/chemistry/qc/lcao/integrals/screening/screen_base.h"
#include "mpqc/chemistry/qc/lcao/integrals/task_integrals_common.h"
#include "mpqc/math/groups/petite_list.h"

namespace mpqc {
namespace lcao {
namespace gaussian {
namespace detail {

extern double integral_engine_precision;

using Engine = libint2::Engine;
using data_pointer = TA::TensorD::pointer;

inline void set_eng_precision(Engine &eng) {
  eng.set_precision(integral_engine_precision);
}

template <typename... Shells>
inline void shell_set(Engine &e, Shells &&... shells) {
  e.compute(std::forward<Shells>(shells)...);
}

TA::TensorD integral_kernel(Engine &eng, TA::Range &&rng,
                            std::array<ShellVec const *, 2> shell_ptrs,
                            Screener &screen, const math::PetiteList &plist);

TA::TensorD integral_kernel(Engine &eng, TA::Range &&rng,
                            std::array<ShellVec const *, 3> shell_ptrs,
                            Screener &screen, const math::PetiteList &plist);

TA::TensorD integral_kernel(Engine &eng, TA::Range &&rng,
                            std::array<ShellVec const *, 4> shell_ptrs,
                            Screener &screen, const math::PetiteList &plist);

template <libint2::Operator libint2_oper>
std::array<TA::TensorD, libint2::operator_traits<libint2_oper>::nopers>
integral_kernel(Engine &eng, TA::Range &&rng,
                std::array<ShellVec const *, 2> shell_ptrs, Screener &screen,
                const math::PetiteList &plist) {
  set_eng_precision(eng);
  auto const &lobound = rng.lobound();

  // get the index of the first basis function in each cluster
  const auto basis_func_offset_cluster0 = lobound[0];
  const auto basis_func_offset_cluster1 = lobound[1];

  // get reference to clusters (ShellVec type)
  auto const &shvec0 = *shell_ptrs[0];
  auto const &shvec1 = *shell_ptrs[1];

  // get the total number of shells in each cluster
  const auto nshells0 = shvec0.size();
  const auto nshells1 = shvec1.size();

  // get the total number of basis functions in the second cluster
  auto ext_c1 = rng.extent_data()[1];

  // create the result tiles
  using ResultType =
      std::array<TA::TensorD, libint2::operator_traits<libint2_oper>::nopers>;
  ResultType result;
  for (auto &tile : result) {
    tile = TA::TensorD(rng, 0.0);
  }

  // grab pointers to the result tiles to make addressing more efficient
  using ResultPtrType =
      std::array<data_pointer,
                 libint2::operator_traits<libint2_oper>::nopers>;
  ResultPtrType array_of_tile_ptr0;
  const auto nopers = libint2::operator_traits<libint2_oper>::nopers;
  for (auto op = 0u; op != nopers; ++op) {
    array_of_tile_ptr0[op] = result[op].data();
  }

  const auto &ints_shell_sets = eng.results();
  // compute integrals
  auto cf_offset0 = 0u;
  auto bf_offset0 = basis_func_offset_cluster0;
  for (auto idx_sh0 = 0ul; idx_sh0 < nshells0; ++idx_sh0) {
    auto const &sh0 = shvec0[idx_sh0];
    const auto nf0 = sh0.size();

    if (plist.is_canonical(bf_offset0)) {

      auto cf_offset1 = 0u;
      auto bf_offset1 = basis_func_offset_cluster1;
      for (auto idx_sh1 = 0ul; idx_sh1 < nshells1; ++idx_sh1) {
        auto const &sh1 = shvec1[idx_sh1];
        const auto nf1 = sh1.size();

        if (plist.is_canonical(bf_offset0, bf_offset1)) {
          shell_set(eng, sh0, sh1);
          const double permutational_multiplicity =
              plist.multiplicity(bf_offset0, bf_offset1);

          for (auto op = 0u; op != nopers; ++op) {
            if (ints_shell_sets[op] != nullptr) {
              const auto *ints_ptr = ints_shell_sets[op];
              auto *tile_ptr = array_of_tile_ptr0[op];

              for (auto f0 = 0ul, f01 = 0ul; f0 != nf0; ++f0) {
                const auto cf1_offset = (f0 + cf_offset0) * ext_c1;
                for (auto f1 = 0ul; f1 != nf1; ++f1, ++f01) {
                  const auto cf01 = cf1_offset + f1 + cf_offset1;
                  const auto mult_int = permutational_multiplicity * ints_ptr[f01];
                  tile_ptr[cf01] = mult_int;
                }
              }
            }
          }
        }

        bf_offset1 += nf1;
        cf_offset1 += nf1;
      }
    }

    bf_offset0 += nf0;
    cf_offset0 += nf0;
  }

  return result;
};

}  // namespace detail
}  // namespace gaussian
}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_TASK_INTEGRAL_KERNELS_H_
