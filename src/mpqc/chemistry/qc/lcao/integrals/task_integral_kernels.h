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
  std::array<std::size_t, 2> lb = {{lobound[0], lobound[1]}};
  std::array<std::size_t, 2> ub = lb;

  // get the index of the first basis function in each cluster
  const auto func_offset_cluster0 = lb[0];
  const auto func_offset_cluster1 = lb[1];

//  auto tile = TA::TensorD(std::move(rng), 0.0);

  // get reference to clusters (ShellVec type)
  auto const &shvec0 = *shell_ptrs[0];
  auto const &shvec1 = *shell_ptrs[1];

  // get the total number of shells in each cluster
  const auto nshells0 = shvec0.size();
  const auto nshells1 = shvec1.size();

  // get the total number of basis functions in the second cluster
  auto ext_c0 = rng.extent_data()[0];
  auto ext_c1 = rng.extent_data()[1];

  // create the result tiles
  using ResultType =
      std::array<TA::TensorD, libint2::operator_traits<libint2_oper>::nopers>;
  ResultType result;
  for (auto &tile : result) {
    tile = TA::TensorD(rng, 0.0);
  }
//  result.fill(tile);

  // grab pointers to the result tiles to make addressing more efficient
  using ResultPtrType =
      std::array<data_pointer,
                 libint2::operator_traits<libint2_oper>::nopers>;
  ResultPtrType array_of_tile_ptr0;
  const auto nopers = libint2::operator_traits<libint2_oper>::nopers;
  for (auto op = 0u; op != nopers; ++op) {
    array_of_tile_ptr0[op] = result[op].data();
  }

  // test
  RowMatrixXd result0_eig(ext_c0, ext_c1);

  const auto &ints_shell_sets = eng.results();
  // compute integrals
  auto cf_offset0 = 0u;
  auto bf_offset0 = func_offset_cluster0;
  for (auto idx_sh0 = 0ul; idx_sh0 < nshells0; ++idx_sh0) {
    auto const &sh0 = shvec0[idx_sh0];
    const auto nf0 = sh0.size();
    const auto lb0 = lb[0];
    ub[0] += nf0;

    if (plist.is_canonical(bf_offset0)) {
      lb[1] = ub[1] = lobound[1];

      // test
//      std::cout << "\nsh0 = " << idx_sh0 << std::endl;
      auto cf_offset1 = 0u;
      auto bf_offset1 = func_offset_cluster1;
      for (auto idx_sh1 = 0ul; idx_sh1 < nshells1; ++idx_sh1) {
        auto const &sh1 = shvec1[idx_sh1];
        const auto nf1 = sh1.size();
        const auto lb1 = lb[1];
        ub[1] += nf1;

        // test
//        std::cout << "\tsh1 = " << idx_sh1 << std::endl;
        if (plist.is_canonical(bf_offset0, bf_offset1)) {
          shell_set(eng, sh0, sh1);
          const double permutational_multiplicity =
              plist.multiplicity(bf_offset0, bf_offset1);

          for (auto op = 0u; op != nopers; ++op) {
            if (ints_shell_sets[op] != nullptr) {
              const auto *ints_ptr = ints_shell_sets[op];
              auto *tile_ptr = array_of_tile_ptr0[op];

              // test
//              std::cout << "\t\toper = " << op << std::endl;
              for (auto f0 = 0ul, f01 = 0ul; f0 != nf0; ++f0) {
                const auto cf1_offset = (f0 + cf_offset0) * ext_c1;
                for (auto f1 = 0ul; f1 != nf1; ++f1, ++f01) {
                  const auto cf01 = cf1_offset + f1 + cf_offset1;
                  const auto mult_int = permutational_multiplicity * ints_ptr[f01];
                  tile_ptr[cf01] = mult_int;

                  // test
//                  if (op == 0) {
//                    std::cout << "\t\t\tcf01 = " << cf01
//                              << ", multiplicity_times_int = " << mult_int
//                              << ", SpheMM = " << tile_ptr[cf01]
//                              << ", array_of_tile_ptr0 = " << array_of_tile_ptr0[op][cf01]
//                              << std::endl;
//
//                    result0_eig(f0 + cf_offset0, f1 + cf_offset1) = mult_int;
//                  }
                }
              }

              // test
//              for (auto f0 = 0ul; f0 != nf0; ++f0) {
//                const auto cf1_offset = (f0 + cf_offset0) * ext_c1;
//                for (auto f1 = 0ul; f1 != nf1; ++f1) {
//                  const auto cf01 = cf1_offset + f1 + cf_offset1;
//                  if (op == 0) {
//                    std::cout << "\n\t\t\tcf01 = " << cf01
//                              << ", reprint array_of_tile_ptr0 = " << array_of_tile_ptr0[op][cf01]
//                              << std::endl;
//                  }
//
//                }
//              }

//              auto shell_ord = 0ul;
//              const auto lb0 = lb[0] - func_offset_cluster0;
//              const auto ub0 = ub[0] - func_offset_cluster0;
//              const auto lb1 = lb[1] - func_offset_cluster1;
//              const auto ub1 = ub[1] - func_offset_cluster1;
//              for (auto el0 = lb0; el0 < ub0; ++el0) {
//                const auto el0_pr e = ext_c1 * el0;
//                data_pointer MADNESS_RESTRICT tile_ptr =
//                    tile_ptr_first + el0_pre;
//
//                for (auto el1 = lb1; el1 < ub1; ++el1, ++shell_ord) {
//                  tile_ptr[el1] =
//                      permutational_multiplicity * ints_ptr[shell_ord];
//                }
//              }

            }
          }
        }
        // test
//        for (auto f0 = 0ul; f0 != nf0; ++f0) {
//          const auto cf1_offset = (f0 + cf_offset0) * ext_c1;
//          for (auto f1 = 0ul; f1 != nf1; ++f1) {
//            const auto cf01 = cf1_offset + f1 + cf_offset1;
//            std::cout << "\n\t\t\tcf01 = " << cf01
//                      << ", re-re-print array_of_tile_ptr0 = " << array_of_tile_ptr0[0][cf01]
//                      << std::endl;
//          }
//        }

        lb[1] = ub[1];
        bf_offset1 += nf1;
        cf_offset1 += nf1;
      }
    }

    lb[0] = ub[0];
    bf_offset0 += nf0;
    cf_offset0 += nf0;
  }

  // test
//  std::cout << "\nresult[oper=0] = \n" << result0_eig << std::endl;

  // test
//  auto volume = ext_c0 * ext_c1;
//  for (auto op = 0; op != nopers; ++op) {
//    if (op == 0) {
//      for (auto ord = 0u; ord != volume; ++ord) {
//        std::cout << "ord = " << ord
//                  << ", array_of_tile_ptr0 = " << array_of_tile_ptr0[op][ord]
//                  << ", result = " << result[op][ord] << std::endl;
//      }
//    }
//  }

  return result;
};

}  // namespace detail
}  // namespace gaussian
}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_TASK_INTEGRAL_KERNELS_H_
