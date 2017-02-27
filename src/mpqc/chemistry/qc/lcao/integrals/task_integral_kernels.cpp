#include "mpqc/chemistry/qc/lcao/integrals/task_integral_kernels.h"
#include "mpqc/math/groups/petite_list.h"

#include <TiledArray/tensor/tensor_map.h>

namespace mpqc {
namespace lcao {
namespace gaussian {

namespace detail {

double integral_engine_precision = 0.0;

TA::TensorD integral_kernel(Engine &eng, TA::Range &&rng,
                            std::array<ShellVec const *, 2> shell_ptrs,
                            Screener &, const math::PetiteList& plist) {
  set_eng_precision(eng);

  auto const &lobound = rng.lobound();
  std::array<std::size_t, 2> lb = {{lobound[0], lobound[1]}};
  std::array<std::size_t, 2> ub = lb;

  const auto tile_start0 = lb[0];
  const auto tile_start1 = lb[1];

  auto tile = TA::TensorD(std::move(rng), 0.0);

  const auto &ints_shell_sets = eng.results();

  auto const &sh0 = *shell_ptrs[0];
  auto const &sh1 = *shell_ptrs[1];
  const auto end0 = sh0.size();
  const auto end1 = sh1.size();

  auto ext1 = tile.range().extent_data()[1];
  data_pointer restrict tile_ptr_first = tile.data();

  for (auto idx0 = 0ul; idx0 < end0; ++idx0) {
    auto const &s0 = sh0[idx0];
    const auto ns0 = s0.size();
    const auto lb0 = lb[0];
    ub[0] += ns0;

    if (plist.is_canonical(lb0)) {
      lb[1] = ub[1] = lobound[1];
      for (auto idx1 = 0ul; idx1 < end1; ++idx1) {
        auto const &s1 = sh1[idx1];
        const auto ns1 = s1.size();
        const auto lb1 = lb[1];
        ub[1] += ns1;

        if (plist.is_canonical(lb0, lb1)) {
          shell_set(eng, s0, s1);
          const auto ns01 = ns0 * ns1;
          assert(ints_shell_sets.size() == 1 &&
             "integral_kernel can't handle multi-shell-set engines");
          if (ints_shell_sets[0] != nullptr) {
            const double permutational_multiplicity = plist.multiplicity(lb0,lb1);
            const auto ints_ptr = ints_shell_sets[0];
            auto shell_ord = 0ul;
            const auto lb0 = lb[0] - tile_start0;
            const auto ub0 = ub[0] - tile_start0;
            const auto lb1 = lb[1] - tile_start1;
            const auto ub1 = ub[1] - tile_start1;
            for (auto el0 = lb0; el0 < ub0; ++el0) {
              const auto el0_pre = ext1 * el0;
              data_pointer restrict tile_ptr = tile_ptr_first + el0_pre;

              for (auto el1 = lb1; el1 < ub1; ++el1, ++shell_ord) {
                tile_ptr[el1] = permutational_multiplicity * ints_ptr[shell_ord];
              }
            }
          }
        }
        lb[1] = ub[1];
      }
    }

    lb[0] = ub[0];
  }

  return tile;
}

TA::TensorD integral_kernel(Engine &eng, TA::Range &&rng,
                            std::array<ShellVec const *, 3> shell_ptrs,
                            Screener &screen, const math::PetiteList& plist) {
  eng.set_precision(integral_engine_precision);

  auto const &lobound = rng.lobound();
  std::array<std::size_t, 3> lb = {{lobound[0], lobound[1], lobound[2]}};
  std::array<std::size_t, 3> ub = lb;

  const auto tile_start0 = lb[0];
  const auto tile_start1 = lb[1];
  const auto tile_start2 = lb[2];

  auto tile = TA::TensorD(std::move(rng), 0.0);

  const auto &ints_shell_sets = eng.results();

  auto const &sh0 = *shell_ptrs[0];
  auto const &sh1 = *shell_ptrs[1];
  auto const &sh2 = *shell_ptrs[2];
  const auto end0 = sh0.size();
  const auto end1 = sh1.size();
  const auto end2 = sh2.size();

  auto ext = tile.range().extent_data();
  auto ext2 = ext[2];
  auto ext12 = ext[1] * ext[2];
  data_pointer restrict tile_ptr_first = tile.data();

  for (auto idx0 = 0ul; idx0 < end0; ++idx0) {
    auto const &s0 = sh0[idx0];
    const auto ns0 = s0.size();
    const auto lb0 = lb[0];
    ub[0] += ns0;

    if (plist.is_canonical(lb0) && !screen.skip(lb0)) {
      lb[1] = ub[1] = lobound[1];
      for (auto idx1 = 0ul; idx1 < end1; ++idx1) {
        auto const &s1 = sh1[idx1];
        const auto ns1 = s1.size();
        const auto lb1 = lb[1];
        ub[1] += ns1;

        if (plist.is_canonical(lb0, lb1) && !screen.skip(lb0, lb1)) {
          lb[2] = ub[2] = lobound[2];
          const auto ns01 = ns0 * ns1;
          for (auto idx2 = 0ul; idx2 < end2; ++idx2) {
            auto const &s2 = sh2[idx2];
            const auto ns2 = s2.size();
            const auto lb2 = lb[2];
            ub[2] += ns2;

            if (plist.is_canonical(lb0, lb1, lb2) && !screen.skip(lb0, lb1, lb2)) {
              shell_set(eng, s0, s1, s2);
              const auto ns012 = ns01 * ns2;
              assert(ints_shell_sets.size() == 1 &&
                     "integral_kernel can't handle multi-shell-set engines");
              if (ints_shell_sets[0] != nullptr) {
                const double permutational_multiplicity = plist.multiplicity(lb0,lb1,lb2);
                const auto ints_ptr = ints_shell_sets[0];
                auto shell_ord = 0ul;
                const auto lb0 = lb[0] - tile_start0;
                const auto ub0 = ub[0] - tile_start0;
                const auto lb1 = lb[1] - tile_start1;
                const auto ub1 = ub[1] - tile_start1;
                const auto lb2 = lb[2] - tile_start2;
                const auto ub2 = ub[2] - tile_start2;
                for (auto el0 = lb0; el0 < ub0; ++el0) {
                  const auto el0_pre = ext12 * el0;
                  data_pointer restrict tile_ptr0 = tile_ptr_first + el0_pre;
                  for (auto el1 = lb1; el1 < ub1; ++el1) {
                    const auto el1_pre = ext2 * el1;
                    data_pointer restrict tile_ptr1 = tile_ptr0 + el1_pre;
                    for (auto el2 = lb2; el2 < ub2; ++el2, ++shell_ord) {
                      tile_ptr1[el2] = permutational_multiplicity * ints_ptr[shell_ord];
                    }
                  }
                }
              }
            }

            lb[2] = ub[2];
          }  // end sh2 for
        }    // end 2 shell screen

        lb[1] = ub[1];
      }  // end sh1 for
    }    // end 1 shell screen

    lb[0] = ub[0];
  }

  return tile;
}

TA::TensorD integral_kernel(Engine &eng, TA::Range &&rng,
                            std::array<ShellVec const *, 4> shell_ptrs,
                            Screener &screen, const math::PetiteList& plist) {
  eng.set_precision(integral_engine_precision);

  auto const &lobound = rng.lobound();
  std::array<std::size_t, 4> lb = {
      {lobound[0], lobound[1], lobound[2], lobound[3]}};
  std::array<std::size_t, 4> ub = lb;

  const auto tile_start0 = lb[0];
  const auto tile_start1 = lb[1];
  const auto tile_start2 = lb[2];
  const auto tile_start3 = lb[3];

  auto tile = TA::TensorD(std::move(rng), 0.0);

  const auto &ints_shell_sets = eng.results();

  auto const &sh0 = *shell_ptrs[0];
  auto const &sh1 = *shell_ptrs[1];
  auto const &sh2 = *shell_ptrs[2];
  auto const &sh3 = *shell_ptrs[3];
  const auto end0 = sh0.size();
  const auto end1 = sh1.size();
  const auto end2 = sh2.size();
  const auto end3 = sh3.size();

  auto ext = tile.range().extent_data();
  auto ext3 = ext[3];
  auto ext23 = ext[2] * ext[3];
  auto ext123 = ext[1] * ext[2] * ext[3];
  data_pointer restrict tile_ptr_first = tile.data();

  for (auto idx0 = 0ul; idx0 < end0; ++idx0) {
    auto const &s0 = sh0[idx0];
    const auto ns0 = s0.size();
    const auto lb0 = lb[0];
    ub[0] += ns0;

    if (plist.is_canonical(lb0) && !screen.skip(lb0)) {
      lb[1] = ub[1] = lobound[1];
      for (auto idx1 = 0ul; idx1 < end1; ++idx1) {
        auto const &s1 = sh1[idx1];
        const auto ns1 = s1.size();
        const auto lb1 = lb[1];
        ub[1] += ns1;

        if (plist.is_canonical(lb0, lb1) && !screen.skip(lb0, lb1)) {
          lb[2] = ub[2] = lobound[2];
          const auto ns01 = ns0 * ns1;
          for (auto idx2 = 0ul; idx2 < end2; ++idx2) {
            auto const &s2 = sh2[idx2];
            const auto ns2 = s2.size();
            const auto lb2 = lb[2];
            ub[2] += ns2;

            if (plist.is_canonical(lb0, lb1, lb2) && !screen.skip(lb0, lb1, lb2)) {
              lb[3] = ub[3] = lobound[3];
              const auto ns012 = ns01 * ns2;
              for (auto idx3 = 0ul; idx3 < end3; ++idx3) {
                auto const &s3 = sh3[idx3];
                const auto ns3 = s3.size();
                const auto lb3 = lb[3];
                ub[3] += ns3;

                if (plist.is_canonical(lb0, lb1, lb2, lb3) && !screen.skip(lb0, lb1, lb2, lb3)) {
                  shell_set(eng, s0, s1, s2, s3);
                  const auto ns0123 = ns012 * ns3;
                  assert(
                      ints_shell_sets.size() == 1 &&
                      "integral_kernel can't handle multi-shell-set engines");
                  if (ints_shell_sets[0] != nullptr) {
                    const double permutational_multiplicity = plist.multiplicity(lb0,lb1,lb2,lb3);
                    const auto ints_ptr = ints_shell_sets[0];
                    auto shell_ord = 0ul;
                    const auto lb0 = lb[0] - tile_start0;
                    const auto ub0 = ub[0] - tile_start0;
                    const auto lb1 = lb[1] - tile_start1;
                    const auto ub1 = ub[1] - tile_start1;
                    const auto lb2 = lb[2] - tile_start2;
                    const auto ub2 = ub[2] - tile_start2;
                    const auto lb3 = lb[3] - tile_start3;
                    const auto ub3 = ub[3] - tile_start3;
                    for (auto el0 = lb0; el0 < ub0; ++el0) {
                      const auto el0_pre = ext123 * el0;
                      data_pointer restrict tile_ptr0 = tile_ptr_first + el0_pre;
                      for (auto el1 = lb1; el1 < ub1; ++el1) {
                        const auto el1_pre = ext23 * el1;
                        data_pointer restrict tile_ptr1 = tile_ptr0 + el1_pre;
                        for (auto el2 = lb2; el2 < ub2; ++el2) {
                          const auto el2_pre = ext3 * el2;
                          data_pointer restrict tile_ptr2 = tile_ptr1 + el2_pre;
                          for (auto el3 = lb3; el3 < ub3; ++el3, ++shell_ord) {
                            tile_ptr2[el3] = permutational_multiplicity * ints_ptr[shell_ord];
                          }
                        }
                      }
                    }
                  }
                }

                lb[3] = ub[3];
              }  // for all sh3
            }    // Skip from estimation with shells 1, 2, and 3

            lb[2] = ub[2];
          }  // end sh2 for
        }    // Skip from estimation with shell 0 and shell 1

        lb[1] = ub[1];
      }  // end sh1 for
    }    // Skip from estimation with shell 0

    lb[0] = ub[0];
  }

  return tile;
}

}  // namespace detail
}  // namespace gaussian
}  // namespace lcao
}  // namespace mpqc
