#include "mpqc/chemistry/qc/integrals/task_integral_kernels.h"

#include <TiledArray/tensor/tensor_map.h>


namespace mpqc {
namespace lcao {
namespace gaussian {

namespace detail {

double integral_engine_precision = 0.0;

TA::TensorD
integral_kernel(Engine &eng, TA::Range &&rng,
                std::array<ShellVec const *, 2> shell_ptrs, Screener &) {

    set_eng_precision(eng);

    auto const &lobound = rng.lobound();
    std::array<std::size_t, 2> lb = {{lobound[0], lobound[1]}};
    std::array<std::size_t, 2> ub = lb;

    auto tile = TA::TensorD(std::move(rng), 0.0);

    // this makes a map we can resize later.
    double dummy = 0.0;
    auto map = TA::make_const_map(&dummy, {0, 0}, {1, 1});

    const auto& ints_shell_sets = eng.results();

    auto const &sh0 = *shell_ptrs[0];
    auto const &sh1 = *shell_ptrs[1];
    const auto end0 = sh0.size();
    const auto end1 = sh1.size();

    for (auto idx0 = 0ul; idx0 < end0; ++idx0) {
        auto const &s0 = sh0[idx0];
        const auto ns0 = s0.size();
        ub[0] += ns0;

        lb[1] = ub[1] = lobound[1];
        for (auto idx1 = 0ul; idx1 < end1; ++idx1) {
            auto const &s1 = sh1[idx1];
            const auto ns1 = s1.size();
            ub[1] += ns1;

            shell_set(eng, s0, s1);
            assert(ints_shell_sets.size() == 1 &&
                   "integral_kernel can't handle multi-shell-set engines");
            if (ints_shell_sets[0] != nullptr) {
              TA::remap(map, ints_shell_sets[0], lb, ub);
              tile.block(lb, ub) = map;
            }

            lb[1] = ub[1];
        }
        lb[0] = ub[0];
    }

    return tile;
}

TA::TensorD
integral_kernel(Engine &eng, TA::Range &&rng,
                std::array<ShellVec const *, 3> shell_ptrs, Screener &screen) {

    eng.set_precision(integral_engine_precision);

    auto const &lobound = rng.lobound();
    std::array<std::size_t, 3> lb = {{lobound[0], lobound[1], lobound[2]}};
    std::array<std::size_t, 3> ub = lb;

    auto tile = TA::TensorD(std::move(rng), 0.0);

    // this makes a map we can resize later.
    double dummy = 0.0;
    auto map = TA::make_const_map(&dummy, {0, 0, 0}, {1, 1, 1});

    const auto& ints_shell_sets = eng.results();

    auto const &sh0 = *shell_ptrs[0];
    auto const &sh1 = *shell_ptrs[1];
    auto const &sh2 = *shell_ptrs[2];
    const auto end0 = sh0.size();
    const auto end1 = sh1.size();
    const auto end2 = sh2.size();

    for (auto idx0 = 0ul; idx0 < end0; ++idx0) {
        auto const &s0 = sh0[idx0];
        const auto ns0 = s0.size();
        const auto lb0 = lb[0];
        ub[0] += ns0;

        if (!screen.skip(lb0)) {
            lb[1] = ub[1] = lobound[1];
            for (auto idx1 = 0ul; idx1 < end1; ++idx1) {
                auto const &s1 = sh1[idx1];
                const auto ns1 = s1.size();
                const auto lb1 = lb[1];
                ub[1] += ns1;

                if (!screen.skip(lb0, lb1)) {
                    lb[2] = ub[2] = lobound[2];
                    for (auto idx2 = 0ul; idx2 < end2; ++idx2) {
                        auto const &s2 = sh2[idx2];
                        const auto ns2 = s2.size();
                        const auto lb2 = lb[2];
                        ub[2] += ns2;

                        if (!screen.skip(lb0, lb1, lb2)) {
                          shell_set(eng, s0, s1, s2);
                          assert(ints_shell_sets.size() == 1 &&
                                 "integral_kernel can't handle multi-shell-set engines");
                          if (ints_shell_sets[0] != nullptr) {
                            TA::remap(map, ints_shell_sets[0], lb, ub);
                            tile.block(lb, ub) = map;
                          }
                        }

                        lb[2] = ub[2];
                    } // end sh2 for
                }     // end 2 shell screen

                lb[1] = ub[1];
            } // end sh1 for
        }     // end 1 shell screen

        lb[0] = ub[0];
    }

    return tile;
}

TA::TensorD
integral_kernel(Engine &eng, TA::Range &&rng,
                std::array<ShellVec const *, 4> shell_ptrs, Screener &screen) {

    eng.set_precision(integral_engine_precision);

    auto const &lobound = rng.lobound();
    std::array<std::size_t, 4> lb
          = {{lobound[0], lobound[1], lobound[2], lobound[3]}};
    std::array<std::size_t, 4> ub = lb;

    auto tile = TA::TensorD(std::move(rng), 0.0);

    // this makes a map we can resize later.
    double dummy = 0.0;
    auto map = TA::make_const_map(&dummy, {0, 0, 0, 0}, {1, 1, 1, 1});

    const auto& ints_shell_sets = eng.results();

    auto const &sh0 = *shell_ptrs[0];
    auto const &sh1 = *shell_ptrs[1];
    auto const &sh2 = *shell_ptrs[2];
    auto const &sh3 = *shell_ptrs[3];
    const auto end0 = sh0.size();
    const auto end1 = sh1.size();
    const auto end2 = sh2.size();
    const auto end3 = sh3.size();

    for (auto idx0 = 0ul; idx0 < end0; ++idx0) {
        auto const &s0 = sh0[idx0];
        const auto ns0 = s0.size();
        const auto lb0 = lb[0];
        ub[0] += ns0;

        if (!screen.skip(lb0)) {

            lb[1] = ub[1] = lobound[1];
            for (auto idx1 = 0ul; idx1 < end1; ++idx1) {
                auto const &s1 = sh1[idx1];
                const auto ns1 = s1.size();
                const auto lb1 = lb[1];
                ub[1] += ns1;

                if (!screen.skip(lb0, lb1)) {

                    lb[2] = ub[2] = lobound[2];
                    for (auto idx2 = 0ul; idx2 < end2; ++idx2) {
                        auto const &s2 = sh2[idx2];
                        const auto ns2 = s2.size();
                        const auto lb2 = lb[2];
                        ub[2] += ns2;


                        if (!screen.skip(lb0, lb1, lb2)) {
                            lb[3] = ub[3] = lobound[3];
                            for (auto idx3 = 0ul; idx3 < end3; ++idx3) {
                                auto const &s3 = sh3[idx3];
                                const auto ns3 = s3.size();
                                const auto lb3 = lb[3];
                                ub[3] += ns3;

                                if (!screen.skip(lb0, lb1, lb2, lb3)) {
                                  shell_set(eng, s0, s1, s2, s3);
                                  assert(ints_shell_sets.size() == 1 &&
                                         "integral_kernel can't handle multi-shell-set engines");
                                  if (ints_shell_sets[0] != nullptr) {
                                    TA::remap(map, ints_shell_sets[0], lb, ub);
                                    tile.block(lb, ub) = map;
                                  }
                                } // screen all

                                lb[3] = ub[3];
                            } // for all sh3
                        }     // screen s0, s1, s2

                        lb[2] = ub[2];
                    } // end sh2 for
                }     // end 2 shell screen

                lb[1] = ub[1];
            } // end sh1 for
        }     // end 1 shell screen

        lb[0] = ub[0];
    }

    return tile;
}

}  // namespace detail
}  // namespace gaussian
}  // namespace lcao
}  // namespace mpqc
