#pragma once

#include "mpqc/math/external/eigen/eigen.h"
#include <tiledarray.h>
#include "../common/typedefs.h"
#include "mpqc/math/external/eigen/eigen.h"
#include <array>
#include <iostream>

namespace mpqc {
namespace TT {

RowMatrixXd one_over_delta(Eigen::VectorXd const &e_vals, int64_t occ) {
    const auto size = e_vals.size();
    const auto vir = size - occ;
    Eigen::MatrixXd delta(occ * vir, occ * vir);
    delta.setZero();

    for (auto i = 0; i < occ; ++i) {
        const auto iv = -e_vals[i];

        for (auto a = 0; a < vir; ++a) {
            const auto ia = i * vir + a;
            const auto av = e_vals[a + occ];

            for (auto j = 0; j < occ; ++j) {
                const auto jv = -e_vals[j];

                for (auto b = 0; b < vir; ++b) {
                    const auto jb = j * vir + b;

                    const auto val = iv + av + jv + e_vals[b + occ];
                    delta(ia, jb) = 1.0 / val;
                }
            }
        }
    }

    return delta;
}

std::array<RowMatrixXd, 4>
TT_TA_ARRAY(TA::DistArray<TA::TensorD, TA::SparsePolicy> const &A) {
    std::cout << "Starting Amplitude TT decomp\n";
    auto A_extent = A.trange().elements_range().extent();
    RowMatrixXd A_eig(A_extent[0] * A_extent[1], A_extent[2] * A_extent[3]);
    A_eig.setZero();

    auto occ = A_extent[0];
    auto vir = A_extent[1];

    std::cout << "Occ = \n" << occ << std::endl;
    std::cout << "Vir = \n" << vir << std::endl;

    for (auto it : A) {
        TA::TensorD const &t = it.get();

        auto lo = t.range().lobound();
        auto up = t.range().upbound();
        auto ext = t.range().extent();

        for (auto i = lo[0]; i < up[0]; ++i) {
            auto i_ord = i * A_extent[1] * A_extent[2] * A_extent[3];

            for (auto j = lo[1]; j < up[1]; ++j) {
                auto ij_ord = i_ord + j * A_extent[2] * A_extent[3];

                for (auto k = lo[2]; k < up[2]; ++k) {
                    auto ijk_ord = ij_ord + A_extent[3] * k;

                    for (auto l = lo[3]; l < up[3]; ++l) {
                        auto ijkl_ord = ijk_ord + l;

                        *(A_eig.data() + ijkl_ord) = t(i,j,k,l);
                    }
                }
            }
        }
    }

    RowMatrixXd G0, G1, G2, G3;

    // Resize D (i, ajb)
    A_eig.resize(occ, vir * occ * vir);

    Eigen::JacobiSVD<RowMatrixXd> svd(A_eig, Eigen::ComputeThinU | Eigen::ComputeThinV);
    svd.setThreshold(1e-3);


    // G0 = (i, r0)
    const auto rank0 = svd.rank();
    std::cout << "Rank0 = " << rank0 << std::endl;
    G0 = svd.matrixU().leftCols(rank0);

    // R = (r0, ajb)
    RowMatrixXd sig = RowMatrixXd(svd.singularValues().asDiagonal())
                        .block(0, 0, rank0, rank0);
    RowMatrixXd R0 = sig * svd.matrixV().leftCols(rank0).transpose();

    // Resize R(r0a, jb)
    R0.resize(rank0 * vir, occ * vir);
    svd.compute(R0);
    const auto rank1 = svd.rank();
    std::cout << "Rank1 = " << rank1 << std::endl;

    // G1 = (r0a, r1)
    G1 = svd.matrixU().leftCols(rank1);

    // R1 = (r1, jb)
    sig = RowMatrixXd(svd.singularValues().asDiagonal()).block(0, 0, rank1, rank1);
    RowMatrixXd R1 = sig * svd.matrixV().leftCols(rank1).transpose();

    // Resize R1(r1j, b);
    R1.resize(rank1 * occ, vir);
    svd.compute(R1);
    const auto rank2 = svd.rank();
    std::cout << "Rank2 = " << rank2 << std::endl;

    // G2 = (r1j,r2)
    G2 = svd.matrixU().leftCols(rank2);

    // G3 = (r2, b)
    sig = RowMatrixXd(svd.singularValues().asDiagonal()).block(0, 0, rank2, rank2);
    G3 = sig * svd.matrixV().leftCols(rank2).transpose();

    return {{G0, G1, G2, G3}};
}

std::array<RowMatrixXd, 4>
TT_one_over_delta(Eigen::VectorXd const &e_vals, int64_t occ) {
    const auto size = e_vals.size();
    const auto vir = size - occ;

    // D(ia,jb)
    auto D = one_over_delta(e_vals, occ);

    RowMatrixXd G0, G1, G2, G3;

    // Resize D (i, ajb)
    D.resize(occ, vir * occ * vir);

    Eigen::JacobiSVD<RowMatrixXd> svd(D, Eigen::ComputeThinU | Eigen::ComputeThinV);
    svd.setThreshold(1e-3);


    // G0 = (i, r0)
    const auto rank0 = svd.rank();
    std::cout << "Rank0 = " << rank0 << std::endl;
    G0 = svd.matrixU().leftCols(rank0);

    // R = (r0, ajb)
    RowMatrixXd sig = RowMatrixXd(svd.singularValues().asDiagonal())
                        .block(0, 0, rank0, rank0);
    RowMatrixXd R0 = sig * svd.matrixV().leftCols(rank0).transpose();

    // Resize R(r0a, jb)
    R0.resize(rank0 * vir, occ * vir);
    svd.compute(R0);
    const auto rank1 = svd.rank();
    std::cout << "Rank1 = " << rank1 << std::endl;

    // G1 = (r0a, r1)
    G1 = svd.matrixU().leftCols(rank1);

    // R1 = (r1, jb)
    sig = RowMatrixXd(svd.singularValues().asDiagonal()).block(0, 0, rank1, rank1);
    RowMatrixXd R1 = sig * svd.matrixV().leftCols(rank1).transpose();

    // Resize R1(r1j, b);
    R1.resize(rank1 * occ, vir);
    svd.compute(R1);
    const auto rank2 = svd.rank();
    std::cout << "Rank2 = " << rank2 << std::endl;

    // G2 = (r1j,r2)
    G2 = svd.matrixU().leftCols(rank2);

    // G3 = (r2, b)
    sig = RowMatrixXd(svd.singularValues().asDiagonal()).block(0, 0, rank2, rank2);
    G3 = sig * svd.matrixV().leftCols(rank2).transpose();

    return {{G0, G1, G2, G3}};
}

TA::DistArray<TA::TensorD, TA::SparsePolicy>
make_Wir0k(madness::World &world, RowMatrixXd const &G0, RowMatrixXd const &W,
           TA::TiledRange1 const &occ_trange, int64_t occ, int64_t rank0) {

    std::vector<int64_t> rvec = {0, rank0};
    TA::TiledRange1 r_trange1(rvec.begin(), rvec.end());
    TA::TiledRange trange({occ_trange, r_trange1, occ_trange});

    TA::TensorD norms(trange.tiles_range(), std::numeric_limits<float>::max());
    TA::SparseShape<float> shape(norms, trange);

    TA::DistArray<TA::TensorD, TA::SparsePolicy> Wir0k(world, trange, shape);

    auto task = [&](TA::Range const &range, int64_t ord) {
        TA::TensorD tile(range, 0.0);
        auto lo_bound = range.lobound();
        auto hi_bound = range.upbound();

        for (auto i = lo_bound[0]; i < hi_bound[0]; ++i) {
            for (auto r = lo_bound[1]; r < hi_bound[1]; ++r) {
                for (auto k = lo_bound[2]; k < hi_bound[2]; ++k) {
                    for (auto z = 0; z < occ; ++z) {
                        tile(i, r, k) += W(z, i) * G0(z, r) * W(z, k);
                    }
                }
            }
        }

        Wir0k.set(ord, tile);
    };

    const auto vol = trange.tiles_range().volume();
    for (auto i = 0; i < vol; ++i) {
        world.taskq.add(task, trange.make_tile_range(i), i);
    }
    world.gop.fence();
    Wir0k.truncate();

    return Wir0k;
}

TA::DistArray<TA::TensorD, TA::SparsePolicy>
make_Ycr2b(madness::World &world, RowMatrixXd const &G3, RowMatrixXd const &Cv,
           RowMatrixXd const &Q, TA::TiledRange1 const &pao_trange, int64_t vir,
           int64_t rank2) {

    std::vector<int64_t> rvec = {0, rank2};
    TA::TiledRange1 r_trange1(rvec.begin(), rvec.end());
    TA::TiledRange trange({pao_trange, r_trange1, pao_trange});

    TA::TensorD norms(trange.tiles_range(), std::numeric_limits<float>::max());
    TA::SparseShape<float> shape(norms, trange);

    TA::DistArray<TA::TensorD, TA::SparsePolicy> Ycr2b(world, trange, shape);

    auto task = [&](TA::Range const &range, int64_t ord) {
        TA::TensorD tile(range, 0.0);
        auto lo_bound = range.lobound();
        auto hi_bound = range.upbound();

        for (auto c = lo_bound[0]; c < hi_bound[0]; ++c) {
            for (auto r = lo_bound[1]; r < hi_bound[1]; ++r) {
                for (auto b = lo_bound[2]; b < hi_bound[2]; ++b) {
                    for (auto z = 0; z < vir; ++z) {
                        tile(c, r, b) += Cv(c, z) * G3(r, z) * Q(z, b);
                    }
                }
            }
        }

        Ycr2b.set(ord, tile);
    };

    const auto vol = trange.tiles_range().volume();
    for (auto i = 0; i < vol; ++i) {
        world.taskq.add(task, trange.make_tile_range(i), i);
    }
    world.gop.fence();
    Ycr2b.truncate();

    return Ycr2b;
}

TA::DistArray<TA::TensorD, TA::SparsePolicy>
make_Wr1lr2j(madness::World &world, RowMatrixXd const &G2, RowMatrixXd const &W,
             TA::TiledRange1 const &occ_trange, int64_t occ, int64_t rank1,
             int64_t rank2) {

    std::vector<int64_t> rvec1 = {0, rank1};
    std::vector<int64_t> rvec2 = {0, rank2};
    TA::TiledRange1 r1_trange1(rvec1.begin(), rvec1.end());
    TA::TiledRange1 r2_trange1(rvec2.begin(), rvec2.end());
    TA::TiledRange trange({r1_trange1, occ_trange, r2_trange1, occ_trange});

    TA::TensorD norms(trange.tiles_range(), std::numeric_limits<float>::max());
    TA::SparseShape<float> shape(norms, trange);

    TA::DistArray<TA::TensorD, TA::SparsePolicy> Wr1lr2j(world, trange, shape);

    auto task = [&](TA::Range const &range, int64_t ord) {
        TA::TensorD tile(range, 0.0);
        auto lo_bound = range.lobound();
        auto hi_bound = range.upbound();

        for (auto r1 = lo_bound[0]; r1 < hi_bound[0]; ++r1) {
            for (auto l = lo_bound[1]; l < hi_bound[1]; ++l) {
                for (auto r2 = lo_bound[2]; r2 < hi_bound[2]; ++r2) {
                    for (auto j = lo_bound[3]; j < hi_bound[3]; ++j) {
                        for (auto z = 0; z < occ; ++z) {
                            auto r1z = r1 * occ + z;
                            tile(r1, l, r2, j) += W(z, l) * G2(r1z, r2)
                                                  * W(z, j);
                        }
                    }
                }
            }
        }

        Wr1lr2j.set(ord, tile);
    };

    const auto vol = trange.tiles_range().volume();
    for (auto i = 0; i < vol; ++i) {
        world.taskq.add(task, trange.make_tile_range(i), i);
    }
    world.gop.fence();
    Wr1lr2j.truncate();

    return Wr1lr2j;
}

TA::DistArray<TA::TensorD, TA::SparsePolicy>
make_Yr0ar1d(madness::World &world, RowMatrixXd const &G1, RowMatrixXd const &Cv,
             RowMatrixXd const &Q, TA::TiledRange1 const &pao_trange, int64_t vir,
             int64_t rank0, int64_t rank1) {

    std::vector<int64_t> rvec0 = {0, rank0};
    std::vector<int64_t> rvec1 = {0, rank1};
    TA::TiledRange1 r0_trange1(rvec0.begin(), rvec0.end());
    TA::TiledRange1 r1_trange1(rvec1.begin(), rvec1.end());
    TA::TiledRange trange({r0_trange1, pao_trange, r1_trange1, pao_trange});

    TA::TensorD norms(trange.tiles_range(), std::numeric_limits<float>::max());
    TA::SparseShape<float> shape(norms, trange);

    TA::DistArray<TA::TensorD, TA::SparsePolicy> Yr0ar1d(world, trange, shape);

    auto task = [&](TA::Range const &range, int64_t ord) {
        TA::TensorD tile(range, 0.0);
        auto lo_bound = range.lobound();
        auto hi_bound = range.upbound();

        for (auto r0 = lo_bound[0]; r0 < hi_bound[0]; ++r0) {
            for (auto d = lo_bound[1]; d < hi_bound[1]; ++d) {
                for (auto r1 = lo_bound[2]; r1 < hi_bound[2]; ++r1) {
                    for (auto a = lo_bound[3]; a < hi_bound[3]; ++a) {
                        for (auto z = 0; z < vir; ++z) {
                            auto r0z = r0 * vir + z;
                            tile(r0, d, r1, a) += Cv(d, z) * G1(r0z, r1)
                                                  * Q(z, a);
                        }
                    }
                }
            }
        }

        Yr0ar1d.set(ord, tile);
    };

    const auto vol = trange.tiles_range().volume();
    for (auto i = 0; i < vol; ++i) {
        world.taskq.add(task, trange.make_tile_range(i), i);
    }
    world.gop.fence();
    Yr0ar1d.truncate();

    return Yr0ar1d;
}

std::array<TA::DistArray<TA::TensorD, TA::SparsePolicy>, 4>
transform_tensors(Eigen::VectorXd const &e_vals, int64_t occ,
                  madness::World &world, TA::TiledRange1 const &occ_trange,
                  TA::TiledRange1 const &vir_trange, RowMatrixXd const &Uocc,
                  RowMatrixXd const &Cv, RowMatrixXd const &Q) {
    const auto vir = e_vals.size() - occ;

    auto mats = TT_one_over_delta(e_vals, occ);

    const auto rank0 = mats[0].cols();
    auto Wir0k = make_Wir0k(world, mats[0], Uocc, occ_trange, occ, rank0);

    const auto rank1 = mats[1].cols();
    auto Yr0ar1d
          = make_Yr0ar1d(world, mats[1], Cv, Q, vir_trange, vir, rank0, rank1);

    const auto rank2 = mats[2].cols();
    auto Wr1lr2j
          = make_Wr1lr2j(world, mats[2], Uocc, occ_trange, occ, rank1, rank2);

    auto Ycr2b = make_Ycr2b(world, mats[3], Cv, Q, vir_trange, vir, rank2);

    return {{Wir0k, Yr0ar1d, Wr1lr2j, Ycr2b}};
}

} // namespace TT
} // namespace mpqc
