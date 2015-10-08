#include "integral_screening_matrices.h"

// #include <memory>
// #include <algorithm>

namespace mpqc {
namespace integrals {

namespace detail {

// Depends on the integrals being 1. DF, 2 obs, 3 obs
ScreeningMatrices
screening_matrix_X(mad::World &world, ShrPool<TwoE_Engine> &engines,
                   basis::Basis const &basis) {

    // Number of clusters
    const auto nclusters = basis.nclusters();

    VectorD cluster_screen_vec = VectorD(nclusters);
    auto sh_screen_vectors = std::vector<std::vector<MatrixD>>(nclusters);

    // Get ptrs to pass to task
    auto const cluster_shells_ptr = &(basis.cluster_shells());
    auto cluster_screen_ptr = &cluster_screen_vec;

    // Capture ptrs by value
    auto cluster_task =
          [=](int64_t ord, std::vector<MatrixD> *shell_screen_ptr) {

        auto const &cluster_shells = *cluster_shells_ptr;
        auto const &cl_shells = cluster_shells[ord];
        const auto nshells = cl_shells.size();

        auto &cvec = *cluster_screen_ptr;
        auto &svec = (*shell_screen_ptr)[0];
        svec = VectorD(nshells);

        auto &eng = engines->local();
        eng.set_precision(0.0);

        for (auto s = 0ul; s < nshells; ++s) {
            auto const &sh = cl_shells[s];
            auto nsh = sh.size();
            const auto *buf
                  = engines->local().compute(sh, unit_shell, sh, unit_shell);

            const auto bmap = Eig::Map<const MatrixD>(buf, nsh, nsh);
            // For each shell get the sqrt of the F norm of the quartet
            svec(s) = std::sqrt(bmap.lpNorm<2>());
        }

        // Q_X(cluster) = |shells|_F
        cvec(ord) = svec.norm();
    };

    // Loop over clusters
    for (auto c = 0; c < nclusters; ++c) {
        sh_screen_vectors[c] = std::vector<MatrixD>(1);
        world.taskq.add(cluster_task, c, &sh_screen_vectors[c]);
    }
    world.gop.fence();

    ScreeningMatrices sc_mats;
    sc_mats.shell_screenings = std::move(sh_screen_vectors);
    sc_mats.cluster_screening = std::move(cluster_screen_vec);

    return sc_mats;
}

ScreeningMatrices
screening_matrix_ab(mad::World &world, ShrPool<TwoE_Engine> &engines,
                    basis::Basis const &basis) {

    const auto nclusters = basis.nclusters();

    MatrixD cluster_screen_mat = MatrixD(nclusters, nclusters);
    auto sh_screening_mats = std::vector<std::vector<MatrixD>>(nclusters);

    // Cluster task
    auto const cluster_shells_ptr = &(basis.cluster_shells());
    auto cl_screen_ptr = &cluster_screen_mat;
    auto cluster_task =
          [=](int64_t c0, std::vector<MatrixD> *c0_sh_screen_mats_ptr) {
        auto &eng = engines->local();
        eng.set_precision(0.);

        auto const &cluster_shells = *cluster_shells_ptr;
        auto &cmat = *cl_screen_ptr;
        auto &sh_mats = *c0_sh_screen_mats_ptr;

        auto const &shells0 = cluster_shells[c0];
        auto const &nshells0 = shells0.size();
        for (auto c1 = 0; c1 < nclusters; ++c1) {
            auto const &shells1 = cluster_shells[c1];
            const auto nshells1 = shells1.size();

            auto &sh_mat = sh_mats[c1];
            sh_mat = MatrixD(nshells0, nshells1);

            for (auto s0 = 0ul; s0 < nshells0; ++s0) {
                auto const &sh0 = shells0[s0];
                const auto nsh0 = sh0.size();

                for (auto s1 = 0ul; s1 < nshells1; ++s1) {
                    auto const &sh1 = shells1[s1];
                    const auto nsh1 = sh1.size();

                    const auto *buf = eng.compute(sh0, sh1, sh0, sh1);
                    const auto bmap = Eig::Map<const MatrixD>(buf, nsh0 * nsh1,
                                                              nsh0 * nsh1);

                    sh_mat(s0, s1) = std::sqrt(bmap.lpNorm<2>());
                }
            }
            cmat(c0, c1) = sh_mat.lpNorm<2>();
        }
    };

    for (auto c0 = 0; c0 < nclusters; ++c0) {
        sh_screening_mats[c0] = std::vector<MatrixD>(nclusters);
        auto current_shell_screen_ptr = &sh_screening_mats[c0];
        world.taskq.add(cluster_task, c0, current_shell_screen_ptr);
    }
    world.gop.fence();

    ScreeningMatrices sc_mats;
    sc_mats.shell_screenings = std::move(sh_screening_mats);
    sc_mats.cluster_screening = std::move(cluster_screen_mat);

    return sc_mats;
}

} // namespace detail
} // namespace integrals
} // namespace mpqc
