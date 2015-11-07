#include "../include/eigen.h"
#include "../common/typedefs.h"

#include "../include/libint.h"
#include "../include/tiledarray.h"

#include "../utility/make_array.h"
#include "../clustering/kmeans.h"

#include "../molecule/atom.h"
#include "../molecule/atom_based_cluster.h"
#include "../molecule/molecule.h"
#include "../molecule/clustering_functions.h"
#include "../molecule/make_clusters.h"

#include "../basis/atom_basisset.h"
#include "../basis/basis_set.h"
#include "../basis/basis.h"

#include "../integrals/integrals.h"

#include "../ta_routines/array_to_eigen.h"

#include <tbb.h>
#include <libint2/diis.h>
#include <iostream>
#include <array>
#include <iomanip>
#include <string>

using namespace mpqc;
namespace bas = mpqc::basis;
namespace ints = mpqc::integrals;
namespace L2 = libint2;

using Matrix = Eig::Matrix<double, Eig::Dynamic, Eig::Dynamic, Eig::RowMajor>;

Matrix compute_soad(molecule::Molecule const &);

template <typename ShrPool>
Matrix
compute_2body_fock_from_soad(ShrPool &engs, molecule::Molecule const &mol,
                             basis::Basis const &obs, Matrix const &H,
                             Matrix const &D,
                             double precision
                             = std::numeric_limits<double>::epsilon());

template <typename ShrPool>
Matrix
compute_2body_G(ShrPool &engs, basis::Basis const &obs, Matrix const &Q,
                Matrix const &D,
                double precision = std::numeric_limits<double>::epsilon());

Matrix update_density(Matrix const &F, Matrix const &S, int64_t occ) {
    Eig::GeneralizedSelfAdjointEigenSolver<Matrix> es(F, S);
    Matrix C = es.eigenvectors().leftCols(occ);
    return C * C.transpose();
}

int main(int argc, char *argv[]) {
    auto &world = madness::initialize(argc, argv);
    std::string mol_file = "";
    std::string basis_name = "";
    int nclusters = 1;
    std::cout << std::setprecision(15);
    double threshold = 1e-10;
    if (argc == 3) {
        mol_file = argv[1];
        basis_name = argv[2];
    } else {
        std::cout << "input is $./program mol_file basis_file";
        return 0;
    }
    TiledArray::SparseShape<float>::threshold(threshold);

    auto clustered_mol = molecule::attach_hydrogens_and_kmeans(
          molecule::read_xyz(mol_file).clusterables(), nclusters);

    auto repulsion_energy = clustered_mol.nuclear_repulsion();
    std::cout << "Nuclear Repulsion Energy: " << repulsion_energy << std::endl;
    auto occ = clustered_mol.occupation(0);

    basis::BasisSet bs(basis_name);
    basis::Basis basis(bs.get_cluster_shells(clustered_mol));

    libint2::init();

    Matrix S, H;
    {
        const auto bs_array = tcc::utility::make_array(basis, basis);

        auto overlap_e
              = ints::make_1body_shr_pool("overlap", basis, clustered_mol);
        auto S_TA = ints::sparse_integrals(world, overlap_e, bs_array);

        auto kinetic_e
              = ints::make_1body_shr_pool("kinetic", basis, clustered_mol);
        auto T = ints::sparse_integrals(world, kinetic_e, bs_array);

        auto nuclear_e
              = ints::make_1body_shr_pool("nuclear", basis, clustered_mol);
        auto V = ints::sparse_integrals(world, nuclear_e, bs_array);

        decltype(T) H_TA;
        H_TA("i,j") = T("i,j") + V("i,j");

        S = tcc::array_ops::array_to_eigen(S_TA);
        H = tcc::array_ops::array_to_eigen(H_TA);
    }

    auto eri_e = ints::make_2body_shr_pool(basis);
    Matrix D, F;
    {
        auto D_min = compute_soad(clustered_mol);
        F = compute_2body_fock_from_soad(eri_e, clustered_mol, basis, H, D_min);
        D = update_density(F, S, occ / 2);
    }

    Matrix Q
          = ints::detail::four_center_Q(world, eri_e, basis.flattened_shells());

    const auto maxiter = 100;
    const auto conv = 1e-12;
    auto iter = 0;
    auto rms_error = 1.0;
    auto ediff_rel = 0.0;
    auto ehf = 0.0;
    auto n2 = D.cols() * D.rows();
    libint2::DIIS<Matrix> diis(2); // start DIIS on second iteration

    std::cout << "Initial D = \n" << D << std::endl;

#if 0
    // prepare for incremental Fock build ...
    Matrix D_diff = D;
    F = H;
    bool reset_incremental_fock_formation = false;
    bool incremental_Fbuild_started = false;
    double start_incremental_F_threshold = 1e-5;
    double next_reset_threshold = 0.0;
    size_t last_reset_iteration = 0;

    do {
        const auto tstart = std::chrono::high_resolution_clock::now();
        ++iter;

        // Last iteration's energy and density
        auto ehf_last = ehf;
        Matrix D_last = D;

        if (not incremental_Fbuild_started
            && rms_error < start_incremental_F_threshold) {
            incremental_Fbuild_started = true;
            reset_incremental_fock_formation = false;
            last_reset_iteration = iter - 1;
            next_reset_threshold = rms_error / 1e1;
            std::cout << "== started incremental fock build" << std::endl;
        }
        if (reset_incremental_fock_formation
            || not incremental_Fbuild_started) {
            F = H;
            D_diff = D;
        }
        if (reset_incremental_fock_formation && incremental_Fbuild_started) {
            reset_incremental_fock_formation = false;
            last_reset_iteration = iter;
            next_reset_threshold = rms_error / 1e1;
            std::cout << "== reset incremental fock build" << std::endl;
        }

        F += compute_2body_G(eri_e, basis, Q, D_diff);

        // compute HF energy with the non-extrapolated Fock matrix
        ehf = D.cwiseProduct(H + F).sum();
        ediff_rel = std::abs((ehf - ehf_last) / ehf);

        // compute SCF error
        Matrix FD_comm = F * D * S - S * D * F;
        rms_error = FD_comm.norm() / n2;
        if (rms_error < next_reset_threshold
            || iter - last_reset_iteration >= 8)
            reset_incremental_fock_formation = true;

        // DIIS extrapolate F
        Matrix F_diis = F;
        diis.extrapolate(F_diis, FD_comm);

        D = update_density(F_diis, S, occ / 2);
        D_diff = D - D_last;

        const auto tstop = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> time_elapsed = tstop - tstart;

        if (iter == 1)
            std::cout << "\n\nIter         E(HF)                 D(E)/E        "
                         " RMS([F,D])/nn       Time(s)\n";
        printf(" %02d %20.12f %20.12e %20.12e %10.5lf\n", iter, ehf + repulsion_energy,
               ediff_rel, rms_error, time_elapsed.count());

    } while (((ediff_rel > conv) || (rms_error > conv)) && (iter < maxiter));

    printf("** Hartree-Fock energy = %20.12f\n", ehf + repulsion_energy);

    return 0;
#endif 
}


Matrix compute_soad(molecule::Molecule const &mol) {
    const auto atoms = mol.atoms();
    size_t nao = 0;
    for (const auto &atom : atoms) {
        const auto Z = atom.charge();
        if (Z == 1 || Z == 2) // H, He
            nao += 1;
        else if (Z <= 10) // Li - Ne
            nao += 5;
        else
            throw "SOAD with Z > 10 is not yet supported";
    }

    // compute the minimal basis density
    Matrix D = Matrix::Zero(nao, nao);
    size_t ao_offset = 0; // first AO of this atom
    for (const auto &atom : atoms) {
        const auto Z = atom.charge();
        if (Z == 1 || Z == 2) {          // H, He
            D(ao_offset, ao_offset) = Z; // all electrons go to the 1s
            ao_offset += 1;
        } else if (Z <= 10) {
            D(ao_offset, ao_offset) = 2; // 2 electrons go to the 1s
            D(ao_offset + 1, ao_offset + 1)
                  = (Z == 3) ? 1
                             : 2; // Li? only 1 electron in 2s, else 2 electrons
            // smear the remaining electrons in 2p orbitals
            const double num_electrons_per_2p = (Z > 4) ? (double)(Z - 4) / 3
                                                        : 0;
            for (auto xyz = 0; xyz != 3; ++xyz)
                D(ao_offset + 2 + xyz, ao_offset + 2 + xyz)
                      = num_electrons_per_2p;
            ao_offset += 5;
        }
    }

    return D * 0.5; // we use densities normalized to # of electrons/2
}

std::vector<int64_t> shell_to_func(std::vector<Shell> const &shells) {
    std::vector<int64_t> sh2f(1, 0);
    for (auto const &sh : shells) {
        sh2f.push_back(sh2f.back() + sh.size());
    }

    return sh2f;
}

template <typename ShrPool>
Matrix
compute_2body_fock_from_soad(ShrPool &engs, molecule::Molecule const &mol,
                             basis::Basis const &obs, Matrix const &H,
                             Matrix const &D, double precision) {

    basis::BasisSet bs("sto-3g");
    basis::Basis min_bs(bs.get_cluster_shells(mol));

    Matrix F = Matrix::Zero(obs.nfunctions(), obs.nfunctions());

    tbb::enumerable_thread_specific<Matrix> F_pool(F);

    auto obs_shells = obs.flattened_shells();
    auto mbs_shells = min_bs.flattened_shells();

    auto shell2bf = shell_to_func(obs_shells);
    auto shell2bf_D = shell_to_func(mbs_shells);

    auto task_f = [&](tbb::blocked_range<std::size_t> const &range) {
        engs->local().set_precision(precision);
        for (auto s1 = range.begin(); s1 != range.end(); ++s1) {
            auto const &sh1 = obs_shells[s1];
            auto bf1_first = shell2bf[s1];
            auto n1 = sh1.size();

            for (auto s2 = 0; s2 <= s1; ++s2) {
                auto const &sh2 = obs_shells[s2];
                auto bf2_first = shell2bf[s2];
                auto n2 = sh2.size();

                for (auto s3 = 0; s3 < mbs_shells.size(); ++s3) {
                    auto const &sh3 = mbs_shells[s3];
                    auto bf3_first = shell2bf_D[s3];
                    auto n3 = sh3.size();

                    auto s4_begin = s3;
                    auto s4_fence = s3 + 1;

                    for (auto s4 = s4_begin; s4 != s4_fence; ++s4) {
                        auto sh4 = mbs_shells[s4];
                        auto bf4_first = shell2bf_D[s4];
                        auto n4 = sh4.size();

                        // compute the permutational degeneracy (i.e. # of
                        // equivalents) of the given shell set
                        auto s12_deg = (s1 == s2) ? 1.0 : 2.0;

                        if (s3 >= s4) {
                            auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
                            auto s1234_deg = s12_deg * s34_deg;
                            // auto s1234_deg = s12_deg;
                            const auto *buf_J
                                  = engs->local().compute(sh1, sh2, sh3, sh4);

                            for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
                                const auto bf1 = f1 + bf1_first;
                                for (auto f2 = 0; f2 != n2; ++f2) {
                                    const auto bf2 = f2 + bf2_first;
                                    for (auto f3 = 0; f3 != n3; ++f3) {
                                        const auto bf3 = f3 + bf3_first;
                                        for (auto f4 = 0; f4 != n4;
                                             ++f4, ++f1234) {
                                            const auto bf4 = f4 + bf4_first;

                                            auto &g = F_pool.local();

                                            const auto value = buf_J[f1234];
                                            const auto value_scal_by_deg
                                                  = value * s1234_deg;
                                            g(bf1, bf2) += 2.0 * D(bf3, bf4)
                                                           * value_scal_by_deg;
                                        }
                                    }
                                }
                            }
                        }

                        const auto *buf_K
                              = engs->local().compute(sh1, sh3, sh2, sh4);

                        for (auto f1 = 0, f1324 = 0; f1 != n1; ++f1) {
                            const auto bf1 = f1 + bf1_first;
                            for (auto f3 = 0; f3 != n3; ++f3) {
                                const auto bf3 = f3 + bf3_first;
                                for (auto f2 = 0; f2 != n2; ++f2) {
                                    const auto bf2 = f2 + bf2_first;
                                    for (auto f4 = 0; f4 != n4; ++f4, ++f1324) {
                                        const auto bf4 = f4 + bf4_first;

                                        auto &g = F_pool.local();

                                        const auto value = buf_K[f1324];
                                        const auto value_scal_by_deg
                                              = value * s12_deg;
                                        g(bf1, bf2) -= D(bf3, bf4)
                                                       * value_scal_by_deg;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    };

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, obs_shells.size()),
                      task_f);

    F = F_pool.combine(
          [](Matrix const &F1, Matrix const &F2) { return Matrix(F1 + F2); });

    // symmetrize the result and return
    return H + 0.5 * (F + F.transpose());
}

template <typename ShrPool>
Matrix
compute_2body_G(ShrPool &engs, basis::Basis const &obs, Matrix const &Q,
                Matrix const &D,
                double precision) {

    auto obs_shells = obs.flattened_shells();
    auto shell2bf = shell_to_func(obs_shells);

    Matrix F = Matrix::Zero(obs.nfunctions(), obs.nfunctions());
    tbb::enumerable_thread_specific<Matrix> F_pool(F);

    auto task_j = [&](tbb::blocked_range<std::size_t> const &range) {
        engs->local().set_precision(precision);
        for (auto s1 = range.begin(); s1 != range.end(); ++s1) {
            auto const &sh1 = obs_shells[s1];
            auto bf1_first = shell2bf[s1];
            auto n1 = sh1.size();

            for (auto s2 = 0; s2 <= s1; ++s2) {
                auto const &sh2 = obs_shells[s2];
                auto bf2_first = shell2bf[s2];
                auto n2 = sh2.size();

                for (auto s3 = 0; s3 < s1; ++s3) {
                    auto const &sh3 = obs_shells[s3];
                    auto bf3_first = shell2bf[s3];
                    auto n3 = sh3.size();

                    const auto s4_max = (s1 == s3) ? s2 : s3;
                    for (auto s4 = 0; s4 <= s4_max; ++s4) {
                        auto sh4 = obs_shells[s4];
                        auto bf4_first = shell2bf[s4];
                        auto n4 = sh4.size();

                        // compute the permutational degeneracy (i.e. # of
                        // equivalents) of the given shell set
                        auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
                        auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
                        auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0)
                                                     : 2.0;
                        auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

                        const auto *buf = engs->local().compute(sh1, sh2, sh3, sh4);

                        // 1) each shell set of integrals contributes up to 6
                        // shell sets of the Fock matrix:
                        //    F(a,b) += (ab|cd) * D(c,d)
                        //    F(c,d) += (ab|cd) * D(a,b)
                        //    F(b,d) -= 1/4 * (ab|cd) * D(a,c)
                        //    F(b,c) -= 1/4 * (ab|cd) * D(a,d)
                        //    F(a,c) -= 1/4 * (ab|cd) * D(b,d)
                        //    F(a,d) -= 1/4 * (ab|cd) * D(b,c)
                        // 2) each permutationally-unique integral (shell set)
                        // must be scaled by its degeneracy,
                        //    i.e. the number of the integrals/sets equivalent
                        //    to it
                        // 3) the end result must be symmetrized
                        for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
                            const auto bf1 = f1 + bf1_first;
                            for (auto f2 = 0; f2 != n2; ++f2) {
                                const auto bf2 = f2 + bf2_first;
                                for (auto f3 = 0; f3 != n3; ++f3) {
                                    const auto bf3 = f3 + bf3_first;
                                    for (auto f4 = 0; f4 != n4; ++f4, ++f1234) {
                                        const auto bf4 = f4 + bf4_first;

                                        const auto value = buf[f1234];

                                        const auto value_scal_by_deg
                                              = value * s1234_deg;

                                        auto &g = F_pool.local();

                                        g(bf1, bf2) += D(bf3, bf4)
                                                       * value_scal_by_deg;
                                        g(bf3, bf4) += D(bf1, bf2)
                                                       * value_scal_by_deg;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    };

    auto task_k = [&](tbb::blocked_range<std::size_t> const &range) {
        engs->local().set_precision(precision);
        for (auto s1 = range.begin(); s1 != range.end(); ++s1) {
            auto const &sh1 = obs_shells[s1];
            auto bf1_first = shell2bf[s1];
            auto n1 = sh1.size();

            for (auto s2 = 0; s2 <= s1; ++s2) {
                auto const &sh2 = obs_shells[s2];
                auto bf2_first = shell2bf[s2];
                auto n2 = sh2.size();

                for (auto s3 = 0; s3 < s1; ++s3) {
                    auto const &sh3 = obs_shells[s3];
                    auto bf3_first = shell2bf[s3];
                    auto n3 = sh3.size();

                    const auto s4_max = (s1 == s3) ? s2 : s3;
                    for (auto s4 = 0; s4 <= s4_max; ++s4) {
                        auto sh4 = obs_shells[s4];
                        auto bf4_first = shell2bf[s4];
                        auto n4 = sh4.size();

                        // compute the permutational degeneracy (i.e. # of
                        // equivalents) of the given shell set
                        auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
                        auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
                        auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0)
                                                     : 2.0;
                        auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

                        const auto *buf = engs->local().compute(sh1, sh2, sh3, sh4);

                        // 1) each shell set of integrals contributes up to 6
                        // shell sets of the Fock matrix:
                        //    F(a,b) += (ab|cd) * D(c,d)
                        //    F(c,d) += (ab|cd) * D(a,b)
                        //    F(b,d) -= 1/4 * (ab|cd) * D(a,c)
                        //    F(b,c) -= 1/4 * (ab|cd) * D(a,d)
                        //    F(a,c) -= 1/4 * (ab|cd) * D(b,d)
                        //    F(a,d) -= 1/4 * (ab|cd) * D(b,c)
                        // 2) each permutationally-unique integral (shell set)
                        // must be scaled by its degeneracy,
                        //    i.e. the number of the integrals/sets equivalent
                        //    to it
                        // 3) the end result must be symmetrized
                        for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
                            const auto bf1 = f1 + bf1_first;
                            for (auto f2 = 0; f2 != n2; ++f2) {
                                const auto bf2 = f2 + bf2_first;
                                for (auto f3 = 0; f3 != n3; ++f3) {
                                    const auto bf3 = f3 + bf3_first;
                                    for (auto f4 = 0; f4 != n4; ++f4, ++f1234) {
                                        const auto bf4 = f4 + bf4_first;

                                        const auto value = buf[f1234];

                                        const auto value_scal_by_deg
                                              = value * s1234_deg;

                                        auto &g = F_pool.local();

                                        g(bf1, bf3) -= 0.25 * D(bf2, bf4)
                                                       * value_scal_by_deg;
                                        g(bf2, bf4) -= 0.25 * D(bf1, bf3)
                                                       * value_scal_by_deg;
                                        g(bf1, bf4) -= 0.25 * D(bf2, bf3)
                                                       * value_scal_by_deg;
                                        g(bf2, bf3) -= 0.25 * D(bf1, bf4)
                                                       * value_scal_by_deg;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    };

    auto j0 = tcc::utility::time::now();
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, obs_shells.size()),
                      task_j);
    auto j1 = tcc::utility::time::now();

    auto k0 = tcc::utility::time::now();
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, obs_shells.size()),
                      task_k);
    auto k1 = tcc::utility::time::now();

    std::cout << "\tJ time = " << tcc_time::duration_in_s(j0, j1) << std::endl;
    std::cout << "\tK time = " << tcc_time::duration_in_s(k0, k1) << std::endl;

    F = F_pool.combine(
          [](Matrix const &F1, Matrix const &F2) { return Matrix(F1 + F2); });

    return F;
}
