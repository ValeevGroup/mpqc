#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_ZRHF_IMPL_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_ZRHF_IMPL_H_

#include "mpqc/chemistry/qc/lcao/scf/zrhf.h"

#include "mpqc/chemistry/qc/lcao/scf/decomposed_rij.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_cond_ortho.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_df_fock_builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_four_center_fock_builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_four_center_j_cadf_k_fock_builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_ma_four_center_fock_builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_ma_ri_j_cadf_k_fock_builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_ma_ri_j_four_center_k_fock_builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_ma_four_center_j_cadf_k_fock_builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_ri_j_cadf_k_fock_builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_soad.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_two_center_builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/util.h"

#include "mpqc/math/external/tiledarray/util.h"

namespace mpqc {
namespace lcao {

/**
 * zRHF member functions
 */

template <typename Tile, typename Policy>
zRHF<Tile, Policy>::zRHF(const KeyVal& kv)
    : PeriodicAOWavefunction<Tile, Policy>(kv), kv_(kv) {}

template <typename Tile, typename Policy>
void zRHF<Tile, Policy>::init(const KeyVal& kv) {
  maxiter_ = kv.value<int64_t>("max_iter", 30);
  bool soad_guess = kv.value<bool>("soad_guess", true);
  print_detail_ = kv.value<bool>("print_detail", false);
  print_max_item_ = kv.value<int64_t>("print_max_item", 100);
  max_condition_num_ = kv.value<double>("max_condition_num", 1.0e8);
  fmix_ = kv.value<double>("fock_mixing", 0.0);
  level_shift_ = kv.value<double>("level_shift", 0.0);

  diis_ = kv.value<std::string>("diis", "none");
  diis_start_ = kv.value<unsigned int>("diis_start", 1);
  diis_num_vecs_ = kv.value<unsigned int>("diis_num_vecs", 5);
  diis_damping_ = kv.value<double>("diis_damping", 0.0);
  diis_mixing_ = kv.value<double>("diis_mixing", 0.0);
  diis_num_iters_group_ = kv.value<unsigned int>("diis_num_iters_group", 1);
  diis_num_extrap_group_ = kv.value<unsigned int>("diis_num_extrap_group", 1);

  iter_ = 0;

  auto& ao_factory = this->ao_factory();

  // retrieve world from periodic ao_factory
  auto& world = ao_factory.world();
  auto unitcell = ao_factory.unitcell();

  auto init_start = mpqc::fenced_now(world);

  init_fock_builder();

  ExEnv::out0() << ao_factory << std::endl;
  ExEnv::out0() << unitcell << std::endl;

  // the unit cell must be electrically neutral
  const auto charge = 0;
  const auto nelectrons = unitcell.total_atomic_number() - charge;
  if (nelectrons % 2 != 0)
    throw InputError("zRHF requires an even number of electrons", __FILE__,
                     __LINE__, "unitcell");
  docc_ = nelectrons / 2;
  dcell_ = unitcell.dcell();

  // retrieve unitcell info from periodic ao_factory
  R_max_ = ao_factory.R_max();
  RJ_max_ = ao_factory.RJ_max();
  RD_max_ = ao_factory.RD_max();
  R_size_ = ao_factory.R_size();
  RJ_size_ = ao_factory.RJ_size();
  RD_size_ = ao_factory.RD_size();

  // read # kpoints from keyval
  using ::mpqc::detail::k_ord_idx;
  nk_ = decltype(nk_)(kv.value<std::array<int, 3>>("k_points").data());
  MPQC_ASSERT((nk_.array() > 0).all() &&
              "k_points cannot have zero or negative elements");
  // assume odd number of k points in each direction, i.e. uniformly spaced k
  MPQC_ASSERT((nk_(0) % 2 == 1) && (nk_(1) % 2 == 1) && (nk_(2) % 2 == 1));

  k_size_ = 1 + k_ord_idx(nk_(0) - 1, nk_(1) - 1, nk_(2) - 1, nk_);
  ExEnv::out0() << "zRHF computational parameters:" << std::endl;
  ExEnv::out0() << indent << "# of k points in each direction: ["
                << nk_.transpose() << "]" << std::endl;

  eps_.resize(k_size_);
  C_.resize(k_size_);

  T_ = ao_factory.compute(L"<κ|T|λ>");  // Kinetic

  // Nuclear-attraction
  {
    ExEnv::out0()
        << "\nComputing Two Center Integral for Periodic System: < κ |V| λ >"
        << std::endl;
    auto t0 = mpqc::fenced_now(world);
    using Builder = scf::PeriodicTwoCenterBuilder<Tile, Policy>;
    auto two_center_builder = std::make_unique<Builder>(ao_factory);
    V_ = two_center_builder->eval(Operator::Type::Nuclear);
    auto t1 = mpqc::fenced_now(world);
    auto dur = mpqc::duration_in_s(t0, t1);
    ExEnv::out0() << " Time: " << dur << " s" << std::endl;
  }

  S_ = ao_factory.compute(L"<κ|λ>");  // Overlap in real space
  Sk_ = transform_real2recip(S_);     // Overlap in reciprocal space
  H_("mu, nu") =
      T_("mu, nu") + V_("mu, nu");  // One-body hamiltonian in real space

  // compute density matrix using soad/core guess
  array_type F_init;
  if (!soad_guess) {
    ExEnv::out0() << "\nUsing CORE guess for initial Fock ..." << std::endl;
    F_init = H_;
  } else {
    F_init = gaussian::periodic_fock_soad(world, unitcell, H_, ao_factory);
  }

  // transform Fock from real to reciprocal space
  Fk_ = transform_real2recip(F_init);
  // compute orthogonalizer matrix
  X_ = utility::conditioned_orthogonalizer(Sk_, k_size_, max_condition_num_,
                                           print_max_item_);
  // compute guess density
  std::tie(D_, Dk_) = compute_density();

  auto init_end = mpqc::fenced_now(world);
  init_duration_ = mpqc::duration_in_s(init_start, init_end);
  ExEnv::out0() << "Periodic RHF Init Time: " << init_duration_ << " s"
                << std::endl;
}

template <typename Tile, typename Policy>
void zRHF<Tile, Policy>::solve(double thresh) {
  auto& ao_factory = this->ao_factory();
  auto& world = ao_factory.world();

  auto rms = 0.0;
  array_type Ddiff;
  auto converged = false;
  auto ezrhf = 0.0;
  auto ediff = 0.0;

  // compute nuclear-repulsion energy
  const auto erep = ao_factory.unitcell().nuclear_repulsion_energy(RJ_max_);
  ExEnv::out0() << "\nNuclear Repulsion Energy: " << erep << std::endl;

  auto etotal = erep;

  TiledArray::DIIS<array_type_z> diis_gamma_point(
      diis_start_, diis_num_vecs_, diis_damping_, diis_num_iters_group_,
      diis_num_extrap_group_, diis_mixing_);
  TiledArray::DIIS<array_type_z> diis_allk(
      diis_start_, diis_num_vecs_, diis_damping_, diis_num_iters_group_,
      diis_num_extrap_group_, diis_mixing_);

  do {
    auto iter_start = mpqc::fenced_now(world);

    iter_++;

    // Save a copy of energy and density
    auto etotal_old = etotal;
    auto D_old = D_;
    array_type_z Fk_old = Fk_;

    if (print_detail_) {
      ExEnv::out0() << "\nIteration: " << iter_ << "\n";
    }

    // F = H + 2J - K
    auto f_start = mpqc::fenced_now(world);
    F_ = build_F(D_, H_, R_max_);
    auto f_end = mpqc::fenced_now(world);

    // compute zRHF energy
    ezrhf = compute_energy();
    etotal = erep + ezrhf;

    const auto fock_lattice_range = f_builder_->fock_lattice_range();
    // extra updates for Fock and energy
    if (need_extra_update_) {
      MPQC_ASSERT(extra_F_.is_initialized());
      F_ = ::mpqc::pbc::detail::add(F_, extra_F_, fock_lattice_range, R_max_);
      etotal += extra_energy_;
    }

    // transform Fock from real to reciprocal space
    auto trans_start = mpqc::fenced_now(world);
    Fk_ = transform_real2recip(F_, fock_lattice_range, nk_);
    auto trans_end = mpqc::fenced_now(world);
    trans_duration_ += mpqc::duration_in_s(trans_start, trans_end);

    // DIIS
    if (diis_ != "none") {
      if (diis_ == "gamma_point") {
        MPQC_ASSERT(k_size_ >= 1 && k_size_ % 2 == 1);
        using ::mpqc::pbc::detail::slice_array_at_k;
        const auto k_ord_gamma_point = (k_size_ - 1) / 2;
        auto D_gp = slice_array_at_k(Dk_, k_ord_gamma_point, nk_);
        auto F_gp = slice_array_at_k(Fk_, k_ord_gamma_point, nk_);
        auto S_gp = slice_array_at_k(Sk_, k_ord_gamma_point, nk_);

        array_type_z grad_gp;
        grad_gp("i, j") = F_gp("i, k") * D_gp("k, l") * S_gp("l, j") -
                          S_gp("i, k") * D_gp("k, l") * F_gp("l, j");

        diis_gamma_point.compute_extrapolation_parameters(grad_gp, true);

        const auto param_computed = diis_gamma_point.parameters_computed();

        using EigenVectorX =
            typename TiledArray::DIIS<array_type_z>::EigenVectorX;

        auto diis_coeffs = EigenVectorX();
        unsigned int diis_nskip = 0;
        auto F_allk_diis = Fk_;
        if (iter_ == 1) {
          diis_allk.reinitialize(&F_allk_diis);
        }
        if (!param_computed) {
          diis_allk.extrapolate(F_allk_diis, diis_coeffs, diis_nskip, true);
        } else {
          diis_coeffs = diis_gamma_point.get_coeffs();
          diis_nskip = diis_gamma_point.get_nskip();
          diis_allk.extrapolate(F_allk_diis, diis_coeffs, diis_nskip, true);
        }
        Fk_ = F_allk_diis;
      } else {
        throw InputError("Currently only Gamma-point DIIS is implemented",
                         __FILE__, __LINE__, "diis");
      }
    }

    // mixing of Fock matrices in reciprocal space
    if (fmix_ > 0.0) {
      Fk_("mu, nu") = (1.0 - fmix_) * Fk_("mu, nu") + fmix_ * Fk_old("mu, nu");
    }

    // compute new density
    auto d_start = mpqc::fenced_now(world);
    std::tie(D_, Dk_) = compute_density();
    // update density in periodic ao_factory
    ao_factory.set_density(D_);
    auto d_end = mpqc::fenced_now(world);
    d_duration_ += mpqc::duration_in_s(d_start, d_end);

    // compute difference with last iteration
    ediff = etotal - etotal_old;
    Ddiff("mu, nu") = D_("mu, nu") - D_old("mu, nu");
    auto volume = Ddiff.trange().elements_range().volume();
    rms = Ddiff("mu, nu").norm() / volume;

    if ((rms <= thresh) || fabs(ediff) <= thresh) converged = true;

    auto iter_end = mpqc::fenced_now(world);
    auto iter_duration = mpqc::duration_in_s(iter_start, iter_end);
    scf_duration_ += iter_duration;

    // Print out information
    if (print_detail_) {
      ExEnv::out0() << "\nzRHF Energy: " << ezrhf << "\n"
                    << "Total Energy: " << etotal << "\n"
                    << "Delta(E): " << ediff << "\n"
                    << "RMS(D): " << rms << "\n"
                    << "Fock Build Time: "
                    << mpqc::duration_in_s(f_start, f_end) << " s\n"
                    << "Transform Fock (Real->Recip) Time: "
                    << mpqc::duration_in_s(trans_start, trans_end) << " s\n"
                    << "Density Time: " << mpqc::duration_in_s(d_start, d_end)
                    << " s\n"
                    << "Iteration Time: " << iter_duration << " s\n";
      print_band_gaps();
    } else {
      std::string niter = "Iter", nEle = "E(HF)", nTot = "E(tot)",
                  nDel = "Delta(E)", nRMS = "RMS(D)", nT = "Time(s)";
      if (iter_ == 1)
        ExEnv::out0() << mpqc::printf("\n\n %4s %20s %20s %20s %20s %20s\n",
                                      niter.c_str(), nEle.c_str(), nTot.c_str(),
                                      nDel.c_str(), nRMS.c_str(), nT.c_str());
      ExEnv::out0() << mpqc::printf(
          " %4d %20.12f %20.12f %20.12f %20.12f %20.3f\n", iter_, ezrhf, etotal,
          ediff, rms, iter_duration);
    }

  } while ((iter_ <= maxiter_) && (!converged));

  // save total energy to energy no matter if zRHF converges
  energy_ = etotal;

  if (!converged) {
    // TODO read a keyval value to determine
    if (1)
      ExEnv::out0() << "\nzRHF: SCF did not converge!\n\n";
    else
      throw MaxIterExceeded("zRHF: SCF did not converge", __FILE__, __LINE__,
                            maxiter_);
  } else {
    ExEnv::out0() << "\nPeriodic Hartree-Fock iterations have converged!\n";
  }

  // print out band gap information
  print_band_gaps();

  // store fock matrix in registry
  auto& registry = this->ao_factory().registry();
  f_builder_->register_fock(F_, registry);

  if (world.rank() == 0) {
    std::cout << "\nTotal Periodic Hartree-Fock energy = " << energy_
              << std::endl;

    if (print_detail_ && k_size_ < print_max_item_) {
      Eigen::IOFormat fmt(5);
      std::cout << "\n k | orbital energies" << std::endl;
      for (auto k = 0; k < k_size_; ++k) {
        std::cout << k << " | " << eps_[k].transpose().format(fmt) << std::endl;
      }
    }

    // print out timings
    std::cout << mpqc::printf("\nTime(s):\n");
    std::cout << mpqc::printf("\tInit:                %20.3f\n",
                              init_duration_);
    std::cout << mpqc::printf("\tCoulomb term:        %20.3f\n", j_duration_);
    std::cout << mpqc::printf("\tExchange term:       %20.3f\n", k_duration_);
    std::cout << mpqc::printf("\tReal->Recip trans:   %20.3f\n",
                              trans_duration_);
    std::cout << mpqc::printf("\tDiag + Density:      %20.3f\n", d_duration_);
    std::cout << mpqc::printf("\tTotal:               %20.3f\n\n",
                              scf_duration_);
  }

  // test
  //  using MA_Builder = ::mpqc::pbc::ma::PeriodicMA<factory_type>;
  //  auto ma_builder = std::make_unique<MA_Builder>(this->ao_factory());
  //  ExEnv::out0() << "\n*** test multipole after converged scf ***\n";
  //  auto elec_moments = ma_builder->compute_elec_multipole_moments(D_);
  //  ExEnv::out0() << "electronic spherical multipole moments:"
  //                << "\nmonopole: " << elec_moments[0]
  //                << "\ndipole m=-1: " << elec_moments[1]
  //                << "\ndipole m=0:  " << elec_moments[2]
  //                << "\ndipole m=1:  " << elec_moments[3]
  //                << "\nquadrupole m=-2: " << elec_moments[4]
  //                << "\nquadrupole m=-1: " << elec_moments[5]
  //                << "\nquadrupole m=0:  " << elec_moments[6]
  //                << "\nquadrupole m=1:  " << elec_moments[7]
  //                << "\nquadrupole m=2:  " << elec_moments[8] << "\n";
}

template <typename Tile, typename Policy>
std::pair<typename zRHF<Tile, Policy>::array_type,
          typename zRHF<Tile, Policy>::array_type_z>
zRHF<Tile, Policy>::compute_density() {
  auto& ao_factory = this->ao_factory();
  auto& world = ao_factory.world();

  using ::mpqc::detail::extend_trange1;

  auto tr0 = Fk_.trange().data()[0];
  auto tr1_real = extend_trange1(tr0, RD_size_);
  auto tr1_recip = extend_trange1(tr0, k_size_);

  const auto ext0 = tr0.extent();
  RowMatrixXd result_real_eig(ext0, tr1_real.extent());
  MatrixZ result_recip_eig(ext0, tr1_recip.extent());
  result_real_eig.setZero();
  result_recip_eig.setZero();

  MatrixzVec F_recip_vec, D_recip_vec;
  MatrixdVec D_real_vec;
  F_recip_vec.resize(k_size_);
  D_recip_vec.resize(k_size_);
  D_real_vec.resize(RD_size_);

  auto fock_eig = array_ops::array_to_eigen(Fk_);

  // parallel impl for F_k diagonalization and D_k build
  auto compute_recip_density = [ext0, this](MatrixZ *F_ptr, MatrixZ *X_ptr, MatrixZ *C_old_ptr, MatrixZ *C_ptr, VectorD *eps_ptr, MatrixZ *D_ptr, bool is_gamma_point, bool do_level_shift) {
    // get references to matrix pointers
    const auto &F = *F_ptr;
    const auto &X = *X_ptr;
    auto &C = *C_ptr;
    auto &eps = *eps_ptr;
    auto &D = *D_ptr;

    // symmetrize Fock
    MatrixZ F_symm = 0.5 * (F.transpose().conjugate() + F);
    if (is_gamma_point) {
      F_symm = this->reverse_phase_factor(F_symm);
    }

    // diagonalize Fock
    if (do_level_shift) {
      // transform Fock from AO to CO basis
      const auto &C_old = *C_old_ptr;
      MatrixZ FCO = C_old.transpose().conjugate() * F_symm * C_old;
      // add energy level shift to diagonal elements of unoccupied orbitals
      auto nobs = FCO.cols();
      for (auto a = docc_; a != nobs; ++a) {
        FCO(a, a) += level_shift_;
      }
      // diagonalize Fock in CO basis
      Eigen::SelfAdjointEigenSolver<MatrixZ> comp_eig_solver_fco(FCO);
      VectorD eps_temp = comp_eig_solver_fco.eigenvalues();
      MatrixZ C_temp = comp_eig_solver_fco.eigenvectors();
      // when k=0 (gamma point), reverse phase factor of complex eigenvectors
      if (is_gamma_point) {
        C_temp = this->reverse_phase_factor(C_temp);
      }
      // transform eigenvectors back to CO coefficients
      C = C_old * C_temp;
      // remove energy level shift from eigenvalues
      for (auto p = 0; p != nobs; ++p) {
        eps(p) = (p < docc_) ? eps_temp(p) : eps_temp(p) - level_shift_;
      }
    } else {
      // Orthogonalize Fock matrix: F' = Xt * F * X
      MatrixZ Ft = X.transpose().conjugate() * F_symm * X;
      // Diagonalize F'
      Eigen::SelfAdjointEigenSolver<MatrixZ> comp_eig_solver(Ft);
      eps = comp_eig_solver.eigenvalues();
      MatrixZ Ctemp = comp_eig_solver.eigenvectors();
      // When k=0 (gamma point), reverse phase factor of complex eigenvectors
      if (is_gamma_point) {
        Ctemp = this->reverse_phase_factor(Ctemp);
      }
      // transform eigenvectors back to CO coefficients
      C = X * Ctemp;
    }

    // compute_density
    MatrixZ C_occ = C_ptr->leftCols(docc_);
    D = C_occ.conjugate() * C_occ.transpose();
  };

  bool do_level_shift = (level_shift_ > 0.0 && iter_ > 0);
  for (int64_t k = 0; k != k_size_; ++k) {
    bool is_gamma_point = (k_size_ > 1 && k == ((k_size_ - 1) / 2));
    MatrixZ *C_old = do_level_shift ? &C_[k] : nullptr;
    F_recip_vec[k] = fock_eig.block(0, k * ext0, ext0, ext0);

    world.taskq.add(compute_recip_density, &(F_recip_vec[k]), &(X_[k]), C_old, &(C_[k]), &(eps_[k]), &(D_recip_vec[k]), is_gamma_point, do_level_shift);
  }
  world.gop.fence();

  // collect all D_k's to a big rectangular matrix
  for (int64_t k = 0; k != k_size_; ++k) {
    result_recip_eig.block(0, k * ext0, ext0, ext0) = D_recip_vec[k];
  }

  // parallel impl for D_r
  const auto denom_inv = 1.0 / double(nk_(0) * nk_(1) * nk_(2));
  auto transform_recip2real = [ext0, denom_inv, this](MatrixZ *D_recip_ptr, RowMatrixXd *D_real_ptr, int64_t R) {
    using ::mpqc::detail::k_vector;
    using ::mpqc::detail::direct_vector;

    // get references to all matrix pointers
    const auto &D_recip = *D_recip_ptr;
    auto &D_real = *D_real_ptr;

    auto vec_R = direct_vector(R, RD_max_, dcell_);
    D_real.setZero(ext0, ext0);

    for (int64_t k = 0; k != k_size_; ++k) {
      auto vec_k = k_vector(k, nk_, dcell_);
      auto exponent = std::exp(I * vec_k.dot(vec_R)) * denom_inv;
      D_real += (exponent * D_recip.block(0, k * ext0, ext0, ext0)).real();
    }
  };

  for (int64_t Rd = 0; Rd != RD_size_; ++Rd) {
    world.taskq.add(transform_recip2real, &result_recip_eig, &(D_real_vec[Rd]), Rd);
  }
  world.gop.fence();

  // collect all D_r's to a big rectangular matrix
  for (int64_t Rd = 0; Rd != RD_size_; ++Rd) {
    result_real_eig.block(0, Rd * ext0, ext0, ext0) = D_real_vec[Rd];
  }

  // convert rectangular matrices to TA::DistArray
  auto result_real = array_ops::eigen_to_array<Tile, Policy>(
      world, result_real_eig, tr0, tr1_real);
  auto result_recip = array_ops::eigen_to_array<TA::TensorZ, Policy>(
      world, result_recip_eig, tr0, tr1_recip);

  return std::make_pair(result_real, result_recip);
}

template <typename Tile, typename Policy>
typename zRHF<Tile, Policy>::array_type_z
zRHF<Tile, Policy>::transform_real2recip(const array_type& matrix,
                                         const Vector3i& real_lattice_range,
                                         const Vector3i& recip_lattice_range) {
  // Make sure range values are all positive
  MPQC_ASSERT((real_lattice_range.array() >= 0).all() &&
              (recip_lattice_range.array() > 0).all());

  using ::mpqc::detail::direct_ord_idx;
  using ::mpqc::detail::direct_vector;
  using ::mpqc::detail::extend_trange1;
  using ::mpqc::detail::k_ord_idx;
  using ::mpqc::detail::k_vector;

  const auto real_lattice_size =
      1 + direct_ord_idx(real_lattice_range, real_lattice_range);
  const Vector3i k_end_3D_idx = (recip_lattice_range.array() - 1).matrix();
  const auto recip_lattice_size =
      1 + k_ord_idx(k_end_3D_idx, recip_lattice_range);
  const auto tiles_range = matrix.trange().tiles_range();
  MPQC_ASSERT(tiles_range.extent(1) % tiles_range.extent(0) == 0);
  MPQC_ASSERT(real_lattice_size ==
              tiles_range.extent(1) / tiles_range.extent(0));

  array_type_z result;
  auto tr0 = matrix.trange().data()[0];
  auto tr1 = extend_trange1(tr0, recip_lattice_size);
  auto& world = matrix.world();

  // Perform real->reciprocal transformation with Eigen
  // TODO: perform it with TA
  auto matrix_eig = array_ops::array_to_eigen(matrix);
  MatrixZ result_eig(tr0.extent(), tr1.extent());
  result_eig.setZero();

  auto threshold = std::numeric_limits<double>::epsilon();
  for (auto R = 0; R < real_lattice_size; ++R) {
    auto bmat =
        matrix_eig.block(0, R * tr0.extent(), tr0.extent(), tr0.extent());
    if (bmat.norm() < bmat.size() * threshold)
      continue;
    else {
      auto vec_R = direct_vector(R, real_lattice_range, dcell_);
      for (auto k = 0; k < recip_lattice_size; ++k) {
        auto vec_k = k_vector(k, recip_lattice_range, dcell_);
        auto exponent = std::exp(I * vec_k.dot(vec_R));
        result_eig.block(0, k * tr0.extent(), tr0.extent(), tr0.extent()) +=
            bmat * exponent;
      }
    }
  }

  result = array_ops::eigen_to_array<TA::TensorZ, TA::SparsePolicy>(
      world, result_eig, tr0, tr1);

  return result;
}

template <typename Tile, typename Policy>
typename zRHF<Tile, Policy>::array_type_z
zRHF<Tile, Policy>::transform_real2recip(const array_type& matrix) {
  using ::mpqc::detail::direct_ord_idx;
  using ::mpqc::detail::k_ord_idx;

  const Vector3i k_end_3D_idx = (nk_.array() - 1).matrix();
  MPQC_ASSERT(k_size_ == 1 + k_ord_idx(k_end_3D_idx, nk_));
  MPQC_ASSERT(R_size_ == 1 + direct_ord_idx(R_max_, R_max_));
  return transform_real2recip(matrix, R_max_, nk_);
}

template <typename Tile, typename Policy>
MatrixZ zRHF<Tile, Policy>::reverse_phase_factor(const MatrixZ& mat0) {
  MatrixZ result(mat0);

  for (auto row = 0; row < mat0.rows(); ++row) {
    for (auto col = 0; col < mat0.cols(); ++col) {
      std::complex<double> comp0 = mat0(row, col);

      double norm = std::abs(comp0);
      if (norm == 0.0) {
        result(row, col) = comp0;
      } else {
        double real = comp0.real();
        double imag = comp0.imag();

        double phi = std::atan(imag / real);

        double R;
        if (std::cos(phi) != 0.0) {
          R = real / std::cos(phi);
        } else {
          R = imag / std::sin(phi);
        }

        std::complex<double> comp1 = comp0 * std::exp(-1.0 * I * phi);

        result(row, col) = comp1;
      }
    }
  }

  return result;
}

template <typename Tile, typename Policy>
void zRHF<Tile, Policy>::print_band_gaps() {
  MPQC_ASSERT(eps_.size() == k_size_);

  VectorD hoco(k_size_), luco(k_size_), direct_gap(k_size_);

  for (auto k = 0; k != k_size_; ++k) {
    const auto homo = eps_[k](docc_ - 1);
    const auto lumo = eps_[k](docc_);
    hoco(k) = homo;
    luco(k) = lumo;
    direct_gap(k) = lumo - homo;
  }

  using ::mpqc::detail::k_3D_idx;

  VectorD::Index k_ord_direct_max;
  const auto direct_gap_max = direct_gap.maxCoeff(&k_ord_direct_max);
  const auto k_3D_direct_max = k_3D_idx(k_ord_direct_max, nk_);

  VectorD::Index k_ord_direct_min;
  const auto direct_gap_min = direct_gap.minCoeff(&k_ord_direct_min);
  const auto k_3D_direct_min = k_3D_idx(k_ord_direct_min, nk_);

  VectorD::Index k_ord_luco_max;
  const auto luco_max = luco.maxCoeff(&k_ord_luco_max);
  const auto k_3D_luco_max = k_3D_idx(k_ord_luco_max, nk_);

  VectorD::Index k_ord_luco_min;
  const auto luco_min = luco.minCoeff(&k_ord_luco_min);
  const auto k_3D_luco_min = k_3D_idx(k_ord_luco_min, nk_);

  VectorD::Index k_ord_hoco_max;
  const auto hoco_max = hoco.maxCoeff(&k_ord_hoco_max);
  const auto k_3D_hoco_max = k_3D_idx(k_ord_hoco_max, nk_);

  VectorD::Index k_ord_hoco_min;
  const auto hoco_min = hoco.minCoeff(&k_ord_hoco_min);
  const auto k_3D_hoco_min = k_3D_idx(k_ord_hoco_min, nk_);

  const auto indirect_gap = luco_min - hoco_max;

  auto unit_factory = UnitFactory::get_default();
  auto hartree2ev = unit_factory->make_unit("eV").from_atomic_units();

  ExEnv::out0() << "\nMax LUCO: " << luco_max * hartree2ev << " eV at k = ("
                << k_3D_luco_max.transpose() << ")"
                << "\nMin LUCO: " << luco_min * hartree2ev << " eV at k = ("
                << k_3D_luco_min.transpose() << ")"
                << "\nMax HOCO: " << hoco_max * hartree2ev << " eV at k = ("
                << k_3D_hoco_max.transpose() << ")"
                << "\nMin HOCO: " << hoco_min * hartree2ev << " eV at k = ("
                << k_3D_hoco_min.transpose() << ")"
                << "\nIndirect band gap: " << indirect_gap * hartree2ev << " eV"
                << "\nMax direct band gap: " << direct_gap_max * hartree2ev
                << " eV at k = (" << k_3D_direct_max.transpose() << ")"
                << "\nMin direct band gap: " << direct_gap_min * hartree2ev
                << " eV at k = (" << k_3D_direct_min.transpose() << ")"
                << std::endl;
}

template <typename Tile, typename Policy>
void zRHF<Tile, Policy>::obsolete() {
  Wavefunction::obsolete();
}

template <typename Tile, typename Policy>
bool zRHF<Tile, Policy>::can_evaluate(Energy* energy) {
  // can only evaluate the energy
  return energy->order() == 0;
}

template <typename Tile, typename Policy>
void zRHF<Tile, Policy>::evaluate(Energy* result) {
  if (!this->computed()) {
    init(kv_);
    solve(result->target_precision(0));
    this->computed_ = true;
    set_value(result, energy_);
  }
}

template <typename Tile, typename Policy>
double zRHF<Tile, Policy>::compute_energy() {
  array_type H_plus_F;
  const auto fock_lattice_range = f_builder_->fock_lattice_range();
  H_plus_F = ::mpqc::pbc::detail::add(H_, F_, R_max_, fock_lattice_range);
  Vector3i HpF_lattice_range;
  if (R_max_(0) >= fock_lattice_range(0) &&
      R_max_(1) >= fock_lattice_range(1) &&
      R_max_(2) >= fock_lattice_range(2)) {
    HpF_lattice_range = R_max_;
  } else if (R_max_(0) <= fock_lattice_range(0) &&
             R_max_(1) <= fock_lattice_range(1) &&
             R_max_(2) <= fock_lattice_range(2)) {
    HpF_lattice_range = fock_lattice_range;
  } else {
    ExEnv::out0() << "\nLattice range of H: " << R_max_.transpose()
                  << "\nLattice range of Fock: "
                  << fock_lattice_range.transpose() << std::endl;
    throw ProgrammingError("Invalid lattice ranges!", __FILE__, __LINE__);
  }
  auto ezrhf = ::mpqc::pbc::detail::dot_product(H_plus_F, D_, HpF_lattice_range,
                                                RD_max_);
  return ezrhf;
}

template <typename Tile, typename Policy>
void zRHF<Tile, Policy>::init_fock_builder() {
  using Builder = scf::ReferencePeriodicFourCenterFockBuilder<
      Tile, Policy, zRHF<Tile, Policy>::factory_type>;
  this->f_builder_ = std::make_unique<Builder>(this->ao_factory());
}

template <typename Tile, typename Policy>
typename zRHF<Tile, Policy>::array_type zRHF<Tile, Policy>::build_F(
    const array_type& D, const array_type& H, const Vector3i& H_lattice_range) {
  auto G = f_builder_->operator()(D);
  const auto fock_lattice_range = f_builder_->fock_lattice_range();
  return ::mpqc::pbc::detail::add(H, G, H_lattice_range, fock_lattice_range);
}

/**
 *  DFzRHF member functions
 */

template <typename Tile, typename Policy>
DFzRHF<Tile, Policy>::DFzRHF(const KeyVal& kv) : zRHF<Tile, Policy>(kv) {}

template <typename Tile, typename Policy>
void DFzRHF<Tile, Policy>::init_fock_builder() {
  using Builder =
      scf::PeriodicDFFockBuilder<Tile, Policy,
                                 DFzRHF<Tile, Policy>::factory_type>;
  this->f_builder_ = std::make_unique<Builder>(this->ao_factory());
}

/**
 *  FourCenterzRHF member functions
 */

template <typename Tile, typename Policy>
FourCenterzRHF<Tile, Policy>::FourCenterzRHF(const KeyVal& kv)
    : zRHF<Tile, Policy>(kv) {}

template <typename Tile, typename Policy>
void FourCenterzRHF<Tile, Policy>::init_fock_builder() {
  using Builder = scf::PeriodicFourCenterFockBuilder<Tile, Policy>;
  this->f_builder_ = std::make_unique<Builder>(this->ao_factory(), true, true);
}

/**
 *  RIJCADFKzRHF member functions
 */

template <typename Tile, typename Policy>
RIJCADFKzRHF<Tile, Policy>::RIJCADFKzRHF(const KeyVal& kv)
    : zRHF<Tile, Policy>(kv) {
  force_shape_threshold_ = kv.value<double>("force_shape_threshold", 0.0);
}

template <typename Tile, typename Policy>
void RIJCADFKzRHF<Tile, Policy>::init_fock_builder() {
  using Builder = scf::PeriodicRIJCADFKFockBuilder<
      Tile, Policy, RIJCADFKzRHF<Tile, Policy>::factory_type>;
  this->f_builder_ =
      std::make_unique<Builder>(this->ao_factory(), force_shape_threshold_);
}

/**
 *  FourCenterJCADFKzRHF member functions
 */

template <typename Tile, typename Policy>
FourCenterJCADFKzRHF<Tile, Policy>::FourCenterJCADFKzRHF(const KeyVal& kv)
    : zRHF<Tile, Policy>(kv) {
  force_shape_threshold_ = kv.value<double>("force_shape_threshold", 0.0);
}

template <typename Tile, typename Policy>
void FourCenterJCADFKzRHF<Tile, Policy>::init_fock_builder() {
  using Builder = scf::PeriodicFourCenterJCADFKFockBuilder<
      Tile, Policy, FourCenterJCADFKzRHF<Tile, Policy>::factory_type>;
  this->f_builder_ =
      std::make_unique<Builder>(this->ao_factory(), force_shape_threshold_);
}

/**
 *  MARIJCADFKzRHF member functions
 */

template <typename Tile, typename Policy>
MARIJCADFKzRHF<Tile, Policy>::MARIJCADFKzRHF(const KeyVal& kv)
    : zRHF<Tile, Policy>(kv) {
  force_shape_threshold_ = kv.value<double>("force_shape_threshold", 0.0);
  ma_energy_threshold_ = kv.value<double>("ma_energy_threshold", 1e-9);
  ma_ws_ = kv.value<double>("ma_well_separateness", 3.0);
  ma_extent_threshold_ = kv.value<double>("ma_extent_threshold", 1e-6);
  ma_extent_smallval_ = kv.value<double>("ma_extent_small_value", 0.01);
}

template <typename Tile, typename Policy>
void MARIJCADFKzRHF<Tile, Policy>::init_fock_builder() {
  using Builder = scf::PeriodicMARIJCADFKFockBuilder<
      Tile, Policy, MARIJCADFKzRHF<Tile, Policy>::factory_type>;
  this->f_builder_ =
      std::make_unique<Builder>(this->ao_factory(), force_shape_threshold_, ma_energy_threshold_, ma_ws_, ma_extent_threshold_, ma_extent_smallval_);
  this->need_extra_update_ = dynamic_cast<Builder&>(*this->f_builder_)
                                 .coulomb_builder()
                                 .multipole_builder()
                                 .CFF_reached();
}

template <typename Tile, typename Policy>
typename MARIJCADFKzRHF<Tile, Policy>::array_type
MARIJCADFKzRHF<Tile, Policy>::build_F(const array_type& D, const array_type& H,
                                      const Vector3i& H_lattice_range) {
  auto G_cnf = this->f_builder_->operator()(D);
  const auto fock_lattice_range = this->f_builder_->fock_lattice_range();
  auto F_cnf =
      ::mpqc::pbc::detail::add(H, G_cnf, H_lattice_range, fock_lattice_range);

  using Builder = scf::PeriodicMARIJCADFKFockBuilder<
      Tile, Policy, MARIJCADFKzRHF<Tile, Policy>::factory_type>;

  auto& ma_builder = dynamic_cast<Builder&>(*this->f_builder_)
                         .coulomb_builder()
                         .multipole_builder();
  if (this->need_extra_update_) {
    this->extra_F_ = ma_builder.get_fock();
    this->extra_energy_ = ma_builder.get_energy();
  }

  return F_cnf;
}

/**
 *  MARIJFourCenterKzRHF member functions
 */

template <typename Tile, typename Policy>
MARIJFourCenterKzRHF<Tile, Policy>::MARIJFourCenterKzRHF(const KeyVal& kv)
    : zRHF<Tile, Policy>(kv) {
  ma_energy_threshold_ = kv.value<double>("ma_energy_threshold", 1e-9);
  ma_ws_ = kv.value<double>("ma_well_separateness", 3.0);
  ma_extent_threshold_ = kv.value<double>("ma_extent_threshold", 1e-6);
  ma_extent_smallval_ = kv.value<double>("ma_extent_small_value", 0.01);
}

template <typename Tile, typename Policy>
void MARIJFourCenterKzRHF<Tile, Policy>::init_fock_builder() {
  using Builder = scf::PeriodicMARIJFourCenterKFockBuilder<Tile, Policy>;
  this->f_builder_ = std::make_unique<Builder>(this->ao_factory(), ma_energy_threshold_, ma_ws_, ma_extent_threshold_, ma_extent_smallval_);
  this->need_extra_update_ = dynamic_cast<Builder&>(*this->f_builder_)
                                 .coulomb_builder()
                                 .multipole_builder()
                                 .CFF_reached();
}

template <typename Tile, typename Policy>
typename MARIJFourCenterKzRHF<Tile, Policy>::array_type
MARIJFourCenterKzRHF<Tile, Policy>::build_F(const array_type& D,
                                            const array_type& H,
                                            const Vector3i& H_lattice_range) {
  auto G_cnf = this->f_builder_->operator()(D);
  const auto fock_lattice_range = this->f_builder_->fock_lattice_range();
  auto F_cnf =
      ::mpqc::pbc::detail::add(H, G_cnf, H_lattice_range, fock_lattice_range);

  using Builder = scf::PeriodicMARIJFourCenterKFockBuilder<Tile, Policy>;

  auto& ma_builder = dynamic_cast<Builder&>(*this->f_builder_)
                         .coulomb_builder()
                         .multipole_builder();
  if (this->need_extra_update_) {
    this->extra_F_ = ma_builder.get_fock();
    this->extra_energy_ = ma_builder.get_energy();
  }

  return F_cnf;
}

/**
 *  MAFourCenterzRHF member functions
 */

template <typename Tile, typename Policy>
MAFourCenterzRHF<Tile, Policy>::MAFourCenterzRHF(const KeyVal& kv)
    : zRHF<Tile, Policy>(kv) {
  ma_energy_threshold_ = kv.value<double>("ma_energy_threshold", 1e-9);
  ma_ws_ = kv.value<double>("ma_well_separateness", 3.0);
  ma_extent_threshold_ = kv.value<double>("ma_extent_threshold", 1e-6);
  ma_extent_smallval_ = kv.value<double>("ma_extent_small_value", 0.01);
}

template <typename Tile, typename Policy>
void MAFourCenterzRHF<Tile, Policy>::init_fock_builder() {
  using Builder = scf::PeriodicMAFourCenterFockBuilder<Tile, Policy>;
  this->f_builder_ = std::make_unique<Builder>(this->ao_factory(), ma_energy_threshold_, ma_ws_, ma_extent_threshold_, ma_extent_smallval_);
  this->need_extra_update_ = dynamic_cast<Builder&>(*this->f_builder_)
                                 .multipole_builder()
                                 .CFF_reached();
}

template <typename Tile, typename Policy>
typename MAFourCenterzRHF<Tile, Policy>::array_type
MAFourCenterzRHF<Tile, Policy>::build_F(const array_type& D,
                                        const array_type& H,
                                        const Vector3i& H_lattice_range) {
  auto G_cnf = this->f_builder_->operator()(D);
  const auto fock_lattice_range = this->f_builder_->fock_lattice_range();
  auto F_cnf =
      ::mpqc::pbc::detail::add(H, G_cnf, H_lattice_range, fock_lattice_range);

  using Builder = scf::PeriodicMAFourCenterFockBuilder<Tile, Policy>;

  auto& ma_builder =
      dynamic_cast<Builder&>(*this->f_builder_).multipole_builder();
  if (this->need_extra_update_) {
    this->extra_F_ = ma_builder.get_fock();
    this->extra_energy_ = ma_builder.get_energy();
  }

  return F_cnf;
}

/**
 *  MAFourCenterJCADFKzRHF member functions
 */

template <typename Tile, typename Policy>
MAFourCenterJCADFKzRHF<Tile, Policy>::MAFourCenterJCADFKzRHF(const KeyVal& kv)
    : zRHF<Tile, Policy>(kv) {
  force_shape_threshold_ = kv.value<double>("force_shape_threshold", 0.0);
  ma_energy_threshold_ = kv.value<double>("ma_energy_threshold", 1e-9);
  ma_ws_ = kv.value<double>("ma_well_separateness", 3.0);
  ma_extent_threshold_ = kv.value<double>("ma_extent_threshold", 1e-6);
  ma_extent_smallval_ = kv.value<double>("ma_extent_small_value", 0.01);
}

template <typename Tile, typename Policy>
void MAFourCenterJCADFKzRHF<Tile, Policy>::init_fock_builder() {
  using Builder = scf::PeriodicMAFourCenterJCADFKFockBuilder<Tile, Policy>;
  this->f_builder_ = std::make_unique<Builder>(this->ao_factory(), force_shape_threshold_, ma_energy_threshold_, ma_ws_, ma_extent_threshold_, ma_extent_smallval_);
  this->need_extra_update_ = dynamic_cast<Builder&>(*this->f_builder_)
      .multipole_builder()
      .CFF_reached();
}

template <typename Tile, typename Policy>
typename MAFourCenterJCADFKzRHF<Tile, Policy>::array_type
MAFourCenterJCADFKzRHF<Tile, Policy>::build_F(const array_type& D,
                                        const array_type& H,
                                        const Vector3i& H_lattice_range) {
  auto G_cnf = this->f_builder_->operator()(D);
  const auto fock_lattice_range = this->f_builder_->fock_lattice_range();
  auto F_cnf =
      ::mpqc::pbc::detail::add(H, G_cnf, H_lattice_range, fock_lattice_range);

  using Builder = scf::PeriodicMAFourCenterJCADFKFockBuilder<Tile, Policy>;

  auto& ma_builder =
      dynamic_cast<Builder&>(*this->f_builder_).multipole_builder();
  if (this->need_extra_update_) {
    this->extra_F_ = ma_builder.get_fock();
    this->extra_energy_ = ma_builder.get_energy();
  }

  return F_cnf;
}


}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_ZRHF_IMPL_H_
