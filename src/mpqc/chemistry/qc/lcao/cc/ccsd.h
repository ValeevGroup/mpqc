//
// Created by Chong Peng on 7/1/15.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_CCSD_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_CCSD_H_

#include <tiledarray.h>

#include "mpqc/chemistry/qc/cc/diis.h"
#include "mpqc/chemistry/qc/cc/solvers.h"
#include "mpqc/chemistry/qc/lcao/cc/solvers.h"
#include "mpqc/chemistry/qc/lcao/expression/trange1_engine.h"
#include "mpqc/chemistry/qc/lcao/mbpt/denom.h"
#include "mpqc/chemistry/qc/lcao/scf/mo_build.h"
#include "mpqc/chemistry/qc/lcao/wfn/lcao_wfn.h"
#include "mpqc/chemistry/qc/properties/energy.h"
#include "mpqc/mpqc_config.h"

namespace mpqc {
namespace lcao {

namespace detail {

inline void print_ccsd(int iter, double dE, double error, double E1,
                       double time_total) {
  if (iter == 1) {
    ExEnv::out0() << indent << mpqc::printf("%3s \t %10s \t %10s \t %15s \t %10s \n", "iter", "deltaE",
                "residual", "energy", "total time/s");
  }
  ExEnv::out0() << indent << mpqc::printf("%3i \t %10.5e \t %10.5e \t %15.12f \t %10.1f \n", iter, dE,
              error, E1, time_total);
}

inline void print_ccsd_direct(int iter, double dE, double error, double E1,
                              double time_u, double time_total) {
  if (iter == 1) {
    ExEnv::out0() << indent << mpqc::printf("%3s \t %10s \t %10s \t %15s \t %10s \t %10s \n", "iter",
                "deltaE", "residual", "energy", "u time/s", "total time/s");
  }
  ExEnv::out0() << indent << mpqc::printf("%3i \t %10.5e \t %10.5e \t %15.12f \t %10.1f \t %10.1f \n", iter,
              dE, error, E1, time_u, time_total);
}
}

/**
 * CCSD class that computed CCSD energy
 */

template <typename Tile, typename Policy>
class CCSD : public LCAOWavefunction<Tile, Policy>,
             public Provides<Energy>,
             public std::enable_shared_from_this<CCSD<Tile, Policy>> {
 public:
  using TArray = TA::DistArray<Tile, Policy>;
  using AOFactory = gaussian::AOFactory<Tile, Policy>;

  CCSD() = default;


  // clang-format off

  /**
   * KeyVal constructor
   * @param kv
   *
   * keywords : all keywords for LCAOWavefunciton
   *
   * | Keyword | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | @c ref | Wavefunction | @c none | reference Wavefunction, need to be a Energy::Provider RHF for example |
   * | @c method | string | @c df if @c df_basis is provided, @c standard otherwise | method to compute the CCSD residual; valid choices are: @c standard (uses 4-index MO integrals throughout), @c direct (uses 4-index MO integrals with up to 3 unoccupied indices, and 4-center AO integrals), @c df (approximates 4-index MO integrals using density fitting), @c direct_df (hybrid between @c df and @c direct that avoids storing MO integrals with 3 unoccupied indices by using DF, see DOI 10.1021/acs.jpca.6b10150 for details) |
   * | @c max_iter | int | @c 30 | maxmium iteration in CCSD |
   * | @c verbose | bool | determined by factory.verbose() | if print more information in CCSD iteration |
   * | @c reduced_abcd_memory | bool | @c true | if @c method=standard , avoid storing an extra abcd intermediate at the cost of increased FLOPs; if @c method=df , avoid storage of (ab|cd) integral in favor of lazy evaluation in batches |
   */

  // clang-format on

  CCSD(const KeyVal &kv) : LCAOWavefunction<Tile, Policy>(kv), kv_(kv) {
    if (kv.exists("ref")) {
      ref_wfn_ = kv.class_ptr<Wavefunction>("ref");
    } else {
      throw InputError("Default Ref in CCSD is not support! \n", __FILE__,
                       __LINE__, "ref");
    }

    df_ = false;
    auto default_method =
        this->lcao_factory().basis_registry()->have(L"Κ") ? "df" : "direct";
    method_ = kv.value<std::string>("method", default_method);
    if (method_ != "df" && method_ != "direct" && method_ != "standard" &&
        method_ != "direct_df") {
      throw InputError("Invalid CCSD method! \n", __FILE__, __LINE__, "method");
    }
    if (method_ == "df" || method_ == "direct_df") {
      df_ = true;
    }

    solver_str_ = kv.value<std::string>("solver", "jacobi_diis");
    if (solver_str_ != "jacobi_diis" && solver_str_ != "pno" && solver_str_ != "svo" && solver_str_ != "exact_pno")
      throw InputError("invalid value for solver keyword", __FILE__, __LINE__, "solver");

    reduced_abcd_memory_ = kv.value<bool>("reduced_abcd_memory", true);

    max_iter_ = kv.value<int>("max_iter", 30);
    verbose_ = kv.value<bool>("verbose", this->lcao_factory().verbose());
    min_iter_ = kv.value<int>("min_iter", 20);
  }

  virtual ~CCSD() {}

  /// protected members
 protected:
  TArray T1_;
  TArray T2_;

  const KeyVal
      kv_;  // the input keyval is kept to avoid heavy initialization in ctor
  std::string solver_str_;
  std::shared_ptr<::mpqc::cc::Solver<TArray>> solver_;
  std::shared_ptr<Wavefunction> ref_wfn_;
  typename AOFactory::DirectTArray direct_ao_array_;
  bool df_;
  bool reduced_abcd_memory_ = false;
  std::string method_;
  std::size_t max_iter_;
  std::size_t min_iter_;
  double target_precision_;
  double computed_precision_ = std::numeric_limits<double>::max();
  bool verbose_;
  double ccsd_corr_energy_;
  // diagonal elements of the Fock matrix (not necessarily the eigenvalues)
  std::shared_ptr<const Eigen::VectorXd> f_pq_diagonal_;

 protected:
  std::shared_ptr<const Eigen::VectorXd> orbital_energy() {
    return f_pq_diagonal_;
  }

  void set_orbital_energy(const Eigen::VectorXd &orbital_energy) {
    f_pq_diagonal_ = std::make_shared<Eigen::VectorXd>(orbital_energy);
  }

 public:
  void obsolete() override {
    ccsd_corr_energy_ = 0.0;
    f_pq_diagonal_.reset();
    LCAOWavefunction<Tile, Policy>::obsolete();
    ref_wfn_->obsolete();
  }

  // get T1 amplitudes
  TArray t1() const {
    if (T1_.is_initialized()) {
      return T1_;
    } else {
      throw std::runtime_error("CCSD T1 amplitudes have not been computed");
    }
  }
  // get T2 amplitudes
  TArray t2() const {
    if (T2_.is_initialized()) {
      return T2_;
    } else {
      throw std::runtime_error("CCSD T2 amplitudes have not been computed");
    }
  }

  bool is_df() const { return df_; }

  void set_t1(const TArray &t1) { T1_ = t1; }

  void set_t2(const TArray &t2) { T2_ = t2; }

  bool print_detail() const { return verbose_; }

  void set_target_precision(double precision) { target_precision_ = precision; }

  const typename AOFactory::DirectTArray &get_direct_ao_integral() const {
    return direct_ao_array_;
  }

 protected:
  bool can_evaluate(Energy *energy) override {
    // can only evaluate the energy
    return energy->order() == 0;
  }

  void evaluate(Energy *energy) override {
    auto target_precision = energy->target_precision(0);
    // compute only if never computed, or requested with higher precision than
    // before
    if (!this->computed() || computed_precision_ > target_precision) {
      // compute reference to higher precision than required of correlation
      // energy
      auto target_ref_precision = target_precision / 100.;
      auto ref_energy =
          std::make_shared<Energy>(ref_wfn_, target_ref_precision);
      ::mpqc::evaluate(*ref_energy, ref_wfn_);

      auto &world = this->wfn_world()->world();
      auto time0 = mpqc::fenced_now(world);

      this->init_sdref(ref_wfn_, target_ref_precision);

      f_pq_diagonal_ = make_diagonal_fpq(this->lcao_factory(),
                                         this->ao_factory(), this->is_df());

      // set up the solver
      if (solver_str_ == "jacobi_diis" || solver_str_ == "exact_pno") {
        const auto n_occ = this->trange1_engine()->get_occ();
        const auto n_frozen = this->trange1_engine()->get_nfrozen();
        const auto n_act_occ = n_occ - n_frozen;
        const auto n_uocc = f_pq_diagonal_->rows() - n_occ;
        solver_ = std::make_shared<cc::JacobiDIISSolver<TArray>>(
            kv_, f_pq_diagonal_->segment(n_frozen, n_act_occ),
            f_pq_diagonal_->segment(n_occ, n_uocc));
      }
      else if (solver_str_ == "pno")
        solver_ = std::make_shared<cc::PNOSolver<TArray,typename LCAOFactory<Tile, Policy>::DirectTArray>>(kv_, this->lcao_factory());
      else
        throw ProgrammingError("unknown solver string", __FILE__, __LINE__);

      // set the precision
      target_precision_ = energy->target_precision(0);

      TArray t1;
      TArray t2;

      if (method_ == "standard") {
        ccsd_corr_energy_ = compute_ccsd_conventional(t1, t2);
      } else if (method_ == "df") {
        ccsd_corr_energy_ = compute_ccsd_df(t1, t2);
      } else if (method_ == "direct" || method_ == "direct_df") {
        // initialize direct integral class
        direct_ao_array_ =
            this->ao_factory().compute_direct(L"(μ ν| G|κ λ)");
        ccsd_corr_energy_ = compute_ccsd_direct(t1, t2);
      }

      T1_ = t1;
      T2_ = t2;

      this->computed_ = true;
      this->set_value(energy, ref_energy->energy() + ccsd_corr_energy_);

      auto time1 = mpqc::fenced_now(world);
      auto duration0 = mpqc::duration_in_s(time0, time1);
      ExEnv::out0() << "CCSD Time in CCSD: " << duration0 << " S" << std::endl;
    }
  }

 protected:
  // store all the integrals in memory
  // used as reference for development
  double compute_ccsd_conventional(TArray &t1, TArray &t2) {
    auto &world = this->wfn_world()->world();
    bool accurate_time = this->lcao_factory().accurate_time();

    if (world.rank() == 0) {
      std::cout << "Use Conventional CCSD Compute" << std::endl;
    }

    auto tmp_time0 = mpqc::now(world, accurate_time);
    // get all two electron integrals
    TArray g_abcd = this->get_abcd();
    TArray g_ijab = this->get_ijab();
    TArray g_ijkl = this->get_ijkl();
    TArray g_iajb = this->get_iajb();
    TArray g_iabc = this->get_iabc();
    TArray g_aibc = this->get_aibc();
    TArray g_ijak = this->get_ijak();
    TArray g_ijka = this->get_ijka();
    this->lcao_factory().registry().purge_formula(L"(i ν| G |κ λ )");

    auto tmp_time1 = mpqc::now(world, accurate_time);
    auto tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
    mpqc::utility::print_par(world, "Integral Prepare Time: ", tmp_time, "\n");

    TArray f_ai = this->get_fock_ai();
    TArray f_ij = this->get_fock_ij();
    TArray f_ab = this->get_fock_ab();

    // initial guess = 0
    t1 = TArray(f_ai.world(), f_ai.trange(), f_ai.shape(), f_ai.pmap());
    t1.fill(0.0);
    TArray g_abij;
    g_abij("a,b,i,j") = g_ijab("i,j,a,b");
    t2 = TArray(g_abij.world(), g_abij.trange(), g_abij.shape(), g_abij.pmap());
    t2.fill(0.0);

    TArray tau;
    tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");

    double E0;
    double E1 = 0.0;
    double dE;

    // optimize t1 and t2
    std::size_t iter = 0ul;
    double error = 1.0;
    TArray r1(f_ai);
    TArray r2(g_abij);

    if (world.rank() == 0) {
      std::cout << "Start Iteration" << std::endl;
      std::cout << "Max Iteration: " << max_iter_ << std::endl;
      std::cout << "Target precision: " << target_precision_ << std::endl;
      std::cout << "AccurateTime: " << accurate_time << std::endl;
      std::cout << "PrintDetail: " << verbose_ << std::endl;
      std::cout << "Reduced ABCD Memory Approach: "
                << (reduced_abcd_memory_ ? "Yes" : "No") << std::endl;
    }

    // CCSD solver loop
    mpqc::time_point time0;
    while (iter < max_iter_) {
      TArray::wait_for_lazy_cleanup(world);

      // zero out singles if want CCD
      //TArray r1_new(r1.world(), r1.trange(), r1.shape()); r1_new.fill(0.0); r1("a,i") = r1_new("a,i");
      error = solver_->error(r1, r2);

      // recompute energy
      E0 = E1;
      E1 = 2.0 * TA::dot(f_ai("a,i") + r1("a,i"), t1("a,i")) +
           TA::dot(g_abij("a,b,i,j") + r2("a,b,i,j"),
                   2 * tau("a,b,i,j") - tau("b,a,i,j"));
      dE = std::abs(E0 - E1);

      if (dE >= target_precision_ || error >= target_precision_ || iter == 0) {
        tmp_time0 = mpqc::now(world, accurate_time);

        assert(solver_);
        solver_->update(t1, t2, r1, r2, dE);

        if (verbose_) {
          mpqc::detail::print_size_info(r2, "R2");
          mpqc::detail::print_size_info(t2, "T2");
        }

        // recompute tau
        tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");
        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (verbose_) {
          mpqc::utility::print_par(world, "Solver::update time: ", tmp_time, "\n");
        }

        // log the iteration
        auto time1 = mpqc::fenced_now(world);
        if (world.rank() == 0 && iter > 0) {
          auto duration = mpqc::duration_in_s(time0, time1);
          detail::print_ccsd(iter, dE, error, E1, duration);
        }
      } else {  // break out of the solver loop, if converged
        // log the iteration
        auto time1 = mpqc::fenced_now(world);
        if (world.rank() == 0) {
          MPQC_ASSERT(iter > 0);
          auto duration = mpqc::duration_in_s(time0, time1);
          detail::print_ccsd(iter, dE, error, E1, duration);
        }

        // When solver_str = exact_pno, switch to PNO solver after CCSD
        // amplitudes are converged
        if (solver_str_ == "exact_pno"){
          ExEnv::out0() << "Switching now to PNO solver" << std::endl;
          solver_ = std::make_shared<cc::PNOSolver<TArray,typename LCAOFactory<Tile, Policy>::DirectTArray>>(kv_, this->lcao_factory());
          assert(solver_);
          solver_->update(t1, t2, r1, r2, dE);

          // recompute tau
          tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");
          tmp_time1 = mpqc::now(world, accurate_time);
          tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
          if (verbose_) {
            mpqc::utility::print_par(world, "Solver::update time: ", tmp_time, "\n");
          }

          ExEnv::out0() << "Recomputing energy after PNO compression" << std::endl;
          E0 = E1;
          E1 = 2.0 * TA::dot(f_ai("a,i") + r1("a,i"), t1("a,i")) +
               TA::dot(g_abij("a,b,i,j") + r2("a,b,i,j"),
                       2 * tau("a,b,i,j") - tau("b,a,i,j"));
          dE = std::abs(E0 - E1);

          time1 = mpqc::fenced_now(world);
          if (world.rank() == 0) {
            MPQC_ASSERT(iter > 0);
            auto duration = mpqc::duration_in_s(time0, time1);
            detail::print_ccsd(iter+1, dE, error, E1, duration);
          }
        } // end if solver_str_ == "exact_pno"

        break;
      }

      // start iteration timer
      time0 = mpqc::fenced_now(world);

      auto t1_time0 = mpqc::now(world, accurate_time);
      TArray h_ki, h_ac;
      {
        // intermediates for t1
        // external index i and a
        // vir index a b c d
        // occ index i j k l
        TArray h_kc;

        // compute residual r1(n) (at convergence r1 = 0)
        // external index i and a
        tmp_time0 = mpqc::now(world, accurate_time);
        r1("a,i") = f_ai("a,i") - 2.0 * (f_ai("c,k") * t1("c,i")) * t1("a,k");

        {
          h_ac("a,c") =
              f_ab("a,c") -
              (2.0 * g_ijab("k,l,c,d") - g_ijab("l,k,c,d")) * tau("a,d,k,l");
          r1("a,i") += h_ac("a,c") * t1("c,i");
        }

        {
          h_ki("k,i") =
              f_ij("k,i") +
              (2.0 * g_ijab("k,l,c,d") - g_ijab("k,l,d,c")) * tau("c,d,i,l");
          r1("a,i") -= t1("a,k") * h_ki("k,i");
        }

        {
          h_kc("k,c") =
              f_ai("c,k") +
              (2.0 * g_ijab("k,l,c,d") - g_ijab("k,l,d,c")) * t1("d,l");
          r1("a,i") += h_kc("k,c") * (2.0 * t2("c,a,k,i") - t2("c,a,i,k") +
                                      t1("c,i") * t1("a,k"));
        }

        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (verbose_) {
          mpqc::utility::print_par(world, "t1 h term time: ", tmp_time, "\n");
        }

        tmp_time0 = mpqc::now(world, accurate_time);
        r1("a,i") += (2.0 * g_ijab("k,i,c,a") - g_iajb("k,a,i,c")) * t1("c,k");

        r1("a,i") +=
            (2.0 * g_iabc("k,a,c,d") - g_iabc("k,a,d,c")) * tau("c,d,k,i");

        r1("a,i") -=
            (2.0 * g_ijak("k,l,c,i") - g_ijak("l,k,c,i")) * tau("c,a,k,l");

        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (verbose_) {
          mpqc::utility::print_par(world, "t1 other time: ", tmp_time, "\n");
        }
      }
      auto t1_time1 = mpqc::now(world, accurate_time);
      auto t1_time = mpqc::duration_in_s(t1_time0, t1_time1);
      if (verbose_) {
        mpqc::utility::print_par(world, "t1 total time: ", t1_time, "\n");
      }

      // intermediates for t2
      // external index i j a b

      auto t2_time0 = mpqc::now(world, accurate_time);

      // compute residual r2(n) (at convergence r2 = 0)

      // permutation part
      tmp_time0 = mpqc::now(world, accurate_time);

      {
        r2("a,b,i,j") =
            (g_iabc("i,c,a,b") - g_iajb("k,b,i,c") * t1("a,k")) * t1("c,j");

        r2("a,b,i,j") -=
            (g_ijak("i,j,a,k") + g_ijab("i,k,a,c") * t1("c,j")) * t1("b,k");
      }
      tmp_time1 = mpqc::now(world, accurate_time);
      tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
      if (verbose_) {
        mpqc::utility::print_par(world, "t2 other time: ", tmp_time, "\n");
      }

      tmp_time0 = mpqc::now(world, accurate_time);
      {
        // compute g intermediate
        TArray g_ki, g_ac;

        g_ki("k,i") = h_ki("k,i") + f_ai("c,k") * t1("c,i") +
                      (2.0 * g_ijka("k,l,i,c") - g_ijka("l,k,i,c")) * t1("c,l");

        g_ac("a,c") = h_ac("a,c") - f_ai("c,k") * t1("a,k") +
                      (2.0 * g_aibc("a,k,c,d") - g_aibc("a,k,d,c")) * t1("d,k");

        r2("a,b,i,j") +=
            g_ac("a,c") * t2("c,b,i,j") - g_ki("k,i") * t2("a,b,k,j");
      }
      tmp_time1 = mpqc::now(world, accurate_time);
      tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
      if (verbose_) {
        mpqc::utility::print_par(world, "t2 g term time: ", tmp_time, "\n");
      }

      tmp_time0 = mpqc::now(world, accurate_time);
      {
        TArray j_akic;
        TArray k_kaic;
        // compute j and k intermediate
        {
          TArray T;

          T("d,b,i,l") = 0.5 * t2("d,b,i,l") + t1("d,i") * t1("b,l");

          j_akic("a,k,i,c") = g_ijab("i,k,a,c");

          j_akic("a,k,i,c") -= g_ijka("l,k,i,c") * t1("a,l");

          j_akic("a,k,i,c") += g_aibc("a,k,d,c") * t1("d,i");

          j_akic("a,k,i,c") -= g_ijab("k,l,c,d") * T("d,a,i,l");

          j_akic("a,k,i,c") += 0.5 *
                               (2.0 * g_ijab("k,l,c,d") - g_ijab("k,l,d,c")) *
                               t2("a,d,i,l");

          k_kaic("k,a,i,c") = g_iajb("k,a,i,c")

                              - g_ijka("k,l,i,c") * t1("a,l")

                              + g_iabc("k,a,d,c") * t1("d,i")

                              - g_ijab("k,l,d,c") * T("d,a,i,l");
          if (verbose_) {
            mpqc::detail::print_size_info(T, "T");
            mpqc::detail::print_size_info(j_akic, "J_akic");
            mpqc::detail::print_size_info(k_kaic, "K_kaic");
          }
        }

        r2("a,b,i,j") += 0.5 * (2.0 * j_akic("a,k,i,c") - k_kaic("k,a,i,c")) *
                         (2.0 * t2("c,b,k,j") - t2("b,c,k,j"));

        r2("a,b,i,j") += -0.5 * k_kaic("k,a,i,c") * t2("b,c,k,j") -
                         k_kaic("k,b,i,c") * t2("a,c,k,j");
      }
      tmp_time1 = mpqc::now(world, accurate_time);
      tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
      if (verbose_) {
        mpqc::utility::print_par(world, "t2 j,k term time: ", tmp_time, "\n");
      }

      // perform the permutation
      r2("a,b,i,j") = r2("a,b,i,j") + r2("b,a,j,i");

      r2("a,b,i,j") += g_ijab("i,j,a,b");

      tmp_time0 = mpqc::now(world, accurate_time);
      {
        TArray a_klij;
        // compute a intermediate
        a_klij("k,l,i,j") = g_ijkl("k,l,i,j");

        a_klij("k,l,i,j") += g_ijka("k,l,i,c") * t1("c,j");

        a_klij("k,l,i,j") += g_ijak("k,l,c,j") * t1("c,i");

        a_klij("k,l,i,j") += g_ijab("k,l,c,d") * tau("c,d,i,j");

        r2("a,b,i,j") += a_klij("k,l,i,j") * tau("a,b,k,l");

        if (verbose_) {
          mpqc::detail::print_size_info(a_klij, "A_klij");
        }
      }
      tmp_time1 = mpqc::now(world, accurate_time);
      tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
      if (verbose_) {
        mpqc::utility::print_par(world, "t2 a term time: ", tmp_time, "\n");
      }

      tmp_time0 = mpqc::now(world, accurate_time);
      {
        // compute b intermediate
        if (reduced_abcd_memory_) {
          // avoid store b_abcd
          TArray b_abij;
          b_abij("a,b,i,j") = g_abcd("a,b,c,d") * tau("c,d,i,j");

          b_abij("a,b,i,j") -= g_aibc("a,k,c,d") * tau("c,d,i,j") * t1("b,k");

          b_abij("a,b,i,j") -= g_iabc("k,b,c,d") * tau("c,d,i,j") * t1("a,k");

          if (verbose_) {
            mpqc::detail::print_size_info(b_abij, "B_abij");
          }

          r2("a,b,i,j") += b_abij("a,b,i,j");
        } else {
          TArray b_abcd;

          b_abcd("a,b,c,d") = g_abcd("a,b,c,d") -
                              g_aibc("a,k,c,d") * t1("b,k") -
                              g_iabc("k,b,c,d") * t1("a,k");

          if (verbose_) {
            mpqc::detail::print_size_info(b_abcd, "B_abcd");
          }

          r2("a,b,i,j") += b_abcd("a,b,c,d") * tau("c,d,i,j");
        }
      }
      tmp_time1 = mpqc::now(world, accurate_time);
      tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
      if (verbose_) {
        mpqc::utility::print_par(world, "t2 b term time: ", tmp_time, "\n");
      }

      auto t2_time1 = mpqc::now(world, accurate_time);
      auto t2_time = mpqc::duration_in_s(t2_time0, t2_time1);
      if (verbose_) {
        mpqc::utility::print_par(world, "t2 total time: ", t2_time, "\n");
      }

      ++iter;
    }  // CCSD solver loop
    if (iter >= max_iter_) {
      utility::print_par(this->wfn_world()->world(),
                         "\n Warning!! Exceed Max Iteration! \n");
    }
    if (world.rank() == 0) {
      std::cout << "CCSD Energy  " << E1 << std::endl;
    }
    return E1;
  }

 private:
  double compute_ccsd_df(TArray &t1, TArray &t2) {
    auto &world = this->wfn_world()->world();
    bool accurate_time = this->lcao_factory().accurate_time();

    if (world.rank() == 0) {
      std::cout << "Use DF CCSD Compute" << std::endl;
    }

    auto tmp_time0 = mpqc::now(world, accurate_time);
    // get all two electron integrals
    TArray g_ijab = this->get_ijab();
    TArray g_ijkl = this->get_ijkl();
    TArray g_iajb = this->get_iajb();
    TArray g_ijak = this->get_ijak();
    TArray g_ijka = this->get_ijka();

    TArray g_abcd;
    TArray g_iabc;
    if (!reduced_abcd_memory_) {
      g_abcd = this->get_abcd();
      g_iabc = this->get_iabc();
    }

    TArray X_ai = this->get_Xai();
    TArray X_ij = this->get_Xij();
    TArray X_ab = this->get_Xab();

    TArray f_ai = this->get_fock_ai();
    TArray f_ij = this->get_fock_ij();
    TArray f_ab = this->get_fock_ab();

    auto tmp_time1 = mpqc::now(world, accurate_time);
    auto tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
    mpqc::utility::print_par(world, "Integral Prepare Time: ", tmp_time, "\n");

    // initial guess = 0
    t1 = TArray(f_ai.world(), f_ai.trange(), f_ai.shape(), f_ai.pmap());
    t1.fill(0.0);
    TArray g_abij;
    g_abij("a,b,i,j") = g_ijab("i,j,a,b");
    t2 = TArray(g_abij.world(), g_abij.trange(), g_abij.shape(), g_abij.pmap());
    t2.fill(0.0);

    TArray tau;
    tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");

    double E0;
    double E1 = 0.0;
    double dE;

    // optimize t1 and t2
    std::size_t iter = 0ul;
    double error = 1.0;
    TArray r1(f_ai);
    TArray r2(g_abij);

    if (world.rank() == 0) {
      std::cout << "Start Iteration" << std::endl;
      std::cout << "Max Iteration: " << max_iter_ << std::endl;
      std::cout << "Target Precision: " << target_precision_ << std::endl;
      std::cout << "AccurateTime: " << accurate_time << std::endl;
      std::cout << "PrintDetail: " << verbose_ << std::endl;
      std::cout << "Reduced ABCD Memory Approach: "
                << (reduced_abcd_memory_ ? "Yes" : "No") << std::endl;
    }

    // CCSD solver loop
    mpqc::time_point time0;
    while (iter < max_iter_) {
      TArray::wait_for_lazy_cleanup(world);

      // zero out singles if want CCD
      //TArray r1_new(r1.world(), r1.trange(), r1.shape()); r1_new.fill(0.0); r1("a,i") = r1_new("a,i");
      error = solver_->error(r1, r2);

      // recompute energy
      E0 = E1;
      E1 = 2.0 * TA::dot(f_ai("a,i") + r1("a,i"), t1("a,i")) +
           TA::dot(g_abij("a,b,i,j") + r2("a,b,i,j"),
                   2 * tau("a,b,i,j") - tau("b,a,i,j"));
      dE = std::abs(E0 - E1);

      if (dE >= target_precision_ || error >= target_precision_ || iter == 0) {
        tmp_time0 = mpqc::now(world, accurate_time);

        assert(solver_);
        solver_->update(t1, t2, r1, r2, dE);

        if (verbose_) {
          mpqc::detail::print_size_info(r2, "R2");
          mpqc::detail::print_size_info(t2, "T2");
        }

        // recompute tau
        tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");
        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (verbose_) {
          mpqc::utility::print_par(world, "Solver::update time: ", tmp_time, "\n");
        }

        // log the iteration
        auto time1 = mpqc::fenced_now(world);
        if (world.rank() == 0 && iter > 0) {
          auto duration = mpqc::duration_in_s(time0, time1);
          detail::print_ccsd(iter, dE, error, E1, duration);
        }
      } else {  // break out of the solver loop, if converged
        // log the iteration
        auto time1 = mpqc::fenced_now(world);
        if (world.rank() == 0) {
          MPQC_ASSERT(iter > 0);
          auto duration = mpqc::duration_in_s(time0, time1);
          detail::print_ccsd(iter, dE, error, E1, duration);
        }

        // When solver_str = exact_pno, switch to PNO solver after CCSD
        // amplitudes are converged
        if (solver_str_ == "exact_pno"){
          ExEnv::out0() << "Switching now to PNO solver" << std::endl;
          solver_ = std::make_shared<cc::PNOSolver<TArray,typename LCAOFactory<Tile, Policy>::DirectTArray>>(kv_, this->lcao_factory());
          assert(solver_);
          solver_->update(t1, t2, r1, r2, dE);

          // recompute tau
          tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");
          tmp_time1 = mpqc::now(world, accurate_time);
          tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
          if (verbose_) {
            mpqc::utility::print_par(world, "Solver::update time: ", tmp_time, "\n");
          }

          ExEnv::out0() << "Recomputing energy after PNO compression" << std::endl;
          E0 = E1;
          E1 = 2.0 * TA::dot(f_ai("a,i") + r1("a,i"), t1("a,i")) +
               TA::dot(g_abij("a,b,i,j") + r2("a,b,i,j"),
                       2 * tau("a,b,i,j") - tau("b,a,i,j"));
          dE = std::abs(E0 - E1);

          time1 = mpqc::fenced_now(world);
          if (world.rank() == 0) {
            MPQC_ASSERT(iter > 0);
            auto duration = mpqc::duration_in_s(time0, time1);
            detail::print_ccsd(iter+1, dE, error, E1, duration);
          }
        } // end if solver_str_ == "exact_pno"

        break;
      }

      // start iteration timer
      time0 = mpqc::fenced_now(world);

      auto t1_time0 = mpqc::now(world, accurate_time);
      TArray h_ki, h_ac;
      {
        // intermediates for t1
        // external index i and a
        // vir index a b c d
        // occ index i j k l
        TArray h_kc;

        TArray X_ai_tau;

        X_ai_tau("X,a,i") = 2.0*X_ai("X,b,j")*tau("a,b,i,j") - X_ai("X,b,j")*tau("a,b,j,i");

        // compute residual r1(n) (at convergence r1 = 0)
        // external index i and a
        tmp_time0 = mpqc::now(world, accurate_time);
        r1("a,i") = f_ai("a,i") - 2.0 * f_ai("c,k") * t1("c,i") * t1("a,k");

        {
          h_ac("a,c") =
              f_ab("a,c") - X_ai_tau("K,a,l") * X_ai("K,c,l");
          r1("a,i") += h_ac("a,c") * t1("c,i");
        }

        {
          h_ki("k,i") =
              f_ij("k,i") + X_ai_tau("K,c,i")*X_ai("K,c,k");
          r1("a,i") -= t1("a,k") * h_ki("k,i");
        }

        {
          h_kc("k,c") =
              f_ai("c,k") +
              (2.0 * g_ijab("k,l,c,d") - g_ijab("k,l,d,c")) * t1("d,l");
          r1("a,i") += h_kc("k,c") * (2.0 * t2("c,a,k,i") - t2("c,a,i,k") +
                                      t1("c,i") * t1("a,k"));
        }

        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (verbose_) {
          mpqc::utility::print_par(world, "t1 h term time: ", tmp_time, "\n");
        }

        tmp_time0 = mpqc::now(world, accurate_time);
        r1("a,i") += (2.0 * g_ijab("k,i,c,a") - g_iajb("k,a,i,c")) * t1("c,k");

        r1("a,i") += X_ai_tau("K,c,i") * X_ab("K,a,c");

        r1("a,i") -= X_ai_tau("K,a,l") * X_ij("K,l,i");

        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (verbose_) {
          mpqc::utility::print_par(world, "t1 other time: ", tmp_time, "\n");
        }
      }
      auto t1_time1 = mpqc::now(world, accurate_time);
      auto t1_time = mpqc::duration_in_s(t1_time0, t1_time1);
      if (verbose_) {
        mpqc::utility::print_par(world, "t1 total time: ", t1_time, "\n");
      }

      // intermediates for t2
      // external index i j a b

      auto t2_time0 = mpqc::now(world, accurate_time);

      // compute residual r2(n) (at convergence r2 = 0)

      // permutation part
      tmp_time0 = mpqc::now(world, accurate_time);

      TArray X_ab_t1;
      X_ab_t1("K,a,i") = X_ab("K,a,b")*t1("b,i");

      {
        r2("a,b,i,j") =
            X_ab_t1("K,b,j") * (X_ai("K,a,i") - X_ij("K,k,i")*t1("a,k")) -
            (g_ijak("i,j,a,k") + g_ijab("i,k,a,c") * t1("c,j")) * t1("b,k");
      }
      tmp_time1 = mpqc::now(world, accurate_time);
      tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
      if (verbose_) {
        mpqc::utility::print_par(world, "t2 other time: ", tmp_time, "\n");
      }

      tmp_time0 = mpqc::now(world, accurate_time);
      {
        // compute g intermediate
        TArray g_ki, g_ac;

        g_ki("k,i") = h_ki("k,i") + f_ai("c,k") * t1("c,i") +
                      (2.0 * g_ijka("k,l,i,c") - g_ijka("l,k,i,c")) * t1("c,l");

        g_ac("a,c") = h_ac("a,c") - f_ai("c,k") * t1("a,k") +
                      2.0 * X_ai("K,d,k") * t1("d,k") * X_ab("K,a,c") - X_ab_t1("K,a,k") * X_ai("K,c,k");

        r2("a,b,i,j") +=
            g_ac("a,c") * t2("c,b,i,j") - g_ki("k,i") * t2("a,b,k,j");
      }
      tmp_time1 = mpqc::now(world, accurate_time);
      tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
      if (verbose_) {
        mpqc::utility::print_par(world, "t2 g term time: ", tmp_time, "\n");
      }

      tmp_time0 = mpqc::now(world, accurate_time);
      {
        TArray j_akic;
        TArray k_kaic;
        // compute j and k intermediate
        {
          TArray T;

          T("d,b,i,l") = 0.5 * t2("d,b,i,l") + t1("d,i") * t1("b,l");

          j_akic("a,k,i,c") =
              g_ijab("i,k,a,c")

              - g_ijka("l,k,i,c") * t1("a,l")

              + X_ab_t1("K,a,i") * X_ai("K,c,k")

              - X_ai("x,d,l") * T("d,a,i,l") * X_ai("x,c,k")

              +
              0.5 * (2.0 * g_ijab("k,l,c,d") - g_ijab("k,l,d,c")) *
                  t2("a,d,i,l");

          k_kaic("k,a,i,c") = g_iajb("k,a,i,c")

                              - g_ijka("k,l,i,c") * t1("a,l")

                              + X_ai("K,d,k") * t1("d,i") * X_ab("K,a,c")

                              - g_ijab("k,l,d,c") * T("d,a,i,l");
          if (verbose_) {
            mpqc::detail::print_size_info(T, "T");
            mpqc::detail::print_size_info(j_akic, "J_akic");
            mpqc::detail::print_size_info(k_kaic, "K_kaic");
          }
        }

        r2("a,b,i,j") += 0.5 * (2.0 * j_akic("a,k,i,c") - k_kaic("k,a,i,c")) *
                         (2.0 * t2("c,b,k,j") - t2("b,c,k,j"));

        TArray tmp;
        tmp("a,b,i,j") = k_kaic("k,a,i,c") * t2("c,b,j,k");
        r2("a,b,i,j") -= 0.5 * tmp("a,b,i,j") + tmp("b,a,i,j");
      }
      tmp_time1 = mpqc::now(world, accurate_time);
      tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
      if (verbose_) {
        mpqc::utility::print_par(world, "t2 j,k term time: ", tmp_time, "\n");
      }

      // perform the permutation
      r2("a,b,i,j") = r2("a,b,i,j") + r2("b,a,j,i");

      r2("a,b,i,j") += g_ijab("i,j,a,b");

      tmp_time0 = mpqc::now(world, accurate_time);
      {
        TArray a_klij;
        // compute a intermediate
        a_klij("k,l,i,j") = g_ijkl("k,l,i,j")

                            + g_ijka("k,l,i,c") * t1("c,j")

                            + g_ijak("k,l,c,j") * t1("c,i")

                            + g_ijab("k,l,c,d") * tau("c,d,i,j");

        r2("a,b,i,j") += a_klij("k,l,i,j") * tau("a,b,k,l");

        if (verbose_) {
          mpqc::detail::print_size_info(a_klij, "A_klij");
        }
      }
      tmp_time1 = mpqc::now(world, accurate_time);
      tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
      if (verbose_) {
        mpqc::utility::print_par(world, "t2 a term time: ", tmp_time, "\n");
      }

      tmp_time0 = mpqc::now(world, accurate_time);
      {
        // compute b intermediate
        // avoid store b_abcd
        TArray b_abij;
        if (!reduced_abcd_memory_) {

          b_abij("a,b,i,j") = tau("c,d,i,j") * g_abcd("a,b,c,d");
          TArray tmp;
          tmp("k,a,i,j") = g_iabc("k,a,c,d") * tau("c,d,i,j");
          b_abij("a,b,i,j") -= tmp("k,a,j,i") * t1("b,k");

          b_abij("a,b,i,j") -= tmp("k,b,i,j") * t1("a,k");

        } else {

          TArray X_ab_t1;
          X_ab_t1("K,a,b") = 0.5*X_ab("K,a,b") - X_ai("K,b,i")*t1("a,i");

          auto g_abcd_iabc_direct = gaussian::df_direct_integrals(X_ab, X_ab_t1, Formula::Notation::Physical);

          b_abij("a,b,i,j") = tau("c,d,i,j") * g_abcd_iabc_direct("a,b,c,d");

          b_abij("a,b,i,j") += b_abij("b,a,j,i");

        }

        if (verbose_) {
          mpqc::detail::print_size_info(b_abij, "B_abij");
        }

        r2("a,b,i,j") += b_abij("a,b,i,j");
      }
      tmp_time1 = mpqc::now(world, accurate_time);
      tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
      if (verbose_) {
        mpqc::utility::print_par(world, "t2 b term time: ", tmp_time, "\n");
      }

      auto t2_time1 = mpqc::now(world, accurate_time);
      auto t2_time = mpqc::duration_in_s(t2_time0, t2_time1);
      if (verbose_) {
        mpqc::utility::print_par(world, "t2 total time: ", t2_time, "\n");
      }

      ++iter;
    }  // CCSD solver loop

    if (iter >= max_iter_) {
      utility::print_par(this->wfn_world()->world(),
                         "\n Warning!! Exceed Max Iteration! \n");
    }
    if (world.rank() == 0) {
      std::cout << "CCSD Energy  " << E1 << std::endl;
    }
    return E1;
  }

  double compute_ccsd_direct(TArray &t1, TArray &t2) {
    auto &world = this->wfn_world()->world();
    // get mo coefficient
    TArray ca = this->get_Ca();
    TArray ci = this->get_Ci();

    if (this->df_) {
      mpqc::utility::print_par(world, "Use Direct-DF CCSD Compute \n");
    } else {
      mpqc::utility::print_par(world, "Use Direct CCSD Compute \n");
    }

    bool accurate_time = this->lcao_factory().accurate_time();

    auto tmp_time0 = mpqc::now(world, accurate_time);

    // get three center integral
    TArray Xab, Xij, Xai;
    if (this->df_) {
      Xab = this->get_Xab();
      Xij = this->get_Xij();
      Xai = this->get_Xai();
    }
    // get all integrals
    TArray g_ijkl = this->get_ijkl();
    TArray g_ijab = this->get_ijab();
    TArray g_iajb = this->get_iajb();
    TArray g_ijka = this->get_ijka();
    TArray g_ijak = this->get_ijak();
    TArray g_iabc;
    if (!this->df_) {
      g_iabc = this->get_iabc();
    }
    TArray f_ai = this->get_fock_ai();
    TArray f_ij = this->get_fock_ij();
    TArray f_ab = this->get_fock_ab();
    this->lcao_factory().registry().purge_formula(L"(i ν| G |κ λ )");

    auto tmp_time1 = mpqc::now(world, accurate_time);
    auto tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
    mpqc::utility::print_par(world, "Integral Prepare Time: ", tmp_time, "\n");

    // initial guess = 0
    t1 = TArray(f_ai.world(), f_ai.trange(), f_ai.shape(), f_ai.pmap());
    t1.fill(0.0);
    TArray g_abij;
    g_abij("a,b,i,j") = g_ijab("i,j,a,b");
    t2 = TArray(g_abij.world(), g_abij.trange(), g_abij.shape(), g_abij.pmap());
    t2.fill(0.0);

    TArray tau;
    tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");

    double E0;
    double E1 = 0.0;
    double dE;

    // optimize t1 and t2
    std::size_t iter = 0ul;
    double error = 1.0;
    TArray r1(f_ai);
    TArray r2(g_abij);

    if (world.rank() == 0) {
      std::cout << "Start Iteration" << std::endl;
      std::cout << "Max Iteration: " << max_iter_ << std::endl;
      std::cout << "Target Precision: " << target_precision_ << std::endl;
      std::cout << "AccurateTime: " << accurate_time << std::endl;
      std::cout << "PrintDetail: " << verbose_ << std::endl;
    };

    // CCSD solver loop
    mpqc::time_point time0;
    double duration_u;
    while (iter < max_iter_) {
      TArray::wait_for_lazy_cleanup(world);

      // zero out singles if want CCD
      //TArray r1_new(r1.world(), r1.trange(), r1.shape()); r1_new.fill(0.0); r1("a,i") = r1_new("a,i");
      error = solver_->error(r1, r2);

      // recompute energy
      E0 = E1;
      E1 = 2.0 * TA::dot(f_ai("a,i") + r1("a,i"), t1("a,i")) +
           TA::dot(g_abij("a,b,i,j") + r2("a,b,i,j"),
                   2 * tau("a,b,i,j") - tau("b,a,i,j"));
      dE = std::abs(E0 - E1);

      if (iter == 0 || iter < min_iter_ || !solver_->is_converged(target_precision_, error, dE)) {
        tmp_time0 = mpqc::now(world, accurate_time);

//        ExEnv::out0() << "iter: " << iter << ", dE: " << dE << ", DeltaE: " << DeltaE << ", error: " << error << std::endl;

        assert(solver_);
        solver_->update(t1, t2, r1, r2, E1);

        if (verbose_) {
          mpqc::detail::print_size_info(r2, "R2");
          mpqc::detail::print_size_info(t2, "T2");
        }

        // recompute tau
        tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");
        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (verbose_) {
          mpqc::utility::print_par(world, "Solver::update time: ", tmp_time, "\n");
        }

        // log the iteration
        auto time1 = mpqc::fenced_now(world);
        if (world.rank() == 0 && iter > 0) {
          auto duration_t = mpqc::duration_in_s(time0, time1);
          detail::print_ccsd_direct(iter, dE, error, E1, duration_u,
                                    duration_t);
        }
      } else {  // break out of the solver loop, if converged
        // log the iteration
        auto time1 = mpqc::fenced_now(world);
        if (world.rank() == 0) {
          MPQC_ASSERT(iter > 0);
          auto duration_t = mpqc::duration_in_s(time0, time1);
          detail::print_ccsd_direct(iter, dE, error, E1, duration_u,
                                    duration_t);
        }

        // When solver_str = exact_pno, switch to PNO solver after CCSD
        // amplitudes are converged
        if (solver_str_ == "exact_pno"){
          ExEnv::out0() << "Switching now to PNO solver" << std::endl;
          solver_ = std::make_shared<cc::PNOSolver<TArray,typename LCAOFactory<Tile, Policy>::DirectTArray>>(kv_, this->lcao_factory());
          assert(solver_);
          solver_->update(t1, t2, r1, r2, E1);

          // recompute tau
          tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");
          tmp_time1 = mpqc::now(world, accurate_time);
          tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
          if (verbose_) {
            mpqc::utility::print_par(world, "Solver::update time: ", tmp_time, "\n");
          }

          ExEnv::out0() << "Recomputing energy after PNO compression" << std::endl;
          E0 = E1;
          E1 = 2.0 * TA::dot(f_ai("a,i") + r1("a,i"), t1("a,i")) +
               TA::dot(g_abij("a,b,i,j") + r2("a,b,i,j"),
                       2 * tau("a,b,i,j") - tau("b,a,i,j"));
          dE = std::abs(E0 - E1);

          time1 = mpqc::fenced_now(world);
          if (world.rank() == 0) {
            MPQC_ASSERT(iter > 0);
            auto duration = mpqc::duration_in_s(time0, time1);
            detail::print_ccsd(iter+1, dE, error, E1, duration);
          }
        } // end if solver_str_ == "exact_pno"

        break;
      }

      // start iteration timer
      time0 = mpqc::fenced_now(world);

      TArray u2_u11;
      // compute half transformed intermediates
      auto tu0 = mpqc::now(world, accurate_time);
      { u2_u11 = this->compute_u2_u11(t2, t1); }
      auto tu1 = mpqc::now(world, accurate_time);
      duration_u = mpqc::duration_in_s(tu0, tu1);

      if (verbose_) {
        mpqc::detail::print_size_info(u2_u11, "U_aaoo");
        mpqc::utility::print_par(world, "u term time: ", duration_u, "\n");
      } else if (iter == 0) {
        mpqc::detail::print_size_info(u2_u11, "U_aaoo");
      }

      auto t1_time0 = mpqc::now(world, accurate_time);
      TArray h_ac, h_ki;
      {
        // intermediates for t1
        // external index i and a

        tmp_time0 = mpqc::now(world, accurate_time);
        h_ac("a,c") =
            f_ab("a,c") -
            (2.0 * g_ijab("k,l,c,d") - g_ijab("l,k,c,d")) * tau("a,d,k,l");

        h_ki("k,i") =
            f_ij("k,i") +
            (2.0 * g_ijab("k,l,c,d") - g_ijab("k,l,d,c")) * tau("c,d,i,l");

        // compute residual r1(n) (at convergence r1 = 0)
        // external index i and a

        r1("a,i") = f_ai("a,i") - 2.0 * f_ai("c,k") * t1("c,i") * t1("a,k");

        r1("a,i") += h_ac("a,c") * t1("c,i") - t1("a,k") * h_ki("k,i");

        {
          TArray h_kc;
          h_kc("k,c") =
              f_ai("c,k") +
              (-g_ijab("k,l,d,c") + 2.0 * g_ijab("k,l,c,d")) * t1("d,l");

          r1("a,i") += h_kc("k,c") * (2.0 * t2("c,a,k,i") - t2("c,a,i,k") +
                                      t1("a,k") * t1("c,i"));
        }

        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (verbose_) {
          mpqc::utility::print_par(world, "t1 h term time: ", tmp_time, "\n");
        }

        tmp_time0 = mpqc::now(world, accurate_time);

        r1("a,i") += (2.0 * g_ijab("k,i,c,a") - g_iajb("k,a,i,c")) * t1("c,k");

        r1("a,i") += (2.0 * u2_u11("p,r,k,i") - u2_u11("p,r,i,k")) * ci("p,k") *
                     ca("r,a");

        r1("a,i") -=
            (2.0 * g_ijak("k,l,c,i") - g_ijak("l,k,c,i")) * tau("c,a,k,l");

        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (verbose_) {
          mpqc::utility::print_par(world, "t1 other time: ", tmp_time, "\n");
        }
      }
      auto t1_time1 = mpqc::now(world, accurate_time);
      auto t1_time = mpqc::duration_in_s(t1_time0, t1_time1);
      if (verbose_) {
        mpqc::utility::print_par(world, "t1 total time: ", t1_time, "\n");
      }

      // intermediates for t2
      // external index i j a b
      auto t2_time0 = mpqc::now(world, accurate_time);
      {
        tmp_time0 = mpqc::now(world, accurate_time);

        // permutation term
        {
          r2("a,b,i,j") = -g_iajb("k,b,i,c") * t1("c,j") * t1("a,k");

          if (this->df_) {
            r2("a,b,i,j") += Xab("X,b,c") * t1("c,j") * Xai("X,a,i");
          } else {
            r2("a,b,i,j") += g_iabc("i,c,a,b") * t1("c,j");
          }

          r2("a,b,i,j") -=
              (g_ijak("i,j,a,k") + g_ijab("i,k,a,c") * t1("c,j")) * t1("b,k");
        }
        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (verbose_) {
          mpqc::utility::print_par(world, "t2 other time: ", tmp_time, "\n");
        }

        tmp_time0 = mpqc::now(world, accurate_time);
        {
          // intermediates g
          TArray g_ki, g_ac;

          g_ki("k,i") =
              h_ki("k,i") + f_ai("c,k") * t1("c,i") +
              (2.0 * g_ijka("k,l,i,c") - g_ijka("l,k,i,c")) * t1("c,l");

          if (this->df_) {
            g_ac("a,c") = h_ac("a,c") - f_ai("c,k") * t1("a,k")

                          + 2.0 * Xai("X,d,k") * t1("d,k") * Xab("X,a,c")

                          - Xab("X,a,d") * t1("d,k") * Xai("X,c,k");
          } else {
            g_ac("a,c") = h_ac("a,c") - f_ai("c,k") * t1("a,k") +
                          2.0 * g_iabc("k,a,d,c") * t1("d,k") -
                          g_iabc("k,a,c,d") * t1("d,k");
          }

          r2("a,b,i,j") +=
              (g_ac("a,c") * t2("c,b,i,j") - g_ki("k,i") * t2("a,b,k,j"));
        }
        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (verbose_) {
          mpqc::utility::print_par(world, "t2 g term time: ", tmp_time, "\n");
        }

        tmp_time0 = mpqc::now(world, accurate_time);
        {
          TArray j_akic;
          TArray k_kaic;
          // compute j and k intermediate
          {
            TArray T;

            T("d,b,i,l") = 0.5 * t2("d,b,i,l") + t1("d,i") * t1("b,l");

            if (this->df_) {
              j_akic("a,k,i,c") =
                  g_ijab("i,k,a,c")

                  - g_ijka("l,k,i,c") * t1("a,l")

                  + (Xab("X,a,d") * t1("d,i")) * Xai("X,c,k")

                  - (Xai("X,d,l") * T("d,a,i,l")) * Xai("X,c,k")

                  +
                  (g_ijab("k,l,c,d") - 0.5 * g_ijab("k,l,d,c")) * t2("a,d,i,l");

              k_kaic("k,a,i,c") = g_iajb("k,a,i,c")

                                  - g_ijka("k,l,i,c") * t1("a,l")

                                  + (Xai("X,d,k") * t1("d,i")) * Xab("X,a,c")

                                  - g_ijab("k,l,d,c") * T("d,a,i,l");
            } else {
              j_akic("a,k,i,c") =
                  g_ijab("i,k,a,c") - g_ijka("l,k,i,c") * t1("a,l")

                  + g_iabc("k,a,c,d") * t1("d,i")

                  - g_ijab("k,l,c,d") * T("d,a,i,l")

                  +
                  0.5 * (2.0 * g_ijab("k,l,c,d") - g_ijab("k,l,d,c")) *
                      t2("a,d,i,l");

              k_kaic("k,a,i,c") = g_iajb("k,a,i,c")

                                  - g_ijka("k,l,i,c") * t1("a,l")

                                  + g_iabc("k,a,d,c") * t1("d,i")

                                  - g_ijab("k,l,d,c") * T("d,a,i,l");
            }

            if (verbose_) {
              mpqc::detail::print_size_info(T, "T");
              mpqc::detail::print_size_info(j_akic, "J_akic");
              mpqc::detail::print_size_info(k_kaic, "K_kaic");
            }
          }

          r2("a,b,i,j") += 0.5 * (2.0 * j_akic("a,k,i,c") - k_kaic("k,a,i,c")) *
                           (2.0 * t2("c,b,k,j") - t2("b,c,k,j"));

          TArray tmp;
          tmp("a,b,i,j") = k_kaic("k,a,i,c") * t2("c,b,j,k");
          r2("a,b,i,j") -= 0.5 * tmp("a,b,i,j") + tmp("b,a,i,j");
        }
        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (verbose_) {
          mpqc::utility::print_par(world, "t2 j,k term time: ", tmp_time, "\n");
        }

        // perform permutation
        r2("a,b,i,j") = r2("a,b,i,j") + r2("b,a,j,i");

        r2("a,b,i,j") += g_ijab("i,j,a,b");

        tmp_time0 = mpqc::now(world, accurate_time);
        {
          // intermediate a
          TArray a_klij;
          a_klij("k,l,i,j") = g_ijkl("k,l,i,j")

                              + g_ijka("k,l,i,c") * t1("c,j")

                              + g_ijak("k,l,c,j") * t1("c,i")

                              + g_ijab("k,l,c,d") * tau("c,d,i,j");

          r2("a,b,i,j") += a_klij("k,l,i,j") * tau("a,b,k,l");

          if (verbose_) {
            mpqc::detail::print_size_info(a_klij, "A_klij");
          }
        }
        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (verbose_) {
          mpqc::utility::print_par(world, "t2 a term time: ", tmp_time, "\n");
        }

        tmp_time0 = mpqc::now(world, accurate_time);
        {
          TArray b_abij;

          b_abij("a,b,i,j") =
              (u2_u11("p,r,i,j") * ca("r,b") -
               ci("r,k") * t1("b,k") * u2_u11("p,r,i,j")) *
                  ca("p,a")

              - u2_u11("p,r,i,j") * ci("p,k") * ca("r,b") * t1("a,k");

          r2("a,b,i,j") += b_abij("a,b,i,j");

          if (verbose_) {
            mpqc::detail::print_size_info(b_abij, "B_abij");
          }
        }
        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (verbose_) {
          mpqc::utility::print_par(world, "t2 b term time: ", tmp_time, "\n");
        }
      }

      auto t2_time1 = mpqc::now(world, accurate_time);
      auto t2_time = mpqc::duration_in_s(t2_time0, t2_time1);
      if (verbose_) {
        mpqc::utility::print_par(world, "t2 total time: ", t2_time, "\n");
      }

      ++iter;
    }  // CCSD solver loop
    if (iter >= max_iter_) {
      utility::print_par(this->wfn_world()->world(),
                         "\n Warning!! Exceed Max Iteration! \n");
    }
    if (world.rank() == 0) {
      std::cout << "CCSD Energy     " << E1 << std::endl;
    }
    return E1;
  }

 protected:
  /// get mo coefficient
  /// occ part
  const TArray get_Ci() {
    return this->lcao_factory()
        .orbital_registry()
        .retrieve(OrbitalIndex(L"i"))
        .coefs();
  }

  /// vir part
  const TArray get_Ca() {
    return this->lcao_factory()
        .orbital_registry()
        .retrieve(OrbitalIndex(L"a"))
        .coefs();
  }

  /// get three center integral (X|ab)
  const TArray get_Xab() {
    return this->lcao_factory().compute(L"(Κ|G|a b)[inv_sqr]");
  }

  /// get three center integral (X|ij)
  const TArray get_Xij() {
    return this->lcao_factory().compute(L"(Κ|G|i j)[inv_sqr]");
  }

  /// get three center integral (X|ai)
  const TArray get_Xai() {
    return this->lcao_factory().compute(L"(Κ|G|a i)[inv_sqr]");
  }

  // get two electron integrals
  // using physical notation <ab|ij>

  /// <ij|ab>
  virtual const TArray get_ijab() {
    std::wstring postfix = df_ ? L"[df]" : L"";
    return this->lcao_factory().compute(L"<i j|G|a b>" + postfix);
  }

  /// <ij|kl>
  virtual const TArray get_ijkl() {
    std::wstring postfix = df_ ? L"[df]" : L"";
    return this->lcao_factory().compute(L"<i j|G|k l>" + postfix);
  }

  /// <ab|cd>
  virtual const TArray get_abcd() {
    std::wstring postfix = df_ ? L"[df]" : L"";
    return this->lcao_factory().compute(L"<a b|G|c d>" + postfix);
  }

  /// <ia|bc>
  virtual const TArray get_iabc() {
    std::wstring postfix = df_ ? L"[df]" : L"";
    return this->lcao_factory().compute(L"<i a|G|b c>" + postfix);
  }

  /// <ai|bc>
  virtual const TArray get_aibc() {
    std::wstring postfix = df_ ? L"[df]" : L"";
    return this->lcao_factory().compute(L"<a i|G|b c>" + postfix);
  }

  /// <ij|ak>
  virtual const TArray get_ijak() {
    std::wstring postfix = df_ ? L"[df]" : L"";
    return this->lcao_factory().compute(L"<i j|G|a k>" + postfix);
  }

  /// <ia|jb>
  virtual const TArray get_iajb() {
    std::wstring postfix = df_ ? L"[df]" : L"";
    return this->lcao_factory().compute(L"<i a|G|j b>" + postfix);
  }

  /// <ij|ka>
  virtual const TArray get_ijka() {
    std::wstring postfix = df_ ? L"[df]" : L"";
    return this->lcao_factory().compute(L"<i j|G|k a>" + postfix);
  }

  /// <a|f|i>
  virtual const TArray get_fock_ai() {
    std::wstring postfix = df_ ? L"[df]" : L"";
    return this->lcao_factory().compute(L"<a|F|i>" + postfix);
  }

  /// <i|f|j>
  virtual const TArray get_fock_ij() {
    std::wstring postfix = df_ ? L"[df]" : L"";
    return this->lcao_factory().compute(L"<i|F|j>" + postfix);
  }

  /// <a|f|b>
  virtual const TArray get_fock_ab() {
    std::wstring postfix = df_ ? L"[df]" : L"";
    return this->lcao_factory().compute(L"<a|F|b>" + postfix);
  }

  /// <p|f|q>
  virtual const TArray get_fock_pq() {
    std::wstring postfix = df_ ? L"[df]" : L"";
    return this->lcao_factory().compute(L"<p|F|q>" + postfix);
  }

 private:
  /// AO integral-direct computation of (ab|cd) ints contributions to the
  /// doubles residual

  /// computes \f$ U^{ij}_{\rho\sigma} \equiv \left( t^{ij}_{\mu \nu} +
  /// t^{i}_{\mu} t^{j}_{\nu} \right) (\mu \rho| \nu \sigma) \f$
  /// @param t2 doubles amplitudes in MO basis
  /// @param t1 singles amplitudes in MO basis
  /// @return U tensor
  TArray compute_u2_u11(const TArray &t2, const TArray &t1) {
    if (direct_ao_array_.array().is_initialized()) {
      TArray Ca = get_Ca();
      TArray tc;
      TArray u2_u11;
      tc("i,q") = Ca("q,c") * t1("c,i");
      u2_u11("p, r, i, j") =
          ((t2("a,b,i,j") * Ca("q,a")) * Ca("s,b") + tc("i,q") * tc("j,s")) *
          direct_ao_array_("p,q,r,s");
//      u2_u11("p, r, i, j") =
//          0.5 * (u2_u11("p, r, i, j") + u2_u11("r, p, j, i"));
      return u2_u11;
    } else {
      throw ProgrammingError(
          "CCSD: integral-direct implementation used, but direct integral not "
          "initialized",
          __FILE__, __LINE__);
    }
  }
};  // class CCSD

#if TA_DEFAULT_POLICY == 0
extern template class CCSD<TA::TensorD, TA::DensePolicy>;
#elif TA_DEFAULT_POLICY == 1
extern template class CCSD<TA::TensorD, TA::SparsePolicy>;
#endif
}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_CCSD_H_
