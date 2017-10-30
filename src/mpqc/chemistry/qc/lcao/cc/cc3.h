// CC3 by Varun Rishi, September 2017.
// Adapted from CCSD code written by Chong Peng on 7/1/15.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_CC3_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_CC3_H_

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

inline void print_CC3(int iter, double dE, double error, double E1,
                       double time) {
  if (iter == 0) {
    ExEnv::out0() << mpqc::printf("%3s \t %10s \t %10s \t %15s \t %10s \n",
                                  "iter", "deltaE", "residual", "energy",
                                  "total time/s");
  }
  ExEnv::out0() << mpqc::printf(
      "%3i \t %10.5e \t %10.5e \t %15.12f \t %10.1f \n", iter, dE, error, E1,
      time);
}

inline void print_CC3_direct(int iter, double dE, double error, double E1,
                              double time1, double time2) {
  if (iter == 0) {
    ExEnv::out0() << mpqc::printf(
        "%3s \t %10s \t %10s \t %15s \t %10s \t %10s \n", "iter", "deltaE",
        "residual", "energy", "u time/s", "total time/s");
  }
  ExEnv::out0() << mpqc::printf(
      "%3i \t %10.5e \t %10.5e \t %15.12f \t %10.1f \t %10.1f \n", iter, dE,
      error, E1, time1, time2);
}
}

/**
 * CC3 class that computes CC3 energy
 */

template <typename Tile, typename Policy>
class CC3 : public LCAOWavefunction<Tile, Policy>,
             public Provides<Energy>,
             public std::enable_shared_from_this<CC3<Tile, Policy>> {
 public:
  using TArray = TA::DistArray<Tile, Policy>;
  using AOFactory = gaussian::AOFactory<Tile, Policy>;

  CC3() = default;

  // clang-format off

  /**
   * KeyVal constructor
   * @param kv
   *
   * keywords : all keywords for LCAOWavefunciton
   *
   * | Keyword | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | ref | Wavefunction | none | reference Wavefunction, need to be a Energy::Provider RHF for example |
   * | method | string | standard or df | method to compute ccsd (valid choices are: standard, direct, df, direct_df), the default depends on whether \c df_basis is provided |
   * | max_iter | int | 30 | maxmium iteration in CCSD |
   * | verbose | bool | default use factory.verbose() | if print more information in CCSD iteration |
   * | reduced_abcd_memory | bool | false | avoid store another abcd term in standard and df method |
   */

  // clang-format on

  CC3(const KeyVal &kv) : LCAOWavefunction<Tile, Policy>(kv), kv_(kv) {
    if (kv.exists("ref")) {
      ref_wfn_ = kv.class_ptr<Wavefunction>("ref");
    } else {
      throw InputError("Default Ref in CCSD is not support! \n", __FILE__,
                       __LINE__, "ref");
    }

    df_ = false;
    auto default_method =
        this->lcao_factory().basis_registry()->have(L"Κ") ? "df" : "standard";
    method_ = kv.value<std::string>("method", default_method);
    if (method_ != "df" && method_ != "direct" && method_ != "standard" &&
        method_ != "direct_df") {
      throw InputError("Invalid CCSD method! \n", __FILE__, __LINE__, "method");
    }
    if (method_ == "df" || method_ == "direct_df") {
      df_ = true;
    }

    solver_str_ = kv.value<std::string>("solver", "jacobi_diis");
    if (solver_str_ != "jacobi_diis" && solver_str_ != "pno")
      throw InputError("invalid value for solver keyword", __FILE__, __LINE__,
                       "solver");

    reduced_abcd_memory_ = kv.value<bool>("reduced_abcd_memory", false);

    max_iter_ = kv.value<int>("max_iter", 30);
    verbose_ = kv.value<bool>("verbose", this->lcao_factory().verbose());
  }

  virtual ~CC3() {}

  /// protected members
 protected:
  TArray T1_;
  TArray T2_;
  TArray T3_;

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
  double target_precision_;
  double computed_precision_ = std::numeric_limits<double>::max();
  bool verbose_;
  double CC3_corr_energy_;
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
    CC3_corr_energy_ = 0.0;
    f_pq_diagonal_.reset();
    LCAOWavefunction<Tile, Policy>::obsolete();
    ref_wfn_->obsolete();
  }

  // get T1 amplitudes
  TArray t1() const {
    if (T1_.is_initialized()) {
      return T1_;
    } else {
      throw std::runtime_error("CC3 T1 amplitudes have not been computed");
    }
  }
  // get T2 amplitudes
  TArray t2() const {
    if (T2_.is_initialized()) {
      return T2_;
    } else {
      throw std::runtime_error("CC3 T2 amplitudes have not been computed");
    }
  }

  // get T3 amplitudes
  TArray t3() const {
    if (T3_.is_initialized()) {
      return T3_;
    } else {
      throw std::runtime_error("CC3 T3 amplitudes have not been computed");
    }
  }

  bool is_df() const { return df_; }

  void set_t1(const TArray &t1) { T1_ = t1; }

  void set_t2(const TArray &t2) { T2_ = t2; }

  void set_t3(const TArray &t3) { T3_ = t3; }

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
      if (solver_str_ == "jacobi_diis") {
        const auto n_occ = this->trange1_engine()->get_occ();
        const auto n_frozen = this->trange1_engine()->get_nfrozen();
        const auto n_act_occ = n_occ - n_frozen;
        const auto n_uocc = f_pq_diagonal_->rows() - n_occ;
        solver_ = std::make_shared<cc::JacobiDIISSolver<TArray>>(
            kv_, f_pq_diagonal_->segment(n_frozen, n_act_occ),
            f_pq_diagonal_->segment(n_occ, n_uocc));
      } else if (solver_str_ == "pno")
        solver_ = std::make_shared<cc::PNOSolver<
            TArray, typename LCAOFactory<Tile, Policy>::DirectTArray>>(
            kv_, this->lcao_factory());
      else
        throw ProgrammingError("unknown solver string", __FILE__, __LINE__);

      // set the precision
      target_precision_ = energy->target_precision(0);

      TArray t1;
      TArray t2;
      TArray t3;

      if (method_ == "standard") {
        CC3_corr_energy_ = compute_CC3_conventional(t1, t2, t3);
      } /*else if (method_ == "df") {
        CC3_corr_energy_ = compute_cc3_df(t1, t2);
      } else if (method_ == "direct" || method_ == "direct_df") {
        // initialize direct integral class
        direct_ao_array_ =
            this->ao_factory().compute_direct(L"(μ ν| G|κ λ)[ab_ab]");
        CC3_corr_energy_ = compute_ccsd_direct(t1, t2);
      }*/

      T1_ = t1;
      T2_ = t2;
      T3_ = t3;

      this->computed_ = true;
      this->set_value(energy, ref_energy->energy() + CC3_corr_energy_);

      auto time1 = mpqc::fenced_now(world);
      auto duration0 = mpqc::duration_in_s(time0, time1);
      ExEnv::out0() << "CC3 Time in CC3: " << duration0 << " S" << std::endl;
    }
  }

 protected:
  // store all the integrals in memory
  // used as reference for development
  double compute_CC3_conventional(TArray &t1, TArray &t2, TArray &t3) {

  //VR ... Initialize
    auto &world = this->wfn_world()->world();
    bool accurate_time = this->lcao_factory().accurate_time();

    auto n_occ = this->trange1_engine()->get_occ();
    auto n_frozen = this->trange1_engine()->get_nfrozen();

    if (world.rank() == 0) {
      std::cout << "Use Conventional CC3 Compute" << std::endl;
    }

    auto tmp_time0 = mpqc::now(world, accurate_time);
    // get all two electron integrals
    TArray g_abcd = this->get_abcd();
    TArray g_ijab = this->get_ijab();
    TArray g_ijkl = this->get_ijkl();
    TArray g_iajb = this->get_iajb();
    TArray g_iabj = this->get_iabj();
    TArray g_iabc = this->get_iabc();
    TArray g_aibc = this->get_aibc();
    TArray g_ijak = this->get_ijak();
    TArray g_ijka = this->get_ijka();
    //some more integrals  _VR
    TArray g_abci = this->get_abci();
    TArray g_cjkl = this->get_aijk();
    TArray g_dabi = this->get_abci();
    TArray g_abij = this->get_abij();
    TArray g_jkbc = this->get_ijab();
    TArray g_aijk = this->get_aijk();
    TArray g_iajk = this->get_iajk() ;
    //_eVR

    this->lcao_factory().registry().purge_formula(L"(i ν| G |κ λ )");

    auto tmp_time1 = mpqc::now(world, accurate_time);
    auto tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
    mpqc::utility::print_par(world, "Integral Prepare Time: ", tmp_time, "\n");

    TArray f_ai = this->get_fock_ai();
    TArray f_ij = this->get_fock_ij();
    TArray f_ab = this->get_fock_ab();

    // store d1 to local
    TArray d1 = create_d_ai<Tile, Policy>(f_ai.world(), f_ai.trange(),
                                          *orbital_energy(), n_occ, n_frozen);

    t1("a,i") = f_ai("a,i") * d1("a,i");
    t1.truncate();

    {
      TArray g_abij;
      g_abij("a,b,i,j") = g_ijab("i,j,a,b");
      t2 = d_abij(g_abij, *orbital_energy(), n_occ, n_frozen);
    }


    // T3

    // Initialize T3 and set all elements to be zero
    //first assign it a size indirectly
    t3("a,b,c,i,j,k") = t2("a,b,i,j") * t1("c,k");

    TArray tau;
    tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");

    double E0 = 0.0;
    double E1 = 2.0 * TA::dot(f_ai("a,i"), t1("a,i")) +
                TA::dot(g_ijab("i,j,a,b"), 2 * tau("a,b,i,j") - tau("b,a,i,j"));
    double mp2 = E1;
    double dE = std::abs(E1 - E0);

    mpqc::utility::print_par(world, "MP2 Energy      ", mp2, "\n");

    // optimize t1 and t2
    std::size_t iter = 0ul;
    double error = 1.0;
    TArray r1;
    TArray r2;
    TArray r3;

    if (world.rank() == 0) {
      std::cout << "Start Iteration" << std::endl;
      std::cout << "Max Iteration: " << max_iter_ << std::endl;
      std::cout << "Target precision: " << target_precision_ << std::endl;
      std::cout << "AccurateTime: " << accurate_time << std::endl;
      std::cout << "PrintDetail: " << verbose_ << std::endl;
      std::cout << "Reduced ABCD Memory Approach: "
                << (reduced_abcd_memory_ ? "Yes" : "No") << std::endl;
    }

    while (iter < max_iter_) {
      // start timer
      auto time0 = mpqc::fenced_now(world);
      TArray::wait_for_lazy_cleanup(world);
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

      //Add T3 contribution to T1 _VR
      {

        TArray g_jkbc_AS;

        g_jkbc_AS("j,k,b,c") = (2.0 * g_jkbc("j,k,b,c") ) - g_jkbc("j,k,c,b");

       // seems like a sign error in eqn 15, Scuseria & Schaefer, CPL (1988)
       // Instead of '+', it should be '-'
        r1("a,i") -= g_jkbc_AS("j,k,b,c") * (t3("b,a,c,j,k,i") - t3("b,c,a,j,k,i"));

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

      // Add contribution of T3 into T2 _VR
      {
         TArray g_akcd=this->get_aibc();
         TArray g_akcd_AS;
         TArray g_klic=this->get_ijka();
         TArray g_klic_AS;
         TArray n2_abij;
         TArray o2_abij;
         TArray g_ijab_AS;

         g_akcd_AS("a,k,c,d") = (2 * g_akcd("a,k,c,d")) - g_akcd("a,k,d,c");

         n2_abij("a,b,i,j") = g_akcd_AS("a,k,c,d") * t3("c,b,d,i,j,k") -
                              g_akcd   ("a,k,c,d") * t3("c,d,b,i,j,k");

         g_klic_AS("k,l,i,c") = (2 * g_klic("k,l,i,c")) - g_klic("l,k,i,c");

         n2_abij("a,b,i,j") -= g_klic_AS("k,l,i,c") * t3("a,b,c,k,j,l") -
                               g_klic   ("k,l,i,c") * t3("a,c,b,k,j,l");

         // f_ai ..

         n2_abij("a,b,i,j") -= f_ai("c,k") * (t3("a,b,c,i,k,j")-t3("a,c,b,i,k,j"));

         // permute
         n2_abij("a,b,i,j") = n2_abij("a,b,i,j") + n2_abij("b,a,j,i");

         //add to the r2 residual
         r2("a,b,i,j") += n2_abij("a,b,i,j");

         //add contribution of T1T3 to the r2 residual
         g_ijab_AS("k,l,c,d")  = (2.0 * g_ijab("k,l,c,d")) - g_ijab("l,k,c,d") ;

         o2_abij("a,b,i,j")  = (g_ijab_AS("k,l,c,d") * (t3("a,b,c,i,j,k") - t3("a,c,b,i,j,k")))* t1("d,l") ;

         o2_abij("a,b,i,j") -= (g_ijab_AS("k,l,c,d") * t3("a,c,b,i,k,l") - g_ijab("k,l,c,d") * t3("c,a,b,i,k,l"))
                               * t1("d,j") ;


         o2_abij("a,b,i,j") -= (g_ijab_AS("k,l,c,d") * t3("a,d,c,i,j,k") - g_ijab("k,l,c,d") * t3("c,d,a,i,j,k"))
                               * t1("b,l") ;

         //permute
         o2_abij("a,b,i,j") = o2_abij("a,b,i,j") + o2_abij("b,a,j,i");

         //add to r2 residual
         r2("a,b,i,j") += o2_abij("a,b,i,j");

      }

      auto t2_time1 = mpqc::now(world, accurate_time);
      auto t2_time = mpqc::duration_in_s(t2_time0, t2_time1);
      if (verbose_) {
        mpqc::utility::print_par(world, "t2 total time: ", t2_time, "\n");
      }

      // compute r3 residual  _VR
      {
          // equations below are from spin-adapted CCSDT implementation
          // by Noga and Bartlett, JCP, 86, 7041 (1987)
          // by Scuseria and Schaefer, CPL 146, 23 (1988), Eq. 9 for T-1a
          // Extension from CCSDT-1b to CC3 requires contribution of T1T2, T1T1T2, and T1T1T1T2 terms to t3/r3 residual
          // For CC3, refer to papers below
          // by Christiansen, Koch, and Jørgensen, JCP, 103, 7429 (1995)
          // by Koch, Christiansen, Jørgensen, de Merás, and Helgaker, JCP, 106, 1808 (1997)
          // Let's construct intermediates as done by Noga and Bartlett

           TArray Chi_dabi;
           TArray Chi_cjkl;
           TArray Chi_im;
           TArray Chi_ae;
           TArray Chi_jkmn;
           TArray Chi_bcef;
           TArray Chi_amei;
           TArray Chi_amie;
           TArray f_aijk;
           TArray f_aibc;
           TArray f_aijb;
           TArray f_aibj;
           TArray f_iabj;
           TArray f_iajb;
           TArray g_aijb;
           TArray g_aibj;
           TArray g_ijab_AS;
           TArray f_int_me;

           // some two-electron integrals

           g_aijb("a,i,j,b")  = g_iabj("i,a,b,j") ;
           g_aibj("a,i,b,j")  = g_iajb("i,a,j,b") ;
           g_ijab_AS("i,j,a,b") = (2 * g_ijab("i,j,a,b")) - g_ijab("i,j,b,a");

           // first define the intermediates

           f_aijk("e,i,m,n")  = g_aijk("e,i,m,n") ;
           f_aijk("e,i,m,n") += ( g_ijab("m,n,e,f") * t1("f,i")) ;

           f_aibc("a,m,e,f")  = g_aibc("a,m,e,f") ;
           f_aibc("a,m,e,f") -= ( g_ijab("n,m,e,f") * t1("a,n")) ;

           f_aijb("a,m,i,e")  = g_aijb("a,m,i,e") ;
           f_aijb("a,m,i,e") += ( g_aibc("a,m,f,e") * t1("f,i")) ;

           f_aibj("a,m,e,i")  = g_aibj("a,m,e,i") ;
           f_aibj("a,m,e,i") += ( g_aibc("a,m,e,f") * t1("f,i")) ;

           f_iabj("i,e,a,m")  = g_iabj("i,e,a,m");
           f_iabj("i,e,a,m") -= ( g_iajk("i,e,n,m") * t1("a,n")) ;

           f_iajb("i,e,m,a")  = g_iajb("i,e,m,a")  ;
           f_iajb("i,e,m,a") -= ( g_iajk("i,e,m,n") * t1("a,n")) ;

           f_int_me("m,e")    = g_ijab_AS ("m,n,e,f") * t1("f,n") ;
           //For non-HF reference with non-zero occupied-virtual block, add this part
           //f_int_me("m,e")   += f_ia("m,e") ;

           // Now the super-intermediates,


           Chi_dabi("b,a,e,i") =  g_dabi("b,a,e,i") + ((f_aijk("e,i,m,n") * t1("a,n")) * t1("b,m"))
                                - f_aijb("a,m,i,e") * t1("b,m") - f_aibj("b,m,e,i") * t1("a,m")
                                + g_abcd("a,b,f,e") * t1("f,i") ;

           Chi_cjkl("a,m,i,j") =  g_cjkl("a,m,i,j") + ((f_aibc("a,m,e,f") * t1("e,i")) * t1("f,j"))
                                + f_iabj("i,e,a,m") * t1("e,j") + f_iajb("j,e,m,a") * t1("e,i")
                                - g_ijkl("i,j,n,m") * t1("a,n") ;


           // compute the residual r3
           TArray r3_1;
           TArray r3_2;


           r3_1("a,b,c,i,j,k")  =   f_ab("a,e") * t3("e,b,c,i,j,k")
                                  - f_ij("i,l") * t3("a,b,c,l,j,k");


           // permute [(i,a) with (j,b) and (i,a) with (k,c)]

           r3_1("a,b,c,i,j,k") = r3_1("a,b,c,i,j,k") + r3_1("b,a,c,j,i,k") + r3_1("c,b,a,k,j,i") ;


           //  the T1 and T2 dependent terms in the t3/r3 equation

           r3_2("a,b,c,i,j,k") =    Chi_dabi("b,a,e,i") * t2("c,e,k,j")
                                  - Chi_cjkl("a,m,i,j") * t2("b,c,m,k") ;

           //permute
           r3_2("a,b,c,i,j,k") = r3_2("a,b,c,i,j,k") + r3_2("a,c,b,i,k,j") +
                                 r3_2("c,a,b,k,i,j") + r3_2("c,b,a,k,j,i") +
                                 r3_2("b,c,a,j,k,i") + r3_2("b,a,c,j,i,k");


           // Add the contribution of T3 dependent terms

           r3("a,b,c,i,j,k")  = r3_1("a,b,c,i,j,k") ;
           r3("a,b,c,i,j,k") += r3_2("a,b,c,i,j,k") ;


      //denominator
      //t3 = d_abcijk(t3_temp, *orbital_energy(), n_occ, n_frozen);

      }
       // end of triples part  _VR

      // error = residual norm per element _edited VR
      error = std::sqrt((std::pow(norm2(r1), 2) + std::pow(norm2(r2), 2)
               + std::pow(norm2(r3), 2)))/
              (size(r1) + size(r2) + size(r3));

      // recompute energy
      E0 = E1;
      E1 = 2.0 * TA::dot(f_ai("a,i") + r1("a,i"), t1("a,i")) +
           TA::dot(g_ijab("i,j,a,b") + r2("a,b,i,j"),
                   2 * tau("a,b,i,j") - tau("b,a,i,j"));
      dE = std::abs(E0 - E1);

      // update the amplitudes, if not converged
      if (dE >= target_precision_ || error >= target_precision_) {
        tmp_time0 = mpqc::now(world, accurate_time);


        assert(solver_);
        solver_->update(t1, t2, t3, r1, r2, r3);

        if (verbose_) {
          mpqc::detail::print_size_info(r2, "R2");
          mpqc::detail::print_size_info(t2, "T2");
        }

        // recompute tau as well
        tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");
        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (verbose_) {
          mpqc::utility::print_par(world, "solver time: ", tmp_time, "\n");
        }

        auto time1 = mpqc::fenced_now(world);
        auto duration = mpqc::duration_in_s(time0, time1);

        if (world.rank() == 0) {
          detail::print_CC3(iter, dE, error, E1, duration);
        }

        iter += 1ul;

      } else {
        auto time1 = mpqc::fenced_now(world);
        auto duration = mpqc::duration_in_s(time0, time1);

        if (world.rank() == 0) {
          detail::print_CC3(iter, dE, error, E1, duration);
        }

        break;
      }
    }
    if (iter >= max_iter_) {
      utility::print_par(this->wfn_world()->world(),
                         "\n Warning!! Exceed Max Iteration! \n");
    }
    if (world.rank() == 0) {
      std::cout << "CC3 Energy  " << E1 << std::endl;
    }
    return E1;
  }

 private:
  /*
  double compute_cc3_df(TArray &t1, TArray &t2) {
    auto &world = this->wfn_world()->world();
    bool accurate_time = this->lcao_factory().accurate_time();

    auto n_occ = this->trange1_engine()->get_occ();
    auto n_frozen = this->trange1_engine()->get_nfrozen();

    if (world.rank() == 0) {
      std::cout << "Use DF CCSD Compute" << std::endl;
    }

    auto tmp_time0 = mpqc::now(world, accurate_time);
    // get all two electron integrals
    TArray g_ijab = this->get_ijab();
    TArray g_ijkl = this->get_ijkl();
    TArray g_abcd = this->get_abcd();

    TArray X_ai = this->get_Xai();
    TArray g_iajb = this->get_iajb();
    TArray g_iabc = this->get_iabc();
    TArray g_ijak = this->get_ijak();
    TArray g_ijka = this->get_ijka();
    auto tmp_time1 = mpqc::now(world, accurate_time);
    auto tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
    mpqc::utility::print_par(world, "Integral Prepare Time: ", tmp_time, "\n");

    TArray f_ai = this->get_fock_ai();
    TArray f_ij = this->get_fock_ij();
    TArray f_ab = this->get_fock_ab();

    // store d1 to local
    TArray d1 = create_d_ai<Tile, Policy>(f_ai.world(), f_ai.trange(),
                                          *orbital_energy(), n_occ, n_frozen);

    t1("a,i") = f_ai("a,i") * d1("a,i");
    t1.truncate();

    {
      TArray g_abij;
      g_abij("a,b,i,j") = g_ijab("i,j,a,b");
      t2 = d_abij(g_abij, *orbital_energy(), n_occ, n_frozen);
    }

    TArray tau;
    tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");

    double E0 = 0.0;
    double E1 = 2.0 * TA::dot(f_ai("a,i"), t1("a,i")) +
                TA::dot(g_ijab("i,j,a,b"), 2 * tau("a,b,i,j") - tau("b,a,i,j"));
    double mp2 = E1;
    double dE = std::abs(E1 - E0);

    mpqc::utility::print_par(world, "MP2 Energy      ", mp2, "\n");

    // optimize t1 and t2
    std::size_t iter = 0ul;
    double error = 1.0;
    TArray r1;
    TArray r2;

    if (world.rank() == 0) {
      std::cout << "Start Iteration" << std::endl;
      std::cout << "Max Iteration: " << max_iter_ << std::endl;
      std::cout << "Target Precision: " << target_precision_ << std::endl;
      std::cout << "AccurateTime: " << accurate_time << std::endl;
      std::cout << "PrintDetail: " << verbose_ << std::endl;
      std::cout << "Reduced ABCD Memory Approach: "
                << (reduced_abcd_memory_ ? "Yes" : "No") << std::endl;
    }

    while (iter < max_iter_) {
      // start timer
      auto time0 = mpqc::fenced_now(world);
      TArray::wait_for_lazy_cleanup(world);

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
        r1("a,i") = f_ai("a,i") - 2.0 * f_ai("c,k") * t1("c,i") * t1("a,k");

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
            (g_iabc("i,c,a,b") - g_iajb("k,b,i,c") * t1("a,k")) * t1("c,j") -
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
                      (2.0 * g_iabc("k,a,d,c") - g_iabc("k,a,c,d")) * t1("d,k");

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

              + g_iabc("k,a,c,d") * t1("d,i")

              - X_ai("x,d,l") * T("d,a,i,l") * X_ai("x,c,k")

              +
              0.5 * (2.0 * g_ijab("k,l,c,d") - g_ijab("k,l,d,c")) *
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
        b_abij("a,b,i,j") = tau("c,d,i,j") * g_abcd("a,b,c,d");

        b_abij("a,b,i,j") -= g_iabc("k,a,d,c") * tau("c,d,i,j") * t1("b,k");

        b_abij("a,b,i,j") -= g_iabc("k,b,c,d") * tau("c,d,i,j") * t1("a,k");

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

      // error = residual norm per element
      error = std::sqrt((std::pow(norm2(r1), 2) + std::pow(norm2(r2), 2))) /
              (size(r1) + size(r2));

      // recompute energy
      E0 = E1;
      E1 = 2.0 * TA::dot(f_ai("a,i") + r1("a,i"), t1("a,i")) +
           TA::dot(g_ijab("i,j,a,b") + r2("a,b,i,j"),
                   2 * tau("a,b,i,j") - tau("b,a,i,j"));
      dE = std::abs(E0 - E1);

      if (dE >= target_precision_ || error >= target_precision_) {
        tmp_time0 = mpqc::now(world, accurate_time);

        assert(solver_);
        solver_->update(t1, t2, r1, r2);

        if (verbose_) {
          mpqc::detail::print_size_info(r2, "R2");
          mpqc::detail::print_size_info(t2, "T2");
        }

        // recompute tau
        tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");
        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (verbose_) {
          mpqc::utility::print_par(world, "diis time: ", tmp_time, "\n");
        }

        auto time1 = mpqc::fenced_now(world);
        auto duration = mpqc::duration_in_s(time0, time1);

        if (world.rank() == 0) {
          detail::print_CC3(iter, dE, error, E1, duration);
        }

        iter += 1ul;
      } else {
        auto time1 = mpqc::fenced_now(world);
        auto duration = mpqc::duration_in_s(time0, time1);

        if (world.rank() == 0) {
          detail::print_CC3(iter, dE, error, E1, duration);
        }

        break;
      }
    }
    if (iter >= max_iter_) {
      utility::print_par(this->wfn_world()->world(),
                         "\n Warning!! Exceed Max Iteration! \n");
    }
    if (world.rank() == 0) {
      std::cout << "CCSD Energy  " << E1 << std::endl;
    }
    return E1;
  }
    */

  /*
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

    auto n_occ = this->trange1_engine()->get_occ();
    auto n_frozen = this->trange1_engine()->get_nfrozen();

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

    TArray d1 = create_d_ai<Tile, Policy>(f_ai.world(), f_ai.trange(),
                                          *orbital_energy(), n_occ, n_frozen);

    t1("a,i") = f_ai("a,i") * d1("a,i");
    t1.truncate();

    {
      TArray g_abij;
      g_abij("a,b,i,j") = g_ijab("i,j,a,b");
      t2 = d_abij(g_abij, *orbital_energy(), n_occ, n_frozen);
    }

    auto tmp_time1 = mpqc::now(world, accurate_time);
    auto tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
    mpqc::utility::print_par(world, "Integral Prepare Time: ", tmp_time, "\n");

    TArray tau;
    tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");

    double E0 = 0.0;
    double E1 = 2.0 * TA::dot(f_ai("a,i"), t1("a,i")) +
                TA::dot(g_ijab("i,j,a,b"), 2 * tau("a,b,i,j") - tau("b,a,i,j"));
    double dE = std::abs(E1 - E0);
    double mp2 = E1;

    mpqc::utility::print_par(world, "MP2 Energy      ", mp2, "\n");

    // optimize t1 and t2
    std::size_t iter = 0ul;
    double error = 1.0;
    TArray r1;
    TArray r2;

    if (world.rank() == 0) {
      std::cout << "Start Iteration" << std::endl;
      std::cout << "Max Iteration: " << max_iter_ << std::endl;
      std::cout << "Target Precision: " << target_precision_ << std::endl;
      std::cout << "AccurateTime: " << accurate_time << std::endl;
      std::cout << "PrintDetail: " << verbose_ << std::endl;
    };

    while (iter < max_iter_) {
      // start timer
      auto time0 = mpqc::fenced_now(world);
      TArray::wait_for_lazy_cleanup(world);

      TArray u2_u11;
      // compute half transformed intermediates
      auto tu0 = mpqc::now(world, accurate_time);
      { u2_u11 = this->compute_u2_u11(t2, t1); }
      auto tu1 = mpqc::now(world, accurate_time);
      auto duration_u = mpqc::duration_in_s(tu0, tu1);

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

          r2("a,b,i,j") += -0.5 * k_kaic("k,a,i,c") * t2("b,c,k,j") -
                           k_kaic("k,b,i,c") * t2("a,c,k,j");
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

      // error = residual norm per element
      error = std::sqrt((std::pow(norm2(r1), 2) + std::pow(norm2(r2), 2))) /
              (size(r1) + size(r2));

      // recompute energy
      E0 = E1;
      E1 = 2.0 * TA::dot(f_ai("a,i") + r1("a,i"), t1("a,i")) +
           TA::dot(g_ijab("i,j,a,b") + r2("a,b,i,j"),
                   2 * tau("a,b,i,j") - tau("b,a,i,j"));
      dE = std::abs(E0 - E1);

      if (dE >= target_precision_ || error >= target_precision_) {
        tmp_time0 = mpqc::now(world, accurate_time);

        assert(solver_);
        solver_->update(t1, t2, r1, r2);

        if (verbose_) {
          mpqc::detail::print_size_info(r2, "R2");
          mpqc::detail::print_size_info(t2, "T2");
        }

        // recompute tau
        tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");
        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (verbose_) {
          mpqc::utility::print_par(world, "diis time: ", tmp_time, "\n");
        }

        auto time1 = mpqc::fenced_now(world);
        auto duration_t = mpqc::duration_in_s(time0, time1);

        if (world.rank() == 0) {
          detail::print_CC3_direct(iter, dE, error, E1, duration_u,
                                    duration_t);
        }

        iter += 1ul;
      } else {
        auto time1 = mpqc::fenced_now(world);
        auto duration_t = mpqc::duration_in_s(time0, time1);

        if (world.rank() == 0) {
          detail::print_CC3_direct(iter, dE, error, E1, duration_u,
                                    duration_t);
        }

        break;
      }
    }
    if (iter >= max_iter_) {
      utility::print_par(this->wfn_world()->world(),
                         "\n Warning!! Exceed Max Iteration! \n");
    }
    if (world.rank() == 0) {
      std::cout << "CCSD Energy     " << E1 << std::endl;
    }
    return E1;
  }
  */
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
    TArray result;
    TArray sqrt = this->ao_factory().compute(L"(Κ|G| Λ)[inv_sqr]");
    TArray three_center = this->lcao_factory().compute(L"(Κ|G|a b)");
    result("K,a,b") = sqrt("K,Q") * three_center("Q,a,b");
    return result;
  }

  /// get three center integral (X|ij)
  const TArray get_Xij() {
    TArray result;
    TArray sqrt = this->ao_factory().compute(L"(Κ|G| Λ)[inv_sqr]");
    TArray three_center = this->lcao_factory().compute(L"(Κ|G|i j)");
    result("K,i,j") = sqrt("K,Q") * three_center("Q,i,j");
    return result;
  }

  /// get three center integral (X|ai)
  const TArray get_Xai() {
    TArray result;
    TArray sqrt = this->ao_factory().compute(L"(Κ|G| Λ)[inv_sqr]");
    TArray three_center = this->lcao_factory().compute(L"(Κ|G|a i)");
    result("K,a,i") = sqrt("K,Q") * three_center("Q,a,i");
    return result;
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

  /// <ia|bj>
  virtual const TArray get_iabj() {
    std::wstring postfix = df_ ? L"[df]" : L"";
    return this->lcao_factory().compute(L"<i a|G|b j>" + postfix);
  }

  /// <ij|ka>
  virtual const TArray get_ijka() {
    std::wstring postfix = df_ ? L"[df]" : L"";
    return this->lcao_factory().compute(L"<i j|G|k a>" + postfix);
  }

  /// <ai|jk>
  const TArray get_aijk() {
    std::wstring post_fix = this->is_df() ? L"[df]" : L"";
    TArray result;
    result = this->lcao_factory().compute(L"<a i|G|j k>" + post_fix);

    //result("a,i,j,k") = result("i,j,k,a");
    return result;
  }
  /// <ia|jk>
  const TArray get_iajk() {
    std::wstring post_fix = this->is_df() ? L"[df]" : L"";
    TArray result;
    result = this->lcao_factory().compute(L"<i a|G|j k>" + post_fix);

    //result("a,i,j,k") = result("i,j,k,a");
    return result;
  }


  /// <ab|ci>
  const TArray get_abci() {
    std::wstring post_fix = this->is_df() ? L"[df]" : L"";
    TArray result;
     result = this->lcao_factory().compute(L"<i a|G|b c>" + post_fix);

    result("a,b,c,i") = result("i,a,b,c");
    return result;
  }


  /// <ab|ij>
  virtual const TArray get_abij() {
    std::wstring postfix = df_ ? L"[df]" : L"";
    return this->lcao_factory().compute(L"<a b|G|i j>" + postfix);
  }
  /// <ab|di>
  virtual const TArray get_abdi() {
    std::wstring postfix = df_ ? L"[df]" : L"";
    return this->lcao_factory().compute(L"<a b|G|d i>" + postfix);
  }
  /// <am|ij>
  virtual const TArray get_amij() {
    std::wstring postfix = df_ ? L"[df]" : L"";
    return this->lcao_factory().compute(L"<a m|G|i j>" + postfix);
  }

  /// <a|f|i>
  virtual const TArray get_fock_ai() {
    std::wstring postfix = df_ ? L"[df]" : L"";
    return this->lcao_factory().compute(L"<a|F|i>" + postfix);
  }

  /// <i|f|j>
  const TArray get_fock_ij() {
    std::wstring postfix = df_ ? L"[df]" : L"";
    return this->lcao_factory().compute(L"<i|F|j>" + postfix);
  }

  /// <a|f|b>
  const TArray get_fock_ab() {
    std::wstring postfix = df_ ? L"[df]" : L"";
    return this->lcao_factory().compute(L"<a|F|b>" + postfix);
  }

  /// <p|f|q>
  const TArray get_fock_pq() {
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
      u2_u11("p, r, i, j") =
          0.5 * (u2_u11("p, r, i, j") + u2_u11("r, p, j, i"));
      return u2_u11;
    } else {
      throw ProgrammingError(
          "CCSD: integral-direct implementation used, but direct integral not "
          "initialized",
          __FILE__, __LINE__);
    }
  }
};  // class CC3

#if TA_DEFAULT_POLICY == 0
extern template class CC3<TA::TensorD, TA::DensePolicy>;
#elif TA_DEFAULT_POLICY == 1
extern template class CC3<TA::TensorD, TA::SparsePolicy>;
#endif
}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_CC3_H_
