//
// Created by Chong Peng on 7/1/15.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_CCSD_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_CCSD_H_

#include <tiledarray.h>

#include "mpqc/chemistry/qc/cc/tpack.h"
#include "mpqc/chemistry/qc/cc/solvers.h"
#include "mpqc/chemistry/qc/lcao/cc/ccsd_r1_r2.h"
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
    std::printf("%3s \t %10s \t %10s \t %15s \t %10s \n", "iter", "deltaE",
                "residual", "energy", "total time/s");
  }
  std::printf("%3i \t %10.5e \t %10.5e \t %15.12f \t %10.1f \n", iter, dE,
              error, E1, time_total);
}

inline void print_ccsd_direct(int iter, double dE, double error, double E1,
                              double time_u, double time_total) {
  if (iter == 1) {
    std::printf("%3s \t %10s \t %10s \t %15s \t %10s \t %10s \n", "iter",
                "deltaE", "residual", "energy", "u time/s", "total time/s");
  }
  std::printf("%3i \t %10.5e \t %10.5e \t %15.12f \t %10.1f \t %10.1f \n", iter,
              dE, error, E1, time_u, time_total);
}
}  // namespace detail

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
   * | @c method | string | @c df if @c df_basis is provided, @c direct otherwise | method to compute the CCSD residual; valid choices are: @c standard (uses 4-index MO integrals throughout), @c direct (uses 4-index MO integrals with up to 3 unoccupied indices, and 4-center AO integrals), @c df (approximates 4-index MO integrals using density fitting), @c direct_df (hybrid between @c df and @c direct that avoids storing MO integrals with 3 unoccupied indices by using DF, see DOI 10.1021/acs.jpca.6b10150 for details) |
   * | @c max_iter | int | @c 30 | maxmium iteration in CCSD |
   * | @c verbose | bool | false | if print more information in CCSD iteration |
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
    if (solver_str_ != "jacobi_diis" && solver_str_ != "pno" &&
        solver_str_ != "svo")
      throw InputError("invalid value for solver keyword", __FILE__, __LINE__,
                       "solver");

    reduced_abcd_memory_ = kv.value<bool>("reduced_abcd_memory", true);

    max_iter_ = kv.value<int>("max_iter", 30);
    verbose_ = kv.value<bool>("verbose", false);
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
  double target_precision_;
  double computed_precision_ = std::numeric_limits<double>::max();
  bool verbose_;
  double ccsd_corr_energy_;
  // diagonal elements of the Fock matrix (not necessarily the eigenvalues)
  std::shared_ptr<const EigenVector<typename Tile::numeric_type>>
      f_pq_diagonal_;

 protected:
  std::shared_ptr<const EigenVector<typename Tile::numeric_type>>
  orbital_energy() {
    return f_pq_diagonal_;
  }

  void set_orbital_energy(
      const EigenVector<typename Tile::numeric_type> &orbital_energy) {
    f_pq_diagonal_ = std::make_shared<EigenVector<typename Tile::numeric_type>>(
        orbital_energy);
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

  bool verbose() const { return verbose_; }

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
      else if (solver_str_ == "svo")
        solver_ = std::make_shared<cc::SVOSolver<
            TArray, typename LCAOFactory<Tile, Policy>::DirectTArray>>(
            kv_, this->lcao_factory());
      else
        throw ProgrammingError("unknown solver string", __FILE__, __LINE__);

      // set the precision
      target_precision_ = energy->target_precision(0);

      TArray t1;
      TArray t2;

      if (method_ == "direct" || method_ == "direct_df") {
        direct_ao_array_ = this->ao_factory().compute_direct(L"(μ ν| G|κ λ)");
      }

      ccsd_corr_energy_ = compute_ccsd(t1, t2);

      T1_ = t1;
      T2_ = t2;

      // delete current solver
      solver_.reset();

      this->computed_ = true;
      this->set_value(energy, ref_energy->energy() + ccsd_corr_energy_);

      auto time1 = mpqc::fenced_now(world);
      auto duration0 = mpqc::duration_in_s(time0, time1);
      ExEnv::out0() << "CCSD Time in CCSD: " << duration0 << " S" << std::endl;
    }
  }

 protected:
  double compute_ccsd(TArray &t1, TArray &t2) {
    auto &world = this->wfn_world()->world();
    bool accurate_time = this->lcao_factory().accurate_time();

    if (method_ == "standard") {
      mpqc::utility::print_par(world, "Use Conventional CCSD Compute \n");
    } else if (method_ == "df") {
      mpqc::utility::print_par(world, "Use DF CCSD Compute \n");
    } else if (method_ == "direct") {
      mpqc::utility::print_par(world, "Use Direct CCSD Compute \n");
    }

    auto tmp_time0 = mpqc::now(world, accurate_time);
    cc::Integrals<TArray> ints;
    // fock matrix
    ints.Fia = this->get_fock_ia();
    ints.Fij = this->get_fock_ij();
    ints.Fab = this->get_fock_ab();
    // get all two electron integrals
    ints.Gijab = this->get_ijab();
    ints.Gijkl = this->get_ijkl();
    ints.Giajb = this->get_iajb();
    ints.Gijka = this->get_ijka();

    if (method_ == "standard" || (method_ == "df" && !reduced_abcd_memory_)) {
      ints.Gabcd = this->get_abcd();
      ints.Giabc = this->get_iabc();
    } else if (method_ == "direct") {
      ints.Giabc = this->get_iabc();
    }

    if (df_) {
      ints.Xai = this->get_Xai();
      ints.Xij = this->get_Xij();
      ints.Xab = this->get_Xab();
    }

    if (method_ == "direct" || method_ == "direct_df") {
      ints.Ci = this->lcao_factory()
                    .orbital_registry()
                    .retrieve(OrbitalIndex(L"i"))
                    .coefs();
      ints.Ca = this->lcao_factory()
                    .orbital_registry()
                    .retrieve(OrbitalIndex(L"a"))
                    .coefs();
    }

    this->lcao_factory().registry().purge_formula(L"(i ν| G |κ λ )");
    this->lcao_factory().registry().purge_formula(L"(a ν| G |κ λ )");

    auto tmp_time1 = mpqc::now(world, accurate_time);
    auto tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
    mpqc::utility::print_par(world, "Integral Prepare Time: ", tmp_time, "\n");

    TArray f_ai;
    f_ai("a,i") = ints.Fia("i,a");
    // initial guess = 0
    t1 = TArray(f_ai.world(), f_ai.trange(), f_ai.shape(), f_ai.pmap());
    t1.fill(0.0);
    TArray g_abij;
    g_abij("a,b,i,j") = ints.Gijab("i,j,a,b");
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
      std::cout << "Verbose: " << verbose_ << std::endl;
      std::cout << "Reduced ABCD Memory Approach: "
                << (reduced_abcd_memory_ ? "Yes" : "No") << std::endl;
    }

    // CCSD solver loop
    mpqc::time_point time0;
    double duration_u;
    while (iter < max_iter_) {
      TArray::wait_for_lazy_cleanup(world);

      // zero out singles if want CCD
      // TArray r1_new(r1.world(), r1.trange(), r1.shape()); r1_new.fill(0.0);
      // r1("a,i") = r1_new("a,i");
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
        solver_->update(t1, t2, r1, r2);

        // recompute tau
        tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");
        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (verbose_) {
          mpqc::utility::print_par(world, "Solver::update time: ", tmp_time,
                                   "\n");
        }

        // log the iteration
        auto time1 = mpqc::fenced_now(world);
        if (world.rank() == 0 && iter > 0) {
          auto duration = mpqc::duration_in_s(time0, time1);
          if(method_ == "direct" || method_ == "direct_df"){
            detail::print_ccsd_direct(iter, dE, error, E1, duration_u, duration);
          }
          else{
            detail::print_ccsd(iter, dE, error, E1, duration);
          }
        }
      } else {  // break out of the solver loop, if converged
        // log the iteration
        auto time1 = mpqc::fenced_now(world);
        MPQC_ASSERT(iter > 0);
        if (world.rank() == 0) {
          auto duration = mpqc::duration_in_s(time0, time1);
          detail::print_ccsd(iter, dE, error, E1, duration);
        }
        break;
      }

      // start iteration timer
      time0 = mpqc::fenced_now(world);

      TArray U;
      auto tu0 = mpqc::now(world, accurate_time);
      if (method_ == "direct" || method_ == "direct_df") {
        U = compute_u2_u11(t2, t1);
      }
      auto tu1 = mpqc::now(world, accurate_time);
      duration_u = mpqc::duration_in_s(tu0, tu1);

      auto t1_time0 = mpqc::now(world, accurate_time);

      if (method_ == "standard") {
        r1 = cc::compute_cs_ccsd_r1(t1, t2, tau, ints);
      } else if (method_ == "df") {
        r1 = cc::compute_cs_ccsd_r1_df(t1, t2, tau, ints);
      } else if (method_ == "direct" || method_ == "direct_df") {
        r1 = cc::compute_cs_ccsd_r1(t1, t2, tau, ints, U);
      }

      auto t1_time1 = mpqc::now(world, accurate_time);
      auto t1_time = mpqc::duration_in_s(t1_time0, t1_time1);
      if (verbose_) {
        mpqc::utility::print_par(world, "t1 total time: ", t1_time, "\n");
      }

      // intermediates for t2
      // external index i j a b

      auto t2_time0 = mpqc::now(world, accurate_time);

      if (method_ == "standard") {
        r2 = cc::compute_cs_ccsd_r2(t1, t2, tau, ints);
      } else if (method_ == "df") {
        r2 = cc::compute_cs_ccsd_r2_df(t1, t2, tau, ints);
      } else if (method_ == "direct") {
        r2 = cc::compute_cs_ccsd_r2(t1, t2, tau, ints, U);
      } else {
        r2 = cc::compute_cs_ccsd_r2_df(t1, t2, tau, ints, U);
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

 protected:
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

  /// <i|f|a>
  virtual const TArray get_fock_ia() {
    std::wstring postfix = df_ ? L"[df]" : L"";
    return this->lcao_factory().compute(L"<i|F|a>" + postfix);
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
      TArray Ca = this->lcao_factory()
                      .orbital_registry()
                      .retrieve(OrbitalIndex(L"a"))
                      .coefs();
      ;
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
