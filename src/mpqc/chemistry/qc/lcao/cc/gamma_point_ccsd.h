#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_GAMMA_POINT_CCSD_H
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_GAMMA_POINT_CCSD_H

#include "mpqc/chemistry/qc/cc/diis.h"
#include "mpqc/chemistry/qc/lcao/mbpt/gamma_point_mp2.h"
#include "mpqc/chemistry/qc/lcao/scf/zrhf.h"

namespace mpqc {
namespace lcao {

namespace detail {

inline void print_gamma_point_ccsd(int iter, double dE, double error, double E1,
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

}  // namespace detail

template <typename Tile, typename Policy>
class GammaPointCCSD : public PeriodicLCAOWavefunction<Tile, Policy>,
                       public Provides<Energy> {
 public:
  using TArray = TA::DistArray<Tile, Policy>;

  GammaPointCCSD() = default;

  GammaPointCCSD(const KeyVal &kv)
      : PeriodicLCAOWavefunction<Tile, Policy>(kv), kv_(kv) {
    if (kv.exists("ref")) {
      ref_wfn_ =
          kv.keyval("ref").class_ptr<PeriodicAOWavefunction<Tile, Policy>>();
    } else {
      throw std::invalid_argument(
          "Default ref wfn in GammaPointCCSD is not supported!");
    }

    max_iter_ = kv.value<int64_t>("max_iter", 30);
    converge_ = kv.value<double>("converge", 1.0E-7);
    print_detail_ = kv.value<bool>("print_detail", false);
  }

  ~GammaPointCCSD() {}

  void obsolete() override {
    gp_ccsd_corr_energy_ = 0.0;
    PeriodicLCAOWavefunction<Tile, Policy>::obsolete();
    ref_wfn_->obsolete();
  }

 private:
  const KeyVal kv_;
  std::shared_ptr<PeriodicAOWavefunction<Tile, Policy>> ref_wfn_;
  int64_t max_iter_;
  double converge_;
  bool print_detail_;
  double gp_ccsd_corr_energy_;
  Matrixz C_;
  Vectorz eps_;

 private:
  /*!
   * \brief This initializes gamma-point MP2
   */
  void init() {
    auto nk = ref_wfn_->nk();
    auto k_size = 1 + detail::k_ord_idx(nk(0) - 1, nk(1) - 1, nk(2) - 1, nk);
    auto unitcell = this->lcao_factory().pao_factory().molecule();

    int64_t gamma_point;
    if (k_size % 2 == 1)
      gamma_point = (k_size - 1) / 2;
    else
      throw std::invalid_argument(
          "# of k points must be odd in order to run gamma-point methods");

    C_ = ref_wfn_->co_coeff()[gamma_point];
    eps_ = ref_wfn_->co_energy()[gamma_point];

    mo_insert_gamma_point(this->lcao_factory(), C_, unitcell, this->occ_block(),
                          this->unocc_block());
  }

  /*!
   * \brief This computes gamma-point CCSD correlation energy
   * \param t1 T1 amplitudes
   * \param t2 T2 amplitudes
   * \return gamma_point CCSD correlation energy
   */
  double compute_gamma_point_ccsd(TArray &t1, TArray &t2) {
    ExEnv::out0() << "Computing conventional gamma-point CCSD ..." << std::endl;

    auto &world = this->wfn_world()->world();
    auto accurate_time = this->lcao_factory().accurate_time();
    const auto charge = 0;
    const auto nelectrons =
        this->lcao_factory().pao_factory().unitcell().total_atomic_number() -
        charge;
    auto n_occ = nelectrons / 2;
    std::size_t n_frozen = 0;

    auto tmp_time0 = mpqc::now(world, accurate_time);
    // get all two electron integrals
    auto g_abij = this->get_abij();
    auto g_ijkl = this->get_ijkl();
    auto g_abcd = this->get_abcd();
    auto g_iajb = this->get_iajb();
    auto g_iabc = this->get_iabc();
    auto g_aibc = this->get_aibc();
    auto g_ijak = this->get_ijak();
    auto g_ijka = this->get_ijka();
    auto tmp_time1 = mpqc::now(world, accurate_time);
    auto tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
    if (print_detail_)
      ExEnv::out0() << "Integral Prepare Time: " << tmp_time << std::endl;

    auto f_ai = this->get_fock_ai();

    // store d1 to local
    auto d1 = create_d_ai<Tile, Policy>(f_ai.world(), f_ai.trange(), eps_,
                                         n_occ, n_frozen);

    // initial guess of t1 and t2
    t1("a, i") = f_ai("a, i") * d1("a, i");
    t1.truncate();
    t2 = d_abij(g_abij, eps_, n_occ, n_frozen);

    TArray tau;
    tau("a, b, i, j") = t2("a, b, i, j") + t1("a, i") * t1("b, j");

    // gamma-point MP2 energy
    double E0 = 0.0;
    std::complex<double> E1_tmp1 = TA::dot(2.0 * f_ai("a,i"), t1("a,i"));
    std::complex<double> E1_tmp2 =
        TA::dot(2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j"), tau("a,b,i,j"));
    double E1 = (E1_tmp1 + E1_tmp2).real();

    double mp2 = E1;
    double dE = std::abs(E1 - E0);
    ExEnv::out0() << "Gamma-Point MP2 Energy = " << mp2 << std::endl;

    // optimize t1 and t2
    std::size_t iter = 0ul;
    double error = 1.0;
    TArray r1;
    TArray r2;

    bool less = kv_.value<bool>("less_memory", false);

    ExEnv::out0() << "Start Iterations" << std::endl;
    ExEnv::out0() << "Max Iteration: " << max_iter_ << std::endl;
    ExEnv::out0() << "Convergence: " << converge_ << std::endl;
    ExEnv::out0() << "PrintDetail: " << print_detail_ << std::endl;
    if (less)
      ExEnv::out0() << "Less Memory Approach: Yes" << std::endl;
    else
      ExEnv::out0() << "Less Memory Approach: No" << std::endl;

    auto diis = this->get_diis(world);

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

        // compute residual r1(n) = t1(n+1) - t1(n)
        // external index i and a
        tmp_time0 = mpqc::now(world, accurate_time);
        r1("a,i") = f_ai("a,i") - 2.0 * f_ai("c,k") * t1("c,i") * t1("a,k");

        {
          h_ac("a,c") =
              -(2.0 * g_abij("c,d,k,l") - g_abij("c,d,l,k")) * tau("a,d,k,l");
          r1("a,i") += h_ac("a,c") * t1("c,i");
        }

        {
          h_ki("k,i") =
              (2.0 * g_abij("c,d,k,l") - g_abij("d,c,k,l")) * tau("c,d,i,l");
          r1("a,i") -= t1("a,k") * h_ki("k,i");
        }

        {
          h_kc("k,c") =
              f_ai("c,k") +
              (2.0 * g_abij("c,d,k,l") - g_abij("d,c,k,l")) * t1("d,l");
          r1("a,i") += h_kc("k,c") * (2.0 * t2("c,a,k,i") - t2("c,a,i,k") +
                                      t1("c,i") * t1("a,k"));
        }

        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (print_detail_) {
          mpqc::utility::print_par(world, "t1 h term time: ", tmp_time, "\n");
        }

        tmp_time0 = mpqc::now(world, accurate_time);
        r1("a,i") += (2.0 * g_abij("c,a,k,i") - g_iajb("k,a,i,c")) * t1("c,k");

        r1("a,i") +=
            (2.0 * g_iabc("k,a,c,d") - g_iabc("k,a,d,c")) * tau("c,d,k,i");

        r1("a,i") -=
            (2.0 * g_ijak("k,l,c,i") - g_ijak("l,k,c,i")) * tau("c,a,k,l");

        r1("a,i") *= d1("a,i");

        r1("a,i") -= t1("a,i");

        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (print_detail_) {
          mpqc::utility::print_par(world, "t1 other time: ", tmp_time, "\n");
        }
      }
      auto t1_time1 = mpqc::now(world, accurate_time);
      auto t1_time = mpqc::duration_in_s(t1_time0, t1_time1);
      if (print_detail_) {
        mpqc::utility::print_par(world, "t1 total time: ", t1_time, "\n");
      }

      // intermediates for t2
      // external index i j a b

      auto t2_time0 = mpqc::now(world, accurate_time);

      // compute residual r2(n) = t2(n+1) - t2(n)

      // permutation part
      tmp_time0 = mpqc::now(world, accurate_time);

      {
        TArray tmp_abic;
        tmp_abic("a,b,i,c") = g_iabc("i,c,a,b") - g_iajb("k,b,i,c") * t1("a,k");
        r2("a,b,i,j") = tmp_abic("a,b,i,c") * t1("c,j");

        TArray tmp_akij;
        tmp_akij("a,k,i,j") = g_ijak("i,j,a,k") + g_abij("a,c,i,k") * t1("c,j");
        r2("a,b,i,j") -= tmp_akij("a,k,i,j") * t1("b,k");

        // TODO: fix it later. Could be a linker problem.
        //        r2("a,b,i,j") =
        //            (g_iabc("i,c,a,b") - g_iajb("k,b,i,c") * t1("a,k")) *
        //            t1("c,j");

        //        r2("a,b,i,j") -=
        //            (g_ijak("i,j,a,k") + g_abij("a,c,i,k") * t1("c,j")) *
        //            t1("b,k");
      }
      tmp_time1 = mpqc::now(world, accurate_time);
      tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
      if (print_detail_) {
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
      if (print_detail_) {
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

          j_akic("a,k,i,c") = g_abij("a,c,i,k");

          j_akic("a,k,i,c") -= g_ijka("l,k,i,c") * t1("a,l");

          j_akic("a,k,i,c") += g_aibc("a,k,d,c") * t1("d,i");

          j_akic("a,k,i,c") -= g_abij("c,d,k,l") * T("d,a,i,l");

          j_akic("a,k,i,c") += 0.5 *
                               (2.0 * g_abij("c,d,k,l") - g_abij("d,c,k,l")) *
                               t2("a,d,i,l");

          k_kaic("k,a,i,c") = g_iajb("k,a,i,c")

                              - g_ijka("k,l,i,c") * t1("a,l")

                              + g_iabc("k,a,d,c") * t1("d,i")

                              - g_abij("d,c,k,l") * T("d,a,i,l");
          if (print_detail_) {
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
      if (print_detail_) {
        mpqc::utility::print_par(world, "t2 j,k term time: ", tmp_time, "\n");
      }

      // perform the permutation
      r2("a,b,i,j") = r2("a,b,i,j") + r2("b,a,j,i");

      r2("a,b,i,j") += g_abij("a,b,i,j");

      tmp_time0 = mpqc::now(world, accurate_time);
      {
        TArray a_klij;
        // compute a intermediate
        a_klij("k,l,i,j") = g_ijkl("k,l,i,j");

        a_klij("k,l,i,j") += g_ijka("k,l,i,c") * t1("c,j");

        a_klij("k,l,i,j") += g_ijak("k,l,c,j") * t1("c,i");

        a_klij("k,l,i,j") += g_abij("c,d,k,l") * tau("c,d,i,j");

        r2("a,b,i,j") += a_klij("k,l,i,j") * tau("a,b,k,l");

        if (print_detail_) {
          mpqc::detail::print_size_info(a_klij, "A_klij");
        }
      }
      tmp_time1 = mpqc::now(world, accurate_time);
      tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
      if (print_detail_) {
        mpqc::utility::print_par(world, "t2 a term time: ", tmp_time, "\n");
      }

      tmp_time0 = mpqc::now(world, accurate_time);
      {
        // compute b intermediate
        if (less) {
          // avoid store b_abcd
          TArray b_abij;
          b_abij("a,b,i,j") = g_abcd("a,b,c,d") * tau("c,d,i,j");

          b_abij("a,b,i,j") -= g_aibc("a,k,c,d") * tau("c,d,i,j") * t1("b,k");

          b_abij("a,b,i,j") -= g_iabc("k,b,c,d") * tau("c,d,i,j") * t1("a,k");

          if (print_detail_) {
            mpqc::detail::print_size_info(b_abij, "B_abij");
          }

          r2("a,b,i,j") += b_abij("a,b,i,j");
        } else {
          TArray b_abcd;

          b_abcd("a,b,c,d") = g_abcd("a,b,c,d") -
                              g_aibc("a,k,c,d") * t1("b,k") -
                              g_iabc("k,b,c,d") * t1("a,k");

          if (print_detail_) {
            mpqc::detail::print_size_info(b_abcd, "B_abcd");
          }

          r2("a,b,i,j") += b_abcd("a,b,c,d") * tau("c,d,i,j");
        }
      }
      tmp_time1 = mpqc::now(world, accurate_time);
      tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
      if (print_detail_) {
        mpqc::utility::print_par(world, "t2 b term time: ", tmp_time, "\n");
      }

      r2 = d_abij(r2, eps_, n_occ, n_frozen);

      r2("a,b,i,j") -= t2("a,b,i,j");

      t1("a,i") = t1("a,i") + r1("a,i");
      t2("a,b,i,j") = t2("a,b,i,j") + r2("a,b,i,j");
      t1.truncate();
      t2.truncate();

      auto t2_time1 = mpqc::now(world, accurate_time);
      auto t2_time = mpqc::duration_in_s(t2_time0, t2_time1);
      if (print_detail_) {
        mpqc::utility::print_par(world, "t2 total time: ", t2_time, "\n");
      }
      // recompute energy
      E0 = E1;

      E1_tmp1 = TA::dot(2.0 * f_ai("a,i"), t1("a,i"));
      E1_tmp2 = TA::dot((2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j")),
                        tau("a,b,i,j"));
      E1 = (E1_tmp1 + E1_tmp2).real();
      dE = std::abs(E0 - E1);

      if (dE >= converge_ || error >= converge_) {
        tmp_time0 = mpqc::now(world, accurate_time);
        cc::T1T2<TArray, TArray> t(t1, t2);
        cc::T1T2<TArray, TArray> r(r1, r2);
        error = r.norm() / (size(t1) + size(t2));  // error = residual norm per element
        diis.extrapolate(t, r);

        // update t1 and t2
        t1("a,i") = t.t1("a,i");
        t2("a,b,i,j") = t.t2("a,b,i,j");

        if (print_detail_) {
          mpqc::detail::print_size_info(r2, "R2");
          mpqc::detail::print_size_info(t2, "T2");
        }

        tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");
        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (print_detail_) {
          mpqc::utility::print_par(world, "diis time: ", tmp_time, "\n");
        }

        auto time1 = mpqc::fenced_now(world);
        auto duration = mpqc::duration_in_s(time0, time1);

        detail::print_gamma_point_ccsd(iter, dE, error, E1, duration);

        iter += 1ul;
      } else {
        auto time1 = mpqc::fenced_now(world);
        auto duration = mpqc::duration_in_s(time0, time1);

        detail::print_gamma_point_ccsd(iter, dE, error, E1, duration);

        break;
      }
    }
    if (iter >= max_iter_) {
      utility::print_par(this->wfn_world()->world(),
                         "\n Warning!! Exceed Max Iteration! \n");
    }

    ExEnv::out0() << "Gamma-Point CCSD Energy = " << E1 << std::endl;

    return E1;
  }

  /// <ab|ij>
  const TArray get_abij() {
    return this->lcao_factory().compute(L"<a b|G|i j>");
  }

  /// <ij|kl>
  const TArray get_ijkl() {
    return this->lcao_factory().compute(L"<i j|G|k l>");
  }

  /// <ab|cd>
  const TArray get_abcd() {
    return this->lcao_factory().compute(L"<a b|G|c d>");
  }

  /// <ia|jb>
  const TArray get_iajb() {
    return this->lcao_factory().compute(L"<i a|G|j b>");
  }

  /// <ia|bc>
  const TArray get_iabc() {
    return this->lcao_factory().compute(L"<i a|G|b c>");
  }

  /// <ai|bc>
  const TArray get_aibc() {
    return this->lcao_factory().compute(L"<a i|G|b c>");
  }

  /// <ij|ak>
  const TArray get_ijak() {
    return this->lcao_factory().compute(L"<i j|G|a k>");
  }

  /// <ij|ka>
  const TArray get_ijka() {
    return this->lcao_factory().compute(L"<i j|G|k a>");
  }

  /// <a|F|i>
  const TArray get_fock_ai() {
    return this->lcao_factory().compute(L"<a|F|i>");
  }

  /// returns a DIIS object
  TA::DIIS<cc::T1T2<TA::DistArray<Tile, Policy>, TA::DistArray<Tile, Policy>>>
  get_diis(const madness::World &world) {
    int n_diis, strt, ngr, ngrdiis;
    double dmp, mf;

    strt = kv_.value<int>("diis_strt", 1);
    n_diis = kv_.value<int>("n_diis", 8);
    ngr = kv_.value<int>("diis_ngr", 2);
    ngrdiis = kv_.value<int>("ngrdiis", 1);
    dmp = kv_.value<double>("diis_dmp", 0.0);
    mf = kv_.value<double>("diis_mf", 0.0);

    if (world.rank() == 0) {
      std::cout << "DIIS Starting Iteration:  " << strt << std::endl;
      std::cout << "DIIS Storing Size:  " << n_diis << std::endl;
      std::cout << "DIIS ngr:  " << ngr << std::endl;
      std::cout << "DIIS ngrdiis:  " << ngrdiis << std::endl;
      std::cout << "DIIS dmp:  " << dmp << std::endl;
      std::cout << "DIIS mf:  " << mf << std::endl;
    }
    TA::DIIS<cc::T1T2<TA::DistArray<Tile, Policy>, TA::DistArray<Tile, Policy>>>
        diis(strt, n_diis, 0.0, ngr, ngrdiis);

    return diis;
  }

  bool can_evaluate(Energy *energy) override {
    // can only evaluate the energy
    return energy->order() == 0;
  }

  void evaluate(Energy *result) override {
    if (!this->computed()) {
      /// cast ref_wfn to Energy::Provider
      auto ref_evaluator =
          std::dynamic_pointer_cast<typename Energy::Provider>(ref_wfn_);
      if (ref_evaluator == nullptr) {
        std::ostringstream oss;
        oss << "RefWavefunction in GammaPointCCSD" << ref_wfn_->class_key()
            << " cannot compute Energy" << std::endl;
        throw InputError(oss.str().c_str(), __FILE__, __LINE__);
      }

      ref_evaluator->evaluate(result);

      double ref_energy = result->energy();

      // initialize
      init();

      // set the precision
      if (converge_ == 0.0) {
        // if no user provided converge limit, use the default one from Energy
        converge_ = result->target_precision(0);
      }

      TArray t1, t2;

      gp_ccsd_corr_energy_ = compute_gamma_point_ccsd(t1, t2);

      this->computed_ = true;
      this->set_value(result, ref_energy + gp_ccsd_corr_energy_);
    }
  }
};

#if TA_DEFAULT_POLICY == 0
extern template class GammaPointCCSD<TA::TensorZ, TA::DensePolicy>;
#elif TA_DEFAULT_POLICY == 1
extern template class GammaPointCCSD<TA::TensorZ, TA::SparsePolicy>;
#endif

}  // namespace lcao
}  // namespace mpqc
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_GAMMA_POINT_CCSD_H
