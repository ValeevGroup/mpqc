//
// Created by Chong Peng on 7/1/15.
//

#ifndef MPQC_CHEMISTRY_QC_CC_CCSD_H
#define MPQC_CHEMISTRY_QC_CC_CCSD_H

#include <tiledarray.h>

#include "mpqc/chemistry/qc/wfn/trange1_engine.h"
#include "mpqc/chemistry/qc/mbpt/denom.h"
#include <mpqc/chemistry/qc/cc/ccsd_intermediates.h>
#include <mpqc/chemistry/qc/cc/diis_ccsd.h>
#include <mpqc/chemistry/qc/integrals/direct_atomic_integral.h>
#include <mpqc/chemistry/qc/scf/mo_build.h>
#include <mpqc/chemistry/qc/wfn/lcao_wfn.h>

namespace mpqc {
namespace cc {

/*
 * CCSD class that computed CCSD energy
 *   Options
 *   BlockSize = int, control the block size in MO, default 16
 *   OccBlockSize = int, control the block size in Occ, overide BlockSize
 *   VirBlockSize = int, control the block size in Vir, overide BlockSize
 *   FrozenCore = bool, control if use frozen core, default False
 *   Direct = bool , control if use direct approach, default True
 *   DIIS_ndi = int, control the number of data sets to retain, default is 5
 *   DIIS_dmp = float, see DIIS in TA
 *   DIIS_strt = int, see DIIS in TA
 *   DIIS_ngr = int, see DIIS in TA
 *   DIIS_ngrdiis = int, see DIIS in TA
 *   DIIS_mf = float, see DIIS in TA
 *   LessMemory = bool, control if store large intermediate in straight
 * approach, default is false
 *   PrintDetail = bool, if do detail printing, default is false
 *   Converge = double, convergence of CCSD energy, default is 1.0e-07
 */

template <typename Tile, typename Policy>
class CCSD : public qc::LCAOWavefunction<Tile, Policy> {
 public:
  using TArray = TA::DistArray<Tile, Policy>;
  using DirectAOIntegral = integrals::DirectAtomicIntegral<Tile, Policy>;

  CCSD() = default;
  ~CCSD();

  /**
   * KeyVal constructor
   * @param kv
   *
   * keywords
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | method | string | standard | method to compute ccsd (standard, df,
   * direct) |
   */
  CCSD(const KeyVal &kv) : qc::LCAOWavefunction<Tile, Policy>(kv), kv_(kv) {

    if (kv.exists("ref")) {
      ref_wfn_ = kv.keyval("ref").class_ptr<qc::Wavefunction>();
    } else {
      throw std::invalid_argument(
          "Default Ref Wfn in CCSD is not support! \n");
    }

    method_ = kv_.value<std::string>("method", "df");
    if (method_ == "df" || method_ == "direct") {
      df_ = true;
    }

    // initialize direct integral class
    if (method_ == "direct") {
      auto direct_ao_integral =
          integrals::detail::construct_direct_atomic_integral<Tile, Policy>(kv);
      direct_ao_ints_ = direct_ao_integral->compute(L"(μ ν| G|κ λ)");
    }

    max_iter_ = kv.value<int>("max_iter", 15);
    converge_ = kv.value<double>("converge", 1.0e-7);
    print_detail_ = kv.value<bool>("print_detail_", false);
  }

 protected:
  TArray T1_;
  TArray T2_;

 private:
  const KeyVal kv_;
  std::shared_ptr<qc::Wavefunction> ref_wfn_;
  typename DirectAOIntegral::DirectTArray direct_ao_ints_;
  bool df_;
  std::string method_;
  std::size_t max_iter_;
  double converge_;
  bool print_detail_;
  double ccsd_corr_energy_;

 public:
  void compute(qc::PropertyBase* pb) override {
    throw std::runtime_error("Not Implemented!!");
  }

  void obsolete() override {
    this->energy_ = 0.0;
    qc::LCAOWavefunction<Tile, Policy>::obsolete();
    ref_wfn_->obsolete();
  }

  /// compute function
  double value() override {

    if(this->energy_ == 0.0){

      double ref_energy = ref_wfn_->value();

      // initialize
      init();

      TArray t1;
      TArray t2;

      if (method_ == "standard") {
        ccsd_corr_energy_ = compute_ccsd_conventional(t1, t2);
      } else if (method_ == "df") {
        ccsd_corr_energy_ = compute_ccsd_df(t1, t2);
      } else if (method_ == "direct") {
        ccsd_corr_energy_ = compute_ccsd_direct(t1, t2);
      }

      T1_ = t1;
      T2_ = t2;

      this->energy_ = ref_energy + ccsd_corr_energy_;
    }

    return this->energy_;
  }

  void set_trange1_engine(const std::shared_ptr<TRange1Engine> &tr1) {
    this->trange1_engine_ = tr1;
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

  bool is_df() const {
    return df_;
  }

  void set_t1(const TArray &t1) { T1_ = t1; }

  void set_t2(const TArray &t2) { T2_ = t2; }

 protected:
  // store all the integrals in memory
  // used as reference for development
  double compute_ccsd_conventional(TArray &t1, TArray &t2) {
    auto &world = this->wfn_world()->world();
    bool accurate_time = this->lcao_factory().accurate_time();

    auto n_occ = this->trange1_engine()->get_occ();
    auto n_frozen = this->trange1_engine()->get_nfrozen();

    if (world.rank() == 0) {
      std::cout << "Use Conventional CCSD Compute" << std::endl;
    }

    auto tmp_time0 = mpqc::now(world, accurate_time);
    // get all two electron integrals
    TArray g_abij = this->get_abij();
    TArray g_ijkl = this->get_ijkl();
    TArray g_abcd = this->get_abcd();
    TArray g_iajb = this->get_iajb();
    TArray g_iabc = this->get_iabc();
    TArray g_aibc = this->get_aibc();
    TArray g_ijak = this->get_ijak();
    TArray g_ijka = this->get_ijka();
    auto tmp_time1 = mpqc::now(world, accurate_time);
    auto tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
    if (print_detail_) {
      mpqc::utility::print_par(world, "Integral Prepare Time: ", tmp_time,
                               "\n");
    }

    TArray f_ai = this->get_fock_ai();

    TArray d1(f_ai.world(), f_ai.trange(), f_ai.shape(), f_ai.pmap());
    // store d1 to local
    create_d_ai(d1, *this->orbital_energy(), n_occ, n_frozen);

    t1("a,i") = f_ai("a,i") * d1("a,i");

    t2 = d_abij(g_abij, *this->orbital_energy(), n_occ, n_frozen);

    TArray tau;
    tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");

    double E0 = 0.0;
    double E1 =
        2.0 * TA::dot(f_ai("a,i"), t1("a,i")) +
        TA::dot(2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j"), tau("a,b,i,j"));
    double mp2 = E1;
    double dE = std::abs(E1 - E0);

    mpqc::utility::print_par(world, "MP2 Energy      ", mp2, "\n");

    // optimize t1 and t2
    std::size_t iter = 0ul;
    double error = 1.0;
    TArray r1;
    TArray r2;

    bool less = kv_.value<bool>("less_memory", false);

    if (world.rank() == 0) {
      std::cout << "Start Iteration" << std::endl;
      std::cout << "Max Iteration" << max_iter_ << std::endl;
      std::cout << "Convergence " << converge_ << std::endl;
      std::cout << "AccurateTime" << accurate_time << std::endl;
      std::cout << "PrintDetail" << print_detail_ << std::endl;
      if (less) {
        std::cout << "Less Memory Approach: Yes" << std::endl;
      } else {
        std::cout << "Less Memory Approach: No" << std::endl;
      }
    }

    auto diis = get_diis(world);

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
        r2("a,b,i,j") =
            (g_iabc("i,c,a,b") - g_iajb("k,b,i,c") * t1("a,k")) * t1("c,j");

        r2("a,b,i,j") -=
            (g_ijak("i,j,a,k") + g_abij("a,c,i,k") * t1("c,j")) * t1("b,k");
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
            detail::print_size_info(T, "T");
            detail::print_size_info(j_akic, "J_akic");
            detail::print_size_info(k_kaic, "K_kaic");
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
          detail::print_size_info(a_klij, "A_klij");
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
            detail::print_size_info(b_abij, "B_abij");
          }

          r2("a,b,i,j") += b_abij("a,b,i,j");
        } else {
          TArray b_abcd;

          b_abcd("a,b,c,d") = g_abcd("a,b,c,d") -
                              g_aibc("a,k,c,d") * t1("b,k") -
                              g_iabc("k,b,c,d") * t1("a,k");

          if (print_detail_) {
            detail::print_size_info(b_abcd, "B_abcd");
          }

          r2("a,b,i,j") += b_abcd("a,b,c,d") * tau("c,d,i,j");
        }
      }
      tmp_time1 = mpqc::now(world, accurate_time);
      tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
      if (print_detail_) {
        mpqc::utility::print_par(world, "t2 b term time: ", tmp_time, "\n");
      }

      d_abij_inplace(r2, *this->orbital_energy(), n_occ, n_frozen);

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
      E1 = 2.0 * TA::dot(f_ai("a,i"), t1("a,i")) +
           TA::dot((2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j")),
                   tau("a,b,i,j"));
      dE = std::abs(E0 - E1);

      if (iter == 0 && world.rank() == 0) {
        std::cout << "iter "
                  << "    deltaE    "
                  << "            residual       "
                  << "      energy     "
                  << " total/second " << std::endl;
      }

      if (dE >= converge_ || error >= converge_) {
        tmp_time0 = mpqc::now(world, accurate_time);
        mpqc::cc::T1T2<double, Tile, Policy> t(t1, t2);
        mpqc::cc::T1T2<double, Tile, Policy> r(r1, r2);
        error = r.norm() / size(t);
        diis.extrapolate(t, r);

        // update t1 and t2
        t1("a,i") = t.first("a,i");
        t2("a,b,i,j") = t.second("a,b,i,j");

        if (print_detail_) {
          detail::print_size_info(r2, "R2");
          detail::print_size_info(t2, "T2");
        }

        tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");
        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (print_detail_) {
          mpqc::utility::print_par(world, "diis time: ", tmp_time, "\n");
        }

        auto time1 = mpqc::fenced_now(world);
        auto duration = mpqc::duration_in_s(time0, time1);

        if (world.rank() == 0) {
          std::cout << iter << "  " << dE << "  " << error << "  " << E1 << "  "
                    << duration << std::endl;
        }

        iter += 1ul;
      } else {
        auto time1 = mpqc::fenced_now(world);
        auto duration = mpqc::duration_in_s(time0, time1);

        if (world.rank() == 0) {
          std::cout << iter << "  " << dE << "  " << error << "  " << E1 << "  "
                    << duration << std::endl;
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

  double compute_ccsd_df(TArray &t1, TArray &t2) {
    auto &world = this->wfn_world()->world();
    bool accurate_time = this->lcao_factory().accurate_time();

    auto n_occ = this->trange1_engine()->get_occ();
    auto n_frozen = this->trange1_engine()->get_nfrozen();

    if (world.rank() == 0) {
      std::cout << "Use DF CCSD Compute" << std::endl;
    }

    auto tmp_time0 = mpqc::now(world, accurate_time);
    // get all two electron integrals
    TArray g_abij = this->get_abij();
    TArray g_ijkl = this->get_ijkl();
    TArray g_abcd = this->get_abcd();
    TArray X_ai = this->get_Xai();
    TArray g_iajb = this->get_iajb();
    TArray g_iabc = this->get_iabc();
    TArray g_aibc = this->get_aibc();
    TArray g_ijak = this->get_ijak();
    TArray g_ijka = this->get_ijka();
    auto tmp_time1 = mpqc::now(world, accurate_time);
    auto tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
    if (print_detail_) {
      mpqc::utility::print_par(world, "Integral Prepare Time: ", tmp_time,
                               "\n");
    }

    TArray f_ai = this->get_fock_ai();

    TArray d1(f_ai.world(), f_ai.trange(), f_ai.shape(), f_ai.pmap());
    // store d1 to local
    create_d_ai(d1, *this->orbital_energy(), n_occ, n_frozen);

    t1("a,i") = f_ai("a,i") * d1("a,i");

    t2 = d_abij(g_abij, *this->orbital_energy(), n_occ, n_frozen);

    TArray tau;
    tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");

    double E0 = 0.0;
    double E1 =
        2.0 * TA::dot(f_ai("a,i"), t1("a,i")) +
        TA::dot(2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j"), tau("a,b,i,j"));
    double mp2 = E1;
    double dE = std::abs(E1 - E0);

    mpqc::utility::print_par(world, "MP2 Energy      ", mp2, "\n");

    // optimize t1 and t2
    std::size_t iter = 0ul;
    double error = 1.0;
    TArray r1;
    TArray r2;

    bool less = kv_.value<bool>("less_memory", false);

    if (world.rank() == 0) {
      std::cout << "Start Iteration" << std::endl;
      std::cout << "Max Iteration" << max_iter_ << std::endl;
      std::cout << "Convergence " << converge_ << std::endl;
      std::cout << "AccurateTime" << accurate_time << std::endl;
      std::cout << "PrintDetail" << print_detail_ << std::endl;
      if (less) {
        std::cout << "Less Memory Approach: Yes" << std::endl;
      } else {
        std::cout << "Less Memory Approach: No" << std::endl;
      }
    }

    auto diis = get_diis(world);

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
        r2("a,b,i,j") =
            (g_iabc("i,c,a,b") - g_iajb("k,b,i,c") * t1("a,k")) * t1("c,j");

        r2("a,b,i,j") -=
            (g_ijak("i,j,a,k") + g_abij("a,c,i,k") * t1("c,j")) * t1("b,k");
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

          j_akic("a,k,i,c") -= X_ai("x,d,l") * T("d,a,i,l") * X_ai("x,c,k");

          j_akic("a,k,i,c") += 0.5 *
                               (2.0 * g_abij("c,d,k,l") - g_abij("d,c,k,l")) *
                               t2("a,d,i,l");

          k_kaic("k,a,i,c") = g_iajb("k,a,i,c")

                              - g_ijka("k,l,i,c") * t1("a,l")

                              + g_iabc("k,a,d,c") * t1("d,i")

                              - g_abij("d,c,k,l") * T("d,a,i,l");
          if (print_detail_) {
            detail::print_size_info(T, "T");
            detail::print_size_info(j_akic, "J_akic");
            detail::print_size_info(k_kaic, "K_kaic");
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
          detail::print_size_info(a_klij, "A_klij");
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
            detail::print_size_info(b_abij, "B_abij");
          }

          r2("a,b,i,j") += b_abij("a,b,i,j");
        } else {
          TArray b_abcd;

          b_abcd("a,b,c,d") = g_abcd("a,b,c,d") -
                              g_aibc("a,k,c,d") * t1("b,k") -
                              g_iabc("k,b,c,d") * t1("a,k");

          if (print_detail_) {
            detail::print_size_info(b_abcd, "B_abcd");
          }

          r2("a,b,i,j") += b_abcd("a,b,c,d") * tau("c,d,i,j");
        }
      }
      tmp_time1 = mpqc::now(world, accurate_time);
      tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
      if (print_detail_) {
        mpqc::utility::print_par(world, "t2 b term time: ", tmp_time, "\n");
      }

      d_abij_inplace(r2, *this->orbital_energy(), n_occ, n_frozen);

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
      E1 = 2.0 * TA::dot(f_ai("a,i"), t1("a,i")) +
           TA::dot((2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j")),
                   tau("a,b,i,j"));
      dE = std::abs(E0 - E1);

      if (iter == 0 && world.rank() == 0) {
        std::cout << "iter "
                  << "    deltaE    "
                  << "            residual       "
                  << "      energy     "
                  << " total/second " << std::endl;
      }

      if (dE >= converge_ || error >= converge_) {
        tmp_time0 = mpqc::now(world, accurate_time);
        mpqc::cc::T1T2<double, Tile, Policy> t(t1, t2);
        mpqc::cc::T1T2<double, Tile, Policy> r(r1, r2);
        error = r.norm() / size(t);
        diis.extrapolate(t, r);

        // update t1 and t2
        t1("a,i") = t.first("a,i");
        t2("a,b,i,j") = t.second("a,b,i,j");

        if (print_detail_) {
          detail::print_size_info(r2, "R2");
          detail::print_size_info(t2, "T2");
        }

        tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");
        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (print_detail_) {
          mpqc::utility::print_par(world, "diis time: ", tmp_time, "\n");
        }

        auto time1 = mpqc::fenced_now(world);
        auto duration = mpqc::duration_in_s(time0, time1);

        if (world.rank() == 0) {
          std::cout << iter << "  " << dE << "  " << error << "  " << E1 << "  "
                    << duration << std::endl;
        }

        iter += 1ul;
      } else {
        auto time1 = mpqc::fenced_now(world);
        auto duration = mpqc::duration_in_s(time0, time1);

        if (world.rank() == 0) {
          std::cout << iter << "  " << dE << "  " << error << "  " << E1 << "  "
                    << duration << std::endl;
        }

        break;
      }
      //        std::cout << indent << scprintf("%-5.0f", iter) <<
      //        scprintf("%-20.10f", Delta_E)
      //        << scprintf("%-15.10f", E_1) << std::endl;
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

  double compute_ccsd_direct(TArray &t1, TArray &t2) {
    auto &world = this->wfn_world()->world();
    // get mo coefficient
    TArray ca = this->get_Ca();
    TArray ci = this->get_Ci();

    // get three center integral
    TArray Xab = this->get_Xab();
    TArray Xij = this->get_Xij();
    TArray Xai = this->get_Xai();

    mpqc::utility::print_par(world, "Use Direct CCSD Compute \n");

    bool accurate_time = this->lcao_factory().accurate_time();

    auto n_occ = this->trange1_engine()->get_occ();
    auto n_frozen = this->trange1_engine()->get_nfrozen();

    auto tmp_time0 = mpqc::now(world, accurate_time);

    TArray g_abij = this->get_abij();

    TArray f_ai = this->get_fock_ai();

    TArray d1(f_ai.world(), f_ai.trange(), f_ai.shape(), f_ai.pmap());

    create_d_ai(d1, *this->orbital_energy(), n_occ, n_frozen);

    t1("a,i") = f_ai("a,i") * d1("a,i");

    t2 = d_abij(g_abij, *this->orbital_energy(), n_occ, n_frozen);

    //      std::cout << t1 << std::endl;
    //      std::cout << t2 << std::endl;

    // get all two electron integrals
    TArray g_ijkl = this->get_ijkl();
    TArray g_iajb = this->get_iajb();
    TArray g_ijak = this->get_ijak();
    TArray g_ijka = this->get_ijka();

    auto tmp_time1 = mpqc::now(world, accurate_time);
    auto tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
    if (print_detail_) {
      mpqc::utility::print_par(world, "Integral Prepare Time: ", tmp_time,
                               "\n");
    }

    TArray tau;
    tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");

    double E0 = 0.0;
    double E1 =
        2.0 * TA::dot(f_ai("a,i"), t1("a,i")) +
        TA::dot(2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j"), tau("a,b,i,j"));
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
      std::cout << "Max Iteration" << max_iter_ << std::endl;
      std::cout << "Convergence " << converge_ << std::endl;
      std::cout << "AccurateTime" << accurate_time << std::endl;
      std::cout << "PrintDetail" << print_detail_ << std::endl;
    };

    auto diis = get_diis(world);

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

      if (print_detail_) {
        detail::print_size_info(u2_u11, "U_aaoo");
        mpqc::utility::print_par(world, "u term time: ", duration_u, "\n");
      } else if (iter == 0) {
        detail::print_size_info(u2_u11, "U_aaoo");
      }

      auto t1_time0 = mpqc::now(world, accurate_time);
      TArray h_ac, h_ki;
      {
        // intermediates for t1
        // external index i and a

        tmp_time0 = mpqc::now(world, accurate_time);
        h_ac("a,c") =
            -(2.0 * g_abij("c,d,k,l") - g_abij("c,d,l,k")) * tau("a,d,k,l");

        h_ki("k,i") =
            (2.0 * g_abij("c,d,k,l") - g_abij("d,c,k,l")) * tau("c,d,i,l");

        // compute residual r1(n) = t1(n+1) - t1(n)
        // external index i and a

        r1("a,i") = f_ai("a,i") - 2.0 * f_ai("c,k") * t1("c,i") * t1("a,k");

        r1("a,i") += h_ac("a,c") * t1("c,i") - t1("a,k") * h_ki("k,i");

        {
          TArray h_kc;
          h_kc("k,c") =
              f_ai("c,k") +
              (-g_abij("d,c,k,l") + 2.0 * g_abij("c,d,k,l")) * t1("d,l");

          r1("a,i") += h_kc("k,c") * (2.0 * t2("c,a,k,i") - t2("c,a,i,k") +
                                      t1("a,k") * t1("c,i"));
        }

        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (print_detail_) {
          mpqc::utility::print_par(world, "t1 h term time: ", tmp_time, "\n");
        }

        tmp_time0 = mpqc::now(world, accurate_time);

        r1("a,i") += (2.0 * g_abij("c,a,k,i") - g_iajb("k,a,i,c")) * t1("c,k");

        r1("a,i") += (2.0 * u2_u11("p,r,k,i") - u2_u11("p,r,i,k")) * ci("p,k") *
                     ca("r,a");

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
      {
        tmp_time0 = mpqc::now(world, accurate_time);
        // permutation term
        {
          r2("a,b,i,j") = Xab("X,b,c") * t1("c,j") * Xai("X,a,i");

          r2("a,b,i,j") -= g_iajb("k,b,i,c") * t1("c,j") * t1("a,k");

          r2("a,b,i,j") -=
              (g_ijak("i,j,a,k") + g_abij("a,c,i,k") * t1("c,j")) * t1("b,k");
        }
        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (print_detail_) {
          mpqc::utility::print_par(world, "t2 other time: ", tmp_time, "\n");
        }

        tmp_time0 = mpqc::now(world, accurate_time);
        {
          // intermediates g
          TArray g_ki, g_ac;

          g_ki("k,i") =
              h_ki("k,i") + f_ai("c,k") * t1("c,i") +
              (2.0 * g_ijka("k,l,i,c") - g_ijka("l,k,i,c")) * t1("c,l");

          g_ac("a,c") = h_ac("a,c") - f_ai("c,k") * t1("a,k")

                        + 2.0 * Xai("X,d,k") * t1("d,k") * Xab("X,a,c")

                        - Xab("X,a,d") * t1("d,k") * Xai("X,c,k");

          r2("a,b,i,j") +=
              (g_ac("a,c") * t2("c,b,i,j") - g_ki("k,i") * t2("a,b,k,j"));
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

            j_akic("a,k,i,c") =
                g_abij("a,c,i,k")

                - g_ijka("l,k,i,c") * t1("a,l")

                + (Xab("X,a,d") * t1("d,i")) * Xai("X,c,k")

                - (Xai("X,d,l") * T("d,a,i,l")) * Xai("X,c,k")

                + (g_abij("c,d,k,l") - 0.5 * g_abij("d,c,k,l")) * t2("a,d,i,l");

            k_kaic("k,a,i,c") = g_iajb("k,a,i,c")

                                - g_ijka("k,l,i,c") * t1("a,l")

                                + (Xai("X,d,k") * t1("d,i")) * Xab("X,a,c")

                                - g_abij("d,c,k,l") * T("d,a,i,l");

            if (print_detail_) {
              detail::print_size_info(T, "T");
              detail::print_size_info(j_akic, "J_akic");
              detail::print_size_info(k_kaic, "K_kaic");
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

        // perform permutation
        r2("a,b,i,j") = r2("a,b,i,j") + r2("b,a,j,i");

        r2("a,b,i,j") += g_abij("a,b,i,j");

        tmp_time0 = mpqc::now(world, accurate_time);
        {
          // intermediate a
          TArray a_klij;
          a_klij("k,l,i,j") = g_ijkl("k,l,i,j")

                              + g_ijka("k,l,i,c") * t1("c,j")

                              + g_ijak("k,l,c,j") * t1("c,i")

                              + g_abij("c,d,k,l") * tau("c,d,i,j");

          r2("a,b,i,j") += a_klij("k,l,i,j") * tau("a,b,k,l");

          if (print_detail_) {
            detail::print_size_info(a_klij, "A_klij");
          }
        }
        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (print_detail_) {
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

          if (print_detail_) {
            detail::print_size_info(b_abij, "B_abij");
          }
        }
        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (print_detail_) {
          mpqc::utility::print_par(world, "t2 b term time: ", tmp_time, "\n");
        }

        d_abij_inplace(r2, *this->orbital_energy(), n_occ, n_frozen);

        r2("a,b,i,j") -= t2("a,b,i,j");
      }

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
      E1 = 2.0 * TA::dot(f_ai("a,i"), t1("a,i")) +
           TA::dot((2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j")),
                   tau("a,b,i,j"));
      dE = std::abs(E0 - E1);

      if (iter == 0 && world.rank() == 0) {
        std::cout << "iter "
                  << "    deltaE    "
                  << "            residual       "
                  << "      energy     "
                  << "    U/second  "
                  << " total/second " << std::endl;
      }

      if (dE >= converge_ || error >= converge_) {
        tmp_time0 = mpqc::now(world, accurate_time);
        mpqc::cc::T1T2<double, Tile, Policy> t(t1, t2);
        mpqc::cc::T1T2<double, Tile, Policy> r(r1, r2);
        error = r.norm() / size(t);
        diis.extrapolate(t, r);

        // update t1 and t2
        t1("a,i") = t.first("a,i");
        t2("a,b,i,j") = t.second("a,b,i,j");

        if (print_detail_) {
          detail::print_size_info(r2, "R2");
          detail::print_size_info(t2, "T2");
        }

        tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");
        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (print_detail_) {
          mpqc::utility::print_par(world, "diis time: ", tmp_time, "\n");
        }

        auto time1 = mpqc::fenced_now(world);
        auto duration_t = mpqc::duration_in_s(time0, time1);

        if (world.rank() == 0) {
          std::cout.precision(15);
          std::cout << iter << "  " << dE << "  " << error << "  " << E1 << "  "
                    << duration_u << " " << duration_t << std::endl;
        }

        iter += 1ul;
      } else {
        auto time1 = mpqc::fenced_now(world);
        auto duration_t = mpqc::duration_in_s(time0, time1);

        if (world.rank() == 0) {
          std::cout.precision(15);
          //                            std::cout.width(20);
          std::cout << iter << "  " << dE << "  " << error << "  " << E1 << "  "
                    << duration_u << " " << duration_t << std::endl;
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

 private:
  virtual void init() {
    if (this->orbital_energy() == nullptr || this->trange1_engine() == nullptr) {
      auto mol = this->lcao_factory().atomic_integral().molecule();
      Eigen::VectorXd orbital_energy;
      this->trange1_engine_ = closed_shell_obs_mo_build_eigen_solve(
          this->lcao_factory(), orbital_energy, mol, this->is_frozen_core(),
          this->occ_block(), this->unocc_block());
      this->orbital_energy_ = std::make_shared<Eigen::VectorXd>(orbital_energy);
    }
  }

  TA::DIIS<mpqc::cc::T1T2<double, Tile, Policy>> get_diis(
      const madness::World &world) {
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
    TA::DIIS<mpqc::cc::T1T2<double, Tile, Policy>> diis(strt, n_diis, 0.0, ngr,
                                                        ngrdiis);

    return diis;
  };

  /// get mo coefficient
  /// occ part
  const TArray get_Ci() {
    return this->lcao_factory()
        .orbital_space()
        .retrieve(OrbitalIndex(L"i"))
        .array();
  }

  /// vir part
  const TArray get_Ca() {
    return this->lcao_factory()
        .orbital_space()
        .retrieve(OrbitalIndex(L"a"))
        .array();
  }

  /// get three center integral (X|ab)
  const TArray get_Xab() {
    TArray result;
    TArray sqrt =
        this->lcao_factory().atomic_integral().compute(L"(Κ|G| Λ)[inv_sqr]");
    TArray three_center = this->lcao_factory().compute(L"(Κ|G|a b)");
    result("K,a,b") = sqrt("K,Q") * three_center("Q,a,b");
    return result;
  }

  /// get three center integral (X|ij)
  const TArray get_Xij() {
    TArray result;
    TArray sqrt =
        this->lcao_factory().atomic_integral().compute(L"(Κ|G| Λ)[inv_sqr]");
    TArray three_center = this->lcao_factory().compute(L"(Κ|G|i j)");
    result("K,i,j") = sqrt("K,Q") * three_center("Q,i,j");
    return result;
  }

  /// get three center integral (X|ai)
  const TArray get_Xai() {
    TArray result;
    TArray sqrt =
        this->lcao_factory().atomic_integral().compute(L"(Κ|G| Λ)[inv_sqr]");
    TArray three_center = this->lcao_factory().compute(L"(Κ|G|a i)");
    result("K,a,i") = sqrt("K,Q") * three_center("Q,a,i");
    return result;
  }

  // get two electron integrals
  // using physical notation <ab|ij>

  /// <ab|ij>
  const TArray get_abij() {
    if (df_) {
      return this->lcao_factory().compute(L"<a b|G|i j>[df]");
    } else {
      return this->lcao_factory().compute(L"<a b|G|i j>");
    }
  }

  /// <ij|kl>
  const TArray get_ijkl() {
    if (df_) {
      return this->lcao_factory().compute(L"<i j|G|k l>[df]");
    } else {
      return this->lcao_factory().compute(L"<i j|G|k l>");
    }
  }

  /// <ab|cd>
  const TArray get_abcd() {
    if (df_) {
      return this->lcao_factory().compute(L"<a b|G|c d>[df]");
    } else {
      return this->lcao_factory().compute(L"<a b|G|c d>");
    }
  }

  /// <ia|bc>
  const TArray get_iabc() {
    if (df_) {
      return this->lcao_factory().compute(L"<i a|G|b c>[df]");
    } else {
      return this->lcao_factory().compute(L"<i a|G|b c>");
    }
  }

  /// <ai|bc>
  const TArray get_aibc() {
    if (df_) {
      return this->lcao_factory().compute(L"<a i|G|b c>[df]");
    } else {
      return this->lcao_factory().compute(L"<a i|G|b c>");
    }
  }

  /// <ij|ak>
  const TArray get_ijak() {
    if (df_) {
      return this->lcao_factory().compute(L"<i j|G|a k>[df]");
    } else {
      return this->lcao_factory().compute(L"<i j|G|a k>");
    }
  }

  /// <ai|jk>
  const TArray get_aijk() {
    if (df_) {
      return this->lcao_factory().compute(L"<a i|G|j k>[df]");
    } else {
      return this->lcao_factory().compute(L"<a i|G|j k>");
    }
  }

  /// <ia|jb>
  const TArray get_iajb() {
    if (df_) {
      return this->lcao_factory().compute(L"<i a|G|j b>[df]");
    } else {
      return this->lcao_factory().compute(L"<i a|G|j b>");
    }
  }

  /// <ij|ka>
  const TArray get_ijka() {
    if (df_) {
      return this->lcao_factory().compute(L"<i j|G|k a>[df]");
    } else {
      return this->lcao_factory().compute(L"<i j|G|k a>");
    }
  }

  /// <a|f|i>
  const TArray get_fock_ai() {
    if (df_) {
      return this->lcao_factory().compute(L"<a|F|i>[df]");
    } else {
      return this->lcao_factory().compute(L"<a|F|i>");
    }
  }

  /// AO integral-direct computation of (ab|cd) ints contributions to the
  /// doubles resudual

  /// computes \f$ U^{ij}_{\rho\sigma} \equiv \left( t^{ij}_{\mu \nu} +
  /// t^{i}_{\mu} t^{j}_{\nu} \right) (\mu \rho| \nu \sigma) \f$
  /// @param t2 doubles amplitudes in MO basis
  /// @param t1 singles amplitudes in MO basis
  /// @return U tensor
  TArray compute_u2_u11(const TArray &t2, const TArray &t1) {
    if (direct_ao_ints_.array().is_initialized()) {
      TArray Ca = get_Ca();
      TArray tc;
      TArray u2_u11;
      tc("i,q") = Ca("q,c") * t1("c,i");
      u2_u11("p, r, i, j") =
          ((t2("a,b,i,j") * Ca("q,a")) * Ca("s,b") + tc("i,q") * tc("j,s")) *
          direct_ao_ints_("p,q,r,s");
      return u2_u11;
    } else {
      throw std::runtime_error(
          "CCSD: integral-direct implementation used, but direct integral not "
          "initialized");
    }
  }
};  // class CCSD

}  // namespace cc
}  // namespace mpqc

#endif  // MPQC_CHEMISTRY_QC_CC_CCSD_H
