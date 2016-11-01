//
// Created by Chong Peng on 7/1/15.
//

#ifndef MPQC_CCSD_H
#define MPQC_CCSD_H

#include <tiledarray.h>



#include "../../../../../utility/cc_utility.h"
#include "../../../../../utility/trange1_engine.h"
#include <mpqc/chemistry/qc/cc/ccsd_intermediates.h>
#include <mpqc/chemistry/qc/cc/diis_ccsd.h>
#include <mpqc/chemistry/qc/cc/mo_block.h>
#include <mpqc/chemistry/qc/scf/mo_build.h>

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
class CCSD {
 public:
  using TArray = TA::DistArray<Tile, Policy>;
  CCSD() = default;

  /// constructor with already computed components
  CCSD(const std::shared_ptr<CCSDIntermediate<Tile, Policy>> &inter,
       const std::shared_ptr<Eigen::VectorXd> &ens,
       const std::shared_ptr<TRange1Engine> &tre, rapidjson::Document &options)
      : ccsd_intermediate_(inter),
        orbital_energy_(ens),
        trange1_engine_(tre),
        options_(std::move(options)) {}

  CCSD(integrals::LCAOFactory<Tile, Policy> &lcao_factory,
       rapidjson::Document &options)
      : options_(std::move(options)) {
    auto &world = lcao_factory.world();

    bool df;
    std::string method =
        options_.HasMember("Method") ? options_["Method"].GetString() : "df";
    if (method == "four center") {
      df = false;
    } else if (method == "df") {
      df = true;
    } else {
      throw std::runtime_error("Wrong CCSD Method");
    }

    std::string screen =
        options_.HasMember("Screen") ? options_["Screen"].GetString() : "";
    int screen_option = 0;
    if (screen == "schwarz") {
      screen_option = 1;
    } else if (screen == "qqr") {
      screen_option = 2;
    }

    typename cc::CCSDIntermediate<Tile,Policy>::DirectTwoElectronArray direct_two_electron_int;
    auto direct =
        options_.HasMember("Direct") ? options_["Direct"].GetBool() : true;
    if (direct) {
      // find the basis
      basis::Basis basis;
      // if VBS
      if (lcao_factory.atomic_integral().orbital_basis_registry().have(
              OrbitalIndex(L"Α"))) {
        basis =
            lcao_factory.atomic_integral().orbital_basis_registry().retrieve(
                OrbitalIndex(L"Α"));
      } else {
        basis =
            lcao_factory.atomic_integral().orbital_basis_registry().retrieve(
                OrbitalIndex(L"μ"));
      }

      std::vector<TA::TiledRange1> tr_04(4, basis.create_trange1());
      TA::TiledRange trange_4(tr_04.begin(), tr_04.end());
      auto time0 = mpqc::fenced_now(world);

      direct_two_electron_int =  cc::make_direct_two_electron_sparse_array(
          world, basis, trange_4, screen_option);

      auto time1 = mpqc::fenced_now(world);
      auto duration = mpqc::duration_in_s(time0, time1);

      if (world.rank() == 0) {
        std::cout << "Time to initialize direct two electron sparse "
                     "integral: "
                  << duration << std::endl;
      }
    }

    ccsd_intermediate_ =
        std::make_shared<mpqc::cc::CCSDIntermediate<Tile, Policy>>(
            lcao_factory, direct_two_electron_int, df);
  }

  /// compute function
  virtual double compute() {
    // initialize
    init(options_);

    TArray t1;
    TArray t2;

    auto direct =
        options_.HasMember("Direct") ? options_["Direct"].GetBool() : true;
    auto df = ccsd_intermediate_->is_df();
    double ccsd_corr = 0.0;
    if (!df) {
      ccsd_corr = compute_ccsd_conventional(t1, t2);
    } else {
      if (direct) {
        ccsd_corr = compute_ccsd_direct(t1, t2);
      } else {
        ccsd_corr = compute_ccsd_df(t1, t2);
      }
    }

    T1_ = t1;
    T2_ = t2;

    //                ccsd_intermediate_->clean_two_electron();

    return ccsd_corr;
  }

  const std::shared_ptr<TRange1Engine> &trange1_engine() const {
    return trange1_engine_;
  }

  void set_trange1_engine(const std::shared_ptr<TRange1Engine>& tr1){
    trange1_engine_ = tr1;
  }

  const rapidjson::Document &options() const { return options_; }

  const std::shared_ptr<CCSDIntermediate<Tile, Policy>> &intermediate() const {
    return ccsd_intermediate_;
  }

  const std::shared_ptr<Eigen::VectorXd> &orbital_energy() const {
    return orbital_energy_;
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

  void set_t1(const TArray& t1){
    T1_ = t1;
  }

  void set_t2(const TArray& t2){
    T2_ = t2;
  }

 protected:
  // store all the integrals in memory
  // used as reference for development
  double compute_ccsd_conventional(TArray &t1, TArray &t2) {
    auto &world = ccsd_intermediate_->get_Ca().world();

    bool print_detail = options_.HasMember("PrintDetail")
                            ? options_["PrintDetail"].GetBool()
                            : false;
    bool accurate_time = options_.HasMember("AccurateTime")
                             ? options_["AccurateTime"].GetBool()
                             : false;

    auto n_occ = trange1_engine_->get_occ();
    auto n_frozen = trange1_engine()->get_nfrozen();

    if (world.rank() == 0) {
      std::cout << "Use Conventional CCSD Compute" << std::endl;
    }

    auto tmp_time0 = mpqc::now(world, accurate_time);
    // get all two electron integrals
    TArray g_abij = ccsd_intermediate_->get_abij();
    TArray g_ijkl = ccsd_intermediate_->get_ijkl();
    TArray g_abcd = ccsd_intermediate_->get_abcd();
    TArray g_iajb = ccsd_intermediate_->get_iajb();
    TArray g_iabc = ccsd_intermediate_->get_iabc();
    TArray g_aibc = ccsd_intermediate_->get_aibc();
    TArray g_ijak = ccsd_intermediate_->get_ijak();
    TArray g_ijka = ccsd_intermediate_->get_ijka();
    auto tmp_time1 = mpqc::now(world, accurate_time);
    auto tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
    if (print_detail) {
      mpqc::utility::print_par(world, "Integral Prepare Time: ", tmp_time,
                               "\n");
    }

    TArray f_ai = ccsd_intermediate_->get_fock_ai();
    world.gop.fence();

    TArray d1(f_ai.world(), f_ai.trange(), f_ai.shape(),
              f_ai.pmap());
    // store d1 to local
    mpqc::cc::create_d_ai(d1, *orbital_energy_, n_occ, n_frozen);

    t1("a,i") = f_ai("a,i") * d1("a,i");

    t2 = d_abij(g_abij, *orbital_energy_, n_occ, n_frozen);

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

    std::size_t max_iter =
        options_.HasMember("MaxIter") ? options_["MaxIter"].GetInt() : 30;
    bool less = options_.HasMember("LessMemory")
                    ? options_["LessMemory"].GetBool()
                    : false;

    double converge = options_.HasMember("Converge")
                          ? options_["Converge"].GetDouble()
                          : 1.0e-7;

    if (world.rank() == 0) {
      std::cout << "Start Iteration" << std::endl;
      std::cout << "Max Iteration" << max_iter << std::endl;
      std::cout << "Convergence " << converge << std::endl;
      std::cout << "AccurateTime" << accurate_time << std::endl;
      std::cout << "PrintDetail" << print_detail << std::endl;
      if (less) {
        std::cout << "Less Memory Approach: Yes" << std::endl;
      } else {
        std::cout << "Less Memory Approach: No" << std::endl;
      }
    }

    auto diis = get_diis(world);

    while (iter < max_iter) {
      // start timer
      auto time0 = mpqc::now();

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
        if (print_detail) {
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
        if (print_detail) {
          mpqc::utility::print_par(world, "t1 other time: ", tmp_time, "\n");
        }
      }
      auto t1_time1 = mpqc::now(world, accurate_time);
      auto t1_time = mpqc::duration_in_s(t1_time0, t1_time1);
      if (print_detail) {
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
      if (print_detail) {
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
      if (print_detail) {
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
          if (print_detail) {
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
      if (print_detail) {
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

        if (print_detail) {
          detail::print_size_info(a_klij, "A_klij");
        }
      }
      tmp_time1 = mpqc::now(world, accurate_time);
      tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
      if (print_detail) {
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

          if (print_detail) {
            detail::print_size_info(b_abij, "B_abij");
          }

          r2("a,b,i,j") += b_abij("a,b,i,j");
        } else {
          TArray b_abcd;

          b_abcd("a,b,c,d") = g_abcd("a,b,c,d") -
                              g_aibc("a,k,c,d") * t1("b,k") -
                              g_iabc("k,b,c,d") * t1("a,k");

          if (print_detail) {
            detail::print_size_info(b_abcd, "B_abcd");
          }

          r2("a,b,i,j") += b_abcd("a,b,c,d") * tau("c,d,i,j");
        }
      }
      tmp_time1 = mpqc::now(world, accurate_time);
      tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
      if (print_detail) {
        mpqc::utility::print_par(world, "t2 b term time: ", tmp_time, "\n");
      }

      d_abij_inplace(r2, *orbital_energy_, n_occ, n_frozen);

      r2("a,b,i,j") -= t2("a,b,i,j");

      t1("a,i") = t1("a,i") + r1("a,i");
      t2("a,b,i,j") = t2("a,b,i,j") + r2("a,b,i,j");
      t1.truncate();
      t2.truncate();

      auto t2_time1 = mpqc::now(world, accurate_time);
      auto t2_time = mpqc::duration_in_s(t2_time0, t2_time1);
      if (print_detail) {
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

      if (dE >= converge || error >= converge) {
        tmp_time0 = mpqc::now(world, accurate_time);
        mpqc::cc::T1T2<double, Tile, Policy> t(t1, t2);
        mpqc::cc::T1T2<double, Tile, Policy> r(r1, r2);
        error = r.norm() / size(t);
        diis.extrapolate(t, r);

        // update t1 and t2
        t1("a,i") = t.first("a,i");
        t2("a,b,i,j") = t.second("a,b,i,j");

        if (print_detail) {
          detail::print_size_info(r2, "R2");
          detail::print_size_info(t2, "T2");
        }

        tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");
        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (print_detail) {
          mpqc::utility::print_par(world, "diis time: ", tmp_time, "\n");
        }

        world.gop.fence();
        auto time1 = mpqc::now();
        auto duration = mpqc::duration_in_s(time0, time1);

        if (world.rank() == 0) {
          std::cout << iter << "  " << dE << "  " << error << "  " << E1 << "  "
                    << duration << std::endl;
        }

        iter += 1ul;
      } else {
        world.gop.fence();
        auto time1 = mpqc::now();
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
    if (iter >= max_iter) {
      throw std::runtime_error("Exceed Max Iteration!");
    }
    if (world.rank() == 0) {
      std::cout << "CCSD Energy  " << E1 << std::endl;
    }
    return E1;
  }

  double compute_ccsd_df(TArray &t1, TArray &t2) {
    auto &world = ccsd_intermediate_->get_Ca().world();

    bool print_detail = options_.HasMember("PrintDetail")
                            ? options_["PrintDetail"].GetBool()
                            : false;
    bool accurate_time = options_.HasMember("AccurateTime")
                             ? options_["AccurateTime"].GetBool()
                             : false;

    auto n_occ = trange1_engine_->get_occ();
    auto n_frozen = trange1_engine()->get_nfrozen();

    if (world.rank() == 0) {
      std::cout << "Use DF CCSD Compute" << std::endl;
    }

    auto tmp_time0 = mpqc::now(world, accurate_time);
    // get all two electron integrals
    TArray g_abij = ccsd_intermediate_->get_abij();
    TArray g_ijkl = ccsd_intermediate_->get_ijkl();
    TArray g_abcd = ccsd_intermediate_->get_abcd();
    TArray X_ai = ccsd_intermediate_->get_Xai();
    TArray g_iajb = ccsd_intermediate_->get_iajb();
    TArray g_iabc = ccsd_intermediate_->get_iabc();
    TArray g_aibc = ccsd_intermediate_->get_aibc();
    TArray g_ijak = ccsd_intermediate_->get_ijak();
    TArray g_ijka = ccsd_intermediate_->get_ijka();
    auto tmp_time1 = mpqc::now(world, accurate_time);
    auto tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
    if (print_detail) {
      mpqc::utility::print_par(world, "Integral Prepare Time: ", tmp_time,
                               "\n");
    }

    TArray f_ai = ccsd_intermediate_->get_fock_ai();
    world.gop.fence();

    TArray d1(f_ai.world(), f_ai.trange(), f_ai.shape(),
              f_ai.pmap());
    // store d1 to local
    mpqc::cc::create_d_ai(d1, *orbital_energy_, n_occ, n_frozen);

    t1("a,i") = f_ai("a,i") * d1("a,i");

    t2 = d_abij(g_abij, *orbital_energy_, n_occ, n_frozen);

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

    std::size_t max_iter =
        options_.HasMember("MaxIter") ? options_["MaxIter"].GetInt() : 30;
    bool less = options_.HasMember("LessMemory")
                    ? options_["LessMemory"].GetBool()
                    : false;

    double converge = options_.HasMember("Converge")
                          ? options_["Converge"].GetDouble()
                          : 1.0e-7;

    if (world.rank() == 0) {
      std::cout << "Start Iteration" << std::endl;
      std::cout << "Max Iteration" << max_iter << std::endl;
      std::cout << "Convergence " << converge << std::endl;
      std::cout << "AccurateTime" << accurate_time << std::endl;
      std::cout << "PrintDetail" << print_detail << std::endl;
      if (less) {
        std::cout << "Less Memory Approach: Yes" << std::endl;
      } else {
        std::cout << "Less Memory Approach: No" << std::endl;
      }
    }

    auto diis = get_diis(world);

    while (iter < max_iter) {
      // start timer
      auto time0 = mpqc::now();

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
        if (print_detail) {
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
        if (print_detail) {
          mpqc::utility::print_par(world, "t1 other time: ", tmp_time, "\n");
        }
      }
      auto t1_time1 = mpqc::now(world, accurate_time);
      auto t1_time = mpqc::duration_in_s(t1_time0, t1_time1);
      if (print_detail) {
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
      if (print_detail) {
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
      if (print_detail) {
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
          if (print_detail) {
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
      if (print_detail) {
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

        if (print_detail) {
          detail::print_size_info(a_klij, "A_klij");
        }
      }
      tmp_time1 = mpqc::now(world, accurate_time);
      tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
      if (print_detail) {
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

          if (print_detail) {
            detail::print_size_info(b_abij, "B_abij");
          }

          r2("a,b,i,j") += b_abij("a,b,i,j");
        } else {
          TArray b_abcd;

          b_abcd("a,b,c,d") = g_abcd("a,b,c,d") -
                              g_aibc("a,k,c,d") * t1("b,k") -
                              g_iabc("k,b,c,d") * t1("a,k");

          if (print_detail) {
            detail::print_size_info(b_abcd, "B_abcd");
          }

          r2("a,b,i,j") += b_abcd("a,b,c,d") * tau("c,d,i,j");
        }
      }
      tmp_time1 = mpqc::now(world, accurate_time);
      tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
      if (print_detail) {
        mpqc::utility::print_par(world, "t2 b term time: ", tmp_time, "\n");
      }

      d_abij_inplace(r2, *orbital_energy_, n_occ, n_frozen);

      r2("a,b,i,j") -= t2("a,b,i,j");

      t1("a,i") = t1("a,i") + r1("a,i");
      t2("a,b,i,j") = t2("a,b,i,j") + r2("a,b,i,j");
      t1.truncate();
      t2.truncate();

      auto t2_time1 = mpqc::now(world, accurate_time);
      auto t2_time = mpqc::duration_in_s(t2_time0, t2_time1);
      if (print_detail) {
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

      if (dE >= converge || error >= converge) {
        tmp_time0 = mpqc::now(world, accurate_time);
        mpqc::cc::T1T2<double, Tile, Policy> t(t1, t2);
        mpqc::cc::T1T2<double, Tile, Policy> r(r1, r2);
        error = r.norm() / size(t);
        diis.extrapolate(t, r);

        // update t1 and t2
        t1("a,i") = t.first("a,i");
        t2("a,b,i,j") = t.second("a,b,i,j");

        if (print_detail) {
          detail::print_size_info(r2, "R2");
          detail::print_size_info(t2, "T2");
        }

        tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");
        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (print_detail) {
          mpqc::utility::print_par(world, "diis time: ", tmp_time, "\n");
        }

        world.gop.fence();
        auto time1 = mpqc::now();
        auto duration = mpqc::duration_in_s(time0, time1);

        if (world.rank() == 0) {
          std::cout << iter << "  " << dE << "  " << error << "  " << E1 << "  "
                    << duration << std::endl;
        }

        iter += 1ul;
      } else {
        world.gop.fence();
        auto time1 = mpqc::now();
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
    if (iter >= max_iter) {
      throw std::runtime_error("Exceed Max Iteration!");
    }
    if (world.rank() == 0) {
      std::cout << "CCSD Energy  " << E1 << std::endl;
    }
    return E1;
  }

  double compute_ccsd_direct(TArray &t1, TArray &t2) {
    // get mo coefficient
    TArray ca = ccsd_intermediate_->get_Ca();
    TArray ci = ccsd_intermediate_->get_Ci();

    // get three center integral
    TArray Xab = ccsd_intermediate_->get_Xab();
    TArray Xij = ccsd_intermediate_->get_Xij();
    TArray Xai = ccsd_intermediate_->get_Xai();

    auto &world = ca.world();

    mpqc::utility::print_par(world, "Use Direct CCSD Compute \n");

    bool print_detail = options_.HasMember("PrintDetail")
                            ? options_["PrintDetail"].GetBool()
                            : false;
    bool accurate_time = options_.HasMember("AccurateTime")
                             ? options_["AccurateTime"].GetBool()
                             : false;

    auto n_occ = trange1_engine_->get_occ();
    auto n_frozen = trange1_engine()->get_nfrozen();

    auto tmp_time0 = mpqc::now(world, accurate_time);

    TArray g_abij = ccsd_intermediate_->get_abij();

    TArray f_ai = ccsd_intermediate_->get_fock_ai();

    world.gop.fence();

    //      std::cout << g_abij << std::endl;

    TArray d1(f_ai.world(), f_ai.trange(), f_ai.shape(),
              f_ai.pmap());

    create_d_ai(d1, *orbital_energy_, n_occ, n_frozen);

    t1("a,i") = f_ai("a,i") * d1("a,i");

    t2 = d_abij(g_abij, *orbital_energy_, n_occ, n_frozen);

    //      std::cout << t1 << std::endl;
    //      std::cout << t2 << std::endl;

    // get all two electron integrals
    TArray g_ijkl = ccsd_intermediate_->get_ijkl();
    TArray g_iajb = ccsd_intermediate_->get_iajb();
    TArray g_ijak = ccsd_intermediate_->get_ijak();
    TArray g_ijka = ccsd_intermediate_->get_ijka();

    auto tmp_time1 = mpqc::now(world, accurate_time);
    auto tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
    if (print_detail) {
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

    std::size_t max_iter =
        options_.HasMember("MaxIter") ? options_["MaxIter"].GetInt() : 30;
    double converge = options_.HasMember("Converge")
                          ? options_["Converge"].GetDouble()
                          : 1.0e-7;

    if (world.rank() == 0) {
      std::cout << "Start Iteration" << std::endl;
      std::cout << "Max Iteration" << max_iter << std::endl;
      std::cout << "Convergence " << converge << std::endl;
      std::cout << "AccurateTime" << accurate_time << std::endl;
      std::cout << "PrintDetail" << print_detail << std::endl;
    };

    auto diis = get_diis(world);

    while (iter < max_iter) {
      // start timer
      auto time0 = mpqc::now();

      TArray::wait_for_lazy_cleanup(world);

      TArray u2_u11;
      // compute half transformed intermediates
      auto tu0 = mpqc::now();
      { u2_u11 = ccsd_intermediate_->compute_u2_u11(t2, t1); }
      world.gop.fence();
      auto tu1 = mpqc::now();
      auto duration_u = mpqc::duration_in_s(tu0, tu1);

      if (print_detail) {
        detail::print_size_info(u2_u11, "U_aaoo");
        mpqc::utility::print_par(world, "u term time: ", duration_u, "\n");
      } else if (iter == 0) {
        detail::print_size_info(u2_u11, "U_aaoo");
      }

      //                    if (g_abij.world().rank() == 0) {
      //                        std::cout << "Time to compute U intermediates
      //                        " << duration << std::endl;
      //                    }

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
        if (print_detail) {
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
        if (print_detail) {
          mpqc::utility::print_par(world, "t1 other time: ", tmp_time, "\n");
        }
      }
      auto t1_time1 = mpqc::now(world, accurate_time);
      auto t1_time = mpqc::duration_in_s(t1_time0, t1_time1);
      if (print_detail) {
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
        if (print_detail) {
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
        if (print_detail) {
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

            if (print_detail) {
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
        if (print_detail) {
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

          if (print_detail) {
            detail::print_size_info(a_klij, "A_klij");
          }
        }
        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (print_detail) {
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

          if (print_detail) {
            detail::print_size_info(b_abij, "B_abij");
          }
        }
        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (print_detail) {
          mpqc::utility::print_par(world, "t2 b term time: ", tmp_time, "\n");
        }

        d_abij_inplace(r2, *orbital_energy_, n_occ, n_frozen);

        r2("a,b,i,j") -= t2("a,b,i,j");
      }

      t1("a,i") = t1("a,i") + r1("a,i");
      t2("a,b,i,j") = t2("a,b,i,j") + r2("a,b,i,j");
      t1.truncate();
      t2.truncate();

      auto t2_time1 = mpqc::now(world, accurate_time);
      auto t2_time = mpqc::duration_in_s(t2_time0, t2_time1);
      if (print_detail) {
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

      if (dE >= converge || error >= converge) {
        tmp_time0 = mpqc::now(world, accurate_time);
        mpqc::cc::T1T2<double, Tile, Policy> t(t1, t2);
        mpqc::cc::T1T2<double, Tile, Policy> r(r1, r2);
        error = r.norm() / size(t);
        diis.extrapolate(t, r);

        // update t1 and t2
        t1("a,i") = t.first("a,i");
        t2("a,b,i,j") = t.second("a,b,i,j");

        if (print_detail) {
          detail::print_size_info(r2, "R2");
          detail::print_size_info(t2, "T2");
        }

        tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");
        tmp_time1 = mpqc::now(world, accurate_time);
        tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
        if (print_detail) {
          mpqc::utility::print_par(world, "diis time: ", tmp_time, "\n");
        }

        world.gop.fence();
        auto time1 = mpqc::now();
        auto duration_t = mpqc::duration_in_s(time0, time1);

        if (world.rank() == 0) {
          std::cout.precision(15);
          //                            std::cout.width(20);
          std::cout << iter << "  " << dE << "  " << error << "  " << E1 << "  "
                    << duration_u << " " << duration_t << std::endl;
        }

        iter += 1ul;
      } else {
        world.gop.fence();
        auto time1 = mpqc::now();
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
    if (iter >= max_iter) {
      throw std::runtime_error("Exceed Max Iteration!");
    }
    if (world.rank() == 0) {
      std::cout << "CCSD Energy     " << E1 << std::endl;
    }
    return E1;
  }

 private:
  void init(const rapidjson::Document &in) {
    if (orbital_energy_ == nullptr || trange1_engine_ == nullptr) {
      auto &lcao_factory = ccsd_intermediate_->lcao_factory();
      auto mol = lcao_factory.atomic_integral().molecule();
      Eigen::VectorXd orbital_energy;
      trange1_engine_ = closed_shell_obs_mo_build_eigen_solve(
          lcao_factory, orbital_energy, in, mol);
      orbital_energy_ = std::make_shared<Eigen::VectorXd>(orbital_energy);
    }
  }

  TA::DIIS<mpqc::cc::T1T2<double, Tile, Policy>> get_diis(
      const madness::World &world) {
    int n_diis, strt, ngr, ngrdiis;
    double dmp, mf;

    strt = options_.HasMember("DIIS_strt") ? options_["DIIS_strt"].GetInt() : 5;
    n_diis = options_.HasMember("DIIS_ndi") ? options_["DIIS_ndi"].GetInt() : 8;
    ngr = options_.HasMember("DIIS_ngr") ? options_["DIIS_ngr"].GetInt() : 2;
    ngrdiis = options_.HasMember("DIIS_ngrdiis")
                  ? options_["DIIS_ngrdiis"].GetInt()
                  : 1;
    dmp =
        options_.HasMember("DIIS_dmp") ? options_["DIIS_dmp"].GetDouble() : 0.0;
    mf = options_.HasMember("DIIS_mf") ? options_["DIIS_mf"].GetDouble() : 0.0;

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

 protected:
  // CCSD intermediate
  std::shared_ptr<mpqc::cc::CCSDIntermediate<Tile, Policy>> ccsd_intermediate_;

  // orbital energy
  std::shared_ptr<Eigen::VectorXd> orbital_energy_;

  // TA::TiledRange1 Engine class
  std::shared_ptr<mpqc::TRange1Engine> trange1_engine_;

  // option member
  rapidjson::Document options_;

  TArray T1_;
  TArray T2_;
};  // class CCSD

}  // namespace cc
}  // namespace mpqc

#endif  // MPQC_CCSD_H
