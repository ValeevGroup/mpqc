//
// Created by Chong Peng on 3/21/17.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_LCAO_CC_CCSD_HBAR_H_
#define SRC_MPQC_CHEMISTRY_QC_LCAO_CC_CCSD_HBAR_H_

/**
 * @brief This file contains functions to compute effective Hamiltonian of CCSD
 *
 * The spin orbital equations for intermediates can be found in
 * Gauss, J., & Stanton, J. F. (1995). The Journal of Chemical Physics, 103(9),
 * 3561–3577.
 *
 */

#include <tuple>

#include <tiledarray.h>

#include "mpqc/chemistry/qc/lcao/factory/lcao_factory.h"

namespace mpqc {
namespace lcao {
namespace cc {

/**
 *
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_W_KlIj(
    LCAOFactoryBase<Tile, Policy>& lcao_factory,
    const TA::DistArray<Tile, Policy>& t1,
    const TA::DistArray<Tile, Policy>& tau, bool df) {
  bool verbose = lcao_factory.verbose();
  auto& world = lcao_factory.world();
  bool accurate_time = lcao_factory.accurate_time();
  auto time0 = mpqc::now(world, accurate_time);

  std::wstring postfix = df ? L"[df]" : L"";
  TA::DistArray<Tile, Policy> W_KlIj;
  auto g_ijka = lcao_factory.compute(L"<i j|G|k a>" + postfix);
  auto g_ijkl = lcao_factory.compute(L"<i j|G|k l>" + postfix);
  auto g_ijab = lcao_factory.compute(L"<i j|G|a b>" + postfix);

  W_KlIj("k,l,i,j") = t1("c,j") * g_ijka("k,l,i,c");
  W_KlIj("k,l,i,j") = W_KlIj("k,l,i,j") + W_KlIj("l,k,j,i") +
                      g_ijkl("k,l,i,j") + g_ijab("k,l,c,d") * tau("c,d,i,j");

  auto time1 = mpqc::now(world, accurate_time);
  auto time = mpqc::duration_in_s(time0, time1);
  if (verbose) {
    ExEnv::out0() << indent
                  << "\nTime to Initialize W_KlIj Intermediates: " << time
                  << " S\n";
  }
  return W_KlIj;
};

template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_W_IbAj(
    LCAOFactoryBase<Tile, Policy>& lcao_factory,
    const TA::DistArray<Tile, Policy>& t1,
    const TA::DistArray<Tile, Policy>& t2, bool df) {
  bool verbose = lcao_factory.verbose();
  auto& world = lcao_factory.world();
  bool accurate_time = lcao_factory.accurate_time();
  auto time0 = mpqc::now(world, accurate_time);

  std::wstring postfix = df ? L"[df]" : L"";

  TA::DistArray<Tile, Policy> W_IbAj;

  auto g_ijab = lcao_factory.compute(L"<i j|G|a b>" + postfix);
  auto g_ijka = lcao_factory.compute(L"<i j|G|k a>" + postfix);

  W_IbAj("i,b,a,j") =
      g_ijab("i,j,a,b") - t1("b,l") * g_ijka("l,i,j,a") +
      t2("b,d,j,l") * (2.0 * g_ijab("i,l,a,d") - g_ijab("l,i,a,d")) -
      t2("d,b,j,l") * g_ijab("i,l,a,d") -
      t1("d,j") * (t1("b,l") * g_ijab("i,l,a,d"));

  if (df) {
    auto Xia = lcao_factory.compute(L"( Λ |G|i a)[inv_sqr]");
    auto Xab = lcao_factory.compute(L"( Λ |G|a b)[inv_sqr]");

    W_IbAj("i,b,a,j") += t1("d,j") * Xab("X,b,d") * Xia("X,i,a");

  } else {
    auto g_iabc = lcao_factory.compute(L"<i a|G|b c>" + postfix);
    W_IbAj("i,b,a,j") += t1("d,j") * g_iabc("i,b,a,d");
  }

  auto time1 = mpqc::now(world, accurate_time);
  auto time = mpqc::duration_in_s(time0, time1);
  if (verbose) {
    ExEnv::out0() << indent
                  << "\nTime to Initialize W_IbAj Intermediates: " << time
                  << " S\n";
  }

  return W_IbAj;
};

template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_W_IbaJ(
    LCAOFactoryBase<Tile, Policy>& lcao_factory,
    const TA::DistArray<Tile, Policy>& t1,
    const TA::DistArray<Tile, Policy>& t2, bool df) {
  bool verbose = lcao_factory.verbose();
  auto& world = lcao_factory.world();
  bool accurate_time = lcao_factory.accurate_time();
  auto time0 = mpqc::now(world, accurate_time);

  std::wstring postfix = df ? L"[df]" : L"";

  TA::DistArray<Tile, Policy> W_IbaJ;

  auto g_ijab = lcao_factory.compute(L"<i j|G|a b>" + postfix);
  auto g_iajb = lcao_factory.compute(L"<i a|G|j b>" + postfix);
  auto g_ijka = lcao_factory.compute(L"<i j|G|k a>" + postfix);

  W_IbaJ("i,b,a,j") = -g_iajb("j,a,i,b") + t1("b,l") * g_ijka("i,l,j,a") +
                      t2("d,b,j,l") * g_ijab("i,l,d,a") +
                      t1("d,j") * (t1("b,l") * g_ijab("i,l,d,a"));

  if (df) {
    auto Xia = lcao_factory.compute(L"( Λ |G|i a)[inv_sqr]");
    auto Xab = lcao_factory.compute(L"( Λ |G|a b)[inv_sqr]");

    W_IbaJ("i,b,a,j") -= t1("d,j") * Xia("X,i,d") * Xab("X,b,a");

  } else {
    auto g_iabc = lcao_factory.compute(L"<i a|G|b c>" + postfix);

    W_IbaJ("i,b,a,j") -= t1("d,j") * g_iabc("i,b,d,a");
  }

  auto time1 = mpqc::now(world, accurate_time);
  auto time = mpqc::duration_in_s(time0, time1);
  if (verbose) {
    ExEnv::out0() << indent
                  << "\nTime to Initialize W_IbaJ Intermediates: " << time
                  << " S\n";
  }
  return W_IbaJ;
};

template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_W_AkCd(
    LCAOFactoryBase<Tile, Policy>& lcao_factory,
    const TA::DistArray<Tile, Policy>& t1, bool df) {
  bool verbose = lcao_factory.verbose();
  auto& world = lcao_factory.world();
  bool accurate_time = lcao_factory.accurate_time();
  auto time0 = mpqc::now(world, accurate_time);

  std::wstring postfix = df ? L"[df]" : L"";

  TA::DistArray<Tile, Policy> W_AkCd;

  auto g_ijab = lcao_factory.compute(L"<i j|G|a b>" + postfix);
  auto g_iabc = lcao_factory.compute(L"<i a|G|b c>" + postfix);

  W_AkCd("a,k,c,d") = g_iabc("k,a,d,c") - t1("a,l") * g_ijab("l,k,c,d");

  auto time1 = mpqc::now(world, accurate_time);
  auto time = mpqc::duration_in_s(time0, time1);
  if (verbose) {
    ExEnv::out0() << indent
                  << "\nTime to Initialize W_AkCd Intermediates: " << time
                  << " S\n";
  }
  return W_AkCd;
};

template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_W_KlIc(
    LCAOFactoryBase<Tile, Policy>& lcao_factory,
    const TA::DistArray<Tile, Policy>& t1, bool df) {
  bool verbose = lcao_factory.verbose();
  auto& world = lcao_factory.world();
  bool accurate_time = lcao_factory.accurate_time();
  auto time0 = mpqc::now(world, accurate_time);
  std::wstring postfix = df ? L"[df]" : L"";

  TA::DistArray<Tile, Policy> W_KlIc;

  auto g_ijab = lcao_factory.compute(L"<i j|G|a b>" + postfix);
  auto g_ijka = lcao_factory.compute(L"<i j|G|k a>" + postfix);

  W_KlIc("k,l,i,c") = g_ijka("k,l,i,c") + t1("d,i") * g_ijab("k,l,d,c");

  auto time1 = mpqc::now(world, accurate_time);
  auto time = mpqc::duration_in_s(time0, time1);
  if (verbose) {
    ExEnv::out0() << indent
                  << "\nTime to Initialize W_KlIc Intermediates: " << time
                  << " S\n";
  }
  return W_KlIc;
};

template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_W_KaIj(
    LCAOFactoryBase<Tile, Policy>& lcao_factory,
    const TA::DistArray<Tile, Policy>& t1,
    const TA::DistArray<Tile, Policy>& t2,
    const TA::DistArray<Tile, Policy>& tau,
    const TA::DistArray<Tile, Policy>& FIA,
    const TA::DistArray<Tile, Policy>& W_KlIj, bool df) {
  bool verbose = lcao_factory.verbose();
  auto& world = lcao_factory.world();
  bool accurate_time = lcao_factory.accurate_time();
  auto time0 = mpqc::now(world, accurate_time);
  std::wstring postfix = df ? L"[df]" : L"";

  TA::DistArray<Tile, Policy> W_KaIj;

  auto g_ijab = lcao_factory.compute(L"<i j|G|a b>" + postfix);
  auto g_iajb = lcao_factory.compute(L"<i a|G|j b>" + postfix);
  auto g_iabc = lcao_factory.compute(L"<i a|G|b c>" + postfix);
  auto g_ijka = lcao_factory.compute(L"<i j|G|k a>" + postfix);

  TA::DistArray<Tile, Policy> t22;
  t22("a,b,i,j") = 2.0 * t2("a,b,i,j") - t2("a,b,j,i");

  W_KaIj("k,a,i,j") =

      g_ijka("i,j,k,a") + FIA("k,c") * t2("c,a,i,j") -

      t1("a,l") * W_KlIj("k,l,i,j") + g_iabc("k,a,c,d") * tau("c,d,i,j") +

      g_ijka("k,l,i,c") * t22("a,c,j,l") -

      g_ijka("l,k,i,c") * t2("a,c,j,l") - g_ijka("l,k,j,c") * t2("c,a,i,l") +

      t1("c,i") * (g_ijab("j,k,a,c") + t22("a,d,j,l") * g_ijab("k,l,c,d") -
                   t2("a,d,j,l") * g_ijab("k,l,d,c")) -

      t1("c,j") * (-g_iajb("i,c,k,a") + t2("a,d,l,i") * g_ijab("k,l,d,c"));

  auto time1 = mpqc::now(world, accurate_time);
  auto time = mpqc::duration_in_s(time0, time1);
  if (verbose) {
    ExEnv::out0() << indent
                  << "\nTime to Initialize W_KaIj Intermediates: " << time
                  << " S\n";
  }
  return W_KaIj;
};

template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_W_AbCi(
    LCAOFactoryBase<Tile, Policy>& lcao_factory,
    gaussian::AOFactoryBase<Tile, Policy>& ao_factory,
    const TA::DistArray<Tile, Policy>& t1,
    const TA::DistArray<Tile, Policy>& t2,
    const TA::DistArray<Tile, Policy>& tau,
    const TA::DistArray<Tile, Policy>& FIA,
    const TA::DistArray<Tile, Policy>& W_AbCd, bool df) {
  bool verbose = lcao_factory.verbose();
  auto& world = lcao_factory.world();
  bool accurate_time = lcao_factory.accurate_time();
  auto time0 = mpqc::now(world, accurate_time);
  std::wstring postfix = df ? L"[df]" : L"";

  auto g_ijab = lcao_factory.compute(L"<i j|G|a b>" + postfix);
  auto g_iajb = lcao_factory.compute(L"<i a|G|j b>" + postfix);
  auto g_iabc = lcao_factory.compute(L"<i a|G|b c>" + postfix);
  auto g_ijka = lcao_factory.compute(L"<i j|G|k a>" + postfix);

  TA::DistArray<Tile, Policy> t22;
  t22("a,b,i,j") = 2.0 * t2("a,b,i,j") - t2("a,b,j,i");

  TA::DistArray<Tile, Policy> W_AbCi;

  W_AbCi("a,b,c,i") =

      g_iabc("i,a,b,c")

      - FIA("k,c") * t2("a,b,k,i")

      + g_ijka("l,k,i,c") * tau("a,b,k,l")

      - g_iabc("k,b,c,d") * t2("a,d,k,i") + g_iabc("k,a,d,c") * t22("b,d,i,k") -
      g_iabc("k,a,c,d") * t2("b,d,i,k")

      - t1("a,k") * (g_ijab("i,k,b,c") + (t22("d,b,l,i") * g_ijab("k,l,c,d") -
                                          t2("d,b,l,i") * g_ijab("k,l,d,c"))) +
      t1("b,k") * (-g_iajb("k,a,i,c") + t2("a,d,l,i") * g_ijab("l,k,c,d"));

  // if W_AbCd term is computed and stored
  if (W_AbCd.is_initialized()) {
    W_AbCi("a,b,c,i") += t1("d,i") * W_AbCd("a,b,c,d");
  } else {
    if (df) {
      auto Xia = lcao_factory.compute(L"( Λ |G|i a)[inv_sqr]");
      auto Xab = lcao_factory.compute(L"( Λ |G|a b)[inv_sqr]");

      TA::DistArray<Tile, Policy> Xt;
      Xt("K,b,i") = (0.5 * Xab("K,b,d") - t1("b,k") * Xia("K,k,d")) * t1("d,i");

      W_AbCi("a,b,c,i") += Xab("K,a,c") * Xt("K,b,i") +
                           Xab("K,b,c") * Xt("K,a,i") +
                           t1("d,i") * g_ijab("k,l,c,d") * tau("a,b,k,l");

    } else {
      W_AbCi("a,b,c,i") += -t1("d,i") * g_iabc("k,a,d,c") * t1("b,k") -
                           t1("d,i") * g_iabc("k,b,c,d") * t1("a,k") +
                           t1("d,i") * g_ijab("k,l,c,d") * tau("a,b,k,l");

      // integral direct term
      auto direct_integral = ao_factory.compute_direct(L"(μ ν| G|κ λ)");
      auto Ca = lcao_factory.orbital_registry().retrieve(OrbitalIndex(L"a"));

      W_AbCi("a,b,c,i") +=
          (t1("d,i") * Ca("s,d") * direct_integral("p,q,r,s")) * Ca("r,b") *
          Ca("q,c") * Ca("p,a");
    }
  }

  auto time1 = mpqc::now(world, accurate_time);
  auto time = mpqc::duration_in_s(time0, time1);
  if (verbose) {
    ExEnv::out0() << indent
                  << "\nTime to Initialize W_AbCi Intermediates: " << time
                  << " S\n";
  }
  return W_AbCi;
};

/*
 *
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_W_AbCd(
    LCAOFactoryBase<Tile, Policy>& lcao_factory,
    const TA::DistArray<Tile, Policy>& t1,
    const TA::DistArray<Tile, Policy>& tau, bool df) {
  bool verbose = lcao_factory.verbose();
  auto& world = lcao_factory.world();
  bool accurate_time = lcao_factory.accurate_time();
  auto time0 = mpqc::now(world, accurate_time);
  std::wstring postfix = df ? L"[df]" : L"";

  TA::DistArray<Tile, Policy> W_AbCd;

  auto g_abcd = lcao_factory.compute(L"<a b|G|c d>" + postfix);
  auto g_ijab = lcao_factory.compute(L"<i j|G|a b>" + postfix);
  auto g_iabc = lcao_factory.compute(L"<i a|G|b c>" + postfix);

  W_AbCd("a,b,c,d") = g_abcd("a,b,c,d") - t1("b,k") * g_iabc("k,a,d,c") -
                      t1("a,k") * g_iabc("k,b,c,d") +
                      tau("a,b,k,l") * g_ijab("k,l,c,d");

  auto time1 = mpqc::now(world, accurate_time);
  auto time = mpqc::duration_in_s(time0, time1);
  if (verbose) {
    ExEnv::out0() << indent
                  << "\nTime to Initialize W_AbCd Intermediates: " << time
                  << " S\n";
  }
  return W_AbCd;
};

/**
 * compute closed shell one particle matrix of effective Hamiltonian of CCSD
 *
 * @return std::tuple of FAB, FIA, FIJ
 *
 */
template <typename Tile, typename Policy>
std::tuple<TA::DistArray<Tile, Policy>, TA::DistArray<Tile, Policy>,
           TA::DistArray<Tile, Policy>>
compute_cs_ccsd_F(LCAOFactoryBase<Tile, Policy>& lcao_factory,
                  gaussian::AOFactoryBase<Tile, Policy>& ao_factory,
                  const TA::DistArray<Tile, Policy>& t1,
                  const TA::DistArray<Tile, Policy>& tau, bool df) {
  bool verbose = lcao_factory.verbose();
  auto& world = lcao_factory.world();
  bool accurate_time = lcao_factory.accurate_time();
  auto time0 = mpqc::now(world, accurate_time);
  std::wstring postfix = df ? L"[df]" : L"";

  auto g_ijab_bar = lcao_factory.compute(L"<i j|G|a b>" + postfix);
  g_ijab_bar("i,j,a,b") = 2.0 * g_ijab_bar("i,j,a,b") - g_ijab_bar("j,i,a,b");

  auto f_ia = lcao_factory.compute(L"<i|F|a>" + postfix);

  // compute FIA
  TA::DistArray<Tile, Policy> FIA;
  { FIA("i,a") = f_ia("i,a") + g_ijab_bar("i,j,a,b") * t1("b,j"); }

  TA::DistArray<Tile, Policy> FIJ;
  {
    FIJ = lcao_factory.compute(L"<i|F|j>" + postfix);
    auto g_ijka = lcao_factory.compute(L"<i j|G|k a>" + postfix);
    FIJ("i,j") += t1("c,l") * (2.0 * g_ijka("i,l,j,c") - g_ijka("l,i,j,c")) +
                  tau("c,d,j,l") * g_ijab_bar("i,l,c,d") -
                  f_ia("i,a") * t1("a,j");
  }

  TA::DistArray<Tile, Policy> FAB;
  {
    FAB = lcao_factory.compute(L"<a|F|b>" + postfix);
    FAB("a,b") -=
        t1("a,m") * f_ia("m,b") + tau("a,d,k,l") * g_ijab_bar("k,l,b,d");

    if (df) {
      //     refactorize with density fitting
      auto Xia = lcao_factory.compute(L"( Λ |G|i a)[inv_sqr]");
      auto Xab = lcao_factory.compute(L"( Λ |G|a b)[inv_sqr]");

      FAB("a,b") += 2.0 * t1("d,k") * Xia("X,k,d") * Xab("X,a,b") -
                    t1("d,k") * Xab("X,a,d") * Xia("X,k,b");

    } else {
      auto g_iabc = lcao_factory.compute(L"<i a|G|b c>" + postfix);
      FAB("a,b") += t1("d,k") * (2.0 * g_iabc("k,a,d,b") - g_iabc("k,a,b,d"));
    }
  }
  auto time1 = mpqc::now(world, accurate_time);
  auto time = mpqc::duration_in_s(time0, time1);
  if (verbose) {
    ExEnv::out0() << indent << "\nTime to Initialize F Intermediates: " << time
                  << " S\n";
  }
  return std::make_tuple(FAB, FIA, FIJ);
}

}  // namespace cc
}  // namespace lcao
}  // namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_LCAO_CC_CCSD_HBAR_H_
