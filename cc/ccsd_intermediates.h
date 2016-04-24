//
// Created by Chong Peng on 7/6/15.
//

#ifndef MPQC_CCSD_INTERMEDIATES_H
#define MPQC_CCSD_INTERMEDIATES_H

#include "../include/tiledarray.h"
#include "../common/namespaces.h"
#include "../integrals/direct_task_integrals.h"
#include "../integrals/molecular_integral.h"
#include "lazy_tile.h"

namespace mpqc {

namespace cc {

// TODO support four center AO integral
// class to compute the two electron integrals and intermediate needed for ccsd
// to construct this object needs:
//   three center integral Xpq
//   mo coefficient Cpi Cqa
//   direct ao integral array (optional for compute_dummy())
// p q r s stands for AO indices
// a b c d stands for virtual MO indices
// i j k l stands for occupied MO indices
// AO integrals uses chemical notation (pq|rs)
// three center uses chemical notation (X|pq)
// MO integrals used physical noation <ij|ab>
template <typename Tile, typename Policy,
          typename DirectTwoElectronArray = cc::DirectTwoElectronSparseArray>
class CCSDIntermediate {
  public:
    using TArray = TA::DistArray<Tile, Policy>;

    CCSDIntermediate(integrals::MolecularIntegral<Tile,Policy>& mo_int,
                     DirectTwoElectronArray &direct_ao)
            : mo_int_(mo_int), direct_ao_(direct_ao)
    {
            have_four_center_ = direct_ao.is_initialized();
    }

    /// clean the three center ingeral
    void clean_three_center() {
    }

    /// clean all the two electron integral computed
    void clean_two_electron() {
    }

    /// get mo coefficient
    /// occ part
    const TArray get_Ci() const {
        return mo_int_.orbital_space()->retrieve(OrbitalIndex(L"i")).array();
    }

    /// vir part
    const TArray get_Ca() const {
        return mo_int_.orbital_space()->retrieve(OrbitalIndex(L"a")).array();
    }

    /// get three center integral (X|ab)
    const TArray get_Xab() const {
        TArray result;
        TArray sqrt = mo_int_.atomic_integral().compute(L"(Κ|G| Λ)[inv_sqr]");
        TArray three_center = mo_int_.compute(L"(Κ|G|a b)");
        result("K,a,b") = sqrt("K,Q")*three_center("Q,a,b");
        return result;
    }

    /// get three center integral (X|ij)
    const TArray get_Xij() const {
        TArray result;
        TArray sqrt = mo_int_.atomic_integral().compute(L"(Κ|G| Λ)[inv_sqr]");
        TArray three_center = mo_int_.compute(L"(Κ|G|i j)");
        result("K,i,j") = sqrt("K,Q")*three_center("Q,i,j");
        return result;
    }

    /// get three center integral (X|ai)
    const TArray get_Xai() const {
        TArray result;
        TArray sqrt = mo_int_.atomic_integral().compute(L"(Κ|G| Λ)[inv_sqr]");
        TArray three_center = mo_int_.compute(L"(Κ|G|a i)");
        result("K,a,i") = sqrt("K,Q")*three_center("Q,a,i");
        return result;
    }

    // get two electron integrals
    // using physical notation <ab|ij>

    /// <ab|ij>
    const TArray get_abij() {
        return mo_int_.compute(L"<a b|G|i j>[df]");
    }

    /// <ij|kl>
    const TArray get_ijkl() {
        return mo_int_.compute(L"<i j|G|k l>[df]");
    }

    /// <ab|cd>
    const TArray get_abcd() {
        return  mo_int_.compute(L"<a b|G|c d>[df]");
    }

    /// <ia|bc>
    const TArray get_iabc() {
        return  mo_int_.compute(L"<i a|G|b c>[df]");
    }

    /// <ai|bc>
    const TArray get_aibc() {
        return mo_int_.compute(L"<a i|G|b c>[df]");
    }

    /// <ij|ak>
    const TArray get_ijak() {
        return mo_int_.compute(L"<i j|G|a k>[df]");
    }

    /// <ij|ka>
    const TArray get_ijka() {
        return mo_int_.compute(L"<i j|G|k a>[df]");
    }

    /// <ia|jb>
    const TArray get_iajb() {
        return mo_int_.compute(L"<i a|G|j b>[df]");
    }

    /// <a|f|i>
    const TArray get_fock_ai(){
        return mo_int_.compute(L"(a|F|i)[df]");
    }


    /// AO integral-direct computation of (ab|cd) ints contributions to the
    /// doubles resudual

    /// computes \f$ U^{ij}_{\rho\sigma} \equiv \left( t^{ij}_{\mu \nu} +
    /// t^{i}_{\mu} t^{j}_{\nu} \right) (\mu \rho| \nu \sigma) \f$
    /// @param t2 doubles amplitudes in MO basis
    /// @param t1 singles amplitudes in MO basis
    /// @return U tensor
    TArray compute_u2_u11(const TArray &t2, const TArray &t1) {
        if (have_four_center_) {
            TArray Ca = get_Ca();
            TArray tc;
            TArray u2_u11;
            tc("i,q") = Ca("q,c") * t1("c,i");
            u2_u11("p, r, i, j") = ((t2("a,b,i,j") * Ca("q,a")) * Ca("s,b")
                                    + tc("i,q") * tc("j,s"))
                                   * direct_ao_("p,q,r,s");
            return u2_u11;
        } else {
            throw std::runtime_error("CCSDIntermediate: integral-direct "
                                     "implementation used, but not "
                                     "initialized");
        }
    }

    TArray compute_u1a(const TArray &t1) {
        TArray Ci = get_Ci();
        TArray Ca = get_Ca();
        if (have_four_center_) {
            TArray u1a;
            u1a("p,r,i,j") = (Ca("q,a")*t1("a,i")*Ci("s,j")) *direct_ao_("p,q,r,s");
            return u1a;
        } else {
            throw std::runtime_error("CCSDIntermediate: integral-direct "
                                             "implementation used, but not "
                                             "initialized");
        }
    }

    TArray compute_u1b(const TArray &t1) {
        TArray Ci = get_Ci();
        TArray Ca = get_Ca();
        if (have_four_center_) {
            TArray u1b;
            u1b("i,r,j,s") = (Ca("p,a")*t1("a,i")*Ci("q,j")) *direct_ao_("p,q,r,s");
            return  u1b;

        } else {
            throw std::runtime_error("CCSDIntermediate: integral-direct "
                                             "implementation used, but not "
                                             "initialized");
        }
    }

    TArray compute_u_irjs() {
        TArray Ci = get_Ci();
        if (have_four_center_) {
            TArray u_ipjr;
            u_ipjr("i,r,j,s") = (Ci("p,i")*Ci("q,j"))*direct_ao_("p,q,r,s");
            return u_ipjr;
        } else {
            throw std::runtime_error("CCSDIntermediate: integral-direct "
                                             "implementation used, but not "
                                             "initialized");
        }
    }

    TArray compute_u_ijqs() {
        TArray Ci = get_Ci();
        if (have_four_center_) {
            TArray u_ijqs;
            u_ijqs("i,j,q,s") = (Ci("p,i")*Ci("r,j")) *direct_ao_("p,q,r,s");
            return u_ijqs;
        } else {
            throw std::runtime_error("CCSDIntermediate: integral-direct "
                                             "implementation used, but not "
                                             "initialized");
        }
    }
  private:

    /// MO Integral Object
    integrals::MolecularIntegral<Tile,Policy>& mo_int_;

    // direct ao in chemical notation (pq|rs)
    DirectTwoElectronArray direct_ao_;

    // check if have direct AO array
    bool have_four_center_;
};
}
}
#endif // MPQC_CCSD_INTERMEDIATES_H
