//
// Created by Chong Peng on 7/6/15.
//

#ifndef TILECLUSTERCHEM_CCSD_INTERMEDIATES_H
#define TILECLUSTERCHEM_CCSD_INTERMEDIATES_H

#include "../include/tiledarray.h"
#include "../common/namespaces.h"
#include "../integrals/direct_task_integrals.h"
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

    CCSDIntermediate(const TArray &Ci, const TArray &Ca, const TArray &Xpq,
                     DirectTwoElectronArray &direct_ao)
            : direct_ao_(direct_ao) {
        {
            // convert to MO, store this temporary,
            // call clean() to clean Xab_ Xij_ Xai_
            if (Xpq.is_initialized()) {
                Xab_("X,a,b") = Xpq("X,p,q") * Ca("q,a") * Ca("p,b");
                utility::print_size_info(Xab_, "X_xvv");
                Xij_("X,i,j") = Xpq("X,p,q") * Ci("q,i") * Ci("p,j");
                utility::print_size_info(Xij_, "X_xoo");
                Xai_("X,a,i") = Xpq("X,p,q") * Ca("q,a") * Ci("p,i");
                utility::print_size_info(Xai_, "X_xvo");
                have_three_center_ = true;
            } else {
                have_three_center_ = false;
            }
            Ci_ = Ci;
            Ca_ = Ca;

            if (direct_ao.is_initialized()) {
                have_four_center_ = true;
            } else {
                have_four_center_ = false;
            }

            // two electron integral
            //                    TArray abij_;
            //                    TArray ijkl_ = TArray();
            //                    TArray abcd_ = TArray();
            //                    TArray iabc_ = TArray();
            //                    TArray aibc_ = TArray();
            //                    TArray ijak_ = TArray();
            //                    TArray ijka_ = TArray();
            //                    TArray iajb_ = TArray();
            TArray abij_;
            TArray ijkl_;
            TArray abcd_;
            TArray iabc_;
            TArray aibc_;
            TArray ijak_;
            TArray ijka_;
            TArray iajb_;
        }
        TArray::wait_for_lazy_cleanup(Xpq.get_world());
    }

    // clean the three center ingeral
    void clean_three_center() {
        Xab_ = TArray();
        Xai_ = TArray();
        Xij_ = TArray();
        have_three_center_ = false;
    }

    // clean all the two electron integral computed
    void clean_two_electron() {
        abij_ = TArray();
        ijkl_ = TArray();
        abcd_ = TArray();
        iabc_ = TArray();
        aibc_ = TArray();
        ijak_ = TArray();
        ijka_ = TArray();
        iajb_ = TArray();
    }

    // get mo coefficient
    // occ part
    const TArray get_Ci() const { return Ci_; }

    // vir part
    const TArray get_Ca() const { return Ca_; }

    // get three center integral (X|ab)
    const TArray get_Xab() const {
        if (have_three_center_) {
            return Xab_;
        } else {
            throw std::runtime_error(
                  "CCSDIntermediate does not have three center");
        }
    }

    // get three center integral (X|ij)
    const TArray get_Xij() const {
        if (have_three_center_) {
            return Xij_;
        } else {
            throw std::runtime_error(
                  "CCSDIntermediate does not have three center");
        }
    }

    // get three center integral (X|ai)
    const TArray get_Xai() const {
        if (have_three_center_) {
            return Xai_;
        } else {
            throw std::runtime_error(
                  "CCSDIntermediate does not have three center");
        }
    }

    // get two electron integrals
    // using physical notation <ab|ij>

    // <ab|ij>
    const TArray get_abij() {
        if (have_three_center_) {
            if (abij_.is_initialized()) {
                return abij_;
            } else {
                abij_("a,b,i,j") = Xai_("X,a,i") * Xai_("X,b,j");
                utility::print_size_info(abij_, "G_vvoo");
                return abij_;
            }
        }  else {
            throw std::runtime_error(
                   "CCSDIntermediate does not have three center");
            return abij_;
        }
        
        // Removing until warning is fixed. 
        // else if (have_four_center_) {

        // } else {
        //     throw std::runtime_error(
        //           "CCSDIntermediate does not have three center");
        // }
    }

    // <ij|kl>
    const TArray get_ijkl() {
        if (have_three_center_) {
            if (ijkl_.is_initialized()) {
                return ijkl_;
            } else {
                ijkl_("i,j,k,l") = Xij_("X,i,k") * Xij_("X,j,l");
                utility::print_size_info(ijkl_, "G_oooo");
                return ijkl_;
            }
        } else {
            throw std::runtime_error(
                  "CCSDIntermediate does not have three center");
            return ijkl_;
        }
    }

    // <ab|cd>
    const TArray get_abcd() {
        if (have_three_center_) {
            if (abcd_.is_initialized()) {
                return abcd_;
            } else {
                abcd_("a,b,c,d") = Xab_("X,a,c") * Xab_("X,b,d");
                utility::print_size_info(abcd_, "G_vvvv");
                return abcd_;
            }
        } else {
            throw std::runtime_error(
                  "CCSDIntermediate does not have three center");
        }
    }

    // <ia|bc>
    const TArray get_iabc() {
        if (have_three_center_) {
            if (iabc_.is_initialized()) {
                return iabc_;
            } else {
                iabc_("i,a,b,c") = Xab_("X,a,c") * Xai_("X,b,i");
                utility::print_size_info(iabc_, "G_ovvv");
                return iabc_;
            }
        } else {
            throw std::runtime_error(
                  "CCSDIntermediate does not have three center");
        }
    }

    // <ai|bc>
    const TArray get_aibc() {
        if (have_three_center_) {
            if (aibc_.is_initialized()) {
                return aibc_;
            } else {
                aibc_("a,i,b,c") = Xai_("X,c,i") * Xab_("X,a,b");
                utility::print_size_info(aibc_, "G_vovv");
                return aibc_;
            }
        } else {
            throw std::runtime_error(
                  "CCSDIntermediate does not have three center");
        }
    }

    // <ij|ak>
    const TArray get_ijak() {
        if (have_three_center_) {
            if (ijak_.is_initialized()) {
                return ijak_;
            } else {
                ijak_("i,j,a,k") = Xai_("X,a,i") * Xij_("X,j,k");
                utility::print_size_info(ijak_, "G_oovo");
                return ijak_;
            }
        } else {
            throw std::runtime_error(
                  "CCSDIntermediate does not have three center");
        }
    }

    // <ij|ka>
    const TArray get_ijka() {
        if (have_three_center_) {
            if (ijka_.is_initialized()) {
                return ijka_;
            } else {
                ijka_("i,j,k,a") = Xai_("X,a,j") * Xij_("X,i,k");
                utility::print_size_info(ijka_, "G_ooov");
                return ijka_;
            }
        } else {
            throw std::runtime_error(
                  "CCSDIntermediate does not have three center");
        }
    }

    // <ia|jb>
    const TArray get_iajb() {
        if (have_three_center_) {
            if (iajb_.is_initialized()) {
                return iajb_;
            } else {
                iajb_("i,a,j,b") = Xab_("X,a,b") * Xij_("X,i,j");
                utility::print_size_info(iajb_, "G_ovov");
                return iajb_;
            }
        } else {
            throw std::runtime_error(
                  "CCSDIntermediate does not have three center");
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
        if (have_four_center_) {
            TArray tc;
            TArray u2_u11;
            tc("i,q") = Ca_("q,c") * t1("c,i");
            u2_u11("p, r, i, j") = ((t2("a,b,i,j") * Ca_("q,a")) * Ca_("s,b")
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
        if (have_four_center_) {
            TArray u1a;
            u1a("p,r,i,j") = (Ca_("q,a")*t1("a,i")*Ci_("s,j")) *direct_ao_("p,q,r,s");
            return u1a;
        } else {
            throw std::runtime_error("CCSDIntermediate: integral-direct "
                                             "implementation used, but not "
                                             "initialized");
        }
    }

    TArray compute_u1b(const TArray &t1) {
        if (have_four_center_) {
            TArray u1b;
            u1b("i,r,j,s") = (Ca_("p,a")*t1("a,i")*Ci_("q,j")) *direct_ao_("p,q,r,s");
            return  u1b;

        } else {
            throw std::runtime_error("CCSDIntermediate: integral-direct "
                                             "implementation used, but not "
                                             "initialized");
        }
    }

    TArray compute_u_irjs() {
        if (have_four_center_) {
            TArray u_ipjr;
            u_ipjr("i,r,j,s") = (Ci_("p,i")*Ci_("q,j"))*direct_ao_("p,q,r,s");
            return u_ipjr;
        } else {
            throw std::runtime_error("CCSDIntermediate: integral-direct "
                                             "implementation used, but not "
                                             "initialized");
        }
    }

    TArray compute_u_ijqs() {
        if (have_four_center_) {
            TArray u_ijqs;
            u_ijqs("i,j,q,s") = (Ci_("p,i")*Ci_("r,j")) *direct_ao_("p,q,r,s");
            return u_ijqs;
        } else {
            throw std::runtime_error("CCSDIntermediate: integral-direct "
                                             "implementation used, but not "
                                             "initialized");
        }
    }
  private:
    // three center integral, need to be cleaned when not needed
    TArray Xab_;
    TArray Xai_;
    TArray Xij_;

    // two electron integral
    TArray abij_;
    TArray ijkl_;
    TArray abcd_;
    TArray iabc_;
    TArray aibc_;
    TArray ijak_;
    TArray ijka_;
    TArray iajb_;

    // mo coefficient
    TArray Ci_;
    TArray Ca_;

    // direct ao
    // in chemical notation (pq|rs)
    DirectTwoElectronArray direct_ao_;

    // check if Xab, Xai, Xij has been cleaned
    bool have_three_center_;

    // check if have direct AO array
    bool have_four_center_;
};
}
}
#endif // TILECLUSTERCHEM_CCSD_INTERMEDIATES_H
