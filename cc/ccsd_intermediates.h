//
// Created by Chong Peng on 7/6/15.
//

#ifndef TILECLUSTERCHEM_CCSD_INTERMEDIATES_H
#define TILECLUSTERCHEM_CCSD_INTERMEDIATES_H

#include "../include/tiledarray.h"
#include "../common/namespaces.h"
#include "integral_generator.h"
#include "lazy_integral.h"

namespace tcc {

    namespace cc {

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
        template<typename Tile, typename Policy,
                typename IntegralGenerator=tcc::cc::TwoBodyIntGenerator<libint2::Coulomb> >
        class CCSDIntermediate {
        public:

            typedef TA::Array <double, 2, Tile, Policy> TArray2;
            typedef TA::Array <double, 3, Tile, Policy> TArray3;
            typedef TA::Array <double, 4, Tile, Policy> TArray4;

            typedef tcc::cc::LazyIntegral<4, IntegralGenerator> LazyTwoElectronTile;
            typedef TA::Array <double, 4, LazyTwoElectronTile, Policy> DirectTwoElectronArray;

            CCSDIntermediate(const TArray3 &Xpq, const TArray2 &Ci, const TArray2 &Ca,
                             const DirectTwoElectronArray& direct_ao=DirectTwoElectronArray()) {
                {
                    // convert to MO, store this temporary,
                    // call clean() to clean Xab_ Xij_ Xai_
                    Xab_("X,a,b") = Xpq("X,p,q") * Ca("q,a") * Ca("p,b");
                    Xij_("X,i,j") = Xpq("X,p,q") * Ci("q,i") * Ci("p,j");
                    Xai_("X,a,i") = Xpq("X,p,q") * Ca("q,a") * Ci("p,i");
                    Ci_= Ci;
                    Ca_ = Ca;

                    if(direct_ao.is_initialized()){
                        direct_ = true;
                    }else{
                        direct_ = false;
                    }
                    direct_ao_ = direct_ao;
                    cleaned_ = false;
                }
                TArray3::wait_for_lazy_cleanup(Xpq.get_world());
            }

            // clean the three center ingeral
            void clean(){
                Xab_ = TArray3();
                Xai_ = TArray3();
                Xij_ = TArray3();
                cleaned_ = true;
            }

            // get mo coefficient
            // occ part
            const TArray2 get_Ci() const{
                return Ci_;
            }

            // vir part
            const TArray2 get_Ca() const{
                return Ca_;
            }

            // get three center integral (X|ab)
            const TArray3 get_Xab() const{
                if(!cleaned_){
                    return Xab_;
                }else{
                    throw std::runtime_error("CCSDIntermediate has been cleaned");
                }
            }

            // get three center integral (X|ij)
            const TArray3 get_Xij() const{
                if(!cleaned_){
                    return Xij_;
                }else{
                    throw std::runtime_error("CCSDIntermediate has been cleaned");
                }
            }

            // get three center integral (X|ai)
            const TArray3 get_Xai() const{
                if(!cleaned_){
                    return Xai_;
                }else{
                    throw std::runtime_error("CCSDIntermediate has been cleaned");
                }
            }

            // get two electron integrals
            // using physical notation <ab|ij>

            // <ab|ij>
            const TArray4 get_abij() const {
                if(!cleaned_){
                    TArray4 abij;
                    abij("a,b,i,j") = Xai_("X,a,i") * Xai_("X,b,j");
                    return abij;
                }else{
                    throw std::runtime_error("CCSDIntermediate has been cleaned");
                }
            }

            // <ij|kl>
            const TArray4 get_ijkl() const {
                if(!cleaned_){
                    TArray4 ijkl;
                    ijkl("i,j,k,l") = Xij_("X,i,k") * Xij_("X,j,l");
                    return ijkl;
                }else{
                    throw std::runtime_error("CCSDIntermediate has been cleaned");
                }
            }

            // <ab|cd>
            const TArray4 get_abcd() const {
                if(!cleaned_){
                    TArray4 abcd;
                    abcd("a,b,c,d") = Xab_("X,a,c") * Xab_("X,b,d");
                    //std::cout << abcd_ << std::endl;
                    return abcd;
                }else{
                    throw std::runtime_error("CCSDIntermediate has been cleaned");
                }
            }

            // <ia|bc>
            const TArray4 get_iabc() const {
                if(!cleaned_){
                    TArray4 iabc;
                    iabc("i,a,b,c") = Xab_("X,a,c") * Xai_("X,b,i");
                    return iabc;
                }else{
                    throw std::runtime_error("CCSDIntermediate has been cleaned");
                }
            }

            // <ai|bc>
            const TArray4 get_aibc() const {
                if(!cleaned_){
                    TArray4 aibc;
                    aibc("a,i,b,c") = Xai_("X,c,i") * Xab_("X,a,b");
                    return aibc;
                }else{
                    throw std::runtime_error("CCSDIntermediate has been cleaned");
                }
            }

            // <ij|ak>
            const TArray4 get_ijak() const {
                if(!cleaned_){
                    TArray4 ijak;
                    ijak("i,j,a,k") = Xai_("X,a,i") * Xij_("X,j,k");
                    return ijak;
                }else{
                    throw std::runtime_error("CCSDIntermediate has been cleaned");
                }
            }

            // <ij|ka>
            const TArray4 get_ijka() const {
                if(!cleaned_){
                    TArray4 ijka;
                    ijka("i,j,k,a") = Xai_("X,a,j") * Xij_("X,i,k");
                    return ijka;
                } else{
                    throw std::runtime_error("CCSDIntermediate has been cleaned");
                }
            }

            // <ia|jb>
            const TArray4 get_iajb() const {
                if(!cleaned_){
                    TArray4 iajb;
                    iajb("i,a,j,b") = Xab_("X,a,b") * Xij_("X,i,j");
                    return iajb;
                } else{
                    throw std::runtime_error("CCSDIntermediate has been cleaned");
                }
            }

            /// AO integral-direct computation of (ab|cd) ints contributions to the doubles resudual

            /// computes \f$ U^{ij}_{\rho\sigma} \equiv \left( t^{ij}_{\mu \nu} + t^{i}_{\mu} t^{j}_{\nu} \right) (\mu \rho| \nu \sigma) \f$
            /// @param t2 doubles amplitudes in MO basis
            /// @param t1 singles amplitudes in MO basis
            /// @return U tensor
            TArray4 compute_u2_u11(const TArray4& t2, const TArray2& t1){
                if (direct_){
                    TArray2 tc;
                    TArray4 u2_u11;
                    tc("i,q") = Ca_("q,c") * t1("c,i");
                    u2_u11("p, r, i, j") = ((t2("a,b,i,j")*Ca_("q,a"))*Ca_("s,b") +
                                                 tc("i,q") * tc("j,s")) * direct_ao_("p,q,r,s");
                    return u2_u11;
                }else{
                    throw std::runtime_error("CCSDIntermediate: integral-direct implementation used, but not initialized");
                }
            }

        private:

            // three center integral, need to be cleaned when not needed
            TArray3 Xab_;
            TArray3 Xai_;
            TArray3 Xij_;

            // mo coefficient
            TArray2 Ci_;
            TArray2 Ca_;

            // direct ao
            // in chemical notation (pq|rs)
            DirectTwoElectronArray direct_ao_;

            // check if Xab, Xai, Xij has been cleaned
            bool cleaned_;

            // check if have direct AO array
            bool direct_;


        };

    }
}
#endif //TILECLUSTERCHEM_CCSD_INTERMEDIATES_H