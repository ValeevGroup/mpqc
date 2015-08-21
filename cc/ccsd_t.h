//
// Created by Chong Peng on 8/21/15.
//

#ifndef TILECLUSTERCHEM_CCSD_T_H_H
#define TILECLUSTERCHEM_CCSD_T_H_H

#include "ccsd.h"

namespace tcc{
    namespace cc{

        // CCSD_T class that compute CCSD(T) (T) calculation
        template<typename Tile, typename Policy>
        class CCSD_T : public CCSD<Tile,Policy> {

        public:

            typedef TA::Array <double, 2, Tile, Policy> TArray2;
            typedef TA::Array <double, 4, Tile, Policy> TArray4;

            CCSD_T(const TArray2 &fock, const Eigen::VectorXd &ens,
                 const std::shared_ptr<TRange1Engine> &tre,
                 const std::shared_ptr<CCSDIntermediate<Tile, Policy>> &inter): CCSD<Tile,Policy>(fock,ens,tre,inter)
            {}


            void compute(){

                TArray2 t1;
                TArray4 t2;

                // compute ccsd first
                double ccsd_corr = CCSD<Tile,Policy>::compute_ccsd(t1,t2);

                double ccsd_t = compute_ccsd_t(t2);

            }

            double compute_ccsd_t(TArray4& t2){
                return 0;
            }
        };

    } // namespace cc
}  // namespace tcc

#endif //TILECLUSTERCHEM_CCSD_T_H_H
