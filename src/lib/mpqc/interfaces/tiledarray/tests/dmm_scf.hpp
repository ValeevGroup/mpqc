/*
 * dmm_scf.hpp
 *
 *  Created on: Aug 14, 2013
 *      Author: drewlewis
 */

#ifndef mpqc_interfaces_tiledarray_tests_dmm_scf_hpp
#define mpqc_interfaces_tiledarray_tests_dmm_scf_hpp

#include "common.hpp"
#include "tiledarray_fock.hpp"
#include <TiledArray/algebra/diis.h>
#include "curvy_steps.hpp"

namespace mpqc {
namespace tests {

    using Array2 = TA::Array<double,2>;
    using Array3 = TA::Array<double,3>;
    using Array4 = TA::Array<double,4>;

    // Function that controls how the desnity will be updated, may change for
    // Different methods. I should make an attemp to wrap this in class.
    void
    Density_Update(Array2 &R, const Array2 &S, const Array2 &F,
                   const std::size_t iter,
                   const std::size_t rot_order = 4){
        curvy_steps::Density_Update(R, S, F, iter, rot_order);
    }

    double DF_DMM(Array2 &D, const Array2 &S, const Array2 &H,
                  Array2 &F, Array2 &G, const Array3 &Eri3,
                  const Array2 &Inv_Eri2, double nuc_repl){

        // create a DIIS object to extrapolate Fock matrix
        TA::DIIS<Array2> diis;

        // How many purification steps to use in Density code
        int scf_iter = 1;
        double energy = 0;
        double error_norminf = 1.0;

        // Begin SCF iterations
        while(error_norminf > 4e-6){
            double iter0 = madness::wall_time(); // Iteration timer

            //Density Update
            double mad_conj0 = madness::wall_time();
            Density_Update(D, S, F, scf_iter);
            double mad_conjf = madness::wall_time();

            //Fock Build
            double Fock0 = madness::wall_time(); // Fock build timer
            G("i,j") = 2.0 *
                 (Eri3("i,j,X") * Inv_Eri2("X,Y") * (Eri3("n,m,Y") * D("m,n"))) -
                 (Eri3("i,n,X") * Inv_Eri2("X,Y") * (Eri3("j,m,Y") * D("m,n")));
            F("i,j") = H("i,j") + G("i,j");
            S.world().gop.fence();

            // Computing gradient for DIIS error calculation
            Array2 gradient = 8 * ( S("i,q") * D("q,x") * F("x,j") -
                                                 F("i,q") * D("q,x") * S("x,j") );

            //Performing DIIS update of the Fock Matrix
            error_norminf = TA::expressions::norminf(gradient("i,j"));
            diis.extrapolate(F, gradient);

            //Energy for this iteration
            energy = TA::expressions::dot( 2.0 * H("i,j") + G("i,j"), D("i,j") );

            S.world().gop.fence(); // End of iteration work
            double iterf = madness::wall_time(); // End iteration timer

            ++scf_iter;
        }

        return energy + nuc_repl;
    }



} // namespace tests
} // namesapce mpqc




#endif /* mpqc_interfaces_tiledarray_tests_dmm_scf_hpp */
