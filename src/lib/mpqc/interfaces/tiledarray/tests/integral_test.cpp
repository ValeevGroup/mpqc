/*
 * integral_test.cpp
 *
 *  Created on: Aug 14, 2013
 *      Author: drewlewis
 */

#include "common.hpp"

using namespace mpqc;
using namespace mpqc::tests;

namespace TA = TiledArray;
int main(int argc, char** argv) {

    madness::World &world = madness::initialize(argc, argv);

    R<Molecule> mol = get_molecule("H2");

    R<Basis> basis = get_basis("3-21G", mol);
    R<Basis> basis_df = get_basis("cc-pVDZ", mol);

    // Make an integral factory that will generate engines
    R<Integral> int_fac = get_integral_factory(argc, argv);



    // Set basis for the overlap, hcore, and two electron integrals
    int_fac->set_basis(basis,basis,basis,basis);
    // Get the engine pools
    IntPool<R<sc::OneBodyInt> > overlap_pool(int_fac->overlap());
    IntPool<R<sc::OneBodyInt> > hcore_pool(int_fac->hcore());
    IntPool<R<sc::TwoBodyInt> > eri_pool(int_fac->electron_repulsion());

    // Set basis for density fiting three center integrals.
    int_fac->set_basis(basis, basis, basis_df, basis_df);
    // Get three center pool
    sc::Ref<sc::TwoBodyThreeCenterInt> eri3 = int_fac->electron_repulsion3();
    IntPool<R<sc::TwoBodyThreeCenterInt> > eri3_pool(eri3);

    // Set basis for density fitting two center integrals
    int_fac->set_basis(basis_df, basis_df, basis_df, basis_df);
    // Get two center engine pool
    IntPool<R<sc::TwoBodyTwoCenterInt> > eri2_pool(
                                         int_fac->electron_repulsion2());

    // Fill TiledArray's with data
    TA::Array<double, 2> S = Integrals(world, overlap_pool);
    TA::Array<double, 2> H = Integrals(world, hcore_pool);
    TA::Array<double, 4> Eri = Integrals(world, eri_pool);
    TA::Array<double, 3> Eri3 = Integrals(world, eri3_pool);
    TA::Array<double, 2> Eri2 = Integrals(world, eri2_pool);

    // Print TiledArrays
    std::cout << "S = \n" << S << std::endl;
    std::cout << "\nH = \n" << H << std::endl;
    std::cout << "\nEri = \n" << Eri << std::endl;
    std::cout << "\nEri3 = \n" << Eri3 << std::endl;
    std::cout << "\nEri2 = \n" << Eri2 << std::endl;

    world.gop.fence();
    madness::finalize();
    return 0;
}
