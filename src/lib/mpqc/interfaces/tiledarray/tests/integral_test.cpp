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

    R<Integral> int_fac = get_integral_factory(argc, argv);



    int_fac->set_basis(basis);
    IntPool<R<sc::OneBodyInt> > overlap_pool(int_fac->overlap());
    IntPool<R<sc::OneBodyInt> > hcore_pool(int_fac->hcore());
    IntPool<R<sc::TwoBodyInt> > eri_pool(int_fac->electron_repulsion());

    int_fac->set_basis(basis, basis, basis_df);
    IntPool<R<sc::TwoBodyThreeCenterInt> > eri3_pool(
                                           int_fac->electron_repulsion3());
    int_fac->set_basis(basis_df, basis_df);
    IntPool<R<sc::TwoBodyTwoCenterInt> > eri2_pool(
                                         int_fac->electron_repulsion2());
    TA::Array<double, 2> S = Integrals(world, overlap_pool,
                                       tiling::tile_by_atom);
    TA::Array<double, 2> H = Integrals(world, hcore_pool,
                                       tiling::tile_by_atom);
    TA::Array<double, 4> Eri = Integrals(world, eri_pool,
                                       tiling::tile_by_atom);
    TA::Array<double, 3> Eri3 = Integrals(world, eri3_pool,
                                          tiling::tile_by_atom);
    TA::Array<double, 2> Eri2 = Integrals(world, eri2_pool,
                                          tiling::tile_by_atom);
    std::cout << "S = \n" << S << std::endl;
    std::cout << "\nH = \n" << H << std::endl;
    std::cout << "\nEri = \n" << Eri << std::endl;
    std::cout << "\nEri3 = \n" << Eri3 << std::endl;
    std::cout << "\nEri2 = \n" << Eri2 << std::endl;

    world.gop.fence();
    madness::finalize();
    return 0;
}
