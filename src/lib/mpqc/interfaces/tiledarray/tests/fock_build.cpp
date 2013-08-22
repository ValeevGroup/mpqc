/*
 * fock_build.cpp
 *
 *  Created on: Aug 14, 2013
 *      Author: drewlewis
 */

#include "tiledarray_fock.hpp"

using namespace mpqc;
using namespace mpqc::tests;

int main(int argc, char** argv){

    // Get madness world, not the preferred way of doing this in MPQC.  Ok in
    // stand alone programs.
    madness::World &world = madness::initialize(argc, argv);

    // Get molecule
    const std::string mol_name = "H2O";
    R<Molecule> mol = get_molecule(mol_name);

    // Get basis set
    R<Basis> basis = get_basis("3-21G", mol);

    // Get an integral factory for generating integral engine.
    R<Integral> int_fac = get_integral_factory(argc, argv);

    // Get the different types of integrals needed for calculation and
    // also get the density matrix.
    TA::Array<double, 2> S = get_overlap(world, mol, basis, int_fac);
    TA::Array<double, 2> H = get_hcore(world, mol, basis, int_fac);
    TA::Array<double, 4> Eri = get_eri(world, mol, basis, int_fac);
    TA::Array<double, 2> D = get_soad_guess(world, mol, basis);
    world.gop.fence(); // Fence to make sure that all of the Array assignments finish

    // Compute G and F contractions.  This should be essentiually like Szabo
    TA::Array<double, 2> G = D("n,m") * (2.0 * Eri("i,j,m,n") - Eri("i,n,m,j"));
    TA::Array<double, 2> F = H("i,j") + G("i,j");
    // Compute Energy = \sum_{ij} (2.0 * H_{ij} + G_{ij}) * D_{ij}
    double energy = TA::expressions::dot( 2.0 * H("i,j") + G("i,j"), D("i,j"));

    world.gop.fence();

    // If the molecule is small enough print things so we can see them.
    if(mol_name == "H2"){
        std::cout << "D = \n" << D << std::endl;
        std::cout << "\nG = \n" << G << std::endl;
        std::cout << "\nF = \n" << F << std::endl;
    }

    std::cout << "\nEnergy guess = " << energy << std::endl;

    return 0;
}
