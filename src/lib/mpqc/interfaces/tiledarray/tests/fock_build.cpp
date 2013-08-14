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

    madness::World &world = madness::initialize(argc, argv);

    const std::string mol_name = "H2O";
    R<Molecule> mol = get_molecule(mol_name);

    R<Basis> basis = get_basis("3-21G", mol);

    R<Integral> int_fac = get_integral_factory(argc, argv);

    TA::Array<double, 2> S = get_overlap(world, mol, basis, int_fac);
    TA::Array<double, 2> H = get_hcore(world, mol, basis, int_fac);
    TA::Array<double, 4> Eri = get_eri(world, mol, basis, int_fac);
    TA::Array<double, 2> D = get_soad_guess(world, mol, basis);
    world.gop.fence();

    TA::Array<double, 2> G = D("n,m") * (2.0 * Eri("i,j,m,n") - Eri("i,n,m,j"));
    TA::Array<double, 2> F = H("i,j") + G("i,j");
    double energy_guess = TA::expressions::dot(
                    2.0 * H("i,j") + G("i,j"), D("i,j"));

    world.gop.fence();

    if(mol_name == "H2"){
        std::cout << "D = \n" << D << std::endl;
        std::cout << "\nG = \n" << G << std::endl;
        std::cout << "\nF = \n" << F << std::endl;
    }

    std::cout << "\nEnergy guess = " << energy_guess << std::endl;

    return 0;
}
