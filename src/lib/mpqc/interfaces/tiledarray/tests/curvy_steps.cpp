/*
 * dmm_scf.cpp
 *
 *  Created on: Aug 14, 2013
 *      Author: drewlewis
 */

#include "dmm_scf.hpp"
#include <fstream>
#include <iomanip>

using namespace mpqc;
using namespace mpqc::tests;

int main(int argc, char** argv){

    madness::World &world = madness::initialize(argc, argv);

    const std::string mol_name = "H2O";
    R<Molecule> mol = get_molecule(mol_name);

    R<Basis> basis = get_basis("STO-3G", mol);
    R<Basis> basis_df = get_basis("cc-pVDZ-RI", mol);

    R<Integral> int_fac = get_integral_factory(argc, argv);

    TA::Array<double, 2> S = get_overlap(world, mol, basis, int_fac);
    TA::Array<double, 2> H = get_hcore(world, mol, basis, int_fac);
    TA::Array<double, 3> Eri3 = get_eri3(world, mol, basis, basis_df, int_fac);
    TA::Array<double, 2> Eri2 = get_eri2(world, mol, basis_df, int_fac);
    TA::Array<double, 2> D = get_soad_guess(world, mol, basis, int_fac);
    world.gop.fence();

    TA::Array<double, 2> Inv_Eri2 = get_inverse(Eri2);
    world.gop.fence();

    // This is not a smart way to construct G
    TA::Array<double, 2> G = 2.0 *
            (Eri3("i,j,X") * Inv_Eri2("X,Y") * ((Eri3("n,m,Y") * D("m,n")))) -
            (Eri3("i,n,X") * Inv_Eri2("X,Y") * ((Eri3("j,m,Y") * D("m,n"))));

    TA::Array<double, 2> F = H("i,j") + G("i,j");
    double energy_guess = TA::expressions::dot(
           2.0 * H("i,j") + G("i,j"), D("i,j"));

    Eigen::MatrixXd Se = TA::array_to_eigen(S);
    Eigen::MatrixXd Fe = TA::array_to_eigen(F);
    Eigen::MatrixXd De = TA::array_to_eigen(D);
    std::cout << "S = \n" << Se << "\n" << "F = \n" << Fe 
        << "\nD = \n" << De << std::endl;

    double nuc_repl = mol->nuclear_repulsion_energy();
    double final_energy = DF_DMM(D, S, H, F, G, Eri3, Inv_Eri2, nuc_repl);
    world.gop.fence();

    std::cout << "Final Energy = " << final_energy << std::endl;
    world.gop.fence();

    return 0;
}


