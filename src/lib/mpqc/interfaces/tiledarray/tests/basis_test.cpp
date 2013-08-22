/*
 * basis_test.cpp
 *
 *  Created on: Aug 14, 2013
 *      Author: drewlewis
 */

#include "common.hpp"

using namespace mpqc;
using namespace mpqc::tests;

int main(){
    R<Molecule> mol = get_molecule("H2O");

    // Get some basis sets.  Getting a regular basis set and a density fititng one
    // Next try with empty basis and see what happens.
    std::cout << "General initializing Basis set test" << std::endl;
    std::cout << "Getting Basis for Water" << std::endl;
    R<Basis> basis = get_basis("STO-3G", mol);
    R<Basis> basis_df = get_basis("cc-pVDZ-RI", mol);
    R<Basis> basis_empty = get_basis("", mol);

    // Print to check for correctness, this eventually needs to be an actual
    // test and not a print.
    std::cout << "Basis = " << std::endl;
    basis->print();
    std::cout << "\nBasis DF = " << std::endl;
    basis_df->print();
    std::cout << "\nBasis Empty = " << std::endl;
    basis_empty->print();

    return 0;
}


