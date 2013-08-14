/*
 * molecule_test.cpp
 *
 *  Created on: Aug 14, 2013
 *      Author: drewlewis
 */

#include "common.hpp"

using namespace mpqc;
using namespace mpqc::tests;

int main() {

    R<Molecule> mol_h2 = get_molecule("H2");
    R<Molecule> mol_h2o = get_molecule("H2O");
    R<Molecule> mol_empty = get_molecule();
    R<Molecule> mol_not_present = get_molecule("CH4");

    std::cout << "H2 = " << std::endl;
    mol_h2->print();
    std::cout << "\nH2O = " << std::endl;
    mol_h2o->print();
    std::cout << "\nempty = " << std::endl;
    mol_empty->print();
    std::cout << "\nnot in list = " << std::endl;
    mol_not_present->print();

    return 0;
}


