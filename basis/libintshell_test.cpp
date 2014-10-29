#include "../include/libint.h"
#include "basisset.h"
#include "atom_basisset.h"

#include <iostream>
#include <vector>

using namespace tcc;
using namespace tcc::basis;

int main(int argc, char *argv[]) {
    std::string basis_file_name;
    if (argc >= 2) {
        basis_file_name = argv[1];
    }

    tcc::basis::BasisSet bs(basis_file_name);
    std::cout << bs << std::endl;

    return 0;
}
