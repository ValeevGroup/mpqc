#include "atom_basisset.h"
#include <iostream>

namespace tcc {
namespace basis {

std::ostream &operator<<(std::ostream &os, AtomBasisSet const &abs) {

    os << "Atomic Number " << abs.atomic_number() << "\n";
    auto nshells = abs.ang_mos().size();
    
    for(auto j = 0ul; j < nshells; ++j){

        std::cout << "\tShell " << j << " with angular momentum "
                  << abs.ang_mo(j) << "\n";

        auto contraction_length = abs.exponents(j).size();
        for (auto i = 0ul; i < contraction_length; ++i) {
            os << "\t\twith exponent " << abs.exponent(j,i) 
               << " having coefficent " << abs.coeff(j,i) << "\n";
        }

    }
    return os;
}

} // namespace basis
} // namespace tcc
