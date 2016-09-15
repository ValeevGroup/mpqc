#include <iostream>
#include <mpqc/chemistry/qc/basis/atom_basisset.h>

namespace mpqc {
namespace basis {

std::ostream &operator<<(std::ostream &os, AtomBasisShell const &abs) {
    bool is_spher = abs.spherical();
    if (is_spher) {
        os << "\tBasis is spherical ";
    } else {
        os << "\tBasis is not spherical ";
    }

    bool is_SP = abs.am_is_SP();
    if (is_SP) {
        os << "Angular Momentum is SP\n";
    } else {
        os << "Angular Momentum is " << abs.angular_momentum() << "\n";
    }

    for (auto i = 0ul; i < abs.contraction_length(); ++i) {
        os << "\t\tExponent " << i << " is " << abs.exponent(i) << " ";
        auto const &coeffs = abs.coeffs();
        if (is_SP) {
            os << "With coefficents ";
            for(auto const &elem : coeffs[i]){
                std::cout << elem << " ";
            }
            std::cout << "\n";
        } else {
            os << "With coefficent " << coeffs[i][0] << "\n";
        }
    }
    return os;
}

std::ostream &operator<<(std::ostream &os, AtomBasisSet const &abs) {

    os << "Atomic Number " << abs.atomic_number() << "\n";
    for (auto const &shell : abs.shells()) {
        os << shell;
    }
    return os;
}

} // namespace basis
} // namespace mpqc
