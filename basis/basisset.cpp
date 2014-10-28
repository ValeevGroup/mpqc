#include "basisset.h"
#include "basis_set_maps.h"

#include <iostream>
#include <fstream>
#include <sstream>

namespace tcc {
namespace basis {

void read_shell_info(std::ifstream &is, std::string &line,
                     std::vector<double> &exponents,
                     std::vector<std::vector<double>> &coeffs,
                     int contraction_length, std::string const &ang_mo) {


    std::vector<double> local_coeffs;
    std::vector<double> pshell_local_coeffs;

    for (auto i = 0; i < contraction_length; ++i) {
        std::getline(is, line);
        std::stringstream atom_ss(line);

        double exponent = 0;
        atom_ss >> exponent;
        exponents.emplace_back(exponent);

        double coeff = 0;
        atom_ss >> coeff;
        local_coeffs.emplace_back(coeff);

        if (ang_mo == "SP") {
            double p_coeff = 0;
            atom_ss >> p_coeff;
            pshell_local_coeffs.emplace_back(p_coeff);
        }
    }

    coeffs.emplace_back(std::move(local_coeffs));
    if(ang_mo == "SP"){
        coeffs.emplace_back(std::move(pshell_local_coeffs));
    }
}

AtomBasisSet read_atom_basis(std::ifstream &is, std::string &line) {

    // Get name of atom
    std::stringstream atom_info(line);
    std::string atom_name;
    atom_info >> atom_name;

    std::vector<AtomBasisShell> shells;

    // Start collecting atom info
    while (std::getline(is, line) && line != "****") {
        std::stringstream ss(line);

        std::string ang_mo;
        ss >> ang_mo;
        auto contraction_length = 0;
        ss >> contraction_length;
        auto normalization = 0.0;
        ss >> normalization;

        std::vector<double> shell_exponents;
        std::vector<std::vector<double>> shell_coeffs;
        read_shell_info(is, line, shell_exponents, shell_coeffs,
                        contraction_length, ang_mo);
        if (ang_mo == "SP") {
            ang_mo = "S";
        }

        shells.emplace_back(ang_mo_map[ang_mo], false,
                            std::move(shell_exponents),
                            std::move(shell_coeffs));
    }

    return AtomBasisSet{elements[atom_name], std::move(shells)};
}


void BasisSet::read_basis(std::string const &s) {

    std::ifstream basis_file(s);
    std::string line;

    // Skip header information
    while (std::getline(basis_file, line) && line != "****") {
        continue;
    }

    while (std::getline(basis_file, line)) {
        if (line != "") {
            atom_bases_.emplace_back(read_atom_basis(basis_file, line));
        }
    }
}

std::ostream &operator<<(std::ostream &os, BasisSet const &bs) {
    for (auto const &atom_bs : bs.basis()) {
        os << atom_bs << "\n";
    }
    return os;
}

} // namespace basis
} // namespace tcc
