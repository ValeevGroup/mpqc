#include "basis_set.h"
#include "basis_set_maps.h"
#include "atom_basisset.h"
#include "cluster_shells.h"
#include "../molecule/cluster_collapse.h"
#include "../molecule/cluster.h"
#include "../molecule/atom.h"
#include "../include/libint.h"

#include <iostream>
#include <memory>
#include <fstream>
#include <sstream>
#include <algorithm>

namespace tcc {
namespace basis {

void read_shell_info(std::ifstream &is, std::string &line,
                     std::vector<double> &exponents,
                     std::vector<std::vector<double>> &coeffs,
                     int contraction_length, std::string const &ang_mo);

AtomBasisSet read_atom_basis(std::ifstream &is, std::string &line);


BasisSet::BasisSet() = default;
BasisSet::BasisSet(BasisSet const &) = default;
BasisSet::BasisSet(BasisSet &&) = default;
BasisSet &BasisSet::operator=(BasisSet const &) = default;
BasisSet &BasisSet::operator=(BasisSet &&) = default;

BasisSet::BasisSet(std::string const &s) : atom_bases_() { read_basis(s); }

std::vector<libint2::Shell>
BasisSet::atom_basis(molecule::Atom const &a) const {

    auto a_basis = std::find_if(atom_bases_.begin(), atom_bases_.end(),
                                [&](AtomBasisSet const &abs) {
        return abs.atomic_number() == a.charge();
    });

    std::vector<libint2::Shell> atom_shells;
    for (auto const &s : a_basis->shells()) {
        std::vector<libint2::Shell::Contraction> cntrs;
        if (!s.am_is_SP()) {
            cntrs.emplace_back(libint2::Shell::Contraction{
                s.angular_momentum(), s.spherical(), s.coeffs()[0]});
        } else {
            cntrs.emplace_back(
                libint2::Shell::Contraction{0, s.spherical(), s.coeffs()[0]});
            cntrs.emplace_back(
                libint2::Shell::Contraction{1, s.spherical(), s.coeffs()[1]});
        }

        libint2::Shell sh{s.exponents(),
                          std::move(cntrs),
                          {{a.center()[0], a.center()[1], a.center()[2]}}};

        atom_shells.push_back(std::move(sh));
    }

    return atom_shells;
}

std::vector<ClusterShells> BasisSet::create_basis(
    std::vector<std::shared_ptr<molecule::Cluster>> const &clusters) const {
    std::vector<ClusterShells> cs;

    auto atom_bs = atom_basis_set();

    for (auto const &c : clusters) {

        std::vector<molecule::Atom> atoms = molecule::collapse_to_atoms(*c);

        unsigned int max_am = 0;
        for (auto const &atom : atoms) {
            auto abs = std::find_if(atom_bs.begin(), atom_bs.end(),
                                    [&](AtomBasisSet const &abs) {
                return abs.atomic_number() == atom.charge();
            });
            max_am = std::max(max_am, abs->max_am());
        }
        const auto n_am = max_am + 1;

        std::vector<std::vector<libint2::Shell>> binned_shells(n_am);
        for (auto const &atom : atoms) {
            auto shells = atom_basis(atom);

            for (auto shell : shells) {
                shell.renorm(); // takes care of normalization
                const auto ang_mo = shell.contr[0].l;
                binned_shells[ang_mo].push_back(std::move(shell));
            }
        }
        cs.emplace_back(std::move(binned_shells), c);
    }

    return cs;
}

std::vector<AtomBasisSet> const &BasisSet::atom_basis_set() const {
    return atom_bases_;
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
    for (auto const &atom_bs : bs.atom_basis_set()) {
        os << atom_bs << "\n";
    }
    return os;
}

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
    if (ang_mo == "SP") {
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

} // namespace basis
} // namespace tcc
