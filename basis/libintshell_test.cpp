#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <map>
#include <memory>
#include <sstream>
#include <algorithm>
#include "../include/libint.h"

struct atom {
    std::array<double, 3> pos_;
    short charge_;
    short mass_;
};

class basis_set {
  public:
    basis_set(std::string const &file_name) : basis_set_() {
        std::ifstream basis_file(file_name);
        std::string line;
        while (std::getline(basis_file, line) && line != "****") {
        } // Do nothing

        while (std::getline(basis_file, line)) {
            std::stringstream ss(line);
            std::string elem;
            ss >> elem;

            auto elem_pair = elements.find(elem);
            if (elem_pair != elements.end()) {
                basis_set_.push_back(
                    get_atom_basis(elem_pair->second, basis_file, line));
            } else if (line.size() != 0) {
                std::cout << "Could not locate element give by " << line
                          << std::endl;
            }
        }
    }

    std::vector<libint2::Shell> get_basis(atom const &a) const {
        std::vector<libint2::Shell> shells;
        auto basis = std::find_if(basis_set_.begin(), basis_set_.end(),
                                  [&](atom_basis const &ab) {
            return ab.atomic_number == a.charge_;
        });
        if (basis != basis_set_.end()) {
            for (auto shell : basis->functions) {
                libint2::Shell sh;
                sh.alpha.insert(sh.alpha.end(), shell.exp.begin(),
                                shell.exp.end());
                sh.contr.push_back(
                    {shell.ang_mo, false, std::vector<double>(shell.coeff)});
                sh.O = {{a.pos_[0], a.pos_[1], a.pos_[2]}};
                shells.push_back(sh);
            }
        }
        return shells;
    }

  private:
    struct atom_basis {
        struct shell {
            int ang_mo;
            std::vector<double> exp;
            std::vector<double> coeff;
        };

        int atomic_number;
        std::vector<shell> functions;
    };

    atom_basis
    get_atom_basis(int atomic_number, std::ifstream &is, std::string &line) {
        atom_basis a;
        a.atomic_number = atomic_number;
        while (std::getline(is, line) && line != "****") {
            std::stringstream ss(line);
            std::string ang_mo;
            auto contraction_length = 0;
            auto normalization = 0.0;

            ss >> ang_mo;
            ss >> contraction_length;
            ss >> normalization;

            atom_basis::shell sh;
            sh.ang_mo = ang_mo_map[ang_mo];
            for (auto i = 0; i < contraction_length; ++i) {
                std::getline(is, line);
                std::stringstream ss(line);
                double exp = 0;
                ss >> exp;
                double coeff = 0;
                ss >> coeff;
                sh.exp.push_back(exp);
                sh.coeff.push_back(coeff);
            }

            a.functions.push_back(std::move(sh));
            if (line.size() == 0) {
                std::getline(is, line); // not sure why I need to advance here.
            }
        }

        return a;
    }

    std::vector<atom_basis> basis_set_;

    std::map<std::string, int> elements{
        {"H", 1},
        {"He", 2},
        {"Li", 3},
        {"Be", 4},
        {"B", 5},
        {"C", 6},
        {"N", 7},
        {"O", 8},
        {"F", 9},
        {"Ne", 10},
        {"Na", 11},
        {"Mg", 12},
        {"Al", 13},
        {"Si", 14},
        {"P", 15},
        {"S", 16},
        {"Cl", 17},
        {"Ar", 18},
        {"K", 19},
        {"Ca", 20},
        {"Sc", 21},
        {"Ti", 22},
        {"V", 23},
        {"Cr", 24},
        {"Mn", 25},
        {"Fe", 26},
        {"Co", 27},
        {"Ni", 28},
        {"Cu", 29},
        {"Zn", 30},
        {"Ga", 31},
        {"Ge", 32},
        {"As", 33},
        {"Se", 34},
        {"Br", 35},
        {"Kr", 36},
    };
    std::map<std::string, int> ang_mo_map{
        {"S", 0}, {"P", 1}, {"D", 2}, {"F", 3}, {"G", 4}, {"H", 5}, {"I", 6}};
};

libint2::Shell hbasis(std::array<double, 3> pos) {
    return libint2::Shell{{3.425250910, 0.623913730, 0.168855400},
                          {{0, false, {0.15432897, 0.53532814, 0.44463454}}},
                          {{pos[0], pos[1], pos[2]}}};
}

std::vector<libint2::Shell> basis(std::vector<atom> const &atoms, basis_set const &bs) {
    std::vector<libint2::Shell> shells;
    shells.reserve(atoms.size());
    for (auto const &a : atoms) {
        auto atom_shells = bs.get_basis(a);
        shells.insert(shells.end(),atom_shells.begin(), atom_shells.end());
    }

    for(auto &s : shells){
      s.renorm();
    }

    return shells;
}

int main(int argc, char *argv[]) {
    std::string basis_file_name;
    if (argc >= 2) {
        basis_file_name = argv[1];
    }

    basis_set bs(basis_file_name);

    libint2::init();
    std::vector<atom> atoms;
    for (auto i = 0; i < 2; ++i) {
        atoms.emplace_back(atom{{{0.0, 0.0, double(i)}}, 1, 1});
    }

    auto shells = basis(atoms, bs);

    libint2::OneBodyEngine engine(libint2::OneBodyEngine::overlap,
                                  shells[0].nprim(), 3, 0);

    for (auto const &s1 : shells) {
        for (auto const &s2 : shells) {

            const auto buf = engine.compute(s1, s2);
            for (auto i = 0ul; i < (s1.size() * s2.size()); ++i) {
                       std::cout << buf[i] << "\t";
            }
        }
        std::cout << std::endl;
    }

    return 0;
}
