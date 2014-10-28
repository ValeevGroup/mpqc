#ifndef TILECLUSTERCHEM_BASIS_ATOM_BASIS_H
#define TILECLUSTERCHEM_BASIS_ATOM_BASIS_H

#include <cassert>
#include <vector>
#include <iosfwd>

namespace tcc {
namespace basis {

class AtomBasisSet {
  public:
    AtomBasisSet() = default;
    AtomBasisSet(AtomBasisSet const &a) = default;
    AtomBasisSet(AtomBasisSet &&a) = default;
    AtomBasisSet &operator=(AtomBasisSet const &a) = default;
    AtomBasisSet &operator=(AtomBasisSet &&a) = default;

    AtomBasisSet(unsigned int atomic_number, std::vector<int> ang_mo,
                 std::vector<std::vector<double>> exp,
                 std::vector<std::vector<double>> coeff)
        : atomic_number_(atomic_number),
          angular_momentums_(std::move(ang_mo)),
          exponents_(std::move(exp)),
          coeffs_(std::move(coeff)) {}

    int atomic_number() const { return atomic_number_; }

    int ang_mo(int i) const { return angular_momentums_.at(i); }
    std::vector<int> const &ang_mos() const { return angular_momentums_; }

    std::vector<std::vector<double>> const &exponents() const {
        return exponents_;
    }
    std::vector<double> const &exponents(int i) const {
        return exponents_.at(i);
    }
    double exponent(int i, int j) const { return exponents_.at(i).at(j); }

    std::vector<std::vector<double>> const &coeffs() const { return coeffs_; }
    std::vector<double> coeffs(int i) const { return coeffs_.at(i); }
    double coeff(int i, int j) const { return coeffs_.at(i).at(j); }

  private:
    unsigned int atomic_number_;
    std::vector<int> angular_momentums_;
    std::vector<std::vector<double>> exponents_;
    std::vector<std::vector<double>> coeffs_;
};

std::ostream &operator<<(std::ostream &os, AtomBasisSet const &abs);

} // namespace basis
} // namespace tcc

#endif // TILECLUSTERCHEM_BASIS_ATOM_BASIS_H
