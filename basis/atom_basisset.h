#pragma once
#ifndef TILECLUSTERCHEM_BASIS_ATOM_BASIS_H
#define TILECLUSTERCHEM_BASIS_ATOM_BASIS_H

#include <cassert>
#include <vector>
#include <iosfwd>

namespace tcc {
namespace basis {

class AtomBasisShell {
  public:
    AtomBasisShell() = default;
    AtomBasisShell(AtomBasisShell const &) = default;
    AtomBasisShell &operator=(AtomBasisShell const &) = default;
    AtomBasisShell(AtomBasisShell &&) = default;
    AtomBasisShell &operator=(AtomBasisShell &&) = default;

    AtomBasisShell(int am, bool spherical, std::vector<double> exps,
                   std::vector<std::vector<double>> coeffs)
        : angular_momentum_(am),
          spherical_(spherical),
          exponents_(std::move(exps)),
          coeffs_(std::move(coeffs)) {}

    bool am_is_SP() const { return coeffs_.size() == 2; }
    int angular_momentum() const { return angular_momentum_; }
    bool spherical() const { return spherical_; }
    unsigned long contraction_length() const { return exponents_.size(); }

    double exponent(int i) const { return exponents_[i]; }
    std::vector<double> const &exponents() const { return exponents_; }
    double coeff(int i, int j = 0) const { return coeffs_[j][i]; }
    std::vector<std::vector<double>> const &coeffs() const { return coeffs_; }


  private:
    int angular_momentum_;
    bool spherical_;
    std::vector<double> exponents_;
    std::vector<std::vector<double>> coeffs_;
};

std::ostream &operator<<(std::ostream &os, AtomBasisShell const &abs);

class AtomBasisSet {
  public:
    AtomBasisSet() = default;
    AtomBasisSet(AtomBasisSet const &a) = default;
    AtomBasisSet(AtomBasisSet &&a) = default;
    AtomBasisSet &operator=(AtomBasisSet const &a) = default;
    AtomBasisSet &operator=(AtomBasisSet &&a) = default;

    AtomBasisSet(unsigned int atomic_number, std::vector<AtomBasisShell> shells)
        : atomic_number_(atomic_number), shells_(std::move(shells)) {}

    int atomic_number() const { return atomic_number_; }

    int ang_mo(int i) const { return shells_[i].angular_momentum(); }
    bool is_SP(int i) const { return shells_[i].am_is_SP(); }

    double exponent(int i, int j) const { return shells_[i].exponent(j); }

    double coeff(int i, int j, int k = 0) const {
        return shells_[i].coeff(j, k);
    }

    std::vector<AtomBasisShell> const &shells() const { return shells_; }
    AtomBasisShell const &shell(int i) const { return shells_[i]; }
    unsigned long nshells() const { return shells_.size(); }

  private:
    unsigned int atomic_number_;
    std::vector<AtomBasisShell> shells_;
};

std::ostream &operator<<(std::ostream &os, AtomBasisSet const &abs);

} // namespace basis
} // namespace tcc

#endif // TILECLUSTERCHEM_BASIS_ATOM_BASIS_H
