#ifndef TILECLUSTERCHEM_BASIS_ATOM_BASIS_H
#define TILECLUSTERCHEM_BASIS_ATOM_BASIS_H

#include <cassert>
#include <vector>

class AtomBasisSet {
  public:
    AtomBasisSet() = default;
    AtomBasisSet(AtomBasisSet const &a) = default;
    AtomBasisSet(AtomBasisSet &&a) = default;
    AtomBasisSet operator=(AtomBasisSet const &a) = default;
    AtomBasisSet operator=(AtomBasisSet &&a) = default;

    AtomBasisSet(std::vector<int> ang_mo, std::vector<double> exp,
              std::vector<std::vector<double>> coeff)
        : angular_momentums_(std::move(ang_mo)), exponents_(std::move(exp)),
          coeffs_(std::move(coeff)) {
      assert(angular_momentums_.size() == coeffs_.size());
    }

    int ang_mo(int i = 0) const { return angular_momentums_.at(i); }
    std::vector<double> const & exponents() const {return exponents_;}
    std::vector<double> const & coeffs(int i = 0) const {return coeffs_.at(i);}
    int coeffs_per_exp() const {return angular_momentums_.size();}

  private:
    std::vector<int> angular_momentums_;
    std::vector<double> exponents_;
    std::vector<std::vector<double>> coeffs_;
};

#endif // TILECLUSTERCHEM_BASIS_ATOM_BASIS_H
