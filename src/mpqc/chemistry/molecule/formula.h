//
// formula.h
//
// Copyright (C) 2016 Virginia Tech
//
// Author: Ed Valeyev <eduard@valeyev.net>
// Maintainer: EV
//
// This file is part of the MPQC Framework. It is distributed under the terms of
// the GNU General Public License, version 3.

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_MOLECULE_FORMULA_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_MOLECULE_FORMULA_H_

#include <string>
#include <utility>
#include <vector>

#include "./molecule.h"

namespace mpqc {

/// @addtogroup ChemistryMolecule
/// @{

/** The MolecularFormula class is used to calculate the molecular
 formula of a Molecule and return it as a string. The format is "CxHy<rest of
 elements>" where
 the "rest" of elements is given in the order of increasing atomic number.
 */
class MolecularFormula {
 private:
  using formula_type = std::vector<std::pair<int, size_t>>;
  formula_type formula_;
  std::string formula_str_;

  void compute(const Molecule& m);

 public:
  explicit MolecularFormula(const Molecule& m);
  ~MolecularFormula() { }

  /// computes and returns molecular formula
  /// @return std::string containing the molecular formula.
  const std::string& string() const { return formula_str_; }
  /// returns the formula as a vector of pairs {atomic number of the element,
  /// number of occurrences}
  const formula_type& formula() const { return formula_; }
};

/// @}
// end of addtogroup ChemistryMolecule

}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_MOLECULE_FORMULA_H_
