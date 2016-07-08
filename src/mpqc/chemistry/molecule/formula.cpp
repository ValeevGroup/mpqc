//
// formula.cc --- implementation of the MolecularFormula class
//
// Copyright (C) 1997 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#include "./formula.h"

#include <map>
#include <string>
#include <sstream>

#include <libint2/chemistry/elements.h>

namespace mpqc {
namespace molecule{

MolecularFormula::MolecularFormula(const mpqc::molecule::Molecule& mol) {
  compute(mol);
}

void MolecularFormula::compute(const mpqc::molecule::Molecule& mol) {
  std::map<int, size_t> count; // maps atomic number -> number of atoms

  auto atoms = mol.atoms();

  for (const auto& atom: atoms) {
    auto Z = atom.charge();
    if (count.find(Z) == count.end()) count[Z] = 0;
    ++count[Z];
  }

  formula_.reserve(count.size());

  std::ostringstream sstr;
  // pattern: CxHy... where "..." are the rest of the elements in the order of increasing atomic number
  if (count.find(6) != count.end()) {
    sstr << "C";
    if (count[6] > 1) sstr << count[6];
    formula_.emplace_back(std::make_pair(6,count[6]));
    count.erase(6);
  }
  if (count.find(1) != count.end()) {
    sstr << "H";
    if (count[1] > 1) sstr << count[1];
    formula_.emplace_back(std::make_pair(1,count[1]));
    count.erase(1);
  }
  for (const auto& it: count) {
    auto Z = it.first;
    auto c = it.second;
    formula_.emplace_back(std::make_pair(Z,c));
    assert(Z != 0);
    sstr << libint2::chemistry::element_info[Z-1].symbol;
    if (c > 1) sstr << c;
  }
  formula_str_ = sstr.str();
}

}}  // namespace mpqc::molecule
