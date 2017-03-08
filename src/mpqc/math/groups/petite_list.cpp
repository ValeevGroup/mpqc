/*
 *  This file is a part of Massively Parallel Quantum Chemistry package (v4).
 *  Copyright (C) 2017  Virginia Tech
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "mpqc/math/groups/petite_list.h"

#include "mpqc/util/misc/exception.h"

namespace mpqc {
namespace math {

std::shared_ptr<const PetiteList> PetiteList::make_trivial() {
  return make(Symmetry::e);
}

std::shared_ptr<const PetiteList> PetiteList::make(Symmetry symmetry) {
  switch (symmetry) {
    case Symmetry::e:
      return std::make_shared<const SymmPetiteList<Symmetry::e>>();
    case Symmetry::aa:
      return std::make_shared<const SymmPetiteList<Symmetry::aa>>();
    case Symmetry::aa_bb:
      return std::make_shared<const SymmPetiteList<Symmetry::aa_bb>>();
    case Symmetry::ab_ab:
      return std::make_shared<const SymmPetiteList<Symmetry::ab_ab>>();
    case Symmetry::aa_aa:
      return std::make_shared<const SymmPetiteList<Symmetry::aa_aa>>();
    default:
      throw ProgrammingError("unknown PetiteList::Symmetry", __FILE__,
                             __LINE__);
  }
}

// clang-format off
std::map<PetiteList::Symmetry, std::string> PetiteList::symmetry_to_string = {
    {PetiteList::Symmetry::e,           "e"},
    {PetiteList::Symmetry::aa,         "aa"},
    {PetiteList::Symmetry::aa_aa,   "aa_aa"},
    {PetiteList::Symmetry::aa_bb,   "aa_bb"},
    {PetiteList::Symmetry::ab_ab,   "ab_ab"}
};
// clang-format on

std::string to_string(PetiteList::Symmetry symmetry) {
  auto iter = PetiteList::symmetry_to_string.find(symmetry);
  if (iter != PetiteList::symmetry_to_string.end()) return iter->second;
  throw ProgrammingError("to_string: unknown PetiteList::Symmetry", __FILE__,
                         __LINE__);
}

}  // namespace math
}  // namespace mpqc
