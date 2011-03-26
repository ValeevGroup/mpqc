//
// psiqtorder.h
//
// Copyright (C) 2008 Martin Torheyden
//
// Author: Martin Torheyden <mtorhey@vt.edu>
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
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifdef __GNUC__
#pragma interface
#endif

#ifndef PSIQTORDER_H_
#define PSIQTORDER_H_

#include <vector>
#include <math/scmat/matrix.h>
#include <math/scmat/blocked.h>

namespace sc {
  
  /** Returns map from symmetry-blocked orbitals to correlated-order orbitals.
   *  The former blocks orbitals by irrep, and within each block orbitals are ordered by occupation or other correlation
   *  attribute (RAS, etc.), then by energy. The latter orders by the attribute, then by
   *  irrep, then by energy.
   *
   *  The common orderings used by Psi are:
   *  QT ("quantum trio") -- orbitals are classified as follows: frozen docc, docc, socc, uocc, frozen uocc
   *  RAS -- orbitals are classified as follows: frozen docc, RAS1, RAS2, RAS3, frozen uocc (not sure about RAS4)
   *  since both orderings define 5 classes of attributes, we only need 1 function
   *
   *  \param class1 -- specifies the number of orbitals of class 1 (e.g., frozen docc) in each irrep
   *  \param class2 -- specifies the number of orbitals of class 2 (e.g., RAS1) in each irrep
   *  \param class3 -- specifies the number of orbitals of class 3 (e.g., RAS2) in each irrep
   *  \param class4 -- specifies the number of orbitals of class 4 (e.g., RAS3) in each irrep
   *  \param class5 -- specifies the number of orbitals of class 5 (e.g., frozen uocc) in each irrep
   *  */
  std::vector<unsigned int> index_map_symmtocorrorder(const std::vector<unsigned int> &class1,
                                                      const std::vector<unsigned int> &class2,
                                                      const std::vector<unsigned int> &class3,
                                                      const std::vector<unsigned int> &class4,
                                                      const std::vector<unsigned int> &class5);

  /// inverts an (isomorphic) index map that maps [0,n) onto itself
  std::vector<unsigned int> index_map_inverse(const std::vector<unsigned int>& map);

}

#endif /*PSIQTORDER_H_*/
