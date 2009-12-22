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
  
  /** Returns map from symmetry-blocked orbitals to QT-ordered orbitals
   *  \param frozen_docc -- vector of number of frozen doubly occupied orbitals in each irrep
   *  \param docc_act -- vector of number of doubly occupied orbitals in each irrep
   *  \param socc_act -- vector of number of singly occupied orbitals in each irrep
   *  \param uocc_act -- vector of number of unoccupied orbitals in each irrep
   *  \param frozen_uocc -- vector of number of frozen unoccupied orbitals in each irrep */
  std::vector<unsigned int> index_map_symmtoqtorder(const std::vector<unsigned int> &frozen_docc,
                                                    const std::vector<unsigned int> &docc_act,
                                                    const std::vector<unsigned int> &socc_act,
                                                    const std::vector<unsigned int> &uocc_act,
                                                    const std::vector<unsigned int> &frozen_uocc);
  
  /** Returns map from symmetry-blocked orbitals to RAS-ordered orbitals
   *  \param frozen_docc -- vector of number of frozen doubly occupied orbitals in each irrep
   *  \param docc_act -- vector of number of doubly occupied orbitals in each irrep
   *  \param socc_act -- vector of number of singly occupied orbitals in each irrep
   *  \param uocc_act -- vector of number of unoccupied orbitals in each irrep
   *  \param frozen_uocc -- vector of number of frozen unoccupied orbitals in each irrep */
  std::vector<unsigned int> index_map_symmtorasorder(const std::vector<unsigned int> &frozen_docc,
                                                     const std::vector<unsigned int> &ras1,
                                                     const std::vector<unsigned int> &ras2,
                                                     const std::vector<unsigned int> &ras3,
                                                     const std::vector<unsigned int> &frozen_uocc);
  
  /// inverts an (isomorphic) index map that maps [0,n) onto itself
  std::vector<unsigned int> index_map_inverse(const std::vector<unsigned int>& map);

}

#endif /*PSIQTORDER_H_*/
