//
// psiqtorder.cc
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


#include <chemistry/qc/psi/psiqtorder.h>
#include <math/scmat/local.h>
#include <valarray>
#include <numeric>
#include <cassert>

namespace sc {
  
  std::vector<unsigned int> index_map_symmtocorrorder(const std::vector<unsigned int> &class1,
                                                    const std::vector<unsigned int> &class2,
                                                    const std::vector<unsigned int> &class3,
                                                    const std::vector<unsigned int> &class4,
                                                    const std::vector<unsigned int> &class5) {
    
    const int nirrep = class1.size();
    std::vector<unsigned int> mos(nirrep);
    for(int i=0; i<nirrep; i++) {
      mos[i] = class1[i] + class2[i] + class3[i] + class4[i] + class5[i];
    }
    
#if 0
    ExEnv::out0() << "frozen_docc:" << std::endl;
    for(int i=0; i<class1.size(); i++) {
      ExEnv::out0() << std::setw(5) << class1[i];
    }
    ExEnv::out0() << std::endl;
    ExEnv::out0() << "docc_act:" << std::endl;
    for(int i=0; i<class2.size(); i++) {
      ExEnv::out0() << std::setw(5) << class2[i];
    }
    ExEnv::out0() << std::endl;
    ExEnv::out0() << "socc_act:" << std::endl;
    for(int i=0; i<class3.size(); i++) {
      ExEnv::out0() << std::setw(5) << class3[i];
    }
    ExEnv::out0() << std::endl;
    ExEnv::out0() << "uocc_act:" << std::endl;
    for(int i=0; i<class4.size(); i++) {
      ExEnv::out0() << std::setw(5) << class4[i];
    }
    ExEnv::out0() << std::endl;
    ExEnv::out0() << "frozen_uocc:" << std::endl;
    for(int i=0; i<class5.size(); i++) {
      ExEnv::out0() << std::setw(5) << class5[i];
    }
    ExEnv::out0() << std::endl;
    
    ExEnv::out0() << "All molecular orbitals:" << std::endl;
    for(int i=0; i<mos.size(); i++) {
      ExEnv::out0() << std::setw(5) << mos[i];
    }
    ExEnv::out0() << std::endl;
#endif
    
    // compute initial offsets
    unsigned int class1_offset = 0;
    unsigned int class2_offset = std::accumulate(class1.begin(),class1.end(),class1_offset);
    unsigned int class3_offset = std::accumulate(class2.begin(),class2.end(),class2_offset);
    unsigned int class4_offset = std::accumulate(class3.begin(),class3.end(),class3_offset);
    unsigned int class5_offset = std::accumulate(class4.begin(),class4.end(),class4_offset);
    const unsigned int nmo = std::accumulate(mos.begin(), mos.end(), 0);
    
    // generate index_map mapping symmetry ordered data to QT ordered data
    std::vector<unsigned int> index_map(nmo);
    unsigned int ind = 0;
    for(int i=0; i<nirrep; i++) {
      
      // irrep i block of frozen_docc
      for(int j=0; j<class1[i]; j++) {
        index_map[ind] = class1_offset;
        
        class1_offset += 1;
        ind += 1;
      }
      
      // irrep i block of docc_act
      for(int j=0; j<class2[i]; j++) {
        index_map[ind] = class2_offset;
        
        class2_offset += 1;
        ind += 1;
      }
      
      // irrep i block of socc_act
      for(int j=0; j<class3[i]; j++) {
        index_map[ind] = class3_offset;
        
        class3_offset += 1;
        ind +=1;
      }
      
      // irrep i block of uocc_act
      for(int j=0; j<class4[i]; j++) {
        index_map[ind] = class4_offset;
        
        class4_offset += 1;
        ind +=1;
      }
      
      // irrep i block of frozen_uocc
      for(int j=0; j<class5[i]; j++) {
        index_map[ind] = class5_offset;
        
        class5_offset += 1;
        ind += 1;
      }
    }
   
#if 0
    ExEnv::out0() << "index_map from symmetric to correlated order:" << std::endl;
    for(int i=0; i<index_map.size(); i++) {
      ExEnv::out0() << " " << index_map[i];
    }
    ExEnv::out0() << std::endl;
#endif
    
    return(index_map);
  }
  
  std::vector<unsigned int> index_map_inverse(const std::vector<unsigned int>& map) {
    typedef std::vector<unsigned int>::iterator iter;
    const unsigned int min_index = * min_element(map.begin(), map.end());
    MPQC_ASSERT(min_index == 0);
    const unsigned int max_index = * max_element(map.begin(), map.end());

    std::vector<unsigned int> imap(max_index+1, UINT_MAX);
    const unsigned int n = map.size();
    for(unsigned int i=0; i<n; ++i)
      imap[map[i]] = i;

    // make sure the map is isomorphic and can be inverted
    MPQC_ASSERT( *max_element(imap.begin(), imap.end()) != UINT_MAX);

    return imap;
  }

}
