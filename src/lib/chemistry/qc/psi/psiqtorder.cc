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

#ifdef __GNUC__
#pragma implementation
#endif


#include <chemistry/qc/psi/psiqtorder.h>
#include <math/scmat/local.h>
#include <valarray>
#include <numeric>
#include <cassert>

namespace sc {
  
  std::vector<unsigned int> index_map_symmtoqtorder(const std::vector<unsigned int> &frozen_docc,
                                                    const std::vector<unsigned int> &docc_act,
                                                    const std::vector<unsigned int> &socc_act,
                                                    const std::vector<unsigned int> &uocc_act,
                                                    const std::vector<unsigned int> &frozen_uocc) {
    
    const int nirrep = frozen_docc.size();
    std::vector<unsigned int> mos(nirrep);
    for(int i=0; i<nirrep; i++) {
      mos[i] = frozen_docc[i] + docc_act[i] + socc_act[i] + uocc_act[i] + frozen_uocc[i];
    }
    
#if 0
    ExEnv::out0() << "frozen_docc:" << std::endl;
    for(int i=0; i<frozen_docc.size(); i++) {
      ExEnv::out0() << std::setw(5) << frozen_docc[i];
    }
    ExEnv::out0() << std::endl;
    ExEnv::out0() << "docc_act:" << std::endl;
    for(int i=0; i<docc_act.size(); i++) {
      ExEnv::out0() << std::setw(5) << docc_act[i];
    }
    ExEnv::out0() << std::endl;
    ExEnv::out0() << "socc_act:" << std::endl;
    for(int i=0; i<socc_act.size(); i++) {
      ExEnv::out0() << std::setw(5) << socc_act[i];
    }
    ExEnv::out0() << std::endl;
    ExEnv::out0() << "uocc_act:" << std::endl;
    for(int i=0; i<uocc_act.size(); i++) {
      ExEnv::out0() << std::setw(5) << uocc_act[i];
    }
    ExEnv::out0() << std::endl;
    ExEnv::out0() << "frozen_uocc:" << std::endl;
    for(int i=0; i<frozen_uocc.size(); i++) {
      ExEnv::out0() << std::setw(5) << frozen_uocc[i];
    }
    ExEnv::out0() << std::endl;
    
    ExEnv::out0() << "All molecular orbitals:" << std::endl;
    for(int i=0; i<mos.size(); i++) {
      ExEnv::out0() << std::setw(5) << mos[i];
    }
    ExEnv::out0() << std::endl;
#endif
    
    // compute initial offsets
    unsigned int frozen_docc_offset = 0;
    unsigned int docc_act_offset = std::accumulate(frozen_docc.begin(),frozen_docc.end(),frozen_docc_offset);
    unsigned int socc_act_offset = std::accumulate(docc_act.begin(),docc_act.end(),docc_act_offset);
    unsigned int uocc_act_offset = std::accumulate(socc_act.begin(),socc_act.end(),socc_act_offset);
    unsigned int frozen_uocc_offset = std::accumulate(uocc_act.begin(),uocc_act.end(),uocc_act_offset);
    const unsigned int nmo = std::accumulate(mos.begin(), mos.end(), 0);
    
    // generate index_map mapping symmetry ordered data to QT ordered data
    std::vector<unsigned int> index_map(nmo);
    unsigned int ind = 0;
    for(int i=0; i<nirrep; i++) {
      
      // irrep i block of frozen_docc
      for(int j=0; j<frozen_docc[i]; j++) {
        index_map[ind] = frozen_docc_offset;
        
        frozen_docc_offset += 1;
        ind += 1;
      }
      
      // irrep i block of docc_act
      for(int j=0; j<docc_act[i]; j++) {
        index_map[ind] = docc_act_offset;
        
        docc_act_offset += 1;
        ind += 1;
      }
      
      // irrep i block of socc_act
      for(int j=0; j<socc_act[i]; j++) {
        index_map[ind] = socc_act_offset;
        
        socc_act_offset += 1;
        ind +=1;
      }
      
      // irrep i block of uocc_act
      for(int j=0; j<uocc_act[i]; j++) {
        index_map[ind] = uocc_act_offset;
        
        uocc_act_offset += 1;
        ind +=1;
      }
      
      // irrep i block of frozen_uocc
      for(int j=0; j<frozen_uocc[i]; j++) {
        index_map[ind] = frozen_uocc_offset;
        
        frozen_uocc_offset += 1;
        ind += 1;
      }
    }
   
#if 0
    ExEnv::out0() << "index_map from symmetric to QT order:" << std::endl;
    for(int i=0; i<index_map.size(); i++) {
      ExEnv::out0() << " " << index_map[i];
    }
    ExEnv::out0() << std::endl;
#endif
    
    return(index_map);
  }
  
  std::vector<unsigned int> index_map_symmtorasorder(const std::vector<unsigned int> &frozen_docc,
                                                     const std::vector<unsigned int> &ras1,
                                                     const std::vector<unsigned int> &ras2,
                                                     const std::vector<unsigned int> &ras3,
                                                     const std::vector<unsigned int> &frozen_uocc) {
    int nirrep = frozen_docc.size();
    std::vector<unsigned int> mos(nirrep);
    for(int i=0; i<nirrep; i++) {
      mos[i] = frozen_docc[i] + ras1[i] + ras2[i] + ras3[i] + frozen_uocc[i];
    }
    
    // compute initial offsets
    unsigned int frozen_docc_offset = 0;
    unsigned int ras1_offset = std::accumulate(frozen_docc.begin(),frozen_docc.end(),frozen_docc_offset);
    unsigned int ras2_offset = std::accumulate(ras1.begin(),ras1.end(),ras1_offset);
    unsigned int ras3_offset = std::accumulate(ras2.begin(),ras2.end(),ras2_offset);
    unsigned int frozen_uocc_offset = std::accumulate(ras3.begin(),ras3.end(),ras3_offset);
    const unsigned int nmo = std::accumulate(mos.begin(), mos.end(), 0);

    // generate index_map mapping symmetry ordered data to ras ordered data
    std::vector<unsigned int> index_map(nmo);
    unsigned int ind = 0;
    for(int i=0; i<nirrep; i++) {
      
      // irrep i block of frozen_docc
      for(int j=0; j<frozen_docc[i]; j++) {
        index_map[ind] = frozen_docc_offset;
        
        frozen_docc_offset += 1;
        ind += 1;
      }
      
      // irrep i block of ras1
      for(int j=0; j<ras1[i]; j++) {
        index_map[ind] = ras1_offset;
        
        ras1_offset += 1;
        ind += 1;
      }
      
      // irrep i block of ras2
      for(int j=0; j<ras2[i]; j++) {
        index_map[ind] = ras2_offset;
        
        ras2_offset += 1;
        ind +=1;
      }
      
      // irrep i block of ras3
      for(int j=0; j<ras3[i]; j++) {
        index_map[ind] = ras3_offset;
        
        ras3_offset += 1;
        ind +=1;
      }
      
      // irrep i block of frozen_uocc
      for(int j=0; j<frozen_uocc[i]; j++) {
        index_map[ind] = frozen_uocc_offset;
        
        frozen_uocc_offset += 1;
        ind += 1;
      }
    }
    
#if 0
    ExEnv::out0() << "index_map from symmetric to ras order:" << std::endl;
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
    assert(min_index == 0);
    const unsigned int max_index = * max_element(map.begin(), map.end());

    std::vector<unsigned int> imap(max_index+1, UINT_MAX);
    const unsigned int n = map.size();
    for(unsigned int i=0; i<n; ++i)
      imap[map[i]] = i;

    // make sure the map is isomorphic and can be inverted
    assert( *max_element(imap.begin(), imap.end()) != UINT_MAX);

    return imap;
  }

}
