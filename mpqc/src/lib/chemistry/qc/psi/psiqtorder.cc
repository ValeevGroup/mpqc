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

namespace sc {
  
  std::vector<unsigned int> index_map_symmtompqcorder(const RefDiagSCMatrix &eigenvalues_symmorder) {
    int nmo = eigenvalues_symmorder.dim().n();
    std::vector<unsigned int> index_map(nmo);
    
    for(int i=0; i<nmo; i++) {
      index_map[i] = i;
    }
    
    for(int i=0; i<nmo; i++) {
      for(int j=(i+1); j<nmo; j++) {
        if(eigenvalues_symmorder.get_element(i) > eigenvalues_symmorder.get_element(j)) {
          double tmp = eigenvalues_symmorder.get_element(i);
          eigenvalues_symmorder.set_element(i,eigenvalues_symmorder.get_element(j));
          eigenvalues_symmorder.set_element(j,tmp);
          
          unsigned int tmpui = index_map[i];
          index_map[i] = index_map[j];
          index_map[j] = tmpui;
        }
      }
    }
    
    std::vector<unsigned int> index_map_inv(nmo);
    for(int i=0; i<nmo; i++) {
      index_map_inv[index_map[i]] = i;
    }
    
    return(index_map_inv);
  }
  
  std::vector<unsigned int> index_map_symmtoqtorder(int nmo,
                                                    const std::vector<unsigned int> &frozen_docc,
                                                    const std::vector<unsigned int> &docc_act,
                                                    const std::vector<unsigned int> &socc_act,
                                                    const std::vector<unsigned int> &uocc_act,
                                                    const std::vector<unsigned int> &frozen_uocc) {
    
    int nirrep = frozen_docc.size();
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
  
  std::vector<unsigned int> index_map_symmtorasorder(int nmo,
                                                     const std::vector<unsigned int> &frozen_docc,
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

  RefSCMatrix coeffsymmtoqtorder(const RefSCMatrix &coeffsymm,
                                 const std::vector<unsigned int> &frozen_docc,
                                 const std::vector<unsigned int> &docc_act,
                                 const std::vector<unsigned int> &socc_act,
                                 const std::vector<unsigned int> &uocc_act,
                                 const std::vector<unsigned int> &frozen_uocc) {
    
    Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
    int nmo = coeffsymm.coldim().n();
    RefSCDimension aodim = coeffsymm.rowdim();
    // TODO should use blocked dimensions here, but there are too many changes to be made downstream at the moment
#if 1
    RefSCDimension qtordermodim = new SCDimension(nmo); 
#else
    int* nfunc_per_block = new int[1];
    nfunc_per_block[0] = nmo;
    RefSCDimension qtordermodim = new SCDimension(nmo, 1, nfunc_per_block, "MOs in QT order");
    if (nmo)
      qtordermodim->blocks()->set_subdim(0, new SCDimension(nfunc_per_block[0]));
#endif
    RefSCMatrix coeffqtorder = localkit->matrix(aodim,qtordermodim);
    RefSCMatrix coeffsymm_nb = localkit->matrix(aodim,qtordermodim);
    
    //coeffsymm_nb.assign(coeffsymm);
    for(int i=0; i<nmo; i++) {
      for(int j=0; j<aodim.n(); j++) {
        coeffsymm_nb.set_element(i,j,coeffsymm.get_element(i,j));
      }
    }
    
    // generate index_map mapping symmetry ordered data to QT ordered data
    std::vector<unsigned int> index_map = index_map_symmtoqtorder(nmo,
                                                                  frozen_docc,
                                                                  docc_act,
                                                                  socc_act,
                                                                  uocc_act,
                                                                  frozen_uocc);
    
    // generate coeffqtorder from coeffsymm_nb by using  index_map
    coeffqtorder.assign(0.0);
    for(int i=0; i<nmo; i++) {
      coeffqtorder.assign_column(coeffsymm_nb.get_column(i),index_map[i]);
    }
    
    return(coeffqtorder);
  }
  
  RefSCMatrix coeffsymmtorasorder(const RefSCMatrix &coeffsymm,
                                  const std::vector<unsigned int> &frozen_docc,
                                  const std::vector<unsigned int> &ras1,
                                  const std::vector<unsigned int> &ras2,
                                  const std::vector<unsigned int> &ras3,
                                  const std::vector<unsigned int> &frozen_uocc) {
    Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
    int nmo = coeffsymm.coldim().n();
    RefSCDimension aodim = coeffsymm.rowdim();
    // TODO should use blocked dimensions here, but there are too many changes to be made downstream at the moment
#if 1
    RefSCDimension rasordermodim = new SCDimension(nmo); 
#else
    int* nfunc_per_block = new int[1];
    nfunc_per_block[0] = nmo;
    RefSCDimension rasordermodim = new SCDimension(nmo, 1, nfunc_per_block, "MOs in ras order");
    if (nmo)
      rasordermodim->blocks()->set_subdim(0, new SCDimension(nfunc_per_block[0]));
#endif
    RefSCMatrix coeffrasorder = localkit->matrix(aodim,rasordermodim);
    RefSCMatrix coeffsymm_nb = localkit->matrix(aodim,rasordermodim);
    
    //coeffsymm_nb.assign(coeffsymm);
    for(int i=0; i<nmo; i++) {
      for(int j=0; j<aodim.n(); j++) {
        coeffsymm_nb.set_element(i,j,coeffsymm.get_element(i,j));
      }
    }
    
    // generate index_map mapping symmetry ordered data to ras ordered data
    std::vector<unsigned int> index_map = index_map_symmtorasorder(nmo,
                                                                  frozen_docc,
                                                                  ras1,
                                                                  ras2,
                                                                  ras3,
                                                                  frozen_uocc);
   
    // generate coeffrasorder from coeffsymm_nb by using  index_map
    coeffrasorder.assign(0.0);
    for(int i=0; i<nmo; i++) {
      coeffrasorder.assign_column(coeffsymm_nb.get_column(i),index_map[i]);
    }
    
    return(coeffrasorder);
  }
  
}
