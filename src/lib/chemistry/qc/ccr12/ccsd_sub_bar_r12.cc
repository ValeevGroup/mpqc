//
// ccsd_sub_bar_r12.cc
//
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@theochem.uni-stuttgart.de>
// Maintainer: TS
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


#include <algorithm>
#include <chemistry/qc/ccr12/ccsd_sub_bar_r12.h>
#include <chemistry/qc/ccr12/tensor.h>
#include <chemistry/qc/ccr12/mtensor.h>

using namespace sc;
using namespace std;
  
  
void CCSD_Sub_Bar_R12::compute_amp(){ //k_i0_offset,z->in.at(0),z->t2()=>z->vd2()
  smith_0_1(tildeV_); //z->vd2()=>out
  smith_0_2(tildeV_); //z->t2(),z->vd2()=>out
}
  
void CCSD_Sub_Bar_R12::smith_0_1(Ref<Tensor>& out){ 
      
  const size_t singles = z->maxtilesize() * z->maxtilesize();
  const size_t doubles = singles * singles;
  double* k_a0      = z->mem()->malloc_local_double(doubles); 
  for (long h3b=0L;h3b<z->noab();++h3b) { 
   for (long h4b=h3b;h4b<z->noab();++h4b) { 
    for (long h1b=0L;h1b<z->noab();++h1b) { 
     for (long h2b=h1b;h2b<z->noab();++h2b) { 
      long tileoffset; 
      tileoffset=(h2b+z->noab()*(h1b+z->noab()*(h4b+z->noab()*(h3b)))); 
      if (out->is_this_local(tileoffset)) { 
       if (!z->restricted() || z->get_spin(h3b)+z->get_spin(h4b)+z->get_spin(h1b)+z->get_spin(h2b)!=8L) { 
        if (z->get_spin(h3b)+z->get_spin(h4b)==z->get_spin(h1b)+z->get_spin(h2b)) { 
         if ((z->get_sym(h3b)^(z->get_sym(h4b)^(z->get_sym(h1b)^z->get_sym(h2b))))==z->irrep_e()) { 
          long h3b_0,h4b_0,h1b_0,h2b_0; 
          z->restricted_4(h3b,h4b,h1b,h2b,h3b_0,h4b_0,h1b_0,h2b_0); 
          z->vd2()->get_block(h2b_0+(z->nab())*(h1b_0+(z->nab())*(h4b_0+z->noab()*(h3b_0))),k_a0); 
          out->add_block(h2b+z->noab()*(h1b+z->noab()*(h4b+z->noab()*(h3b))),k_a0); 
         } 
        } 
       } 
      } 
     } 
    } 
   } 
  } 
  z->mem()->free_local_double(k_a0); 
  z->mem()->sync(); 
} 
  
void CCSD_Sub_Bar_R12::smith_0_2(Ref<Tensor>& out){ 
      
  const size_t singles = z->maxtilesize() * z->maxtilesize();
  const size_t doubles = singles * singles;
  double* k_a0      = z->mem()->malloc_local_double(doubles); 
  double* k_a0_sort = z->mem()->malloc_local_double(doubles); 
  double* k_a1      = z->mem()->malloc_local_double(doubles); 
  double* k_c_sort  = z->mem()->malloc_local_double(doubles); 
  for (long h3b=0L;h3b<z->noab();++h3b) { 
   for (long h4b=h3b;h4b<z->noab();++h4b) { 
    for (long h1b=0L;h1b<z->noab();++h1b) { 
     for (long h2b=h1b;h2b<z->noab();++h2b) { 
      long tileoffset; 
      tileoffset=(h2b+z->noab()*(h1b+z->noab()*(h4b+z->noab()*(h3b)))); 
      if (out->is_this_local(tileoffset)) { 
       if (!z->restricted() || z->get_spin(h3b)+z->get_spin(h4b)+z->get_spin(h1b)+z->get_spin(h2b)!=8L) { 
        if (z->get_spin(h3b)+z->get_spin(h4b)==z->get_spin(h1b)+z->get_spin(h2b)) { 
         if ((z->get_sym(h3b)^(z->get_sym(h4b)^(z->get_sym(h1b)^z->get_sym(h2b))))==(z->irrep_t()^z->irrep_e())) { 
          long dimc=z->get_range(h3b)*z->get_range(h4b)*z->get_range(h1b)*z->get_range(h2b); 
          std::fill(k_c_sort,k_c_sort+dimc,0.0); 
          for (long p5b=z->noab();p5b<z->noab()+z->nvab();++p5b) { 
           for (long p6b=p5b;p6b<z->noab()+z->nvab();++p6b) { 
            if (z->get_spin(p5b)+z->get_spin(p6b)==z->get_spin(h1b)+z->get_spin(h2b)) { 
             if ((z->get_sym(p5b)^(z->get_sym(p6b)^(z->get_sym(h1b)^z->get_sym(h2b))))==z->irrep_t()) { 
              long p5b_0,p6b_0,h1b_0,h2b_0; 
              z->restricted_4(p5b,p6b,h1b,h2b,p5b_0,p6b_0,h1b_0,h2b_0); 
              long h3b_1,h4b_1,p5b_1,p6b_1; 
              z->restricted_4(h3b,h4b,p5b,p6b,h3b_1,h4b_1,p5b_1,p6b_1); 
              long dim_common=z->get_range(p5b)*z->get_range(p6b); 
              long dima0_sort=z->get_range(h1b)*z->get_range(h2b); 
              long dima1_sort=z->get_range(h3b)*z->get_range(h4b); 

              z->t2()->get_block(h2b_0+z->noab()*(h1b_0+z->noab()*(p6b_0-z->noab()+z->nvab()*(p5b_0-z->noab()))),k_a0); 
              z->sort_indices4(k_a0,k_a0_sort,z->get_range(p5b),z->get_range(p6b),z->get_range(h1b),z->get_range(h2b),2,3,0,1,+1.0); 

              z->vd2()->get_block(p6b_1+(z->nab())*(p5b_1+(z->nab())*(h4b_1+z->noab()*(h3b_1))),k_a1); 
              double factor=1.0; 
              if (p5b==p6b) factor*=0.5; 
              z->smith_dgemm(dima0_sort,dima1_sort,dim_common,factor,k_a0_sort,dim_common,k_a1,dim_common,1.0,k_c_sort,dima0_sort); 
             } 
            } 
           } 
          } 
          out->add_block(h2b+z->noab()*(h1b+z->noab()*(h4b+z->noab()*(h3b))),k_c_sort); 
         } 
        } 
       } 
      } 
     } 
    } 
   } 
  } 
  z->mem()->free_local_double(k_a1); 
  z->mem()->free_local_double(k_a0_sort); 
  z->mem()->free_local_double(k_a0); 
  z->mem()->free_local_double(k_c_sort); 
  z->mem()->sync(); 
}


