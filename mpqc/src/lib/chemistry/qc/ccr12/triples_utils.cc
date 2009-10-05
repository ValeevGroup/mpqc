//
// triples_utils.cc
//
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@qtp.ufl.edu>
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
#include <chemistry/qc/ccr12/ccr12_triples.h> 
#include <chemistry/qc/ccr12/tensor.h>
using namespace sc;

void CCR12_Triples::offset_hhphhh(Ref<Tensor>& t) {
  long size=0L;
  for (long h4b=0L;h4b<z->noab();++h4b) { 
   for (long h5b=h4b;h5b<z->noab();++h5b) { 
    for (long p6b=z->noab();p6b<z->noab()+z->nvab();++p6b) { 
     for (long h1b=0L;h1b<z->noab();++h1b) { 
      for (long h2b=h1b;h2b<z->noab();++h2b) { 
       for (long h3b=h2b;h3b<z->noab();++h3b) { 
        if (!z->restricted() || z->get_spin(h4b)+z->get_spin(h5b)+z->get_spin(p6b)+z->get_spin(h1b)+z->get_spin(h2b)+z->get_spin(h3b)!=12L) { 
         if (z->get_spin(h4b)+z->get_spin(h5b)+z->get_spin(p6b)==z->get_spin(h1b)+z->get_spin(h2b)+z->get_spin(h3b)) { 
          if ((z->get_sym(h4b)^(z->get_sym(h5b)^(z->get_sym(p6b)^(z->get_sym(h1b)^(z->get_sym(h2b)^z->get_sym(h3b))))))==z->irrep_t()) { 
           t->input_offset(h3b+z->noab()*(h2b+z->noab()*(h1b+z->noab()*(p6b-z->noab()+z->nvab()*(h5b+z->noab()*h4b)))),size);
           size+=z->get_range(h1b)*z->get_range(h2b)*z->get_range(h3b)*z->get_range(h4b)*z->get_range(h5b)*z->get_range(p6b);
          }
         }
        }
       }
      }
     }
    }
   }
  }
  t->set_filesize(size);
  t->createfile();
}


void CCR12_Triples::denom_contraction(){ 

Ref<Tensor> denom = new Tensor("Denom", z->mem());
z->offset_gt2(denom, false);
      
// CAUTION : the loops are reordered to minimize the inversion (which is allowed as we don't have a sequential data read)
for (long p6b=z->noab();p6b<z->noab()+z->nvab();++p6b) { 
 for (long h1b=0L;h1b<z->noab();++h1b) { 
  for (long h2b=h1b;h2b<z->noab();++h2b) { 
   for (long h3b=h2b;h3b<z->noab();++h3b) { 

    // Create denominator here!!!!

    for (long h4b=0L;h4b<z->noab();++h4b) { 
     for (long h5b=h4b;h5b<z->noab();++h5b) { 

      long tileoffset; 
      tileoffset=(h3b+z->noab()*(h2b+z->noab()*(h1b+z->noab()*(p6b-z->noab()+z->nvab()*(h5b+z->noab()*(h4b)))))); 
      if (intermediate_->is_this_local(tileoffset)) { 
       if (!z->restricted() || z->get_spin(h4b)+z->get_spin(h5b)+z->get_spin(p6b)+z->get_spin(h1b)+z->get_spin(h2b)+z->get_spin(h3b)!=12L) { 
        if (z->get_spin(h4b)+z->get_spin(h5b)+z->get_spin(p6b)==z->get_spin(h1b)+z->get_spin(h2b)+z->get_spin(h3b)) { 
         if ((z->get_sym(h4b)^(z->get_sym(h5b)^(z->get_sym(p6b)^(z->get_sym(h1b)^(z->get_sym(h2b)^z->get_sym(h3b))))))==(z->irrep_t()^z->irrep_t())) { 
          long dimc=z->get_range(h4b)*z->get_range(h5b)*z->get_range(p6b)*z->get_range(h1b)*z->get_range(h2b)*z->get_range(h3b); 
          double* k_c_sort=z->mem()->malloc_local_double(dimc); 
          std::fill(k_c_sort,k_c_sort+(size_t)dimc,0.0); 
          for (long h7b=0L;h7b<z->noab();++h7b) { 
           for (long h8b=h7b;h8b<z->noab();++h8b) { 
            if (z->get_spin(h4b)+z->get_spin(h5b)==z->get_spin(h7b)+z->get_spin(h8b)) { 
             if ((z->get_sym(h4b)^(z->get_sym(h5b)^(z->get_sym(h7b)^z->get_sym(h8b))))==z->irrep_t()) { 
              long h4b_0,h5b_0,h7b_0,h8b_0; 
              z->restricted_4(h4b,h5b,h7b,h8b,h4b_0,h5b_0,h7b_0,h8b_0); 
              long h7b_1,h8b_1,p6b_1,h1b_1,h2b_1,h3b_1; 
              z->restricted_6(h7b,h8b,p6b,h1b,h2b,h3b,h7b_1,h8b_1,p6b_1,h1b_1,h2b_1,h3b_1); 
              long dim_common=z->get_range(h7b)*z->get_range(h8b); 
              long dima0_sort=z->get_range(h4b)*z->get_range(h5b); 
              long dima0=dim_common*dima0_sort; 
              long dima1_sort=z->get_range(p6b)*z->get_range(h1b)*z->get_range(h2b)*z->get_range(h3b); 
              long dima1=dim_common*dima1_sort; 
              if (dima0>0L && dima1>0L) { 
               double* k_a0_sort=z->mem()->malloc_local_double(dima0); 
               double* k_a0=z->mem()->malloc_local_double(dima0); 

               // read denom tensor
               denom->get_block(h8b_0+z->noab()*(h7b_0+z->noab()*(h5b_0+z->noab()*(h4b_0))),k_a0); 
               z->sort_indices4(k_a0,k_a0_sort,z->get_range(h4b),z->get_range(h5b),z->get_range(h7b),z->get_range(h8b),1,0,3,2,+1.0); 
               z->mem()->free_local_double(k_a0); 
               double* k_a1_sort=z->mem()->malloc_local_double(dima1); 
               double* k_a1=z->mem()->malloc_local_double(dima1); 

               // read doubles intermediate (right-hand-side numerator)
               doubles_intermediate_->get_block(h3b_1+z->noab()*(h2b_1+z->noab()*(h1b_1+z->noab()*(p6b_1-z->noab()+z->nvab()*(h8b_1+z->noab()*(h7b_1))))),k_a1); 
               z->sort_indices6(k_a1,k_a1_sort,z->get_range(h7b),z->get_range(h8b),z->get_range(p6b),z->get_range(h1b),z->get_range(h2b),z->get_range(h3b),5,4,3,2,1,0,+1.0); 
               z->mem()->free_local_double(k_a1); 
               double factor=1.0; 
               if (h7b==h8b) { 
                factor=factor/2.0; 
               } 
               z->smith_dgemm(dima0_sort,dima1_sort,dim_common,factor,k_a0_sort,dim_common,k_a1_sort,dim_common,1.0,k_c_sort,dima0_sort); 
               z->mem()->free_local_double(k_a1_sort); 
               z->mem()->free_local_double(k_a0_sort); 
              } 
             } 
            } 
           } 
          } 
          double* k_c=z->mem()->malloc_local_double(dimc); 
          z->sort_indices6(k_c_sort,k_c,z->get_range(h3b),z->get_range(h2b),z->get_range(h1b),z->get_range(p6b),z->get_range(h5b),z->get_range(h4b),5,4,3,2,1,0,+0.5/0.5); 
          intermediate_->add_block(h3b+z->noab()*(h2b+z->noab()*(h1b+z->noab()*(p6b-z->noab()+z->nvab()*(h5b+z->noab()*(h4b))))),k_c); 
          z->mem()->free_local_double(k_c); 
          z->mem()->free_local_double(k_c_sort); 
         } 
        } 
       } 
      } 
     } 
    } 
   } 
  } 
 } 
} 
z->mem()->sync(); 
} 

