//
// triples_energy.cc
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
#include <chemistry/qc/mbptr12/blas.h>

using namespace sc;
  

double CCR12_Triples::get_energy() { 

 double energy = 0.0;
 long count = 0L;
      
 for (long h1b=0L;h1b<z->noab();++h1b) { 
  for (long h2b=h1b;h2b<z->noab();++h2b) { 
   for (long p3b=z->noab();p3b<z->noab()+z->nvab();++p3b) { 
    for (long h4b=0L;h4b<z->noab();++h4b) { 
     for (long h5b=h4b;h5b<z->noab();++h5b) { 
      for (long h6b=h5b;h6b<z->noab();++h6b) { 
// this will be updated
       if (count%z->mem()->n() == z->mem()->me()){ 
        if (z->get_spin(h1b)+z->get_spin(h2b)+z->get_spin(p3b)==z->get_spin(h4b)+z->get_spin(h5b)+z->get_spin(h6b)) { 
         if ((z->get_sym(h1b)^(z->get_sym(h2b)^(z->get_sym(p3b)^(z->get_sym(h4b)^(z->get_sym(h5b)^z->get_sym(h6b))))))==z->irrep_t()) { 
          long h1b_0,h2b_0,p3b_0,h4b_0,h5b_0,h6b_0; 
          z->restricted_6(h1b,h2b,p3b,h4b,h5b,h6b,h1b_0,h2b_0,p3b_0,h4b_0,h5b_0,h6b_0); 
          long h4b_1,h5b_1,h6b_1,h1b_1,h2b_1,p3b_1; 
          z->restricted_6(h4b,h5b,h6b,h1b,h2b,p3b,h4b_1,h5b_1,h6b_1,h1b_1,h2b_1,p3b_1); 
          const int dim=z->get_range(h1b)*z->get_range(h2b)*z->get_range(h4b)*z->get_range(h5b)*z->get_range(h6b)*z->get_range(p3b); 
          if (dim>0L) { 
           double* k_a0=z->mem()->malloc_local_double(dim); 
           singles_intermediate_->get_block(h6b_0+z->noab()*(h5b_0+z->noab()*(h4b_0+z->noab()*(p3b_0-z->noab()+z->nvab()*(h2b_0+z->noab()*(h1b_0))))),k_a0); 
           double* k_a1=z->mem()->malloc_local_double(dim); 
           intermediate_->get_block(h6b_0+z->noab()*(h5b_0+z->noab()*(h4b_0+z->noab()*(p3b_0-z->noab()+z->nvab()*(h2b_0+z->noab()*(h1b_0))))),k_a1); 
           double factor=1.0; 
           if (h1b==h2b) factor=factor/2.0; 
           if (h4b==h5b && h5b==h6b) factor=factor/6.0; 
           if (h4b==h5b && h4b!=h6b) factor=factor/2.0; 
           if (h5b==h6b && h5b!=h4b) factor=factor/2.0; 
           if (h4b==h6b && h4b!=h5b) { 
            factor=factor/2.0; 
            assert(false);
           } 
           const int unit = 1;
           energy += factor * F77_DDOT(&dim, k_a0, &unit, k_a1, &unit);
           z->mem()->free_local_double(k_a1); 
           z->mem()->free_local_double(k_a0); 
          } 
         } 
        } 
       } 
       ++count;
      } 
     } 
    } 
   } 
  } 
 }

 z->mem()->sync();
 Ref<MessageGrp> msg_=MessageGrp::get_default_messagegrp();
 msg_->sum(energy);

 return energy;
} 

