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
#include <iostream>

using namespace sc;
using namespace std;


double CCR12_Triples::get_energy_ip() {

 double energy = 0.0;

 const size_t maxsize1 = z->maxtilesize();
 const size_t maxsize3 = maxsize1 * maxsize1 * maxsize1;
 const size_t size_alloc = maxsize3 * maxsize3;
 double* work0 = z->mem()->malloc_local_double(size_alloc);
 double* work1 = z->mem()->malloc_local_double(size_alloc);

 const long noab = z->noab();
 const long nvab = z->nvab();

 long count = 0L;
 for (long h1b = 0L; h1b < noab; ++h1b) {
  for (long h2b = h1b; h2b < noab+nvab; ++h2b) {
   for (long p3b = noab; p3b < noab + nvab; ++p3b) {
    for (long h4b = 0L; h4b < noab; ++h4b) {
     for (long h5b = h4b; h5b < noab; ++h5b) {
      for (long h6b = h5b; h6b < noab; ++h6b, ++count) {
       // the most primitive way of parallelizing...
       if (count%z->mem()->n() == z->mem()->me()){
        if (z->get_spin(h1b) + z->get_spin(h2b) + z->get_spin(p3b) ==
            z->get_spin(h4b) + z->get_spin(h5b) + z->get_spin(h6b)) {
         if ((z->get_sym(h1b)^(z->get_sym(h2b)^(z->get_sym(p3b)^
                 (z->get_sym(h4b)^(z->get_sym(h5b)^z->get_sym(h6b)))))) == z->irrep_t()) {

          // For RHF reference: mapping beta indices to alpha
          long h1ba, h2ba, p3ba, h4ba, h5ba, h6ba;
          z->restricted_6(h1b, h2b, p3b, h4b, h5b, h6b,
                          h1ba, h2ba, p3ba, h4ba, h5ba, h6ba);
          const int dim = z->get_range(h1b) * z->get_range(h2b)
                        * z->get_range(h4b) * z->get_range(h5b)
                        * z->get_range(h6b) * z->get_range(p3b);

          // read blocks
          const long tag = h6ba + noab * (h5ba + noab * (h4ba + noab * (p3ba - noab + nvab * (h2ba + (noab+nvab) * h1ba))));
          singles_intermediate_->get_block(tag, work0);
          rhs_intermediate_->get_block(tag, work1);

          // prefactor
          double factor = 1.0;
          if (h1b==h2b) factor *= 0.5;
          if (h4b==h5b && h5b==h6b) factor /= 6.0;
          else if (h4b==h5b || h5b==h6b) factor *= 0.5;

          // adds the contribution from this block
          const int unit = 1;
          energy += factor * F77_DDOT(&dim, work0, &unit, work1, &unit);

         }
        }
       }
      }
     }
    }
   }
  }
 }
 z->mem()->free_local_double(work1);
 z->mem()->free_local_double(work0);

 z->mem()->sync();
 Ref<MessageGrp> msg_ = MessageGrp::get_default_messagegrp();
 msg_->sum(energy);

 return energy;
}
  

double CCR12_Triples::get_energy() { 

 double energy = 0.0;
 
 const size_t maxsize1 = z->maxtilesize();
 const size_t maxsize3 = maxsize1 * maxsize1 * maxsize1;
 const size_t size_alloc = maxsize3 * maxsize3;
 double* work0 = z->mem()->malloc_local_double(size_alloc); 
 double* work1 = z->mem()->malloc_local_double(size_alloc); 

 const long noab = z->noab();
 const long nvab = z->nvab();
      
 long count = 0L;
 for (long h1b = 0L; h1b < noab; ++h1b) {
  for (long h2b = h1b; h2b < noab; ++h2b) {
   for (long p3b = noab; p3b < noab + nvab; ++p3b) {
    for (long h4b = 0L; h4b < noab; ++h4b) {
     for (long h5b = h4b; h5b < noab; ++h5b) {
      for (long h6b = h5b; h6b < noab; ++h6b, ++count) {
       // the most primitive way of parallelizing...
       if (count%z->mem()->n() == z->mem()->me()){ 
        if (z->get_spin(h1b) + z->get_spin(h2b) + z->get_spin(p3b) ==
            z->get_spin(h4b) + z->get_spin(h5b) + z->get_spin(h6b)) { 
         if ((z->get_sym(h1b)^(z->get_sym(h2b)^(z->get_sym(p3b)^
        		 (z->get_sym(h4b)^(z->get_sym(h5b)^z->get_sym(h6b)))))) == z->irrep_t()) {

          // For RHF reference: mapping beta indices to alpha
          long h1ba, h2ba, p3ba, h4ba, h5ba, h6ba; 
          z->restricted_6(h1b, h2b, p3b, h4b, h5b, h6b,
                          h1ba, h2ba, p3ba, h4ba, h5ba, h6ba);
          const int dim = z->get_range(h1b) * z->get_range(h2b)
                        * z->get_range(h4b) * z->get_range(h5b)
                        * z->get_range(h6b) * z->get_range(p3b); 

          // read blocks
          const long tag = h6ba + noab * (h5ba + noab * (h4ba + noab * (p3ba - noab + nvab * (h2ba + noab * h1ba))));
          singles_intermediate_->get_block(tag, work0); 
          rhs_intermediate_->get_block(tag, work1);

          // prefactor
          double factor = 1.0; 
          if (h1b==h2b) factor *= 0.5; 
          if (h4b==h5b && h5b==h6b) factor /= 6.0; 
          else if (h4b==h5b || h5b==h6b) factor *= 0.5;

          // adds the contribution from this block
          const int unit = 1;
          energy += factor * F77_DDOT(&dim, work0, &unit, work1, &unit);
          //cout << setprecision(15) << endl;
          //cout << h1b << h2b << p3b << h4b << h5b << h6b << " " << factor * F77_DDOT(&dim, work0, &unit, work1, &unit) << endl;

         }
        }
       }
      }
     }
    }
   }
  }
 }
 z->mem()->free_local_double(work1); 
 z->mem()->free_local_double(work0);

 z->mem()->sync();
 Ref<MessageGrp> msg_ = MessageGrp::get_default_messagegrp();
 msg_->sum(energy);

 return energy;
} 

