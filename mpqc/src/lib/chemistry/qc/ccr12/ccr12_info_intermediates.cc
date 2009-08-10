//
// ccr12_info_intermediates.cc
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
#include <chemistry/qc/ccr12/ccr12_info.h>
using namespace sc;
  
// Resable intermediate F12^(aA)_(ii) * gt2^(ii)_(ii) 
// which is ubiquitous in equations.
  
void CCR12_Info::update_qy(){ 
      
for (long p3b=noab();p3b<noab()+nvab();++p3b) { 
 for (long q4b=noab()+nvab();q4b<nab();++q4b) { 
  for (long h1b=0L;h1b<noab();++h1b) { 
   for (long h2b=h1b;h2b<noab();++h2b) { 
    long tileoffset; 
    tileoffset=(h2b+noab()*(h1b+noab()*(q4b-noab()-nvab()+ncab()*(p3b-noab())))); 
    if (qy()->is_this_local(tileoffset)) { 
     if (!restricted() || get_spin(p3b)+get_spin(q4b)+get_spin(h1b)+get_spin(h2b)!=8L) { 
      if (get_spin(p3b)+get_spin(q4b)==get_spin(h1b)+get_spin(h2b)) { 
       if ((get_sym(p3b)^(get_sym(q4b)^(get_sym(h1b)^get_sym(h2b))))==(irrep_e()^irrep_t())) { 
        long dimc=get_range(p3b)*get_range(q4b)*get_range(h1b)*get_range(h2b); 
        double* k_c_sort=mem()->malloc_local_double(dimc); 
        std::fill(k_c_sort,k_c_sort+(size_t)dimc,0.0); 
        for (long h5b=0L;h5b<noab();++h5b) { 
         for (long h6b=h5b;h6b<noab();++h6b) { 
          if (get_spin(p3b)+get_spin(q4b)==get_spin(h5b)+get_spin(h6b)) { 
           if ((get_sym(p3b)^(get_sym(q4b)^(get_sym(h5b)^get_sym(h6b))))==irrep_e()) { 
            long p3b_0,q4b_0,h5b_0,h6b_0; 
            restricted_4(p3b,q4b,h5b,h6b,p3b_0,q4b_0,h5b_0,h6b_0); 
            long h5b_1,h6b_1,h1b_1,h2b_1; 
            restricted_4(h5b,h6b,h1b,h2b,h5b_1,h6b_1,h1b_1,h2b_1); 
            long dim_common=get_range(h5b)*get_range(h6b); 
            long dima0_sort=get_range(p3b)*get_range(q4b); 
            long dima0=dim_common*dima0_sort; 
            long dima1_sort=get_range(h1b)*get_range(h2b); 
            long dima1=dim_common*dima1_sort; 
            if (dima0>0L && dima1>0L) { 
             double* k_a0_sort=mem()->malloc_local_double(dima0); 
             double* k_a0=mem()->malloc_local_double(dima0); 
             fr2()->get_block(h6b_0+noab()*(h5b_0+noab()*(q4b_0+(nab())*(p3b_0))),k_a0); 
             sort_indices4(k_a0,k_a0_sort,get_range(p3b),get_range(q4b),get_range(h5b),get_range(h6b),1,0,3,2,+1.0,false); 
             mem()->free_local_double(k_a0); 
             double* k_a1_sort=mem()->malloc_local_double(dima1); 
             double* k_a1=mem()->malloc_local_double(dima1); 
             gt2()->get_block(h2b_1+noab()*(h1b_1+noab()*(h6b_1+noab()*(h5b_1))),k_a1); 
             sort_indices4(k_a1,k_a1_sort,get_range(h5b),get_range(h6b),get_range(h1b),get_range(h2b),3,2,1,0,+1.0,false); 
             mem()->free_local_double(k_a1); 
             double factor=1.0; 
             if (h5b==h6b) { 
              factor=factor/2.0; 
             } 
             smith_dgemm(dima0_sort,dima1_sort,dim_common,factor,k_a0_sort,dim_common,k_a1_sort,dim_common,1.0,k_c_sort,dima0_sort); 
             mem()->free_local_double(k_a1_sort); 
             mem()->free_local_double(k_a0_sort); 
            } 
           } 
          } 
         } 
        } 
        double* k_c=mem()->malloc_local_double(dimc); 
        sort_indices4(k_c_sort,k_c,get_range(h2b),get_range(h1b),get_range(q4b),get_range(p3b),3,2,1,0,+0.5/0.5,false); 
        qy()->add_block(h2b+noab()*(h1b+noab()*(q4b-noab()-nvab()+ncab()*(p3b-noab()))),k_c); 
        mem()->free_local_double(k_c); 
        mem()->free_local_double(k_c_sort); 
       } 
      } 
     } 
    } 
   } 
  } 
 } 
} 
mem()->sync(); 
} 


void CCR12_Info::update_ly(){ 

for (long h3b=0L;h3b<noab();++h3b) { 
 for (long h4b=h3b;h4b<noab();++h4b) { 
  for (long p1b=noab();p1b<noab()+nvab();++p1b) { 
   for (long q2b=noab()+nvab();q2b<nab();++q2b) { 
    long tileoffset; 
    tileoffset=(q2b-noab()-nvab()+ncab()*(p1b-noab()+nvab()*(h4b+noab()*(h3b)))); 
    if (ly()->is_this_local(tileoffset)) { 
     if (!restricted() || get_spin(h3b)+get_spin(h4b)+get_spin(p1b)+get_spin(q2b)!=8L) { 
      if (get_spin(h3b)+get_spin(h4b)==get_spin(p1b)+get_spin(q2b)) { 
       if ((get_sym(h3b)^(get_sym(h4b)^(get_sym(p1b)^get_sym(q2b))))==(irrep_e()^irrep_t())) { 
        long dimc=get_range(h3b)*get_range(h4b)*get_range(p1b)*get_range(q2b); 
        double* k_c_sort=mem()->malloc_local_double(dimc); 
        std::fill(k_c_sort,k_c_sort+(size_t)dimc,0.0); 
        for (long h5b=0L;h5b<noab();++h5b) { 
         for (long h6b=h5b;h6b<noab();++h6b) { 
          if (get_spin(h5b)+get_spin(h6b)==get_spin(p1b)+get_spin(q2b)) { 
           if ((get_sym(h5b)^(get_sym(h6b)^(get_sym(p1b)^get_sym(q2b))))==irrep_e()) { 
            long h5b_0,h6b_0,p1b_0,q2b_0; 
            restricted_4(h5b,h6b,p1b,q2b,h5b_0,h6b_0,p1b_0,q2b_0); 
            long h3b_1,h4b_1,h5b_1,h6b_1; 
            restricted_4(h3b,h4b,h5b,h6b,h3b_1,h4b_1,h5b_1,h6b_1); 
            long dim_common=get_range(h5b)*get_range(h6b); 
            long dima0_sort=get_range(p1b)*get_range(q2b); 
            long dima0=dim_common*dima0_sort; 
            long dima1_sort=get_range(h3b)*get_range(h4b); 
            long dima1=dim_common*dima1_sort; 
            if (dima0>0L && dima1>0L) { 
             double* k_a0_sort=mem()->malloc_local_double(dima0); 
             double* k_a0=mem()->malloc_local_double(dima0); 
             fd2()->get_block(q2b_0+(nab())*(p1b_0+(nab())*(h6b_0+noab()*(h5b_0))),k_a0); 
             sort_indices4(k_a0,k_a0_sort,get_range(h5b),get_range(h6b),get_range(p1b),get_range(q2b),3,2,1,0,+1.0,false); 
             mem()->free_local_double(k_a0); 
             double* k_a1_sort=mem()->malloc_local_double(dima1); 
             double* k_a1=mem()->malloc_local_double(dima1); 
             glambda2()->get_block(h6b_1+noab()*(h5b_1+noab()*(h4b_1+noab()*(h3b_1))),k_a1); 
             sort_indices4(k_a1,k_a1_sort,get_range(h3b),get_range(h4b),get_range(h5b),get_range(h6b),1,0,3,2,+1.0,false); 
             mem()->free_local_double(k_a1); 
             double factor=1.0; 
             if (h5b==h6b) { 
              factor=factor/2.0; 
             } 
             smith_dgemm(dima0_sort,dima1_sort,dim_common,factor,k_a0_sort,dim_common,k_a1_sort,dim_common,1.0,k_c_sort,dima0_sort); 
             mem()->free_local_double(k_a1_sort); 
             mem()->free_local_double(k_a0_sort); 
            } 
           } 
          } 
         } 
        } 
        double* k_c=mem()->malloc_local_double(dimc); 
        sort_indices4(k_c_sort,k_c,get_range(h4b),get_range(h3b),get_range(q2b),get_range(p1b),1,0,3,2,+0.5/0.5,false); 
        ly()->add_block(q2b-noab()-nvab()+ncab()*(p1b-noab()+nvab()*(h4b+noab()*(h3b))),k_c); 
        mem()->free_local_double(k_c); 
        mem()->free_local_double(k_c_sort); 
       } 
      } 
     } 
    } 
   } 
  } 
 } 
} 
mem()->sync(); 
} 


void CCR12_Info::prod_iiii(Ref<Tensor>& a, Ref<Tensor>& b, Ref<Tensor>& out){ 
      
long tileoffset; 
if (out->is_this_local(0L)) { 
  long dimc=1; 
  double* k_c_sort=mem()->malloc_local_double(dimc); 
  std::fill(k_c_sort,k_c_sort+(size_t)dimc,0.0); 
  for (long h1b=0L;h1b<noab();++h1b) { 
   for (long h2b=h1b;h2b<noab();++h2b) { 
    for (long h3b=0L;h3b<noab();++h3b) { 
     for (long h4b=h3b;h4b<noab();++h4b) { 
      if (get_spin(h1b)+get_spin(h2b)==get_spin(h3b)+get_spin(h4b)) { 
       if ((get_sym(h1b)^(get_sym(h2b)^(get_sym(h3b)^get_sym(h4b))))==irrep_e()) { 
        long h1b_0,h2b_0,h3b_0,h4b_0; 
        restricted_4(h1b,h2b,h3b,h4b,h1b_0,h2b_0,h3b_0,h4b_0); 
        long h3b_1,h4b_1,h1b_1,h2b_1; 
        restricted_4(h3b,h4b,h1b,h2b,h3b_1,h4b_1,h1b_1,h2b_1); 
        long dim_common=get_range(h1b)*get_range(h2b)*get_range(h3b)*get_range(h4b); 
        long dima0_sort=1L; 
        long dima0=dim_common*dima0_sort; 
        long dima1_sort=1L; 
        long dima1=dim_common*dima1_sort; 
        if (dima0>0L && dima1>0L) { 
         double* k_a0_sort=mem()->malloc_local_double(dima0); 
         double* k_a0=mem()->malloc_local_double(dima0); 
         a->get_block(h4b_0+noab()*(h3b_0+noab()*(h2b_0+noab()*(h1b_0))),k_a0); 
         sort_indices4(k_a0,k_a0_sort,get_range(h1b),get_range(h2b),get_range(h3b),get_range(h4b),3,2,1,0,+1.0,false); 
         mem()->free_local_double(k_a0); 
         double* k_a1_sort=mem()->malloc_local_double(dima1); 
         double* k_a1=mem()->malloc_local_double(dima1); 
         b->get_block(h2b_1+noab()*(h1b_1+noab()*(h4b_1+noab()*(h3b_1))),k_a1); 
         sort_indices4(k_a1,k_a1_sort,get_range(h3b),get_range(h4b),get_range(h1b),get_range(h2b),1,0,3,2,+1.0,false); 
         mem()->free_local_double(k_a1); 
         double factor=1.0; 
         if (h1b==h2b) { 
          factor=factor/2.0; 
         } 
         if (h3b==h4b) { 
          factor=factor/2.0; 
         } 
         smith_dgemm(dima0_sort,dima1_sort,dim_common,factor,k_a0_sort,dim_common,k_a1_sort,dim_common,1.0,k_c_sort,dima0_sort); 
         mem()->free_local_double(k_a1_sort); 
         mem()->free_local_double(k_a0_sort); 
        } 
       } 
      } 
     } 
    } 
   } 
  } 
  double* k_c=mem()->malloc_local_double(dimc); 
  sort_indices0(k_c_sort,k_c,1.0,false); 
  out->add_block((0),k_c); 
  mem()->free_local_double(k_c); 
  mem()->free_local_double(k_c_sort); 
} 
mem()->sync(); 
} 
