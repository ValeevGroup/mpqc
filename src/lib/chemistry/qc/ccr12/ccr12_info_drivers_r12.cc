//
// ccr12_info_drivers_r12.cc
//
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
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

#include <sstream>
#include <algorithm>
#include <cassert>
#include <chemistry/qc/ccr12/ccr12_info.h>
#include <chemistry/qc/ccr12/tensor.h>
#include <chemistry/qc/ccr12/mtensor.h>

using namespace sc;
using namespace std;

void CCR12_Info::jacobi_t2_and_gt2_(const Ref<Tensor>& d_r2_, Ref<Tensor>& d_t2_,
                                    const Ref<Tensor>& d_gr2_,Ref<Tensor>& d_gt2_){

  MPQC_ASSERT(need_w1());
  
  /// TODO most likely we can reduce the number of temporary tensors in the following.
  /// So far two T2-sized tensor has been created at the same time here, which 
  /// might be prohibitive in really large calculations

  /// obtain d_r2/Diag
  Ref<Tensor> t2_mp2 = new Tensor("t2_mp2_local",mem_);
  offset_t2(t2_mp2, false);
  jacobi_t2_(d_r2_, t2_mp2);

  /// obtain duplicate of d_gr2  
  Ref<Tensor> vv = d_gr2_->copy();

  /// obtain Ad = f^alpha_a R^ij_alpha b (Eq. 13 of Valeev and Janssen 2004 JCP) 
  Ref<Tensor> ad = new Tensor("ad_local",mem_);
  offset_l2(ad);
  // form_ad(ad);

  /// obtain \v{V} (Eq. 18 of Valeev and Janssen 2004 JCP)
  form_adt(ad,t2_mp2,vv);

  /// obtain copt  (Eq. 17 of Valeev and Janssen 2004 JCP)

  Ref<Tensor> copt = d_gr2_->clone();
  denom_contraction(vv, copt);
  d_gt2_->daxpy(copt, 1.0);

  /// obtain effective T2 (The numerator of Eq. 20 of Valeev and Janssen 2004 JCP)
  /// we will reuse ad; every time when we ad->get_block, we need to invert it.
  t2_mp2->zero();
  t2_mp2->daxpy(d_r2_, 1.0); //data copy
  form_ca(copt, ad, t2_mp2);
  
  /// devide by the denominator and add it to d_t2_
  jacobi_t2_(t2_mp2, d_t2_);

}


void CCR12_Info::form_ad(Ref<Tensor>& out) {

 for (long h3b=0L;h3b<noab();++h3b) { 
  for (long h4b=h3b;h4b<noab();++h4b) { 
   for (long p1b=noab();p1b<noab()+nvab();++p1b) { 
    for (long p2b=noab();p2b<noab()+nvab();++p2b) { 
     long tileoffset; 
     if (p1b<p2b) 
      tileoffset=(p2b-noab()+nvab()*(p1b-noab()+nvab()*(h4b+noab()*(h3b)))); 
     else if (p2b<=p1b) 
      tileoffset=(p1b-noab()+nvab()*(p2b-noab()+nvab()*(h4b+noab()*(h3b)))); 
     if (out->is_this_local(tileoffset)) { 
      if (!restricted() || get_spin(h3b)+get_spin(h4b)+get_spin(p1b)+get_spin(p2b)!=8L) { 
       if (get_spin(h3b)+get_spin(h4b)==get_spin(p1b)+get_spin(p2b)) { 
        if ((get_sym(h3b)^(get_sym(h4b)^(get_sym(p1b)^get_sym(p2b))))==(irrep_e()^irrep_f())) { 
         long dimc=get_range(h3b)*get_range(h4b)*get_range(p1b)*get_range(p2b); 
         double* k_c_sort=mem()->malloc_local_double(dimc); 
         std::fill(k_c_sort,k_c_sort+(size_t)dimc,0.0); 
         for (long q5b=noab()+nvab();q5b<nab();++q5b) { 
          if (get_spin(h3b)+get_spin(h4b)==get_spin(p1b)+get_spin(q5b)) { 
           if ((get_sym(h3b)^(get_sym(h4b)^(get_sym(p1b)^get_sym(q5b))))==irrep_e()) { 
            long h3b_0,h4b_0,p1b_0,q5b_0; 
            restricted_4(h3b,h4b,p1b,q5b,h3b_0,h4b_0,p1b_0,q5b_0); 
            long q5b_1,p2b_1; 
            restricted_2(q5b,p2b,q5b_1,p2b_1); 
            long dim_common=get_range(q5b); 
            long dima0_sort=get_range(h3b)*get_range(h4b)*get_range(p1b); 
            long dima0=dim_common*dima0_sort;
            long dima1_sort=get_range(p2b); 
            long dima1=dim_common*dima1_sort; 
            if (dima0>0L && dima1>0L) { 
             double* k_a0_sort=mem()->malloc_local_double(dima0); 
             double* k_a0=mem()->malloc_local_double(dima0); 
             fd2()->get_block(q5b_0+nab()*(p1b_0+nab()*(h4b_0+noab()*h3b_0)),k_a0); 
             sort_indices4(k_a0,k_a0_sort,get_range(h3b),get_range(h4b),get_range(p1b),get_range(q5b),2,1,0,3,+1.0); 
             mem()->free_local_double(k_a0); 
             double* k_a1_sort=mem()->malloc_local_double(dima1); 
             double* k_a1=mem()->malloc_local_double(dima1); 
             f1()->get_block(p2b_1+nab()*q5b_1,k_a1); 
             sort_indices2(k_a1,k_a1_sort,get_range(q5b),get_range(p2b),1,0,+1.0); 
             mem()->free_local_double(k_a1); 
             double factor=1.0; 
             smith_dgemm(dima0_sort,dima1_sort,dim_common,factor,k_a0_sort,dim_common,k_a1_sort,dim_common,1.0,k_c_sort,dima0_sort); 
             mem()->free_local_double(k_a1_sort); 
             mem()->free_local_double(k_a0_sort); 
            } 
           } 
          } 
         } 
         double* k_c=mem()->malloc_local_double(dimc); 
         if (p2b>=p1b) { 
          sort_indices4(k_c_sort,k_c,get_range(p2b),get_range(p1b),get_range(h4b),get_range(h3b),3,2,1,0,+1.0); 
          out->add_block(p2b-noab()+nvab()*(p1b-noab()+nvab()*(h4b+noab()*h3b)),k_c); 
         } 
         if (p1b>=p2b) { 
          sort_indices4(k_c_sort,k_c,get_range(p2b),get_range(p1b),get_range(h4b),get_range(h3b),3,2,0,1,-1.0); 
          out->add_block(p1b-noab()+nvab()*(p2b-noab()+nvab()*(h4b+noab()*h3b)),k_c); 
         } 
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


void CCR12_Info::form_adt(const Ref<Tensor>& inp_l2, const Ref<Tensor>& inp_t2, Ref<Tensor>& out) {
 for (long h3b=0L;h3b<noab();++h3b) { 
  for (long h4b=h3b;h4b<noab();++h4b) { 
   for (long h1b=0L;h1b<noab();++h1b) { 
    for (long h2b=h1b;h2b<noab();++h2b) { 
     long tileoffset; 
     tileoffset=(h2b+noab()*(h1b+noab()*(h4b+noab()*(h3b)))); 
     if (out->is_this_local(tileoffset)) { 
      if (!restricted() || get_spin(h3b)+get_spin(h4b)+get_spin(h1b)+get_spin(h2b)!=8L) { 
       if (get_spin(h3b)+get_spin(h4b)==get_spin(h1b)+get_spin(h2b)) { 
        if ((get_sym(h3b)^(get_sym(h4b)^(get_sym(h1b)^get_sym(h2b))))==(irrep_t()^irrep_e())) { 
         long dimc=get_range(h3b)*get_range(h4b)*get_range(h1b)*get_range(h2b); 
         double* k_c_sort=mem()->malloc_local_double(dimc); 
         std::fill(k_c_sort,k_c_sort+(size_t)dimc,0.0); 
         for (long p5b=noab();p5b<noab()+nvab();++p5b) { 
          for (long p6b=p5b;p6b<noab()+nvab();++p6b) { 
           if (get_spin(p5b)+get_spin(p6b)==get_spin(h1b)+get_spin(h2b)) { 
            if ((get_sym(p5b)^(get_sym(p6b)^(get_sym(h1b)^get_sym(h2b))))==irrep_t()) { 
             long p5b_0,p6b_0,h1b_0,h2b_0; 
             restricted_4(p5b,p6b,h1b,h2b,p5b_0,p6b_0,h1b_0,h2b_0); 
             long h3b_1,h4b_1,p5b_1,p6b_1; 
             restricted_4(h3b,h4b,p5b,p6b,h3b_1,h4b_1,p5b_1,p6b_1); 
             long dim_common=get_range(p5b)*get_range(p6b); 
             long dima0_sort=get_range(h1b)*get_range(h2b); 
             long dima0=dim_common*dima0_sort; 
             long dima1_sort=get_range(h3b)*get_range(h4b); 
             long dima1=dim_common*dima1_sort; 
             if (dima0>0L && dima1>0L) { 
              double* k_a0_sort=mem()->malloc_local_double(dima0); 
              double* k_a0=mem()->malloc_local_double(dima0); 
              inp_t2->get_block(h2b_0+noab()*(h1b_0+noab()*(p6b_0-noab()+nvab()*(p5b_0-noab()))),k_a0); 
              sort_indices4(k_a0,k_a0_sort,get_range(p5b),get_range(p6b),get_range(h1b),get_range(h2b),3,2,1,0,+1.0); 
              mem()->free_local_double(k_a0); 
              double* k_a1_sort=mem()->malloc_local_double(dima1); 
              double* k_a1=mem()->malloc_local_double(dima1); 
              inp_l2->get_block(p6b_1-noab()+nvab()*(p5b_1-noab()+nvab()*(h4b_1+noab()*(h3b_1))),k_a1); 
              sort_indices4(k_a1,k_a1_sort,get_range(h3b),get_range(h4b),get_range(p5b),get_range(p6b),1,0,3,2,+1.0); 
              mem()->free_local_double(k_a1); 
              double factor=1.0; 
              if (p5b==p6b) factor=factor/2.0; 
              smith_dgemm(dima0_sort,dima1_sort,dim_common,factor,k_a0_sort,dim_common,k_a1_sort,dim_common,1.0,k_c_sort,dima0_sort); 
              mem()->free_local_double(k_a1_sort); 
              mem()->free_local_double(k_a0_sort); 
             } 
            } 
           } 
          } 
         } 
         double* k_c=mem()->malloc_local_double(dimc); 
         sort_indices4(k_c_sort,k_c,get_range(h4b),get_range(h3b),get_range(h2b),get_range(h1b),1,0,3,2,+0.5/0.5); 
         out->add_block(h2b+noab()*(h1b+noab()*(h4b+noab()*(h3b))),k_c); 
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


void CCR12_Info::form_ca(const Ref<Tensor>& c_, const Ref<Tensor>& ad_, Ref<Tensor>& out) {
 for (long p3b=noab();p3b<noab()+nvab();++p3b) { 
  for (long p4b=p3b;p4b<noab()+nvab();++p4b) { 
   for (long h1b=0L;h1b<noab();++h1b) { 
    for (long h2b=h1b;h2b<noab();++h2b) { 
     long tileoffset; 
     tileoffset=(h2b+noab()*(h1b+noab()*(p4b-noab()+nvab()*(p3b-noab())))); 
     if (out->is_this_local(tileoffset)) { 
      if (!restricted() || get_spin(p3b)+get_spin(p4b)+get_spin(h1b)+get_spin(h2b)!=8L) { 
       if (get_spin(p3b)+get_spin(p4b)==get_spin(h1b)+get_spin(h2b)) { 
        if ((get_sym(p3b)^(get_sym(p4b)^(get_sym(h1b)^get_sym(h2b))))==(irrep_t()^irrep_t())) { 
         long dimc=get_range(p3b)*get_range(p4b)*get_range(h1b)*get_range(h2b); 
         double* k_c_sort=mem()->malloc_local_double(dimc); 
         std::fill(k_c_sort,k_c_sort+(size_t)dimc,0.0); 
         for (long h5b=0L;h5b<noab();++h5b) { 
          for (long h6b=h5b;h6b<noab();++h6b) { 
           if (get_spin(p3b)+get_spin(p4b)==get_spin(h5b)+get_spin(h6b)) { 
            if ((get_sym(p3b)^(get_sym(p4b)^(get_sym(h5b)^get_sym(h6b))))==irrep_t()) { 
             long p3b_0,p4b_0,h5b_0,h6b_0; 
             restricted_4(p3b,p4b,h5b,h6b,p3b_0,p4b_0,h5b_0,h6b_0); 
             long h5b_1,h6b_1,h1b_1,h2b_1; 
             restricted_4(h5b,h6b,h1b,h2b,h5b_1,h6b_1,h1b_1,h2b_1); 
             long dim_common=get_range(h5b)*get_range(h6b); 
             long dima0_sort=get_range(p3b)*get_range(p4b); 
             long dima0=dim_common*dima0_sort; 
             long dima1_sort=get_range(h1b)*get_range(h2b); 
             long dima1=dim_common*dima1_sort; 
             if (dima0>0L && dima1>0L) { 
              double* k_a0_sort=mem()->malloc_local_double(dima0); 
              double* k_a0=mem()->malloc_local_double(dima0); 

              //a_->get_block(h6b_0+noab()*(h5b_0+noab()*(p4b_0-noab()+nvab()*(p3b_0-noab()))),k_a0); 
              //sort_indices4(k_a0,k_a0_sort,get_range(p3b),get_range(p4b),get_range(h5b),get_range(h6b),1,0,3,2,+1.0); 
              // we are reusing a conjugate tensor
              ad_->get_block(p4b_0-noab()+nvab()*(p3b_0-noab()+nvab()*(h6b_0+noab()*h5b_0)),k_a0); 
              //sort_indices4(k_a0,k_a0_sort,get_range(h5b),get_range(h6b),get_range(p3b),get_range(p4b),2,3,0,1,+1.0); 
              sort_indices4(k_a0,k_a0_sort,get_range(h5b),get_range(h6b),get_range(p3b),get_range(p4b),3,2,1,0,+1.0); 

              mem()->free_local_double(k_a0); 
              double* k_a1_sort=mem()->malloc_local_double(dima1); 
              double* k_a1=mem()->malloc_local_double(dima1); 
              c_->get_block(h2b_1+noab()*(h1b_1+noab()*(h6b_1+noab()*h5b_1)),k_a1); 
              sort_indices4(k_a1,k_a1_sort,get_range(h5b),get_range(h6b),get_range(h1b),get_range(h2b),3,2,1,0,+1.0); 
              mem()->free_local_double(k_a1); 
              double factor=1.0; 
              if (h5b==h6b) factor=factor/2.0; 
              smith_dgemm(dima0_sort,dima1_sort,dim_common,factor,k_a0_sort,dim_common,k_a1_sort,dim_common,1.0,k_c_sort,dima0_sort); 
              mem()->free_local_double(k_a1_sort); 
              mem()->free_local_double(k_a0_sort); 
             } 
            } 
           } 
          } 
         } 
         double* k_c=mem()->malloc_local_double(dimc); 
         sort_indices4(k_c_sort,k_c,get_range(h2b),get_range(h1b),get_range(p4b),get_range(p3b),3,2,1,0,+0.5/0.5); 
         out->add_block(h2b+noab()*(h1b+noab()*(p4b-noab()+nvab()*(p3b-noab()))),k_c); 
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

