//
// parenthesis2q.cc
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
#include <chemistry/qc/ccr12/parenthesis2q.h>

using namespace sc;
using namespace std;
  
static ClassDesc Parenthesis2q_cd(
  typeid(Parenthesis2q),"Parenthesis2q",1,"public RefCount"
  ,0,0,0);
  
Parenthesis2q::Parenthesis2q(CCR12_Info* info): z(info){
}


Parenthesis2q::~Parenthesis2q(){
}


double Parenthesis2q::compute_energy(Ref<Parenthesis2tNum> eval_left, 
                                     Ref<Parenthesis2tNum> eval_right){

 double energy=0.0;
 double dummy = 0.0;

 // precalculating the intermediates
 Ref<RegionTimer> timer_ = new RegionTimer();
 timer_->enter("(2)Q intermediates");
 const double start = timer_->get_wall_time();
 eval_left->compute_amp(&dummy,0L,0L,0L,0L,0L,0L,0L,0L,1L);
 eval_right->compute_amp(&dummy,0L,0L,0L,0L,0L,0L,0L,0L,1L);
 if (z->mem()->me()==0) {
   ExEnv::out0() << indent << fixed    << "Elapsed time [ "<< "(2)Q intermediates" << " ]: "
                           << setw(10) << setprecision(2)  << timer_->get_wall_time() - start << endl;
 }

 long count=0L;

 const size_t maxsize1 = z->maxtilesize();
 const size_t maxsize2 = maxsize1 * maxsize1;
 const size_t maxsize4 = maxsize2 * maxsize2;
 const size_t size_alloc = maxsize4 * maxsize4;

 double* data_right=z->mem()->malloc_local_double(size_alloc);
 double* data_left=z->mem()->malloc_local_double(size_alloc);
 double* data_left_sorted=z->mem()->malloc_local_double(size_alloc);

 for (long t_p5b=z->noab();t_p5b<z->noab()+z->nvab();++t_p5b){
  if (true) ExEnv::out0() << indent << indent << "o Block " << setw(3) << t_p5b << " being processed." << endl;
  for (long t_p6b=t_p5b;    t_p6b<z->noab()+z->nvab();++t_p6b){
   for (long t_p7b=t_p6b;    t_p7b<z->noab()+z->nvab();++t_p7b){
    for (long t_p8b=t_p7b;    t_p8b<z->noab()+z->nvab();++t_p8b){
     for (long t_h1b=0L;   t_h1b<z->noab();++t_h1b){
      for (long t_h2b=t_h1b;t_h2b<z->noab();++t_h2b){
       for (long t_h3b=t_h2b;t_h3b<z->noab();++t_h3b){
        for (long t_h4b=t_h3b;t_h4b<z->noab();++t_h4b){
// this will be updated
         if (count%z->mem()->n()==z->mem()->me()){ 
   
          if((z->get_sym(t_p5b)^(z->get_sym(t_p6b)^(z->get_sym(t_p7b)^(z->get_sym(t_p8b)^(z->get_sym(t_h1b)^(z->get_sym(t_h2b)^(z->get_sym(t_h3b)^z->get_sym(t_h4b))))))))==z->irrep_t()){
          if(z->get_spin(t_p5b)+z->get_spin(t_p6b)+z->get_spin(t_p7b)+z->get_spin(t_p8b)==z->get_spin(t_h1b)+z->get_spin(t_h2b)+z->get_spin(t_h3b)+z->get_spin(t_h4b)){
          if(!z->restricted() || z->get_spin(t_p5b)+z->get_spin(t_p6b)+z->get_spin(t_p7b)+z->get_spin(t_p8b)
                                +z->get_spin(t_h1b)+z->get_spin(t_h2b)+z->get_spin(t_h3b)+z->get_spin(t_h4b)<=12L){
           long size=z->get_range(t_p5b)*z->get_range(t_p6b)*z->get_range(t_p7b)*z->get_range(t_p8b)
                    *z->get_range(t_h1b)*z->get_range(t_h2b)*z->get_range(t_h3b)*z->get_range(t_h4b);
           std::fill(data_right,data_right+size,0.0);
           std::fill(data_left,data_left+size,0.0);
   
           eval_left->compute_amp( data_left, t_h1b,t_h2b,t_h3b,t_h4b,t_p5b,t_p6b,t_p7b,t_p8b,2L);  
           eval_right->compute_amp(data_right,t_p5b,t_p6b,t_p7b,t_p8b,t_h1b,t_h2b,t_h3b,t_h4b,2L);
   
#if 0
           z->sort_indices8(data_left,data_left_sorted,z->get_range(t_h1b),z->get_range(t_h2b),z->get_range(t_h3b),z->get_range(t_h4b),
                                                       z->get_range(t_p5b),z->get_range(t_p6b),z->get_range(t_p7b),z->get_range(t_p8b),
                                                       4,5,6,7,0,1,2,3,1.0);
#else
           z->sort_indices2(data_left,data_left_sorted,z->get_range(t_h1b)*z->get_range(t_h2b)*z->get_range(t_h3b)*z->get_range(t_h4b),
                                                       z->get_range(t_p5b)*z->get_range(t_p6b)*z->get_range(t_p7b)*z->get_range(t_p8b),
                                                       1,0,1.0);
#endif
   
           double factor=1.0;
           if (z->restricted() && z->get_spin(t_p5b)+z->get_spin(t_p6b)+z->get_spin(t_p7b)+z->get_spin(t_p8b)
                                 +z->get_spin(t_h1b)+z->get_spin(t_h2b)+z->get_spin(t_h3b)+z->get_spin(t_h4b)!=12L) factor*=2.0; 
   
           if      (t_p5b==t_p6b && t_p6b==t_p7b && t_p7b==t_p8b) factor *= (1.0/24.0);
           else if (t_p5b==t_p6b && t_p6b==t_p7b) factor *= 0.166666666666666667;
           else if (t_p6b==t_p7b && t_p7b==t_p8b) factor *= 0.166666666666666667;
           else if (t_p5b==t_p6b && t_p7b==t_p8b) factor *= 0.25;
           else if (t_p5b==t_p6b || t_p6b==t_p7b || t_p7b==t_p8b) factor *= 0.5;
   
           if      (t_h1b==t_h2b && t_h2b==t_h3b && t_h3b==t_h4b) factor *= (1.0/24.0);
           else if (t_h1b==t_h2b && t_h2b==t_h3b) factor *= 0.166666666666666667;
           else if (t_h2b==t_h3b && t_h3b==t_h4b) factor *= 0.166666666666666667;
           else if (t_h1b==t_h2b && t_h3b==t_h4b) factor *= 0.25;
           else if (t_h1b==t_h2b || t_h2b==t_h3b || t_h3b==t_h4b) factor *= 0.5;
   
           long iall=0;
           for (long p5=0;p5<z->get_range(t_p5b);++p5) {
            const double ep5=z->get_orb_energy(z->get_offset(t_p5b)+p5);
            for (long p6=0;p6<z->get_range(t_p6b);++p6) {
             const double ep6=z->get_orb_energy(z->get_offset(t_p6b)+p6);
             for (long p7=0;p7<z->get_range(t_p7b);++p7) {
              const double ep7=z->get_orb_energy(z->get_offset(t_p7b)+p7);
              for (long p8=0;p8<z->get_range(t_p8b);++p8) {
               const double ep8=z->get_orb_energy(z->get_offset(t_p8b)+p8);
               const double eps = ep5 + ep6 + ep7 + ep8; 

               for (long h1=0;h1<z->get_range(t_h1b);++h1) {
                const double eh1=z->get_orb_energy(z->get_offset(t_h1b)+h1);
                for (long h2=0;h2<z->get_range(t_h2b);++h2) {
                 const double eh2=z->get_orb_energy(z->get_offset(t_h2b)+h2);
                 for (long h3=0;h3<z->get_range(t_h3b);++h3) {
                  const double eh3=z->get_orb_energy(z->get_offset(t_h3b)+h3);
                  const double eh123 = eh1 + eh2 + eh3 - eps;

                  long h4a = z->get_offset(t_h4b);
                  for (long h4=0;h4<z->get_range(t_h4b);++h4, ++h4a, ++iall) {
                   const double eh4=z->get_orb_energy(h4a);
  
                   const double numerator=data_left_sorted[iall]*data_right[iall]*factor;
                   energy+=numerator/(eh123+eh4);
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
          }
         } // end of parallel loop
         ++count;
        }
       }
      }
     }
    }
   }
  }
 }

 z->mem()->free_local_double(data_left_sorted);
 z->mem()->free_local_double(data_left);
 z->mem()->free_local_double(data_right);


 // deallocating the intermediates
 eval_left->compute_amp(&dummy,0L,0L,0L,0L,0L,0L,0L,0L,3L);
 eval_right->compute_amp(&dummy,0L,0L,0L,0L,0L,0L,0L,0L,3L);

 z->mem()->sync();
 Ref<MessageGrp> msg_=MessageGrp::get_default_messagegrp();
 msg_->sum(energy);

 return energy;
}  

