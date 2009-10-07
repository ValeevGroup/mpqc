//
// ccsd_pt.cc
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
#include <chemistry/qc/ccr12/ccsd_pt.h>

using namespace sc;
  

static ClassDesc CCSD_PT_cd(
  typeid(CCSD_PT),"CCSD_PT",1,"public Parenthesis2t"
  ,0,0,0);
  

double CCSD_PT::compute_energy(Ref<PTNum> eval_left, Ref<PTNum> eval_right){

 double energy = 0.0;

 long count = 0L;

 const size_t maxtile = z->maxtilesize();
 const size_t mem_singles = maxtile * maxtile;
 const size_t mem_doubles = mem_singles * mem_singles;
 const size_t mem_triples = mem_singles * mem_doubles;

 double* k_a0     = z->mem()->malloc_local_double(mem_doubles); 
 double* k_a0_sort= z->mem()->malloc_local_double(mem_doubles); 
 double* k_a1     = z->mem()->malloc_local_double(mem_doubles); 
 double* k_a1_sort= z->mem()->malloc_local_double(mem_doubles); 
 double* k_c_sort = z->mem()->malloc_local_double(mem_triples); 

 double* doubles = z->mem()->malloc_local_double(mem_triples);
 double* singles = z->mem()->malloc_local_double(mem_triples);

 double* work_doubles[6] = {doubles, k_a0, k_a0_sort, k_a1, k_a1_sort, k_c_sort};
 double* work_singles[6] = {singles, k_a0, k_a0_sort, k_a1, k_a1_sort, k_c_sort};

 for (long t_p4b = z->noab(); t_p4b<z->noab() + z->nvab(); ++t_p4b){
  for (long t_p5b = t_p4b;     t_p5b<z->noab() + z->nvab(); ++t_p5b){
   for (long t_p6b = t_p5b;     t_p6b<z->noab() + z->nvab(); ++t_p6b){

    for (long t_h1b = 0L;    t_h1b<z->noab(); ++t_h1b){
     for (long t_h2b = t_h1b; t_h2b<z->noab(); ++t_h2b){
      for (long t_h3b  =t_h2b; t_h3b<z->noab(); ++t_h3b){

       // Naive way of parallelization; will be updated
       if (count%z->mem()->n() == z->mem()->me()){ 

        if(!z->restricted() ||  z->get_spin(t_p4b) + z->get_spin(t_p5b) + z->get_spin(t_p6b)
                              + z->get_spin(t_h1b) + z->get_spin(t_h2b) + z->get_spin(t_h3b) < 9L){
        if(z->get_spin(t_p4b) + z->get_spin(t_p5b) + z->get_spin(t_p6b) == z->get_spin(t_h1b) + z->get_spin(t_h2b) + z->get_spin(t_h3b)){
        if((z->get_sym(t_p4b)^(z->get_sym(t_p5b)^(z->get_sym(t_p6b)^(z->get_sym(t_h1b)^(z->get_sym(t_h2b)^z->get_sym(t_h3b)))))) == 0L){
 
         const long size = z->get_range(t_p4b) * z->get_range(t_p5b) * z->get_range(t_p6b)
                         * z->get_range(t_h1b) * z->get_range(t_h2b) * z->get_range(t_h3b);

         std::fill(doubles, doubles + size, 0.0);
         std::fill(singles, singles + size, 0.0);

         eval_left->compute_amp(work_singles, t_p4b, t_p5b, t_p6b, t_h1b, t_h2b, t_h3b, 2L);  
         eval_right->compute_amp(work_doubles, t_p4b, t_p5b, t_p6b, t_h1b, t_h2b, t_h3b, 2L);  

         double factor = 1.0;
         if (z->restricted()) factor *= 2.0; 
 
         if      (t_p4b == t_p5b && t_p5b == t_p6b) factor *= 0.166666666666666667;
         else if (t_p4b == t_p5b || t_p5b == t_p6b) factor *= 0.5; 
 
         if      (t_h1b == t_h2b && t_h2b == t_h3b) factor *= 0.166666666666666667; 
         else if (t_h1b == t_h2b || t_h2b == t_h3b) factor *= 0.5; 

         long iall = 0L;
         for (long p4 = 0L; p4 < z->get_range(t_p4b); ++p4) {
          const double ep4 = z->get_orb_energy(z->get_offset(t_p4b) + p4);
          for (long p5 = 0L; p5 < z->get_range(t_p5b); ++p5) {
           const double ep5 = z->get_orb_energy(z->get_offset(t_p5b) + p5);
           for (long p6 = 0L; p6 < z->get_range(t_p6b); ++p6) {
            const double ep6 = z->get_orb_energy(z->get_offset(t_p6b) + p6);
            const double eps = ep4 + ep5 + ep6;

            for (long h1 = 0L; h1 < z->get_range(t_h1b); ++h1) {
             const double eh1 = z->get_orb_energy(z->get_offset(t_h1b) + h1);
             for (long h2 = 0L; h2 < z->get_range(t_h2b); ++h2) {
              const double eh2 = z->get_orb_energy(z->get_offset(t_h2b) + h2);
              for (long h3 = 0L; h3 < z->get_range(t_h3b); ++h3, ++iall) {
               const double eh3 = z->get_orb_energy(z->get_offset(t_h3b) + h3);

               const double numerator = factor * (singles[iall] + doubles[iall]) * doubles[iall];
               energy += numerator / (eh1 + eh2 + eh3 - eps);
              }
             }
            }
           }
          }
         }
        }
        }
        }
       } // end parallel loop 
       ++count;
      }
     }
    }
   }
  }
 }
 z->mem()->free_local_double(singles);
 z->mem()->free_local_double(doubles);
 z->mem()->free_local_double(k_c_sort); 
 z->mem()->free_local_double(k_a0); 
 z->mem()->free_local_double(k_a0_sort); 
 z->mem()->free_local_double(k_a1); 
 z->mem()->free_local_double(k_a1_sort); 

 z->mem()->sync();
 Ref<MessageGrp> msg_=MessageGrp::get_default_messagegrp();
 msg_->sum(energy);

 return energy;
}  
  
