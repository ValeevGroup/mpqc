#include <util/misc/exenv.h>
#include <util/state/stateio.h>
#include <math/scmat/blocked.h>
#include <util/class/scexception.h>
#include <algorithm>

#include <chemistry/qc/ccr12/parenthesis2t.h>
#include <chemistry/qc/ccr12/parenthesis2tnum.h>

using namespace sc;
  
static ClassDesc Parenthesis2t_cd(
  typeid(Parenthesis2t),"Parenthesis2t",1,"virtual public RefCount"
  ,0,0,0);
  
Parenthesis2t::Parenthesis2t(CCR12_Info* info): z(info){
}


Parenthesis2t::~Parenthesis2t(){
}


double Parenthesis2t::compute(Parenthesis2tNum* eval_left, 
                              Parenthesis2tNum* eval_right){

 double energy=0.0;
 double* dummy;

 // precalculating the intermediates
 eval_left->compute_amp(dummy,0L,0L,0L,0L,0L,0L,1L);  
 eval_right->compute_amp(dummy,0L,0L,0L,0L,0L,0L,1L);  

 long count=0;

 for (long t_p4b=z->noab();t_p4b<z->noab()+z->nvab();++t_p4b){
  for (long t_p5b=t_p4b;t_p5b<z->noab()+z->nvab();++t_p5b){
   for (long t_p6b=t_p5b;t_p6b<z->noab()+z->nvab();++t_p6b){
    for (long t_h1b=0;t_h1b<z->noab();++t_h1b){
     for (long t_h2b=t_h1b;t_h2b<z->noab();++t_h2b){
      for (long t_h3b=t_h2b;t_h3b<z->noab();++t_h3b){
// this will be updated
       if (count%z->mem()->n()==z->mem()->me()){ 
 
        if(!z->restricted() || z->get_spin(t_p4b)+z->get_spin(t_p5b)+z->get_spin(t_p6b)
                              +z->get_spin(t_h1b)+z->get_spin(t_h2b)+z->get_spin(t_h3b)<9L){ // should not be be 9L
        if(z->get_spin(t_p4b)+z->get_spin(t_p5b)+z->get_spin(t_p6b)==z->get_spin(t_h1b)+z->get_spin(t_h2b)+z->get_spin(t_h3b)){
        if((z->get_sym(t_p4b)^(z->get_sym(t_p5b)^(z->get_sym(t_p6b)^(z->get_sym(t_h1b)^(z->get_sym(t_h2b)^z->get_sym(t_h3b))))))==z->irrep_t()){
         long size=z->get_range(t_p4b)*z->get_range(t_p5b)*z->get_range(t_p6b)
                  *z->get_range(t_h1b)*z->get_range(t_h2b)*z->get_range(t_h3b);
         double* data_right=z->mem()->malloc_local_double(size);
         std::fill(data_right,data_right+size,0.0);
         double* data_left=z->mem()->malloc_local_double(size);
         std::fill(data_left,data_left+size,0.0);
 
         eval_left->compute_amp( data_left, t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2L);  
         eval_right->compute_amp(data_right,t_p4b,t_p5b,t_p6b,t_h1b,t_h2b,t_h3b,2L);  
 
         double* data_left_sorted=z->mem()->malloc_local_double(size);
         z->sort_indices6(data_left,data_left_sorted,z->get_range(t_h1b),z->get_range(t_h2b),z->get_range(t_h3b),
                                                     z->get_range(t_p4b),z->get_range(t_p5b),z->get_range(t_p6b),
                                                     3,4,5,0,1,2,1.0,false);
         z->mem()->free_local_double(data_left);
 
         double factor=1.0;
         if (z->restricted()) factor*=2.0; 
 
         if      (t_p4b==t_p5b && t_p5b==t_p6b) factor/=6.0; 
         else if (t_p4b==t_p5b || t_p5b==t_p6b) factor/=2.0; 
 
         if      (t_h1b==t_h2b && t_h2b==t_h3b) factor/=6.0; 
         else if (t_h1b==t_h2b || t_h2b==t_h3b) factor/=2.0; 
 
         long iall=0;
         for (long p4=0;p4<z->get_range(t_p4b);++p4) {
          const double ep4=z->get_orb_energy(z->get_offset(t_p4b)+p4);
          for (long p5=0;p5<z->get_range(t_p5b);++p5) {
           const double ep5=z->get_orb_energy(z->get_offset(t_p5b)+p5);
           for (long p6=0;p6<z->get_range(t_p6b);++p6) {
            const double ep6=z->get_orb_energy(z->get_offset(t_p6b)+p6);
            for (long h1=0;h1<z->get_range(t_h1b);++h1) {
             const double eh1=z->get_orb_energy(z->get_offset(t_h1b)+h1);
             for (long h2=0;h2<z->get_range(t_h2b);++h2) {
              const double eh2=z->get_orb_energy(z->get_offset(t_h2b)+h2);
              for (long h3=0;h3<z->get_range(t_h3b);++h3,++iall) {
               const double eh3=z->get_orb_energy(z->get_offset(t_h3b)+h3);
 
               const double numerator=data_left_sorted[iall]*data_right[iall]*factor;
               energy+=numerator/(eh1+eh2+eh3-ep4-ep5-ep6);
              }
             }
            }
           }
          }
         }
         z->mem()->free_local_double(data_left_sorted);
         z->mem()->free_local_double(data_right);
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


 // deallocating the intermediates
 eval_left->compute_amp(dummy,0L,0L,0L,0L,0L,0L,3L);  
 eval_right->compute_amp(dummy,0L,0L,0L,0L,0L,0L,3L);  

 z->mem()->sync();
 Ref<MessageGrp> msg_=MessageGrp::get_default_messagegrp();
 msg_->sum(energy);

 return energy;
}  
  
