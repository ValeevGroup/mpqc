//
#ifndef _chemistry_qc_ccr12_ccsd_p2_t_num_h
#define _chemistry_qc_ccr12_ccsd_p2_t_num_h

#ifdef __GNUC__
#pragma interface
#endif

#include <string>
#include <math/scmat/blocked.h>
#include <util/misc/compute.h>
#include <util/group/memory.h>
#include <util/group/message.h>
#include <util/group/thread.h>
#include <chemistry/qc/ccr12/ccr12_info.h>

namespace sc {

class Parenthesis2tNum : virtual public RefCount {

  protected:
   CCR12_Info* z;

   std::vector<Tensor*> in; 
   std::vector<Tensor*> kn; 
   std::vector<Tensor*> i1xn; 

  public:
   Parenthesis2tNum(CCR12_Info* info);
    
   ~Parenthesis2tNum();

   virtual void compute_amp(double*,const long,const long,const long,const long,const long,const long,const long);
   virtual void compute_amp(double*,const long,const long,const long,const long,const long,const long,const long,const long,const long);

};



}

#endif


