//
#ifndef _chemistry_qc_ccr12_ccsd_ccsd_e_h
#define _chemistry_qc_ccr12_ccsd_ccsd_e_h

#ifdef __GNUC__
#pragma interface
#endif

#include <string>
#include <util/misc/compute.h>
#include <util/group/memory.h>
#include <util/group/message.h>
#include <util/group/thread.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>

#include <chemistry/qc/ccr12/tensor.h>
#include <chemistry/qc/ccr12/ccr12.h>
#include <chemistry/qc/ccr12/ccsd.h>
#include <chemistry/qc/ccr12/ccr12_info.h>

namespace sc {

class CCSD_E {

  protected:
   CCR12_Info* z;
 
   std::vector<Tensor*> in;
   std::vector<Tensor*> kn;

   void offset_smith_0_1();
   void smith_0_1_0(); //z->f1()=>z->in[1]
   void smith_1_2(); //z->t1(),z->v2()=>z->in[1]
   void smith_0_1(Ref<Tensor>& out); //z->t1(),z->in[1]=>out
   void smith_0_3(Ref<Tensor>& out); //z->t2(),z->v2()=>out

  public:
   CCSD_E(CCR12_Info* info);
    
   ~CCSD_E();
   void compute_amp(Ref<Tensor>& out);

};



}

#endif


