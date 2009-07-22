//
#ifndef _chemistry_qc_ccr12_ccsdt_ccsdt_t1_h
#define _chemistry_qc_ccr12_ccsdt_ccsdt_t1_h

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
#include <chemistry/qc/ccr12/ccsdt.h>


namespace sc {

class CCSDT_T1 {

  protected:
   CCR12_Info* z;
   
   std::vector<Tensor*> in;
   std::vector<Tensor*> kn;

   void offset_k0();
   void smith_k0(); //z->t1(),z->v2()=>z->kn.at(0)
   void smith_0_1(Ref<Tensor>&); //z->f1()=>out
   void offset_smith_0_2();
   void smith_0_2_0(); //z->f1()=>z->in.at(1)
   void offset_smith_1_4();
   void smith_1_4_0(); //z->f1()=>z->in.at(2)
   void smith_1_4_1(); //z->kn.at(0)=>z->in.at(2)
   void smith_1_4(); //z->t1(),z->in.at(2)=>z->in.at(1)
   void smith_1_7(); //z->t1(),z->v2()=>z->in.at(1)
   void smith_1_14(); //z->t2(),z->v2()=>z->in.at(1)
   void smith_0_2(Ref<Tensor>&); //z->t1(),z->in.at(1)=>out
   void offset_smith_0_3();
   void smith_0_3_0(); //z->f1()=>z->in.at(1)
   void smith_1_8(); //z->t1(),z->v2()=>z->in.at(1)
   void smith_0_3(Ref<Tensor>&); //z->t1(),z->in.at(1)=>out
   void offset_smith_0_5();
   void smith_0_5_0(); //z->f1()=>z->in.at(1)
   void smith_0_5_1(); //z->kn.at(0)=>z->in.at(1)
   void smith_0_5(Ref<Tensor>&); //z->t2(),z->in.at(1)=>out
   void smith_0_6(Ref<Tensor>&); //z->t1(),z->v2()=>out
   void offset_smith_0_9();
   void smith_0_9_0(); //z->v2()=>z->in.at(1)
   void smith_1_13(); //z->t1(),z->v2()=>z->in.at(1)
   void smith_0_9(Ref<Tensor>&); //z->t2(),z->in.at(1)=>out
   void smith_0_10(Ref<Tensor>&); //z->t2(),z->v2()=>out
   void smith_0_15(Ref<Tensor>&); //z->t3(),z->v2()=>out


  public:
   CCSDT_T1(CCR12_Info* info);
   ~CCSDT_T1();
    
   void compute_amp(Ref<Tensor>& out);

};



}

#endif


