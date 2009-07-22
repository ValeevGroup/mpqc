//
#ifndef _chemistry_qc_ccr12_ccsd_ccsd_r12_t1_h
#define _chemistry_qc_ccr12_ccsd_ccsd_r12_t1_h

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


namespace sc {

class CCSD_R12_T1 {

  protected:
   CCR12_Info* z;
   std::vector<Tensor*> in;
   std::vector<Tensor*> kn;

   void offset_k0();
   void smith_k0(); //z->t1(),z->v2()=>kn.at(0)
   void smith_0_1(Ref<Tensor>&); //z->f1()=>Ref<Tensor>&
   void offset_smith_0_2();
   void smith_0_2_0(); //z->f1()=>in.at(1)
   void offset_smith_1_4();
   void smith_1_4_0(); //z->f1()=>in.at(2)
   void smith_1_4_1(); //kn.at(0)=>in.at(2)
   void smith_1_4(); //z->t1(),in.at(2)=>in.at(1)
   void smith_1_8(); //z->t1(),z->v2()=>in.at(1)
   void smith_1_17(); //z->t2(),z->v2()=>in.at(1)
   void smith_1_20(); //z->c2(),z->vr2()=>in.at(1)
   void smith_0_2(Ref<Tensor>&); //z->t1(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_3();
   void smith_0_3_0(); //z->f1()=>in.at(1)
   void smith_1_9(); //z->t1(),z->v2()=>in.at(1)
   void smith_0_3(Ref<Tensor>&); //z->t1(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_5();
   void smith_0_5_0(); //z->f1()=>in.at(1)
   void smith_0_5_1(); //kn.at(0)=>in.at(1)
   void smith_0_5(Ref<Tensor>&); //z->t2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_6();
   void smith_0_6_0(); //z->f1()=>in.at(1)
   void smith_1_18(); //z->t1(),z->v2()=>in.at(1)
   void smith_0_6(Ref<Tensor>&); //z->qy(),in.at(1)=>Ref<Tensor>&
   void smith_0_7(Ref<Tensor>&); //z->t1(),z->v2()=>Ref<Tensor>&
   void offset_smith_0_10();
   void smith_0_10_0(); //z->v2()=>in.at(1)
   void smith_1_16(); //z->t1(),z->v2()=>in.at(1)
   void smith_0_10(Ref<Tensor>&); //z->t2(),in.at(1)=>Ref<Tensor>&
   void smith_0_11(Ref<Tensor>&); //z->t2(),z->v2()=>Ref<Tensor>&
   void offset_smith_0_12();
   void smith_0_12_0(); //z->v2()=>in.at(1)
   void smith_1_19(); //z->t1(),z->v2()=>in.at(1)
   void smith_0_12(Ref<Tensor>&); //z->qy(),in.at(1)=>Ref<Tensor>&
   void smith_0_13(Ref<Tensor>&); //z->c2(),z->vr2()=>Ref<Tensor>&

  public:
   CCSD_R12_T1(CCR12_Info* info);
   ~CCSD_R12_T1();
    
   void compute_amp(Ref<Tensor>& out);

};



}

#endif


