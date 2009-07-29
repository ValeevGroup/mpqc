//
#ifndef _chemistry_qc_ccr12_ccsd_lambda_ccsdp12_t2_h
#define _chemistry_qc_ccr12_ccsd_lambda_ccsdp12_t2_h

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

class LAMBDA_CCSDPR12_T2 {

  protected:
   CCR12_Info* z;

   std::vector<Tensor*> in;
   std::vector<Tensor*> kn;

   void offset_k0();
   void smith_k0(); //z->t1(),z->v2()=>kn.at(0)
   void offset_k1();
   void smith_k1(); //z->t1(),z->v2()=>kn.at(1)
   void offset_k2();
   void smith_k2(); //z->t1(),z->l2()=>kn.at(2)
   void smith_0_1(Ref<Tensor>&); //z->v2()=>Ref<Tensor>&
   void offset_smith_0_2();
   void smith_0_2_0(); //z->v2()=>in.at(1)
   void smith_0_2_1(); //kn.at(1)=>in.at(1)
   void smith_0_2(Ref<Tensor>&); //z->l1(),in.at(1)=>Ref<Tensor>&
   void smith_0_3(Ref<Tensor>&); //z->l1(),z->v2()=>Ref<Tensor>&
   void smith_0_4(Ref<Tensor>&); //z->l1(),kn.at(0)=>Ref<Tensor>&
   void offset_smith_0_6();
   void smith_1_6(); //z->t1(),z->l1()=>in.at(1)
   void smith_1_35(); //z->t2(),z->l2()=>in.at(1)
   void smith_0_6(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_7();
   void smith_0_7_0(); //z->f1()=>in.at(1)
   void smith_1_15(); //z->t1(),z->v2()=>in.at(1)
   void smith_1_24(); //z->t1(),kn.at(0)=>in.at(1)
   void smith_1_29(); //z->t2(),z->v2()=>in.at(1)
   void smith_1_36(); //z->c2(),z->vr2()=>in.at(1)
   void smith_0_7(Ref<Tensor>&); //z->l2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_8();
   void smith_0_8_0(); //z->f1()=>in.at(1)
   void smith_1_17(); //z->t1(),z->v2()=>in.at(1)
   void smith_1_26(); //z->t1(),kn.at(0)=>in.at(1)
   void smith_1_31(); //z->t2(),z->v2()=>in.at(1)
   void smith_1_38(); //z->qy(),z->v2()=>in.at(1)
   void smith_0_8(Ref<Tensor>&); //z->l2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_9();
   void smith_0_9_0(); //z->f1()=>in.at(1)
   void smith_1_21(); //z->t1(),z->v2()=>in.at(1)
   void smith_0_9(Ref<Tensor>&); //z->ly(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_10();
   void smith_0_10_0(); //z->v2()=>in.at(1)
   void offset_smith_1_16();
   void smith_1_16_0(); //z->v2()=>in.at(2)
   void smith_1_16_1(); //kn.at(1)=>in.at(2)
   void smith_1_16(); //z->t1(),in.at(2)=>in.at(1)
   void smith_1_30(); //z->t2(),z->v2()=>in.at(1)
   void smith_1_37(); //z->c2(),z->vr2()=>in.at(1)
   void smith_0_10(Ref<Tensor>&); //z->l2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_11();
   void smith_0_11_0(); //z->v2()=>in.at(1)
   void smith_1_19(); //z->t1(),z->v2()=>in.at(1)
   void smith_1_27(); //z->t1(),kn.at(1)=>in.at(1)
   void smith_1_39(); //z->qy(),z->v2()=>in.at(1)
   void smith_0_11(Ref<Tensor>&); //z->l2(),in.at(1)=>Ref<Tensor>&
   void smith_0_12(Ref<Tensor>&); //z->l2(),z->v2()=>Ref<Tensor>&
   void offset_smith_0_13();
   void smith_0_13_0(); //z->v2()=>in.at(1)
   void smith_1_22(); //z->t1(),z->v2()=>in.at(1)
   void smith_0_13(Ref<Tensor>&); //z->ly(),in.at(1)=>Ref<Tensor>&
   void smith_0_14(Ref<Tensor>&); //z->lc2(),z->vd2()=>Ref<Tensor>&
   void smith_0_18(Ref<Tensor>&); //z->v2(),kn.at(2)=>Ref<Tensor>&
   void smith_0_20(Ref<Tensor>&); //z->v2(),kn.at(2)=>Ref<Tensor>&
   void offset_smith_0_23();
   void smith_1_23(); //z->t1(),z->ly()=>in.at(1)
   void smith_0_23(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_28();
   void smith_1_28(); //z->t1(),kn.at(2)=>in.at(1)
   void smith_1_34(); //z->t2(),z->l2()=>in.at(1)
   void smith_0_28(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_32();
   void smith_1_32(); //z->t2(),z->l2()=>in.at(1)
   void smith_0_32(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_33();
   void smith_1_33(); //z->t2(),z->l2()=>in.at(1)
   void smith_0_33(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_40();
   void smith_1_40(); //z->l2(),z->qy()=>in.at(1)
   void smith_0_40(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&

  public:
   LAMBDA_CCSDPR12_T2(CCR12_Info* info);
   ~LAMBDA_CCSDPR12_T2();
    
   void compute_amp(Ref<Tensor>& out);

};



}

#endif


