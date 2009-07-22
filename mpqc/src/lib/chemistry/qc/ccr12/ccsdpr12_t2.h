//
#ifndef _chemistry_qc_ccr12_ccsd_ccsdpr12_t2_h
#define _chemistry_qc_ccr12_ccsd_ccsdpr12_t2_h

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

class CCSDPR12_T2 {

  protected:
   CCR12_Info* z;
   std::vector<Tensor*> in;
   std::vector<Tensor*> kn;

   void offset_k0();
   void smith_k0(); //z->t1(),z->v2()=>kn.at(0)
   void offset_k1();
   void smith_k1(); //z->t1(),z->v2()=>kn.at(1)
   void offset_k2();
   void smith_k2(); //z->t1(),z->v2()=>kn.at(2)
   void offset_k3();
   void smith_k3(); //z->gt2(),z->vr2()=>kn.at(3)
   void offset_k4();
   void smith_k4(); //z->t2(),z->v2()=>kn.at(4)
   void offset_k5();
   void smith_k5(); //z->t1(),z->v2()=>kn.at(5)
   void offset_k6();
   void smith_k6(); //z->t1(),kn.at(2)=>kn.at(6)
   void offset_smith_0_1();
   void smith_0_1_0(); //z->f1()=>in.at(1)
   void smith_1_25(); //z->t1(),z->v2()=>in.at(1)
   void smith_1_34(); //z->t1(),kn.at(5)=>in.at(1)
   void smith_1_39(); //z->t2(),z->v2()=>in.at(1)
   void smith_1_44(); //z->gt2(),z->vr2()=>in.at(1)
   void smith_0_1(Ref<Tensor>&); //z->t2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_2();
   void smith_0_2_0(); //z->f1()=>in.at(1)
   void smith_1_26(); //z->t1(),z->v2()=>in.at(1)
   void smith_1_37(); //z->t2(),z->v2()=>in.at(1)
   void smith_1_42(); //z->qy(),z->v2()=>in.at(1)
   void smith_0_2(Ref<Tensor>&); //z->t2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_3();
   void smith_0_3_0(); //z->f1()=>in.at(1)
   void smith_1_17(); //z->t1(),z->v2()=>in.at(1)
   void smith_1_40(); //z->t2(),z->v2()=>in.at(1)
   void smith_0_3(Ref<Tensor>&); //z->qy(),in.at(1)=>Ref<Tensor>&
   void smith_0_4(Ref<Tensor>&); //z->v2()=>Ref<Tensor>&
   void offset_smith_0_5();
   void smith_0_5_0(); //z->v2()=>in.at(1)
   void offset_smith_1_7();
   void smith_1_7_0(); //z->v2()=>in.at(2)
   void smith_1_7_1(); //kn.at(0)=>in.at(2)
   void smith_1_7_2(); //kn.at(6)=>in.at(2)
   void smith_1_7_3(); //kn.at(3)=>in.at(2)
   void smith_1_7_4(); //kn.at(4)=>in.at(2)
   void smith_1_7(); //z->t1(),in.at(2)=>in.at(1)
   void offset_smith_1_8();
   void smith_1_8_0(); //z->v2()=>in.at(2)
   void smith_1_8_1(); //kn.at(1)=>in.at(2)
   void smith_1_8(); //z->t1(),in.at(2)=>in.at(1)
   void offset_smith_1_18();
   void smith_1_18_0(); //z->v2()=>in.at(2)
   void smith_2_29(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_18(); //z->qy(),in.at(2)=>in.at(1)
   void smith_1_20(); //z->gt2(),z->vr2()=>in.at(1)
   void offset_smith_1_21();
   void smith_1_21_0(); //z->v2()=>in.at(2)
   void smith_1_21_1(); //kn.at(2)=>in.at(2)
   void smith_1_21(); //z->t2(),in.at(2)=>in.at(1)
   void smith_1_22(); //z->t2(),z->v2()=>in.at(1)
   void offset_smith_1_28();
   void smith_2_28(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_28(); //z->qy(),in.at(2)=>in.at(1)
   void smith_1_35(); //z->t2(),kn.at(5)=>in.at(1)
   void smith_0_5(Ref<Tensor>&); //z->t1(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_6();
   void smith_0_6_0(); //z->v2()=>in.at(1)
   void smith_1_9(); //z->t1(),z->v2()=>in.at(1)
   void smith_0_6(Ref<Tensor>&); //z->t1(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_10();
   void smith_0_10_0(); //z->v2()=>in.at(1)
   void smith_0_10_1(); //kn.at(0)=>in.at(1)
   void smith_0_10_2(); //kn.at(6)=>in.at(1)
   void smith_0_10_3(); //kn.at(4)=>in.at(1)
   void smith_0_10_4(); //kn.at(3)=>in.at(1)
   void smith_0_10(Ref<Tensor>&); //z->t2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_11();
   void smith_0_11_0(); //z->v2()=>in.at(1)
   void smith_0_11_1(); //kn.at(1)=>in.at(1)
   void smith_1_36(); //z->t2(),z->v2()=>in.at(1)
   void smith_1_41(); //z->qy(),z->v2()=>in.at(1)
   void smith_0_11(Ref<Tensor>&); //z->t2(),in.at(1)=>Ref<Tensor>&
   void smith_0_12(Ref<Tensor>&); //z->t2(),z->v2()=>Ref<Tensor>&
   void offset_smith_0_13();
   void smith_0_13_0(); //z->v2()=>in.at(1)
   void smith_1_19(); //z->t1(),z->v2()=>in.at(1)
   void smith_0_13(Ref<Tensor>&); //z->qy(),in.at(1)=>Ref<Tensor>&
   void smith_0_14(Ref<Tensor>&); //z->gt2(),z->vr2()=>Ref<Tensor>&

  public:
   CCSDPR12_T2(CCR12_Info* info);
   ~CCSDPR12_T2();
    
   void compute_amp(Ref<Tensor>& out);

};



}

#endif


