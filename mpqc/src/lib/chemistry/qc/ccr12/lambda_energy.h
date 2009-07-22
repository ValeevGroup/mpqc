//
#ifndef _chemistry_qc_ccr12_lambda_energy_h
#define _chemistry_qc_ccr12_lambda_energy_h

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

class LAMBDA_ENERGY {

  protected:
   CCR12_Info* z;
   std::vector<Tensor*> in;
   std::vector<Tensor*> kn;

   void offset_k0();
   void smith_k0(); //z->t2(),z->ly()=>kn.at(0)
   void offset_k1();
   void smith_k1(); //z->lc2(),z->xs2()=>kn.at(1)
   void offset_k2();
   void smith_k2(); //z->t1(),z->ly()=>kn.at(2)
   void offset_k3();
   void smith_k3(); //z->qx(),z->lx()=>kn.at(3)
   void offset_k4();
   void smith_k4(); //z->qx(),z->ly()=>kn.at(4)
   void offset_k5();
   void smith_k5(); //z->qy(),z->lx()=>kn.at(5)
   void offset_k6();
   void smith_k6(); //z->qy(),z->ly()=>kn.at(6)
   void offset_k7();
   void smith_k7(); //z->qx(),z->lx()=>kn.at(7)
   void offset_k8();
   void smith_k8(); //z->qy(),z->ly()=>kn.at(8)
   void offset_k9();
   void smith_k9(); //z->qx(),z->ly()=>kn.at(9)
   void offset_k10();
   void smith_k10(); //z->qy(),z->lx()=>kn.at(10)
   void offset_k11();
   void smith_k11(); //z->qy(),z->ly()=>kn.at(11)
   void offset_k12();
   void smith_k12(); //z->c2(),kn.at(1)=>kn.at(12)
   void offset_k13();
   void smith_k13(); //z->qy(),kn.at(2)=>kn.at(13)
   void offset_k14();
   void smith_k14(); //z->c2(),kn.at(1)=>kn.at(14)
   void smith_0_1(Ref<Tensor>&); //z->f1(),kn.at(0)=>Ref<Tensor>&
   void smith_0_2(Ref<Tensor>&); //z->f1(),kn.at(12)=>Ref<Tensor>&
   void offset_smith_0_3();
   void offset_smith_1_3();
   void smith_1_3_0(); //z->bs2()=>in.at(2)
   void smith_1_3_1(); //z->ps2()=>in.at(2)
   void smith_1_3(); //z->lc2(),in.at(2)=>in.at(1)
   void offset_smith_1_46();
   void offset_smith_2_46();
   void smith_3_46(); //z->t2(),z->v2()=>in.at(3)
   void smith_2_46(); //z->lc2(),in.at(3)=>in.at(2)
   void offset_smith_2_47();
   void smith_3_47(); //z->t2(),z->v2()=>in.at(3)
   void smith_2_47(); //z->lc2(),in.at(3)=>in.at(2)
   void smith_1_46(); //z->xs2(),in.at(2)=>in.at(1)
   void smith_1_68(); //z->vr2(),kn.at(14)=>in.at(1)
   void smith_1_69(); //z->vr2(),kn.at(12)=>in.at(1)
   void smith_0_3(Ref<Tensor>&); //z->c2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_4();
   void smith_1_4(); //z->f1(),kn.at(12)=>in.at(1)
   void offset_smith_1_5();
   void smith_1_5_0(); //z->v2()=>in.at(2)
   void smith_2_5(); //z->qx(),z->f1()=>in.at(2)
   void offset_smith_2_10();
   void smith_2_10_0(); //z->v2()=>in.at(3)
   void smith_3_21(); //z->t1(),z->v2()=>in.at(3)
   void smith_2_10(); //z->t1(),in.at(3)=>in.at(2)
   void smith_2_24(); //z->t2(),z->v2()=>in.at(2)
   void smith_2_32(); //z->qx(),z->v2()=>in.at(2)
   void smith_2_33(); //z->qy(),z->v2()=>in.at(2)
   void smith_2_39(); //z->c2(),z->vr2()=>in.at(2)
   void smith_1_5(); //z->ly(),in.at(2)=>in.at(1)
   void smith_1_9(); //z->lc2(),z->vd2()=>in.at(1)
   void offset_smith_1_11();
   void smith_2_11(); //z->t1(),z->lc2()=>in.at(2)
   void smith_1_11(); //z->vd2(),in.at(2)=>in.at(1)
   void offset_smith_1_22();
   void smith_1_22_0(); //kn.at(0)=>in.at(2)
   void smith_1_22_1(); //kn.at(10)=>in.at(2)
   void smith_1_22(); //z->v2(),in.at(2)=>in.at(1)
   void offset_smith_1_23();
   void smith_1_23_0(); //kn.at(5)=>in.at(2)
   void smith_2_23(); //z->t2(),z->ly()=>in.at(2)
   void smith_1_23(); //z->v2(),in.at(2)=>in.at(1)
   void smith_1_25(); //z->v2(),kn.at(12)=>in.at(1)
   void smith_1_26(); //z->v2(),kn.at(14)=>in.at(1)
   void offset_smith_1_27();
   void smith_1_27_0(); //kn.at(7)=>in.at(2)
   void smith_1_27_1(); //kn.at(8)=>in.at(2)
   void smith_1_27(); //z->v2(),in.at(2)=>in.at(1)
   void smith_1_29(); //z->v2(),kn.at(9)=>in.at(1)
   void smith_1_31(); //z->v2(),kn.at(11)=>in.at(1)
   void offset_smith_1_34();
   void smith_1_34_0(); //kn.at(3)=>in.at(2)
   void smith_2_35(); //z->qy(),z->ly()=>in.at(2)
   void smith_1_34(); //z->v2(),in.at(2)=>in.at(1)
   void smith_1_36(); //z->v2(),kn.at(4)=>in.at(1)
   void smith_1_38(); //z->v2(),kn.at(6)=>in.at(1)
   void offset_smith_1_40();
   void smith_1_40_0(); //kn.at(13)=>in.at(2)
   void smith_2_40(); //z->t1(),kn.at(12)=>in.at(2)
   void smith_1_40(); //z->v2(),in.at(2)=>in.at(1)
   void offset_smith_1_41();
   void smith_2_41(); //z->t1(),kn.at(14)=>in.at(2)
   void smith_2_45(); //z->qy(),kn.at(2)=>in.at(2)
   void smith_1_41(); //z->v2(),in.at(2)=>in.at(1)
   void offset_smith_1_42();
   void smith_2_42(); //z->qx(),kn.at(2)=>in.at(2)
   void smith_1_42(); //z->v2(),in.at(2)=>in.at(1)
   void offset_smith_1_44();
   void smith_2_44(); //z->qx(),kn.at(2)=>in.at(2)
   void smith_1_44(); //z->v2(),in.at(2)=>in.at(1)
   void smith_0_4(Ref<Tensor>&); //z->t1(),in.at(1)=>Ref<Tensor>&
   void smith_0_6(Ref<Tensor>&); //z->f1(),kn.at(13)=>Ref<Tensor>&
   void offset_smith_0_7();
   void smith_0_7_0(); //z->vd2()=>in.at(1)
   void smith_1_13(); //z->t2(),z->vd2()=>in.at(1)
   void smith_0_7(Ref<Tensor>&); //z->lc2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_12();
   void smith_1_12(); //z->ly(),z->v2()=>in.at(1)
   void smith_1_48(); //z->v2(),kn.at(9)=>in.at(1)
   void smith_1_50(); //z->v2(),kn.at(4)=>in.at(1)
   void offset_smith_1_52();
   void smith_2_52(); //z->qx(),z->v2()=>in.at(2)
   void smith_1_52(); //z->ly(),in.at(2)=>in.at(1)
   void smith_0_12(Ref<Tensor>&); //z->t2(),in.at(1)=>Ref<Tensor>&
   void smith_0_14(Ref<Tensor>&); //z->v2(),kn.at(14)=>Ref<Tensor>&
   void smith_0_15(Ref<Tensor>&); //z->v2(),kn.at(3)=>Ref<Tensor>&
   void offset_smith_0_16();
   void smith_1_16(); //z->ly(),z->v2()=>in.at(1)
   void smith_1_56(); //z->v2(),kn.at(4)=>in.at(1)
   void offset_smith_1_62();
   void smith_2_62(); //z->qx(),z->v2()=>in.at(2)
   void smith_1_62(); //z->ly(),in.at(2)=>in.at(1)
   void smith_1_63(); //z->v2(),kn.at(7)=>in.at(1)
   void smith_0_16(Ref<Tensor>&); //z->qy(),in.at(1)=>Ref<Tensor>&
   void smith_0_17(Ref<Tensor>&); //z->v2(),kn.at(4)=>Ref<Tensor>&
   void smith_0_18(Ref<Tensor>&); //z->v2(),kn.at(5)=>Ref<Tensor>&
   void smith_0_19(Ref<Tensor>&); //z->v2(),kn.at(6)=>Ref<Tensor>&
   void offset_smith_0_49();
   void smith_1_49(); //z->t2(),kn.at(11)=>in.at(1)
   void smith_1_51(); //z->t2(),kn.at(6)=>in.at(1)
   void offset_smith_1_53();
   void smith_1_53_0(); //kn.at(0)=>in.at(2)
   void smith_1_53_1(); //kn.at(10)=>in.at(2)
   void smith_1_53(); //z->qy(),in.at(2)=>in.at(1)
   void smith_1_59(); //z->qy(),kn.at(5)=>in.at(1)
   void smith_0_49(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_54();
   void smith_1_54(); //z->qx(),kn.at(3)=>in.at(1)
   void smith_1_60(); //z->qx(),kn.at(7)=>in.at(1)
   void smith_0_54(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_55();
   void offset_smith_1_55();
   void smith_2_55(); //z->qy(),z->v2()=>in.at(2)
   void smith_1_55(); //z->qx(),in.at(2)=>in.at(1)
   void offset_smith_1_61();
   void smith_2_61(); //z->qy(),z->v2()=>in.at(2)
   void smith_1_61(); //z->qx(),in.at(2)=>in.at(1)
   void smith_0_55(Ref<Tensor>&); //z->ly(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_57();
   void smith_1_57(); //z->qx(),kn.at(5)=>in.at(1)
   void smith_1_58(); //z->qy(),kn.at(6)=>in.at(1)
   void smith_1_64(); //z->qy(),kn.at(8)=>in.at(1)
   void smith_1_65(); //z->qx(),kn.at(10)=>in.at(1)
   void smith_1_66(); //z->qy(),kn.at(11)=>in.at(1)
   void smith_0_57(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&

  public:
   LAMBDA_ENERGY(CCR12_Info* info);
   ~LAMBDA_ENERGY();
    
   void compute_amp(Ref<Tensor>& out);

};



}

#endif


