//
#ifndef _chemistry_qc_ccr12_lambdapenergy_h
#define _chemistry_qc_ccr12_lambdapenergy_h

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

class LAMBDAPENERGY {

  protected:
   CCR12_Info* z;
   std::vector<Tensor*> in;
   std::vector<Tensor*> kn;

   void offset_k0();
   void smith_k0(); //z->t2(),z->ly()=>kn.at(0)
   void smith_0_1(Ref<Tensor>&); //z->f1(),kn.at(0)=>Ref<Tensor>&
   void offset_smith_0_2();
   void offset_smith_1_2();
   void smith_2_2(); //z->lc2(),z->xs2()=>in.at(2)
   void smith_1_2(); //z->c2(),in.at(2)=>in.at(1)
   void smith_0_2(Ref<Tensor>&); //z->f1(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_3();
   void smith_1_3(); //z->lc2(),z->bs2()=>in.at(1)
   void smith_0_3(Ref<Tensor>&); //z->c2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_4();
   void smith_0_4_0(); //z->vd2()=>in.at(1)
   void smith_1_10(); //z->t2(),z->vd2()=>in.at(1)
   void smith_0_4(Ref<Tensor>&); //z->lc2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_5();
   void offset_smith_1_5();
   void smith_1_5_0(); //z->v2()=>in.at(2)
   void offset_smith_2_7();
   void smith_2_7_0(); //z->v2()=>in.at(3)
   void smith_3_11(); //z->t1(),z->v2()=>in.at(3)
   void smith_2_7(); //z->t1(),in.at(3)=>in.at(2)
   void smith_2_14(); //z->t2(),z->v2()=>in.at(2)
   void smith_1_5(); //z->ly(),in.at(2)=>in.at(1)
   void smith_1_6(); //z->lc2(),z->vd2()=>in.at(1)
   void offset_smith_1_8();
   void smith_2_8(); //z->t1(),z->lc2()=>in.at(2)
   void smith_1_8(); //z->vd2(),in.at(2)=>in.at(1)
   void smith_1_12(); //z->v2(),kn.at(0)=>in.at(1)
   void offset_smith_1_13();
   void smith_2_13(); //z->t2(),z->ly()=>in.at(2)
   void smith_1_13(); //z->v2(),in.at(2)=>in.at(1)
   void smith_0_5(Ref<Tensor>&); //z->t1(),in.at(1)=>out
   void offset_smith_0_9();
   void smith_1_9(); //z->ly(),z->v2()=>in.at(1)
   void smith_0_9(Ref<Tensor>&); //z->t2(),in.at(1)=>out

  public:
   LAMBDAPENERGY(CCR12_Info* info);
   ~LAMBDAPENERGY();
    
   void compute_amp(Ref<Tensor>& out);

};



}

#endif


