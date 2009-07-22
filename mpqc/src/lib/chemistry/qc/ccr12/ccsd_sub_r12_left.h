//
#ifndef _chemistry_qc_ccr12_ccsd_sub_r12_left_h
#define _chemistry_qc_ccr12_ccsd_sub_r12_left_h

#ifdef __GNUC__
#pragma interface
#endif

#include <string>
#include <vector>
#include <math/scmat/blocked.h>
#include <util/misc/compute.h>
#include <util/group/memory.h>
#include <util/group/message.h>
#include <util/group/thread.h>

#include <chemistry/qc/ccr12/ccr12_info.h>

namespace sc {

class CCSD_SUB_R12_LEFT { 

  protected:

    CCR12_Info* z;

    std::vector<Tensor*> in; 
    std::vector<Tensor*> kn; 

    void offset_k0();
    void smith_k0(); //z->fx(),z->f1()=>kn.at(0)
    void offset_k1();
    void smith_k1(); //z->l1(),z->fx()=>kn.at(1)
    void offset_k2();
    void smith_k2(); //z->t1(),z->v2()=>kn.at(2)
    void offset_k3();
    void smith_k3(); //z->t1(),z->l2()=>kn.at(3)
    void offset_k4();
    void smith_k4(); //z->fx(),kn.at(2)=>kn.at(4)
    void smith_0_1(Ref<Tensor>&); //z->vr2()=>out
    void offset_smith_0_2();
    void smith_0_2_0(); //kn.at(0)=>in.at(1)
    void smith_0_2_1(); //z->vr2()=>in.at(1)
    void smith_0_2_2(); //kn.at(4)=>in.at(1)
    void smith_0_2(Ref<Tensor>&); //z->l1(),in.at(1)=>out
    void smith_0_3(Ref<Tensor>&); //z->v2(),kn.at(1)=>out
    void offset_smith_0_6();
    void smith_1_6(); //z->t1(),kn.at(1)=>in.at(1)
    void offset_smith_1_21();
    void smith_2_21(); //z->t2(),z->l2()=>in.at(2)
    void smith_1_21(); //z->fx(),in.at(2)=>in.at(1)
    void smith_0_6(Ref<Tensor>&); //z->v2(),in.at(1)=>out
    void offset_smith_0_7();
    void smith_1_7(); //z->t1(),z->l1()=>in.at(1)
    void smith_1_23(); //z->t2(),z->l2()=>in.at(1)
    void smith_0_7(Ref<Tensor>&); //z->vr2(),in.at(1)=>out
    void offset_smith_0_8();
    void smith_0_8_0(); //z->vr2()=>in.at(1)
    void offset_smith_1_8();
    void smith_1_8_0(); //z->f1()=>in.at(2)
    void smith_2_12(); //z->t1(),z->v2()=>in.at(2)
    void smith_2_19(); //z->t2(),z->v2()=>in.at(2)
    void smith_1_8(); //z->fx(),in.at(2)=>in.at(1)
    void offset_smith_1_9();
    void smith_1_9_0(); //kn.at(0)=>in.at(2)
    void smith_1_9_1(); //kn.at(4)=>in.at(2)
    void smith_1_9(); //z->t1(),in.at(2)=>in.at(1)
    void smith_0_8(Ref<Tensor>&); //z->l2(),in.at(1)=>out
    void offset_smith_0_10();
    void offset_smith_1_10();
    void smith_1_10_0(); //z->v2()=>in.at(2)
    void smith_2_14(); //z->t1(),z->v2()=>in.at(2)
    void smith_1_10(); //z->l2(),in.at(2)=>in.at(1)
    void smith_1_13(); //z->v2(),kn.at(3)=>in.at(1)
    void offset_smith_1_17();
    void smith_2_17(); //z->t1(),kn.at(3)=>in.at(2)
    void smith_2_20(); //z->t2(),z->l2()=>in.at(2)
    void smith_1_17(); //z->v2(),in.at(2)=>in.at(1)
    void smith_0_10(Ref<Tensor>&); //z->fx(),in.at(1)=>out
    void smith_0_15(Ref<Tensor>&); //z->vr2(),kn.at(3)=>out
    void offset_smith_0_18();
    void smith_1_18(); //z->t1(),kn.at(3)=>in.at(1)
    void smith_1_22(); //z->t2(),z->l2()=>in.at(1)
    void smith_0_18(Ref<Tensor>&); //z->vr2(),in.at(1)=>out

  public:
    CCSD_SUB_R12_LEFT(CCR12_Info* info);
    ~CCSD_SUB_R12_LEFT();

    void compute_amp(Ref<Tensor>&);

};



}

#endif


