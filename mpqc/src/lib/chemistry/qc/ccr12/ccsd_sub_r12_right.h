//
#ifndef _chemistry_qc_ccr12_ccsd_sub_r12_right_h
#define _chemistry_qc_ccr12_ccsd_sub_r12_right_h

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

class CCSD_SUB_R12_RIGHT { 

  protected:

    CCR12_Info* z;

    std::vector<Tensor*> in; 
    std::vector<Tensor*> kn; 

    void offset_k0();
    void smith_k0(); //z->t1(),z->fy()=>kn.at(0)
    void offset_k1();
    void smith_k1(); //z->t1(),z->v2()=>kn.at(1)
    void offset_smith_0_1();
    void smith_0_1_0(); //z->vd2()=>in.at(1)
    void offset_smith_1_1();
    void smith_1_1_0(); //z->f1()=>in.at(2)
    void smith_2_10(); //z->t1(),z->v2()=>in.at(2)
    void smith_1_1(); //z->fy(),in.at(2)=>in.at(1)
    void smith_1_12(); //z->v2(),kn.at(0)=>in.at(1)
    void smith_0_1(Ref<Tensor>&); //z->t2(),in.at(1)=>out
    void smith_0_2(Ref<Tensor>&); //z->vd2()=>out
    void smith_0_3(Ref<Tensor>&); //z->v2(),kn.at(0)=>out
    void offset_smith_0_4();
    void smith_0_4_0(); //z->vd2()=>in.at(1)
    void smith_1_5(); //z->v2(),kn.at(0)=>in.at(1)
    void smith_1_6(); //z->t1(),z->vd2()=>in.at(1)
    void smith_0_4(Ref<Tensor>&); //z->t1(),in.at(1)=>out
    void offset_smith_0_7();
    void offset_smith_1_7();
    void smith_1_7_0(); //z->v2()=>in.at(2)
    void smith_1_7_1(); //kn.at(1)=>in.at(2)
    void smith_1_7(); //z->t2(),in.at(2)=>in.at(1)
    void offset_smith_1_9();
    void smith_2_9(); //z->t1(),kn.at(1)=>in.at(2)
    void smith_1_9(); //z->t1(),in.at(2)=>in.at(1)
    void smith_0_7(Ref<Tensor>&); //z->fy(),in.at(1)=>out

  public:
    CCSD_SUB_R12_RIGHT(CCR12_Info* info);
    ~CCSD_SUB_R12_RIGHT();

    void compute_amp(Ref<Tensor>&);

};



}

#endif


