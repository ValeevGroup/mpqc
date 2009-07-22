//
#ifndef _chemistry_qc_ccr12_ccsd_2t_right_h
#define _chemistry_qc_ccr12_ccsd_2t_right_h

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
#include <chemistry/qc/ccr12/parenthesis2tnum.h>

namespace sc {

class CCSD_2T_RIGHT : public Parenthesis2tNum {

  protected:

   void offset_k0();
   void smith_k0(); //z->t1(),z->v2()=>kn.at(0)
   void offset_k1();
   void smith_k1(); //z->t1(),z->v2()=>kn.at(1)
   void offset_smith_0_1();
   void smith_0_1_0(); //z->v2()=>in.at(1x0)
   void offset_smith_1_1();
   void smith_1_1_0(); //z->f1()=>in.at(2)
   void smith_2_19(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_1(); //z->t2(),in.at(2)=>in.at(1x0)
   void offset_smith_1_4();
   void smith_1_4_0(); //z->v2()=>in.at(2)
   void offset_smith_2_8();
   void smith_2_8_0(); //z->v2()=>in.at(3)
   void smith_2_8_1(); //kn.at(1)=>in.at(3)
   void smith_2_8(); //z->t1(),in.at(3)=>in.at(2)
   void smith_2_22(); //z->t2(),z->v2()=>in.at(2)
   void smith_1_4(); //z->t1(),in.at(2)=>in.at(1x0)
   void offset_smith_1_5();
   void smith_1_5_0(); //z->v2()=>in.at(2)
   void smith_1_5_1(); //kn.at(0)=>in.at(2)
   void smith_1_5(); //z->t1(),in.at(2)=>in.at(1x0)
   void offset_smith_1_13();
   void smith_1_13_0(); //z->v2()=>in.at(2)
   void smith_1_13_1(); //kn.at(1)=>in.at(2)
   void smith_1_13(); //z->t2(),in.at(2)=>in.at(1x0)
   void smith_1_15(); //z->t2(),z->v2()=>in.at(1x0)
   void smith_0_1(double*,const long,const long,const long,const long,const long,const long);
   void offset_smith_0_3();
   void smith_0_3_0(); //z->v2()=>in.at(1x1)
   void offset_smith_1_6();
   void smith_1_6_0(); //z->v2()=>in.at(2)
   void smith_1_6_1(); //kn.at(0)=>in.at(2)
   void offset_smith_2_10();
   void smith_2_10_0(); //z->v2()=>in.at(3)
   void smith_2_10_1(); //kn.at(1)=>in.at(3)
   void smith_2_10(); //z->t1(),in.at(3)=>in.at(2)
   void smith_2_20(); //z->t2(),z->v2()=>in.at(2)
   void smith_1_6(); //z->t1(),in.at(2)=>in.at(1x1)
   void smith_1_7(); //z->t1(),z->v2()=>in.at(1x1)
   void offset_smith_1_12();
   void smith_1_12_0(); //z->v2()=>in.at(2)
   void smith_1_12_1(); //kn.at(1)=>in.at(2)
   void smith_1_12(); //z->t2(),in.at(2)=>in.at(1x1)
   void smith_1_14(); //z->t2(),z->v2()=>in.at(1x1)
   void smith_0_3(double*,const long,const long,const long,const long,const long,const long);

  public:
   CCSD_2T_RIGHT(CCR12_Info* info);

   void compute_amp(double*,const long,const long,const long,const long,const long,const long,const long);

};



}

#endif


