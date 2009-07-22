//
#ifndef _chemistry_qc_ccr12_ccsd_2q_left_h
#define _chemistry_qc_ccr12_ccsd_2q_left_h

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

class CCSD_2Q_LEFT : public Parenthesis2tNum {

  protected:
  
   void smith_0_1(double*,const long,const long,const long,const long,const long,const long,const long,const long);

  public:
   CCSD_2Q_LEFT(CCR12_Info* info);

   void compute_amp(double*,const long,const long,const long,const long,const long,const long,const long,const long,const long);

};



}

#endif


