//
#ifndef _chemistry_qc_ccr12_ccsd_ccsd_pt_right_h
#define _chemistry_qc_ccr12_ccsd_ccsd_pt_right_h

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
#include <chemistry/qc/ccr12/ccsd.h>
#include <chemistry/qc/ccr12/parenthesis2tnum.h>


namespace sc {

class CCSD_PT_RIGHT : public Parenthesis2tNum {

  protected:

   void smith_0_1(double*, const long,const long,const long,const long,const long,const long);
   void smith_0_2(double*, const long,const long,const long,const long,const long,const long);


  public:
   CCSD_PT_RIGHT(CCR12_Info* info);
    
   void compute_amp(double*,const long,const long,const long,const long,const long,const long,const long);

};



}

#endif


