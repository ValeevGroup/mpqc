// hand written code by Toru Shiozaki
#ifndef _chemistry_qc_ccr12_ccsd_ccsd_pt_h
#define _chemistry_qc_ccr12_ccsd_ccsd_pt_h

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
#include <chemistry/qc/ccr12/ccr12.h>
#include <chemistry/qc/ccr12/ccsd_pt_left.h>
#include <chemistry/qc/ccr12/ccsd_pt_right.h>
#include <chemistry/qc/ccr12/parenthesis2t.h>

namespace sc {

class CCSD_PT : public Parenthesis2t {

  protected:

  public:
   CCSD_PT(CCR12_Info* info);
   
   double compute(Ref<CCSD_PT_LEFT>&, Ref<CCSD_PT_RIGHT>&);

};



}

#endif


