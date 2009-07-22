//
#ifndef _chemistry_qc_ccr12_ccsd_ccsdt_h
#define _chemistry_qc_ccr12_ccsd_ccsdt_h

#ifdef __GNUC__
#pragma interface
#endif

#include <scconfig.h>
#include <scdirlist.h>
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

//class R12IntEval;
class R12IntEvalInfo;

class CCSDT: public CCR12 {

  public:
    CCSDT(StateIn&);
    CCSDT(const Ref<KeyVal>&);
    ~CCSDT();
    void compute();

};

}

#endif


