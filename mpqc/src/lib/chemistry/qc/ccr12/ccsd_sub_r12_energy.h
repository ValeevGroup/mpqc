//
#ifndef _chemistry_qc_ccr12_ccsd_sub_r12_energy_h
#define _chemistry_qc_ccr12_ccsd_sub_r12_energy_h

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

class CCSD_SUB_R12_ENERGY { 

  protected:

    CCR12_Info* z;

  public:
    CCSD_SUB_R12_ENERGY(CCR12_Info* info);
    ~CCSD_SUB_R12_ENERGY();

    void compute_amp(Ref<Tensor>&, Ref<Tensor>&, Ref<Tensor>&);

};



}

#endif


