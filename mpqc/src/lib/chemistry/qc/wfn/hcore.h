
#ifndef _chemistry_qc_wfn_hcore_h
#define _chemistry_qc_wfn_hcore_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/wfn/accum.h>

class AccumHCore: public AccumDIH {
#   define CLASSNAME AccumHCore
#   define HAVE_CTOR
#   define HAVE_STATEIN_CTOR
#   define HAVE_KEYVAL_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    AccumHCore();
    AccumHCore(StateIn&);
    AccumHCore(const RefKeyVal&);
    ~AccumHCore();

    void save_data_state(StateOut&);
    
    void accum(const RefSymmSCMatrix& h);
};
SavableState_REF_dec(AccumHCore);

#endif
