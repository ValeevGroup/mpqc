
#ifndef _chemistry_qc_scf_hcore_h
#define _chemistry_qc_scf_hcore_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/scf/grscf.h>

class AccumHCore: public AccumDIH {
#   define CLASSNAME AccumHCore
#   define HAVE_STATEIN_CTOR
#   define HAVE_KEYVAL_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    AccumHCore(StateIn&);
    AccumHCore(const RefKeyVal&);
    ~AccumHCore();

    void save_data_state(StateOut&);
    
    void accum(const RefSymmSCMatrix& h);
};

#endif
