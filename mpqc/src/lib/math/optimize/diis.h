
#ifndef _math_optimize_diis_h
#define _math_optimize_diis_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/optimize/scextrap.h>

class DIIS: public SelfConsistentExtrapolation {
#   define CLASSNAME DIIS
#   define HAVE_STATEIN_CTOR
#   define HAVE_KEYVAL_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    int start;
    int ndiis;
    int iter;
    double damping_factor;

    double * btemp;
    double ** bold;
    double ** bmat;

    RefSCExtrapData dtemp_data;
    RefSCExtrapError dtemp_error;

    RefSCExtrapData Ldata;

    RefSCExtrapData *diism_data;
    RefSCExtrapError *diism_error;

    void init();
  public:
    DIIS();
    DIIS(StateIn&);
    DIIS(const RefKeyVal&);
    ~DIIS();

    void save_data_state(StateOut&);
    
    int extrapolate(const RefSCExtrapData& data,
                    const RefSCExtrapError& error);

    void reinitialize();
};

#endif
