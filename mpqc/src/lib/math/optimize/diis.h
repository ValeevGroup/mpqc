
#ifndef _math_optimize_diis_h
#define _math_optimize_diis_h

#ifdef __GNUC__
#pragma interface
#endif

#include <stdio.h>
#include <math/array/math_lib.h>
#include <math/optimize/scextrap.h>

class DIIS: public SelfConsistentExtrapolation {
#   define CLASSNAME DIIS
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  protected:
    int start;
    int ndiis;
    int iter;
    double damping_factor;

    double_vector_t btemp;
    double_matrix_t bold;
    double_matrix_t bmat;

    RefSCExtrapData dtemp_data;
    RefSCExtrapError dtemp_error;

    RefSCExtrapData Ldata;

    RefSCExtrapData *diism_data;
    RefSCExtrapError *diism_error;

    void init();
  public:
    DIIS();
    DIIS(const RefKeyVal&);
    ~DIIS();

    int extrapolate(const RefSCExtrapData& data,
                    const RefSCExtrapError& error);
};

#endif
