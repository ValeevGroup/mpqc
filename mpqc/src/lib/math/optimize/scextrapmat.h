
#ifndef _math_optimize_scextrapmat_h
#define _math_optimize_scextrapmat_h

#include <math/optimize/scextrap.h>
#include <math/scmat/matrix.h>

class SymmSCMatrixSCExtrapData: public SCExtrapData {
#   define CLASSNAME SymmSCMatrixSCExtrapData
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    RefSymmSCMatrix m;
  public:
    SymmSCMatrixSCExtrapData(StateIn&);
    SymmSCMatrixSCExtrapData(const RefSymmSCMatrix&);

    void save_data_state(StateOut&);
    
    SCExtrapData* copy();
    void zero();
    void accumulate_scaled(double, const RefSCExtrapData&);
};

class SymmSCMatrix2SCExtrapData: public SCExtrapData {
#   define CLASSNAME SymmSCMatrix2SCExtrapData
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    RefSymmSCMatrix m1;
    RefSymmSCMatrix m2;
  public:
    SymmSCMatrix2SCExtrapData(StateIn&);
    SymmSCMatrix2SCExtrapData(const RefSymmSCMatrix&, const RefSymmSCMatrix&);

    void save_data_state(StateOut&);
    
    SCExtrapData* copy();
    void zero();
    void accumulate_scaled(double, const RefSCExtrapData&);
};

class SymmSCMatrixNSCExtrapData: public SCExtrapData {
#   define CLASSNAME SymmSCMatrixNSCExtrapData
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    int n_;
    RefSymmSCMatrix *m;
  public:
    SymmSCMatrixNSCExtrapData(StateIn&);
    SymmSCMatrixNSCExtrapData(int n, RefSymmSCMatrix*);

    void save_data_state(StateOut&);
    
    SCExtrapData* copy();
    void zero();
    void accumulate_scaled(double, const RefSCExtrapData&);
};

class SymmSCMatrixSCExtrapError: public SCExtrapError {
#   define CLASSNAME SymmSCMatrixSCExtrapError
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    RefSymmSCMatrix m;
  public:
    SymmSCMatrixSCExtrapError(StateIn&);
    SymmSCMatrixSCExtrapError(const RefSymmSCMatrix&);

    void save_data_state(StateOut&);
    
    double error();
    double scalar_product(const RefSCExtrapError&);
};

#endif
