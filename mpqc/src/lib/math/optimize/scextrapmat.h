
#ifndef _math_optimize_scextrapmat_h
#define _math_optimize_scextrapmat_h

#include <math/optimize/scextrap.h>
#include <math/scmat/matrix.h>

class SymmSCMatrixSCExtrapData: public SCExtrapData {
#   define CLASSNAME SymmSCMatrixSCExtrapData
#   include <util/class/classd.h>
  private:
    RefSymmSCMatrix m;
  public:
    SymmSCMatrixSCExtrapData(const RefSymmSCMatrix&);
    SCExtrapData* copy();
    void zero();
    void accumulate_scaled(double, const RefSCExtrapData&);
};

class SymmSCMatrix2SCExtrapData: public SCExtrapData {
#   define CLASSNAME SymmSCMatrix2SCExtrapData
#   include <util/class/classd.h>
  private:
    RefSymmSCMatrix m1;
    RefSymmSCMatrix m2;
  public:
    SymmSCMatrix2SCExtrapData(const RefSymmSCMatrix&, const RefSymmSCMatrix&);
    SCExtrapData* copy();
    void zero();
    void accumulate_scaled(double, const RefSCExtrapData&);
};

class SymmSCMatrixSCExtrapError: public SCExtrapError {
#   define CLASSNAME SymmSCMatrixSCExtrapError
#   include <util/class/classd.h>
  private:
    RefSymmSCMatrix m;
  public:
    SymmSCMatrixSCExtrapError(const RefSymmSCMatrix&);
    double error();
    double scalar_product(const RefSCExtrapError&);
};

#endif
