
#ifndef _math_scmat_local_h
#define _math_scmat_local_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/scmat/block.h>
#include <math/scmat/matrix.h>
#include <math/scmat/abstract.h>

class LocalSCDimension: public SCDimension {
#   define CLASSNAME LocalSCDimension
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    int n_;
  public:
    LocalSCDimension(int n);
    LocalSCDimension(const RefKeyVal&);
    LocalSCDimension(StateIn&);
    ~LocalSCDimension();
    void save_data_state(StateOut&);
    int n();
    SCMatrix* create_matrix(SCDimension*);
    SymmSCMatrix* create_symmmatrix();
    DiagSCMatrix* create_diagmatrix();
    SCVector* create_vector();
};
SavableState_REF_dec(LocalSCDimension);

class LocalSCVector: public SCVector {
    friend class LocalSCMatrix;
    friend class LocalSymmSCMatrix;
    friend class LocalDiagSCMatrix;
#   define CLASSNAME LocalSCVector
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    RefLocalSCDimension d;
    RefSCVectorSimpleBlock block;

    void resize(int);
  public:
    LocalSCVector();
    LocalSCVector(LocalSCDimension*);
    LocalSCVector(const RefKeyVal&);
    LocalSCVector(StateIn&);
    ~LocalSCVector();
    void save_data_state(StateOut&);
    void assign(double);
    void assign(SCVector*);
    void assign(const double*);

    RefSCDimension dim();
    void set_element(int,double);
    double get_element(int);
    void accumulate_product(SymmSCMatrix*,SCVector*);
    void accumulate_product(SCMatrix*,SCVector*);
    void accumulate(SCVector*);
    double scalar_product(SCVector*);
    void element_op(const RefSCVectorElementOp&);
    void print(const char* title=0,ostream& out=cout, int =10);
};

class LocalSCMatrix: public SCMatrix {
    friend class LocalSymmSCMatrix;
    friend class LocalDiagSCMatrix;
    friend LocalSCVector;
#   define CLASSNAME LocalSCMatrix
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    RefLocalSCDimension d1;
    RefLocalSCDimension d2;
    RefSCMatrixRectBlock block;
    double** rows;
  private:
    // utility functions
    int compute_offset(int,int);
    void resize(int,int);
  public:
    LocalSCMatrix();
    LocalSCMatrix(const RefKeyVal&);
    LocalSCMatrix(StateIn&);
    LocalSCMatrix(LocalSCDimension*,LocalSCDimension*);
    ~LocalSCMatrix();

    // implementations and overrides of virtual functions
    void save_data_state(StateOut&);
    RefSCDimension rowdim();
    RefSCDimension coldim();
    double get_element(int,int);
    void set_element(int,int,double);
    void accumulate_outer_product(SCVector*,SCVector*);
    void accumulate_product(SCMatrix*,SCMatrix*);
    void accumulate_product(SCMatrix*,SymmSCMatrix*);
    void accumulate_product(SCMatrix*,DiagSCMatrix*);
    void accumulate(SCMatrix*);
    void transpose_this();
    double invert_this();
    double solve_this(SCVector*);
    double determ_this();
    double trace();
    void gen_invert_this();
    void element_op(const RefSCRectElementOp&);
    void print(const char* title=0,ostream& out=cout, int =10);
};

class LocalSymmSCMatrix: public SymmSCMatrix {
    friend class LocalSCMatrix;
    friend class LocalDiagSCMatrix;
    friend LocalSCVector;
#   define CLASSNAME LocalSymmSCMatrix
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    RefLocalSCDimension d;
    RefSCMatrixLTriBlock block;
    double** rows;
  private:
    // utility functions
    int compute_offset(int,int);
    void resize(int n);
  public:
    LocalSymmSCMatrix();
    LocalSymmSCMatrix(const RefKeyVal&);
    LocalSymmSCMatrix(StateIn&);
    LocalSymmSCMatrix(LocalSCDimension*);
    ~LocalSymmSCMatrix();

    // implementations and overrides of virtual functions
    void save_data_state(StateOut&);
    RefSCDimension dim();
    double get_element(int,int);
    void set_element(int,int,double);
    void accumulate_product(SCMatrix*,SCMatrix*);
    void accumulate(SymmSCMatrix*);
    double invert_this();
    double solve_this(SCVector*);
    double trace();
    double determ_this();
    void gen_invert_this();

    double scalar_product(SCVector*);
    void diagonalize(DiagSCMatrix*,SCMatrix*);
    void accumulate_symmetric_outer_product(SCVector*);
    void accumulate_symmetric_product(SCMatrix*);
    void accumulate_symmetric_sum(SCMatrix*);
    void accumulate_transform(SCMatrix*,SymmSCMatrix*);
    void accumulate_transform(SCMatrix*,DiagSCMatrix*);
    void element_op(const RefSCSymmElementOp&);
    void print(const char* title=0,ostream& out=cout, int =10);
};

class LocalDiagSCMatrix: virtual public DiagSCMatrix {
    friend LocalSCMatrix;
    friend LocalSymmSCMatrix;
    friend LocalSCVector;
#   define CLASSNAME LocalDiagSCMatrix
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    RefLocalSCDimension d;
    RefSCMatrixDiagBlock block;
    void resize(int n);
  public:
    LocalDiagSCMatrix();
    LocalDiagSCMatrix(const RefKeyVal&);
    LocalDiagSCMatrix(StateIn&);
    LocalDiagSCMatrix(LocalSCDimension*);
    ~LocalDiagSCMatrix();

    // implementations and overrides of virtual functions
    void save_data_state(StateOut&);
    RefSCDimension dim();
    double get_element(int);
    void set_element(int,double);
    void accumulate(DiagSCMatrix*);
    double invert_this();
    double determ_this();
    double trace();
    void gen_invert_this();

    void element_op(const RefSCDiagElementOp&);
    void print(const char* title=0,ostream& out=cout, int =10);
};

#endif
