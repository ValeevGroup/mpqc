
#ifndef _math_scmat_matrix_h
#define _math_scmat_matrix_h
#ifdef __GNUC__
#pragma interface
#endif

#include <iostream.h>
#include <util/container/array.h>
#include <util/container/set.h>
#include <util/state/state.h>

#include <math/scmat/abstract.h>

class SCVectordouble;
class SCMatrixdouble;
class SymmSCMatrixdouble;
class DiagSCMatrixdouble;

class SCMatrixBlockIter;
class SCMatrixRectBlock;
class SCMatrixLTriBlock;
class SCMatrixDiagBlock;
class SCVectorSimpleBlock;

class SCElementOp: public SavableState {
#   define CLASSNAME SCElementOp
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SCElementOp();
    SCElementOp(StateIn&s): SavableState(s) {}
    virtual ~SCElementOp();
    // If duplicates of |SCElementOp| are made then if |has_collect| returns
    // true then collect is called in such a way that each duplicated
    // |SCElementOp| is the argument of |collect| once.  The default
    // return value of |has_collect| is 0 and |collect|'s default action
    // is do nothing.
    virtual int has_collect();
    virtual void collect(RefSCElementOp&);
    virtual void process(SCMatrixBlockIter&) = 0;
    virtual void process(SCMatrixRectBlock*);
    virtual void process(SCMatrixLTriBlock*);
    virtual void process(SCMatrixDiagBlock*);
    virtual void process(SCVectorSimpleBlock*);
};
DCRef_declare(SCElementOp);
SSRef_declare(SCElementOp);

class SCElementOp2: public SavableState {
#   define CLASSNAME SCElementOp2
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SCElementOp2();
    SCElementOp2(StateIn&s): SavableState(s) {}
    virtual ~SCElementOp2();
    // If duplicates of |SCElementOp2| are made then if |has_collect| returns
    // true then collect is called in such a way that each duplicated
    // |SCElementOp| is the argument of |collect| once.  The default
    // return value of |has_collect| is 0 and |collect|'s default action
    // is do nothing.
    virtual int has_collect();
    virtual void collect(RefSCElementOp2&);
    virtual void process(SCMatrixBlockIter&,SCMatrixBlockIter&) = 0;
    virtual void process(SCMatrixRectBlock*,SCMatrixRectBlock*);
    virtual void process(SCMatrixLTriBlock*,SCMatrixLTriBlock*);
    virtual void process(SCMatrixDiagBlock*,SCMatrixDiagBlock*);
    virtual void process(SCVectorSimpleBlock*,SCVectorSimpleBlock*);
};
DCRef_declare(SCElementOp2);
SSRef_declare(SCElementOp2);

class SCElementOp3: public SavableState {
#   define CLASSNAME SCElementOp3
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SCElementOp3();
    SCElementOp3(StateIn&s): SavableState(s) {}
    virtual ~SCElementOp3();
    // If duplicates of |SCElementOp3| are made then if |has_collect| returns
    // true then collect is called in such a way that each duplicated
    // |SCElementOp| is the argument of |collect| once.  The default
    // return value of |has_collect| is 0 and |collect|'s default action
    // is do nothing.
    virtual int has_collect();
    virtual void collect(RefSCElementOp3&);
    virtual void process(SCMatrixBlockIter&,
                         SCMatrixBlockIter&,
                         SCMatrixBlockIter&) = 0;
    virtual void process(SCMatrixRectBlock*,
                         SCMatrixRectBlock*,
                         SCMatrixRectBlock*);
    virtual void process(SCMatrixLTriBlock*,
                         SCMatrixLTriBlock*,
                         SCMatrixLTriBlock*);
    virtual void process(SCMatrixDiagBlock*,
                         SCMatrixDiagBlock*,
                         SCMatrixDiagBlock*);
    virtual void process(SCVectorSimpleBlock*,
                         SCVectorSimpleBlock*,
                         SCVectorSimpleBlock*);
};
DCRef_declare(SCElementOp3);
SSRef_declare(SCElementOp3);

class SCElementScalarProduct: public SCElementOp2 {
#   define CLASSNAME SCElementScalarProduct
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    double product;
  public:
    SCElementScalarProduct();
    SCElementScalarProduct(StateIn&);
    ~SCElementScalarProduct();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&,SCMatrixBlockIter&);
    int has_collect();
    void collect(RefSCElementOp&);
    double result();
};
SavableState_REF_dec(SCElementScalarProduct);

class SCDestructiveElementProduct: public SCElementOp2 {
#   define CLASSNAME SCDestructiveElementProduct
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    SCDestructiveElementProduct();
    SCDestructiveElementProduct(StateIn&);
    ~SCDestructiveElementProduct();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&,SCMatrixBlockIter&);
};

class SCElementScale: public SCElementOp {
#   define CLASSNAME SCElementScale
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    double scale;
  public:
    SCElementScale(double a);
    SCElementScale(StateIn&);
    ~SCElementScale();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&);
};

class SCElementAssign: public SCElementOp {
#   define CLASSNAME SCElementAssign
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    double assign;
  public:
    SCElementAssign(double a);
    SCElementAssign(StateIn&);
    ~SCElementAssign();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&);
};

class SCElementSquareRoot: public SCElementOp {
#   define CLASSNAME SCElementSquareRoot
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    SCElementSquareRoot();
    SCElementSquareRoot(double a);
    SCElementSquareRoot(StateIn&);
    ~SCElementSquareRoot();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&);
};

class SCElementInvert: public SCElementOp {
#   define CLASSNAME SCElementInvert
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    SCElementInvert();
    SCElementInvert(double a);
    SCElementInvert(StateIn&);
    ~SCElementInvert();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&);
};

class SCElementShiftDiagonal: public SCElementOp {
#   define CLASSNAME SCElementShiftDiagonal
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    double shift_diagonal;
  public:
    SCElementShiftDiagonal(double a);
    SCElementShiftDiagonal(StateIn&);
    ~SCElementShiftDiagonal();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&);
};

class SCElementMaxAbs: public SCElementOp {
#   define CLASSNAME SCElementMaxAbs
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    double r;
  public:
    SCElementMaxAbs();
    SCElementMaxAbs(StateIn&);
    ~SCElementMaxAbs();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&);
    int has_collect();
    void collect(RefSCElementOp&);
    double result();
};
SavableState_REF_dec(SCElementMaxAbs);

DCRef_declare(SCDimension);
SSRef_declare(SCDimension);
class RefSCDimension: public SSRefSCDimension {
    // standard overrides
  public:
    RefSCDimension();
    RefSCDimension(const RefSCDimension&);
    RefSCDimension(SCDimension *);
    RefSCDimension(const RefDescribedClassBase&);
    ~RefSCDimension();
    RefSCDimension& operator=(SCDimension* cr);
    RefSCDimension& operator=(const RefDescribedClassBase & c);
    RefSCDimension& operator=(const RefSCDimension & c);
    operator int();

    // dimension specific functions
  public:
    int n();
};

class RefSCMatrix;
class RefSymmSCMatrix;
SavableState_named_REF_dec(RefSSSCVector,SCVector);
class RefSCVector: public RefSSSCVector {
    // standard overrides
  public:
    RefSCVector();
    RefSCVector(StateIn&);
    RefSCVector(const RefSCDimension&);
    RefSCVector(const RefSCVector&);
    RefSCVector(SCVector *);
    // don't allow automatic conversion from any reference to a
    // described class
    //RefSCVector(RefDescribedClassBase&);
    ~RefSCVector();
    RefSCVector& operator=(SCVector* cr);
    //RefSCVector& operator=( RefDescribedClassBase & c);
    RefSCVector& operator=(const RefSCVector & c);

    // vector specific members
  public:
    void set_element(int,double) const;
    double get_element(int) const;
    int n() const;
    RefSCDimension dim() const;
    SCVectordouble operator()(int) const;
    SCVectordouble operator[](int) const;
    RefSCVector operator+(const RefSCVector&a) const;
    RefSCVector operator-(const RefSCVector&a) const;
    RefSCVector operator*(double) const;
    RefSCVector clone() const;
    RefSCVector copy() const;
    RefSCMatrix outer_product(const RefSCVector&) const;
    RefSymmSCMatrix symmetric_outer_product() const;
    double maxabs() const;
    double scalar_product(const RefSCVector&) const;
    double dot(const RefSCVector&) const;
    void normalize() const;
    void assign(const RefSCVector&) const;
    void assign(double) const;
    void assign(const double*) const;
    void convert(double*) const;
    void scale(double) const;
    void accumulate(const RefSCVector&) const;
    void element_op(const RefSCElementOp&) const;
    void element_op(const RefSCElementOp2&,
                    const RefSCVector&) const;
    void element_op(const RefSCElementOp3&,
                    const RefSCVector&,
                    const RefSCVector&) const;
    void print(ostream&out) const;
    void print(const char*title=0, ostream&out=cout, int precision=10) const;
};
RefSCVector operator*(double,const RefSCVector&);
ARRAY_dec(RefSCVector);
SET_dec(RefSCVector);

class RefSymmSCMatrix;
class RefDiagSCMatrix;
SavableState_named_REF_dec(RefSSSCMatrix,SCMatrix);
class RefSCMatrix: public RefSSSCMatrix {
    // standard overrides
  public:
    RefSCMatrix();
    RefSCMatrix(StateIn&);
    RefSCMatrix(const RefSCMatrix&);
    RefSCMatrix(SCMatrix *);
    //RefSCMatrix(RefDescribedClassBase&);
    ~RefSCMatrix();
    RefSCMatrix& operator=(SCMatrix* cr);
    //RefSCMatrix& operator=( RefDescribedClassBase & c);
    RefSCMatrix& operator=(const RefSCMatrix & c);

    // matrix specific members
  public:
    RefSCMatrix(const RefSCDimension&,const RefSCDimension&);
    RefSCVector operator*(const RefSCVector&) const;
    RefSCMatrix operator*(const RefSCMatrix&) const;
    RefSCMatrix operator*(const RefSymmSCMatrix&) const;
    RefSCMatrix operator*(const RefDiagSCMatrix&) const;
    RefSCMatrix operator*(double) const;
    RefSCMatrix operator+(const RefSCMatrix&) const;
    RefSCMatrix operator-(const RefSCMatrix&) const;
    RefSCMatrix t() const;
    RefSCMatrix i() const;
    RefSCMatrix gi() const;
    RefSCMatrix clone() const;
    RefSCMatrix copy() const;
    void set_element(int,int,double) const;
    double get_element(int,int) const;
    void accumulate_outer_product(const RefSCVector&,const RefSCVector&) const;
    void accumulate_product(const RefSCMatrix&,const RefSCMatrix&) const;
    void assign(const RefSCMatrix&) const;
    void scale(double) const;
    void assign(double) const;
    void assign(const double*) const;
    void assign(const double**) const;
    void convert(double*) const;
    void convert(double**) const;
    void accumulate(const RefSCMatrix&) const;
    void element_op(const RefSCElementOp&) const;
    void element_op(const RefSCElementOp2&,
                    const RefSCMatrix&) const;
    void element_op(const RefSCElementOp3&,
                    const RefSCMatrix&,
                    const RefSCMatrix&) const;
    int nrow() const;
    int ncol() const;
    double solve_lin(const RefSCVector&) const;
    double determ() const;
    double trace() const;
    RefSCDimension rowdim() const;
    RefSCDimension coldim() const;
    SCMatrixdouble operator()(int i,int j) const;
    void print(ostream&) const;
    void print(const char*title=0,ostream&out=cout, int =10) const;
};
RefSCMatrix operator*(double,RefSCMatrix&);

SavableState_named_REF_dec(RefSSSymmSCMatrix,SymmSCMatrix);
class RefSymmSCMatrix: public RefSSSymmSCMatrix {
    // standard overrides
  public:
    RefSymmSCMatrix();
    RefSymmSCMatrix(StateIn&);
    RefSymmSCMatrix(const RefSymmSCMatrix&);
    RefSymmSCMatrix(SymmSCMatrix *);
    //RefSymmSCMatrix(RefDescribedClassBase&);
    ~RefSymmSCMatrix();
    RefSymmSCMatrix& operator=(SymmSCMatrix* cr);
    //RefSymmSCMatrix& operator=( RefDescribedClassBase & c);
    RefSymmSCMatrix& operator=(const RefSymmSCMatrix & c);

    // matrix specific members
  public:
    RefSymmSCMatrix(const RefSCDimension&);
    RefSCMatrix operator*(const RefSCMatrix&) const;
    RefSCVector operator*(const RefSCVector&a) const;
    RefSymmSCMatrix operator*(double) const;
    RefSymmSCMatrix operator+(const RefSymmSCMatrix&) const;
    RefSymmSCMatrix operator-(const RefSymmSCMatrix&) const;
    RefSymmSCMatrix i() const;
    RefSymmSCMatrix gi() const;
    RefSymmSCMatrix clone() const;
    RefSymmSCMatrix copy() const;
    void set_element(int,int,double) const;
    double get_element(int,int) const;
    void accumulate_symmetric_outer_product(const RefSCVector&) const;
    double scalar_product(const RefSCVector&) const;
    void accumulate_symmetric_product(const RefSCMatrix&) const;
    void accumulate_symmetric_sum(const RefSCMatrix&) const;
    void accumulate_transform(const RefSCMatrix&,const RefSymmSCMatrix&) const;
    void accumulate_transform(const RefSCMatrix&,const RefDiagSCMatrix&) const;
    void assign(const RefSymmSCMatrix&) const;
    void scale(double) const;
    void assign(double) const;
    void assign(const double*) const;
    void assign(const double**) const;
    void convert(double*) const;
    void convert(double**) const;
    void accumulate(const RefSymmSCMatrix&) const;
    void element_op(const RefSCElementOp&) const;
    void element_op(const RefSCElementOp2&,
                    const RefSymmSCMatrix&) const;
    void element_op(const RefSCElementOp3&,
                    const RefSymmSCMatrix&,
                    const RefSymmSCMatrix&) const;
    double solve_lin(const RefSCVector&) const;
    double trace() const;
    double determ() const;
    RefDiagSCMatrix eigvals() const;
    RefSCMatrix eigvecs() const;
    void diagonalize(const RefDiagSCMatrix& eigvals,
                     const RefSCMatrix& eigvecs) const;
    int n() const;
    RefSCDimension dim() const;
    SymmSCMatrixdouble operator()(int i,int j) const;
    void print(ostream&) const;
    void print(const char*title=0,ostream&out=cout, int =10) const;
};
RefSymmSCMatrix operator*(double,const RefSymmSCMatrix&);

SavableState_named_REF_dec(RefSSDiagSCMatrix,DiagSCMatrix);
class RefDiagSCMatrix: public RefSSDiagSCMatrix {
    // standard overrides
  public:
    RefDiagSCMatrix();
    RefDiagSCMatrix(StateIn&);
    RefDiagSCMatrix(const RefDiagSCMatrix&);
    RefDiagSCMatrix(DiagSCMatrix *);
    //RefDiagSCMatrix(RefDescribedClassBase&);
    ~RefDiagSCMatrix();
    RefDiagSCMatrix& operator=(DiagSCMatrix* cr);
    //RefDiagSCMatrix& operator=( RefDescribedClassBase & c);
    RefDiagSCMatrix& operator=(const RefDiagSCMatrix & c);

    // matrix specific members
  public:
    RefDiagSCMatrix(const RefSCDimension&);
    RefSCMatrix operator*(const RefSCMatrix&) const;
    RefDiagSCMatrix operator*(double) const;
    RefDiagSCMatrix operator+(const RefDiagSCMatrix&) const;
    RefDiagSCMatrix operator-(const RefDiagSCMatrix&) const;
    RefDiagSCMatrix i() const;
    RefDiagSCMatrix gi() const;
    RefDiagSCMatrix clone() const;
    RefDiagSCMatrix copy() const;
    void set_element(int,double) const;
    double get_element(int) const;
    void assign(const RefDiagSCMatrix&) const;
    void scale(double) const;
    void assign(double) const;
    void assign(const double*) const;
    void convert(double*) const;
    void accumulate(const RefDiagSCMatrix&) const;
    void element_op(const RefSCElementOp&) const;
    void element_op(const RefSCElementOp2&,
                    const RefDiagSCMatrix&) const;
    void element_op(const RefSCElementOp3&,
                    const RefDiagSCMatrix&,
                    const RefDiagSCMatrix&) const;
    int n() const;
    double trace() const;
    double determ() const;
    RefSCDimension dim() const;
    DiagSCMatrixdouble operator()(int i) const;
    void print(ostream&) const;
    void print(const char*title=0,ostream&out=cout, int =10) const;
};
RefDiagSCMatrix operator*(double,const RefDiagSCMatrix&);

class SCVectordouble {
   friend class RefSCVector;
  private:
    RefSCVector vector;
    int i;
    
    SCVectordouble(SCVector*,int);
  public:
    SCVectordouble(const SCVectordouble&);
    ~SCVectordouble();
    double operator=(double a);
    double operator=(const SCVectordouble&);
    operator double();
    double val() const;
};

class SCMatrixdouble {
   friend class RefSCMatrix;
  private:
    RefSCMatrix matrix;
    int i;
    int j;
    
    SCMatrixdouble(SCMatrix*,int,int);
  public:
    SCMatrixdouble(const SCMatrixdouble&);
    ~SCMatrixdouble();
    double operator=(double a);
    double operator=(const SCMatrixdouble&);
    operator double();
    double val() const;
};

class SymmSCMatrixdouble {
   friend class RefSymmSCMatrix;
  private:
    RefSymmSCMatrix matrix;
    int i;
    int j;
    
    SymmSCMatrixdouble(SymmSCMatrix*,int,int);
  public:
    SymmSCMatrixdouble(const SCMatrixdouble&);
    ~SymmSCMatrixdouble();
    double operator=(double a);
    double operator=(const SymmSCMatrixdouble&);
    operator double();
    double val() const;
};

class DiagSCMatrixdouble {
   friend class RefDiagSCMatrix;
  private:
    RefDiagSCMatrix matrix;
    int i;
    int j;
    
    DiagSCMatrixdouble(DiagSCMatrix*,int,int);
  public:
    DiagSCMatrixdouble(const SCMatrixdouble&);
    ~DiagSCMatrixdouble();
    double operator=(double a);
    double operator=(const DiagSCMatrixdouble&);
    operator double();
    double val() const;
};

#ifdef INLINE_FUNCTIONS
#include <math/scmat/matrix_i.h>
#endif

#endif
