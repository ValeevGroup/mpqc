
#ifndef _math_scmat_matrix_h
#define _math_scmat_matrix_h
#ifdef __GNUC__
#pragma interface
#endif

#include <iostream.h>
#include <util/container/array.h>
#include <util/container/set.h>
#include <util/state/state.h>

class SCDimension;
class SCVector;
class SCMatrix;
class SymmSCMatrix;
class DiagSCMatrix;
class RefDiagSCMatrix;
class SCVectordouble;
class SCMatrixdouble;
class SymmSCMatrixdouble;
class DiagSCMatrixdouble;
class RefSCElementOp;
class SCMatrixBlockIter;
class SCMatrixRectBlock;
class SCMatrixLTriBlock;
class SCMatrixDiagBlock;
class SCVectorSimpleBlock;

class SCElementOp: virtual public SavableState {
#   define CLASSNAME SCElementOp
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SCElementOp();
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
SavableState_REF_dec(SCElementOp);

class SCRectElementOp: virtual public SCElementOp {
#   define CLASSNAME SCRectElementOp
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SCRectElementOp();
    ~SCRectElementOp();
};
SavableState_REF_dec(SCRectElementOp);

class SCDiagElementOp: virtual public SCElementOp {
#   define CLASSNAME SCDiagElementOp
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SCDiagElementOp();
    ~SCDiagElementOp();
};
SavableState_REF_dec(SCDiagElementOp);

class SCSymmElementOp: virtual public SCElementOp {
#   define CLASSNAME SCSymmElementOp
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SCSymmElementOp();
    ~SCSymmElementOp();
};
SavableState_REF_dec(SCSymmElementOp);

class SCVectorElementOp: virtual public SCElementOp {
#   define CLASSNAME SCVectorElementOp
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SCVectorElementOp();
    ~SCVectorElementOp();
};
SavableState_REF_dec(SCVectorElementOp);

class SCElementScale: virtual public SCDiagElementOp,
                      virtual public SCSymmElementOp,
                      virtual public SCRectElementOp,
                      virtual public SCVectorElementOp {
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

class SCElementAssign: virtual public SCDiagElementOp,
                       virtual public SCSymmElementOp,
                       virtual public SCRectElementOp,
                       virtual public SCVectorElementOp {
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

class SCElementSquareRoot: virtual public SCDiagElementOp,
                       virtual public SCSymmElementOp,
                       virtual public SCRectElementOp,
                       virtual public SCVectorElementOp {
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

class SCElementInvert: virtual public SCDiagElementOp,
                       virtual public SCSymmElementOp,
                       virtual public SCRectElementOp,
                       virtual public SCVectorElementOp {
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

class SCElementShiftDiagonal: virtual public SCDiagElementOp,
                       virtual public SCSymmElementOp,
                       virtual public SCRectElementOp,
                       virtual public SCVectorElementOp {
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

class SCElementMaxAbs: virtual public SCDiagElementOp,
                       virtual public SCSymmElementOp,
                       virtual public SCRectElementOp,
                       virtual public SCVectorElementOp {
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

SavableState_named_REF_dec(RefSSSCDimension,SCDimension);
class RefSCDimension: public RefSSSCDimension {
    // standard overrides
  public:
    RefSCDimension();
    RefSCDimension(RefSCDimension&);
    RefSCDimension(SCDimension *);
    RefSCDimension(RefDescribedClassBase&);
    ~RefSCDimension();
    RefSCDimension& operator=(SCDimension* cr);
    RefSCDimension& operator=( RefDescribedClassBase & c);
    RefSCDimension& operator=( RefSCDimension & c);
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
    RefSCVector(RefSCDimension&);
    RefSCVector(RefSCVector&);
    RefSCVector(SCVector *);
    // don't allow automatic conversion from any reference to a
    // described class
    //RefSCVector(RefDescribedClassBase&);
    ~RefSCVector();
    RefSCVector& operator=(SCVector* cr);
    //RefSCVector& operator=( RefDescribedClassBase & c);
    RefSCVector& operator=( RefSCVector & c);

    // vector specific members
  public:
    void set_element(int,double);
    double get_element(int);
    int n();
    RefSCDimension dim();
    SCVectordouble operator()(int);
    SCVectordouble operator[](int);
    RefSCVector operator+(RefSCVector&a);
    RefSCVector operator-(RefSCVector&a);
    RefSCVector operator*(double);
    RefSCVector clone();
    RefSCVector copy();
    RefSCMatrix outer_product(RefSCVector&);
    RefSymmSCMatrix symmetric_outer_product();
    double maxabs();
    double scalar_product(RefSCVector&);
    double dot(RefSCVector&);
    void normalize();
    void assign(RefSCVector&);
    void assign(double);
    void assign(const double*);
    void convert(double*);
    void scale(double);
    void accumulate(RefSCVector&);
    void element_op(RefSCVectorElementOp&);
    void print(ostream&out);
    void print(const char*title=0, ostream&out=cout, int precision=10);
};
RefSCVector operator*(double,RefSCVector&);
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
    RefSCMatrix(RefSCMatrix&);
    RefSCMatrix(SCMatrix *);
    //RefSCMatrix(RefDescribedClassBase&);
    ~RefSCMatrix();
    RefSCMatrix& operator=(SCMatrix* cr);
    //RefSCMatrix& operator=( RefDescribedClassBase & c);
    RefSCMatrix& operator=( RefSCMatrix & c);

    // matrix specific members
  public:
    RefSCMatrix(RefSCDimension&,RefSCDimension&);
    RefSCVector operator*(RefSCVector&);
    RefSCMatrix operator*(RefSCMatrix&);
    RefSCMatrix operator*(RefSymmSCMatrix&);
    RefSCMatrix operator*(RefDiagSCMatrix&);
    RefSCMatrix operator*(double);
    RefSCMatrix operator+(RefSCMatrix&);
    RefSCMatrix operator-(RefSCMatrix&);
    RefSCMatrix t();
    RefSCMatrix i();
    RefSCMatrix clone();
    RefSCMatrix copy();
    void set_element(int,int,double);
    double get_element(int,int);
    void accumulate_outer_product(RefSCVector&,RefSCVector&);
    void accumulate_product(RefSCMatrix&,RefSCMatrix&);
    void assign(RefSCMatrix&);
    void scale(double);
    void assign(double);
    void assign(const double*);
    void assign(const double**);
    void convert(double*);
    void convert(double**);
    void accumulate(RefSCMatrix&);
    void element_op(RefSCRectElementOp&);
    int nrow();
    int ncol();
    RefSCDimension rowdim();
    RefSCDimension coldim();
    SCMatrixdouble operator()(int i,int j);
    void print(ostream&);
    void print(const char*title=0,ostream&out=cout, int =10);
};
RefSCMatrix operator*(double,RefSCMatrix&);

SavableState_named_REF_dec(RefSSSymmSCMatrix,SymmSCMatrix);
class RefSymmSCMatrix: public RefSSSymmSCMatrix {
    // standard overrides
  public:
    RefSymmSCMatrix();
    RefSymmSCMatrix(StateIn&);
    RefSymmSCMatrix(RefSymmSCMatrix&);
    RefSymmSCMatrix(SymmSCMatrix *);
    //RefSymmSCMatrix(RefDescribedClassBase&);
    ~RefSymmSCMatrix();
    RefSymmSCMatrix& operator=(SymmSCMatrix* cr);
    //RefSymmSCMatrix& operator=( RefDescribedClassBase & c);
    RefSymmSCMatrix& operator=( RefSymmSCMatrix & c);

    // matrix specific members
  public:
    RefSymmSCMatrix(RefSCDimension&);
    RefSCMatrix operator*(RefSCMatrix&);
    RefSCVector operator*(RefSCVector&a);
    RefSymmSCMatrix operator*(double);
    RefSymmSCMatrix operator+(RefSymmSCMatrix&);
    RefSymmSCMatrix operator-(RefSymmSCMatrix&);
    RefSymmSCMatrix i();
    RefSymmSCMatrix clone();
    RefSymmSCMatrix copy();
    void set_element(int,int,double);
    double get_element(int,int);
    void accumulate_symmetric_outer_product(RefSCVector&);
    double scalar_product(RefSCVector&);
    void accumulate_symmetric_product(RefSCMatrix&);
    void accumulate_symmetric_sum(RefSCMatrix&);
    void accumulate_transform(RefSCMatrix&,RefSymmSCMatrix&);
    void assign(RefSymmSCMatrix&);
    void scale(double);
    void assign(double);
    void assign(const double*);
    void assign(const double**);
    void convert(double*);
    void convert(double**);
    void accumulate(RefSymmSCMatrix&);
    void element_op(RefSCSymmElementOp&);
    RefDiagSCMatrix eigvals();
    RefSCMatrix eigvecs();
    void diagonalize(RefDiagSCMatrix& eigvals,RefSCMatrix& eigvecs);
    int n();
    RefSCDimension dim();
    SymmSCMatrixdouble operator()(int i,int j);
    void print(ostream&);
    void print(const char*title=0,ostream&out=cout, int =10);
};
RefSymmSCMatrix operator*(double,RefSymmSCMatrix&);

SavableState_named_REF_dec(RefSSDiagSCMatrix,DiagSCMatrix);
class RefDiagSCMatrix: public RefSSDiagSCMatrix {
    // standard overrides
  public:
    RefDiagSCMatrix();
    RefDiagSCMatrix(StateIn&);
    RefDiagSCMatrix(RefDiagSCMatrix&);
    RefDiagSCMatrix(DiagSCMatrix *);
    //RefDiagSCMatrix(RefDescribedClassBase&);
    ~RefDiagSCMatrix();
    RefDiagSCMatrix& operator=(DiagSCMatrix* cr);
    //RefDiagSCMatrix& operator=( RefDescribedClassBase & c);
    RefDiagSCMatrix& operator=( RefDiagSCMatrix & c);

    // matrix specific members
  public:
    RefDiagSCMatrix(RefSCDimension&);
    RefSCMatrix operator*(RefSCMatrix&);
    RefDiagSCMatrix operator*(double);
    RefDiagSCMatrix operator+(RefDiagSCMatrix&);
    RefDiagSCMatrix operator-(RefDiagSCMatrix&);
    RefDiagSCMatrix i();
    RefDiagSCMatrix clone();
    RefDiagSCMatrix copy();
    void set_element(int,double);
    double get_element(int);
    void assign(RefDiagSCMatrix&);
    void scale(double);
    void assign(double);
    void assign(const double*);
    void convert(double*);
    void accumulate(RefDiagSCMatrix&);
    void element_op(RefSCDiagElementOp&);
    int n();
    RefSCDimension dim();
    DiagSCMatrixdouble operator()(int i);
    void print(ostream&);
    void print(const char*title=0,ostream&out=cout, int =10);
};
RefDiagSCMatrix operator*(double,RefDiagSCMatrix&);

class SCVectordouble {
   friend class RefSCVector;
  private:
    RefSCVector vector;
    int i;
    
    SCVectordouble(SCVector*,int);
  public:
    SCVectordouble(SCVectordouble&);
    ~SCVectordouble();
    double operator=(double a);
    operator double();
    double val();
};

class SCMatrixdouble {
   friend class RefSCMatrix;
  private:
    RefSCMatrix matrix;
    int i;
    int j;
    
    SCMatrixdouble(SCMatrix*,int,int);
  public:
    SCMatrixdouble(SCMatrixdouble&);
    ~SCMatrixdouble();
    double operator=(double a);
    operator double();
    double val();
};

class SymmSCMatrixdouble {
   friend class RefSymmSCMatrix;
  private:
    RefSymmSCMatrix matrix;
    int i;
    int j;
    
    SymmSCMatrixdouble(SymmSCMatrix*,int,int);
  public:
    SymmSCMatrixdouble(SCMatrixdouble&);
    ~SymmSCMatrixdouble();
    double operator=(double a);
    operator double();
    double val();
};

class DiagSCMatrixdouble {
   friend class RefDiagSCMatrix;
  private:
    RefDiagSCMatrix matrix;
    int i;
    int j;
    
    DiagSCMatrixdouble(DiagSCMatrix*,int,int);
  public:
    DiagSCMatrixdouble(SCMatrixdouble&);
    ~DiagSCMatrixdouble();
    double operator=(double a);
    operator double();
    double val();
};

#ifdef INLINE_FUNCTIONS
#include <math/scmat/matrix_i.h>
#endif

#endif
