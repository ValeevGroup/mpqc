#ifndef _math_scmat_matrix_h
#define _math_scmat_matrix_h

#include <iostream.h>
#include <util/state/state.h>

class SCDimension;
SavableState_REF_dec(SCDimension);

class SCVectorSimpleBlock: public SCMatrixBlock {
#   define CLASSNAME SCVectorSimpleBlock
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    SCVectorSimpleBlock(int istart,int iend);
    SCVectorSimpleBlock(int istart,int iend);
    SCVectorSimpleBlock(StateIn&);
    virtual ~SCVectorSimpleBlock();
    void save_data_state(StateOut&);
    int istart;
    int iend;
    double* data;
};

SavableState_REF_dec(SCVectorSimpleBlock);

class SCMatrixBlock: virtual public SavableState {
#   define CLASSNAME SCMatrixBlock
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SCMatrixBlock();
    virtual ~SCMatrixBlock();
};

SavableState_REF_dec(SCMatrixBlock);

class SCMatrixRectBlock: public SCMatrixBlock {
#   define CLASSNAME SCMatrixRectBlock
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    SCMatrixRectBlock(int is, int ie, int js, int je);
    SCMatrixRectBlock(StateIn&);
    virtual ~SCMatrixRectBlock();
    void save_data_state(StateOut&);
    int istart;
    int jstart;
    int iend;
    int jend;
    double* data;
};
SavableState_REF_dec(SCMatrixRectBlock);

class SCMatrixLTriBlock: public SCMatrixBlock {
#   define CLASSNAME SCMatrixLTriBlock
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    SCMatrixLTriBlock(int s,int e);
    SCMatrixLTriBlock(StateIn&);
    virtual ~SCMatrixLTriBlock();
    void save_data_state(StateOut&);
    int start;
    int end;
    double* data;
};
SavableState_REF_dec(SCMatrixLTriBlock);

class SCMatrixDiagBlock: public SCMatrixBlock {
#   define CLASSNAME SCMatrixDiagBlock
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    SCMatrixDiagBlock(int istart,int iend,int jstart);
    SCMatrixDiagBlock(int istart,int iend);
    SCMatrixDiagBlock(StateIn&);
    virtual ~SCMatrixDiagBlock();
    void save_data_state(StateOut&);
    int istart;
    int jstart;
    int iend;
    double* data;
};
SavableState_REF_dec(SCMatrixDiagBlock);

class SCVector;
class SCVectordouble;
SavableState_named_REF_dec(RefSSSCVector,SCVector);
class RefSCVector: public RefSSSCVector {
    // standard overrides
  public:
    RefSCVector();
    RefSCVector(RefSCVector&);
    RefSCVector(SCVector *);
    RefSCVector(RefDescribedClassBase&);
    ~RefSCVector();
    RefSCVector& operator=(SCVector* cr);
    RefSCVector& operator=( RefDescribedClassBase & c);
    RefSCVector& operator=( RefSCVector & c);

    // vector specific members
  public:
}

class SCMatrixBlockIter {
  public:
    SCMatrixBlockIter();
    virtual ~SCMatrixBlockIter();
    virtual int i() = 0;
    virtual int j() = 0;
    virtual void set(double) = 0;
    virtual double get() = 0;
    virtual operator int() = 0;
    virtual void operator++() = 0; // prefix ++
    virtual void reset() = 0;
};

class RefSCElementOp;
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
};
SavableState_REF_dec(SCElementOp);

class SCRectElementOp: virtual public SCElementOp {
};
SavableState_REF_dec(SCRectElementOp);

class SCDiagElementOp: virtual public SCElementOp {
};
SavableState_REF_dec(SCDiagElementOp);

class SCSymmElementOp: virtual public SCElementOp {
};
SavableState_REF_dec(SCSymmElementOp);

class SCVectorElementOp: virtual public SCElementOp {
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

class SCMatrix;
class SCMatrixdouble;
SavableState_named_REF_dec(RefSSSCMatrix,SCMatrix);
class RefSCMatrix: public RefSSSCMatrix {
    // standard overrides
  public:
    RefSCMatrix();
    RefSCMatrix(RefSCMatrix&);
    RefSCMatrix(SCMatrix *);
    RefSCMatrix(RefDescribedClassBase&);
    ~RefSCMatrix();
    RefSCMatrix& operator=(SCMatrix* cr);
    RefSCMatrix& operator=( RefDescribedClassBase & c);
    RefSCMatrix& operator=( RefSCMatrix & c);

    // matrix specific members
  public:
    RefSCMatrix(RefSCDimension&,RefSCDimension&);
    RefSCMatrix operator*(RefSCMatrix&);
    RefSCMatrix operator+(RefSCMatrix&);
    RefSCMatrix operator-(RefSCMatrix&);
    RefSCMatrix t();
    RefSCMatrix i();
    RefSCMatrix clone();
    void set_element(int,int,double);
    double get_element(int,int);
    void accumulate_product(RefSCMatrix&,RefSCMatrix&);
    void copy(RefSCMatrix&);
    void scale(double);
    void assign(double);
    void accumulate(RefSCMatrix&);
    int nrow();
    int ncol();
    RefSCDimension rowdim();
    RefSCDimension coldim();
    SCMatrixdouble operator()(int i,int j);
    void print(ostream&);
    void print(const char*title=0,ostream&out=cout, int =10);
};

class SymmSCMatrix;
class RefDiagSCMatrix;
SavableState_named_REF_dec(RefSSSymmSCMatrix,SymmSCMatrix);
class RefSymmSCMatrix: public RefSSSymmSCMatrix {
    // standard overrides
  public:
    RefSymmSCMatrix();
    RefSymmSCMatrix(RefSymmSCMatrix&);
    RefSymmSCMatrix(SymmSCMatrix *);
    RefSymmSCMatrix(RefDescribedClassBase&);
    ~RefSymmSCMatrix();
    RefSymmSCMatrix& operator=(SymmSCMatrix* cr);
    RefSymmSCMatrix& operator=( RefDescribedClassBase & c);
    RefSymmSCMatrix& operator=( RefSymmSCMatrix & c);

    // matrix specific members
  public:
    RefSymmSCMatrix(RefSCDimension&);
    RefSymmSCMatrix operator+(RefSymmSCMatrix&);
    RefSymmSCMatrix operator-(RefSymmSCMatrix&);
    RefSymmSCMatrix i();
    RefSymmSCMatrix clone();
    void set_element(int,int,double);
    double get_element(int,int);
    void accumulate_symmetric_product(RefSymmSCMatrix&);
    void copy(RefSymmSCMatrix&);
    void scale(double);
    void assign(double);
    void accumulate(RefSymmSCMatrix&);
    RefDiagSCMatrix eigvals();
    RefSCMatrix eigvecs();
    void diagonalize(RefDiagSCMatrix& eigvals,RefSCMatrix& eigvecs);
    int n();
    RefSCDimension dim();
    SCMatrixdouble operator()(int i,int j);
    void print(ostream&);
    void print(const char*title=0,ostream&out=cout, int =10);
};

class DiagSCMatrix;
SavableState_named_REF_dec(RefSSDiagSCMatrix,DiagSCMatrix);
class RefDiagSCMatrix: public RefSSDiagSCMatrix {
    // standard overrides
  public:
    RefDiagSCMatrix();
    RefDiagSCMatrix(RefDiagSCMatrix&);
    RefDiagSCMatrix(DiagSCMatrix *);
    RefDiagSCMatrix(RefDescribedClassBase&);
    ~RefDiagSCMatrix();
    RefDiagSCMatrix& operator=(DiagSCMatrix* cr);
    RefDiagSCMatrix& operator=( RefDescribedClassBase & c);
    RefDiagSCMatrix& operator=( RefDiagSCMatrix & c);

    // matrix specific members
  public:
    RefDiagSCMatrix(RefSCDimension&);
    RefDiagSCMatrix operator+(RefDiagSCMatrix&);
    RefDiagSCMatrix operator-(RefDiagSCMatrix&);
    RefDiagSCMatrix i();
    RefDiagSCMatrix clone();
    void set_element(int,double);
    double get_element(int);
    void copy(RefDiagSCMatrix&);
    void scale(double);
    void assign(double);
    void accumulate(RefDiagSCMatrix&);
    int n();
    RefSCDimension dim();
    SCMatrixdouble operator()(int i);
    SCMatrixdouble operator()(int i,int j);
    void print(ostream&);
    void print(const char*title=0,ostream&out=cout, int =10);
};

class SCVectordouble {
   friend class RefSCVector;
   friend class RefSymmSCVector;
   friend class RefDiagSCVector;
  private:
    RefSCVector vector;
    int i;
    
    SCVectordouble(SCVector*,int,int);
  public:
    SCVectordouble(SCVectordouble&);
    ~SCVectordouble();
    double operator=(double a);
    operator double();
    double val();
};

class SCMatrixdouble {
   friend class RefSCMatrix;
   friend class RefSymmSCMatrix;
   friend class RefDiagSCMatrix;
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

#ifdef INLINE_FUNCTIONS
#include <math/scmat/matrix_i.h>
#endif

#endif
