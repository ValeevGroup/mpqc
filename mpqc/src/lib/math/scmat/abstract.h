
#ifndef _math_scmat_abstract_h
#define _math_scmat_abstract_h

#include <math/scmat/matrix.h>

class SCDimension: virtual public SavableState {
#   define CLASSNAME SCDimension
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SCDimension();
    virtual ~SCDimension();
    virtual int n() = 0;
    virtual SCMatrix* create_matrix(SCDimension*) = 0;
    virtual SymmSCMatrix* create_symmmatrix() = 0;
    virtual DiagSCMatrix* create_diagmatrix() = 0;
    virtual SCVector* create_vector() = 0;
};

class SCVector: virtual public SavableState {
#   define CLASSNAME SCMatrix
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SCVector();
    SCVector(StateIn&);

    // concrete functions (some can be overridden)
    virtual ~SCVector();
    void save_data_state(StateOut&);
    virtual int n();
    virtual double maxabs(); // maximum absolute value of the elements
    virtual void assign(double);
    virtual void scale(double);
    virtual void copy(SCVector*);
    virtual void print(ostream&);

    virtual void set_element(int,double) = 0;
    virtual double get_element(int) = 0;
    virtual void accumulate_product(SCMatrix*,SCVector*) = 0;
    virtual void accumulate(SCVector*) = 0;
    virtual double scalar_product(SCVector*) = 0;
    virtual void element_op(RefSCVectorElementOp&) = 0;
    virtual void print(const char* title=0,ostream& out=cout, int =10) = 0;
};


@ This is the base class for all matrices.
@<|SCMatrix| declaration@>=
class RefSCElementOp;
class SCMatrix: virtual public SavableState {
#   define CLASSNAME SCMatrix
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:

    // concrete functions (some can be overridden)
    SCMatrix();
    SCMatrix(StateIn&);
    virtual ~SCMatrix();
    void save_data_state(StateOut&);
    virtual int nrow();
    virtual int ncol();
    virtual double maxabs(); // maximum absolute value of the elements
    virtual void assign(double);
    virtual void scale(double);
    virtual void shift_diagonal(double);
    virtual void unit();
    virtual void copy(SCMatrix*);
    virtual void print(ostream&);

    // pure virtual functions
    virtual RefSCDimension rowdim() = 0;
    virtual RefSCDimension coldim() = 0;
    virtual double get_element(int,int) = 0;
    virtual void set_element(int,int,double) = 0;
    virtual void accumulate_product(SCMatrix*,SCMatrix*) = 0;
    virtual void accumulate(SCMatrix*) = 0;
    virtual void transpose_this() = 0;
    virtual double invert_this() = 0;
    virtual void element_op(RefSCRectElementOp&) = 0;
    virtual void print(const char* title=0,ostream& out=cout, int =10) = 0;
};

class SymmSCMatrix: virtual public SavableState {
#   define CLASSNAME SymmSCMatrix
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SymmSCMatrix();
    SymmSCMatrix(StateIn&);
    virtual double maxabs(); // maximum absolute value of the elements
    virtual void assign(double);
    virtual void scale(double);
    virtual void shift_diagonal(double);
    virtual void unit();
    virtual void copy(SymmSCMatrix*);
    virtual void print(ostream&);

    // pure virtual functions
    virtual RefSCDimension dim() = 0;
    virtual void diagonalize(DiagSCMatrix*,SCMatrix*) = 0;
    virtual void accumulate_symmetric_product(SCMatrix*) = 0;
    virtual void accumulate_transform(SCMatrix*,SymmSCMatrix*) = 0;
    virtual double get_element(int,int) = 0;
    virtual void set_element(int,int,double) = 0;
    virtual void accumulate(SymmSCMatrix*) = 0;
    virtual double invert_this() = 0;
    virtual void element_op(RefSCSymmElementOp&) = 0;
    virtual void print(const char* title=0,ostream& out=cout, int =10) = 0;
};

class DiagSCMatrix: virtual public SavableState {
#   define CLASSNAME DiagSCMatrix
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    DiagSCMatrix();
    DiagSCMatrix(StateIn&);

    virtual double maxabs(); // maximum absolute value of the elements
    virtual void assign(double);
    virtual void scale(double);
    virtual void copy(DiagSCMatrix*);
    virtual void print(ostream&);

    // pure virtual functions
    virtual RefSCDimension dim() = 0;
    virtual double get_element(int,int) = 0;
    virtual void set_element(int,int,double) = 0;
    virtual void accumulate(DiagSCMatrix*) = 0;
    virtual double invert_this() = 0;
    virtual void element_op(RefSCDiagElementOp&) = 0;
    virtual void print(const char* title=0,ostream& out=cout, int =10) = 0;

    // implementations of parent routines
    void save_data_state(StateOut&);
};

#endif
