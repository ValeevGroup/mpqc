
#ifndef _math_scmat_abstract_h
#define _math_scmat_abstract_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/state/state.h>
#include <iostream.h>

class SCMatrix;
class SymmSCMatrix;
class DiagSCMatrix;
class SCVector;

class SSRefSCElementOp;
typedef class SSRefSCElementOp RefSCElementOp;

class SSRefSCElementOp2;
typedef class SSRefSCElementOp2 RefSCElementOp2;

class SSRefSCElementOp3;
typedef class SSRefSCElementOp3 RefSCElementOp3;

class RefSCDimension;

class SCDimension: public SavableState {
#   define CLASSNAME SCDimension
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SCDimension();
    SCDimension(StateIn&s): SavableState(s) {}
    virtual ~SCDimension();
    virtual int n() = 0;
    virtual SCMatrix* create_matrix(SCDimension*) = 0;
    SCMatrix* create_matrix(const RefSCDimension&);
    virtual SymmSCMatrix* create_symmmatrix() = 0;
    virtual DiagSCMatrix* create_diagmatrix() = 0;
    virtual SCVector* create_vector() = 0;
};

class SCVector: public SavableState {
#   define CLASSNAME SCVector
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SCVector();
    SCVector(StateIn&);

    // concrete functions (some can be overridden)
    virtual SCVector* copy();
    virtual SCVector* clone();
    virtual ~SCVector();
    void save_data_state(StateOut&);
    virtual int n();
    virtual double maxabs(); // maximum absolute value of the elements
    virtual void normalize();
    virtual void assign(double);
    virtual void assign(const double*);
    virtual void convert(double*);
    virtual void assign(SCVector*);
    virtual void scale(double);
    virtual void print(ostream&);

    virtual RefSCDimension dim() = 0;
    virtual void set_element(int,double) = 0;
    virtual double get_element(int) = 0;
    virtual void accumulate_product(SymmSCMatrix*,SCVector*) = 0;
    virtual void accumulate_product(SCMatrix*,SCVector*) = 0;
    virtual void accumulate(SCVector*) = 0;
    virtual double scalar_product(SCVector*) = 0;
    virtual void element_op(const RefSCElementOp&) = 0;
    virtual void element_op(const RefSCElementOp2&,
                            SCVector*) = 0;
    virtual void element_op(const RefSCElementOp3&,
                            SCVector*,SCVector*) = 0;
    virtual void print(const char* title=0,ostream& out=cout, int =10) = 0;
};

class SCMatrix: public SavableState {
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
    virtual void assign(const double*);
    virtual void assign(const double**);
    virtual void convert(double*);
    virtual void convert(double**);
    virtual void assign(SCMatrix*);
    virtual void scale(double);
    virtual void shift_diagonal(double);
    virtual void unit();
    virtual void print(ostream&);
    virtual SCMatrix* copy();
    virtual SCMatrix* clone();

    // pure virtual functions
    virtual RefSCDimension rowdim() = 0;
    virtual RefSCDimension coldim() = 0;
    virtual double get_element(int,int) = 0;
    virtual void set_element(int,int,double) = 0;
    virtual void accumulate_outer_product(SCVector*,SCVector*) = 0;
    virtual void accumulate_product(SCMatrix*,SCMatrix*) = 0;
    virtual void accumulate_product(SCMatrix*,SymmSCMatrix*) = 0;
    virtual void accumulate_product(SCMatrix*,DiagSCMatrix*) = 0;
    virtual void accumulate_product(SymmSCMatrix*,SCMatrix*);
    virtual void accumulate_product(DiagSCMatrix*,SCMatrix*);
    virtual void accumulate(SCMatrix*) = 0;
    virtual void transpose_this() = 0;
    virtual double trace() =0;
    virtual double invert_this() = 0;
    virtual double determ_this() = 0;
    virtual double solve_this(SCVector*) = 0;
    virtual void gen_invert_this() = 0;
    virtual void element_op(const RefSCElementOp&) = 0;
    virtual void element_op(const RefSCElementOp2&,
                            SCMatrix*) = 0;
    virtual void element_op(const RefSCElementOp3&,
                            SCMatrix*,SCMatrix*) = 0;
    virtual void print(const char* title=0,ostream& out=cout, int =10) = 0;
};

class SymmSCMatrix: public SavableState {
#   define CLASSNAME SymmSCMatrix
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SymmSCMatrix();
    SymmSCMatrix(StateIn&);
    void save_data_state(StateOut&);
    virtual double maxabs(); // maximum absolute value of the elements
    virtual void assign(double);
    virtual void assign(const double*);
    virtual void assign(const double**);
    virtual void convert(double*);
    virtual void convert(double**);
    virtual void assign(SymmSCMatrix*);
    virtual void scale(double);
    virtual void shift_diagonal(double);
    virtual void unit();
    virtual void print(ostream&);
    virtual int n();
    virtual SymmSCMatrix* copy();
    virtual SymmSCMatrix* clone();

    // pure virtual functions
    virtual RefSCDimension dim() = 0;
    virtual void diagonalize(DiagSCMatrix*,SCMatrix*) = 0;
    virtual void accumulate_symmetric_product(SCMatrix*) = 0;
    virtual void accumulate_symmetric_sum(SCMatrix*) = 0;
    virtual void accumulate_transform(SCMatrix*,SymmSCMatrix*) = 0;
    virtual void accumulate_transform(SCMatrix*,DiagSCMatrix*) = 0;
    virtual double get_element(int,int) = 0;
    virtual void set_element(int,int,double) = 0;
    virtual void accumulate(SymmSCMatrix*) = 0;
    virtual double trace() = 0;
    virtual double invert_this() = 0;
    virtual double determ_this() = 0;
    virtual double solve_this(SCVector*) = 0;
    virtual void gen_invert_this() = 0;
    virtual void element_op(const RefSCElementOp&) = 0;
    virtual void element_op(const RefSCElementOp2&,
                            SymmSCMatrix*) = 0;
    virtual void element_op(const RefSCElementOp3&,
                            SymmSCMatrix*,SymmSCMatrix*) = 0;
    virtual void print(const char* title=0,ostream& out=cout, int =10) = 0;
    virtual void accumulate_symmetric_outer_product(SCVector*) = 0;
    virtual double scalar_product(SCVector*) = 0;
};

class DiagSCMatrix: public SavableState {
#   define CLASSNAME DiagSCMatrix
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    DiagSCMatrix();
    DiagSCMatrix(StateIn&);
    void save_data_state(StateOut&);

    virtual double maxabs(); // maximum absolute value of the elements
    virtual void assign(double);
    virtual void assign(const double*);
    virtual void convert(double*);
    virtual void scale(double);
    virtual void assign(DiagSCMatrix*);
    virtual void print(ostream&);
    virtual int n();
    virtual DiagSCMatrix* copy();
    virtual DiagSCMatrix* clone();

    // pure virtual functions
    virtual RefSCDimension dim() = 0;
    virtual double get_element(int) = 0;
    virtual void set_element(int,double) = 0;
    virtual void accumulate(DiagSCMatrix*) = 0;
    virtual double trace() = 0;
    virtual double determ_this() = 0;
    virtual double invert_this() = 0;
    virtual void gen_invert_this() = 0;
    virtual void element_op(const RefSCElementOp&) = 0;
    virtual void element_op(const RefSCElementOp2&,
                            DiagSCMatrix*) = 0;
    virtual void element_op(const RefSCElementOp3&,
                            DiagSCMatrix*,DiagSCMatrix*) = 0;
    virtual void print(const char* title=0,ostream& out=cout, int =10) = 0;
};

#endif
