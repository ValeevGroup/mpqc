
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

//texi The @code{SCDimension} class is used to determine the size and
// blocking of matrices.
class SCDimension: public SavableState {
#   define CLASSNAME SCDimension
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  private:
    char *name_;
  public:
    //texi Create a dimension with an optional name.  The
    // name is a copy of the @code{'\0'} terminated string @var{name}.
    SCDimension(const char* name = 0);

    SCDimension(StateIn&s);
    virtual ~SCDimension();
    void save_data_state(StateOut&);

    //texi Test to see if two dimensions are equivalent.
    virtual int equiv(SCDimension*) const = 0;
    
    //texi Return the dimension.
    virtual int n() = 0;
    //texi Create matrices or vectors.
    virtual SCMatrix* create_matrix(SCDimension*) = 0;
    SCMatrix* create_matrix(const RefSCDimension&);
    virtual SymmSCMatrix* create_symmmatrix() = 0;
    virtual DiagSCMatrix* create_diagmatrix() = 0;
    virtual SCVector* create_vector() = 0;
    //texi Return the name of the dimension.  If no name was given
    // to the constructor, then return @code{0}.
    const char* name() { return name_; }
};

//texi The @code{SCVector} class is the abstract base class for
// @code{double} valued vectors.
class SCVector: public SavableState {
#   define CLASSNAME SCVector
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SCVector();
    SCVector(StateIn&);

    // concrete functions (some can be overridden)
    //texi Return a vector with the same dimension and same elements.
    virtual SCVector* copy();
    //texi Return a vector with the same dimension but uninitialized memory.
    virtual SCVector* clone();

    virtual ~SCVector();
    void save_data_state(StateOut&);
    //texi Return the length of the vector.
    virtual int n();
    //texi Return the maximum absolute value element of this vector.
    virtual double maxabs();
    //texi Normalize this.
    virtual void normalize();
    //texi Assign each element to a random number between -1 and 1
    virtual void randomize();
    //texi Assign all elements of this to @var{val}.
    virtual void assign(double val);
    //texi Assign element @var{i} to @var{v[i]} for all @var{i}.
    virtual void assign(const double* v);
    //texi Assign @var{v[i]} to element @var{i} for all @var{i}.
    virtual void convert(double* v);
    //texi Make @code{this} have the same elements as @var{v}.
    // The dimensions must match.
    virtual void assign(SCVector* v);
    //texi Multiply each element by @var{val}.
    virtual void scale(double val);

    //texi Return the @code{RefSCDimension} corresponding to this vector.
    virtual RefSCDimension dim() = 0;
    //texi Set element @var{i} to @var{val}.
    virtual void set_element(int,double) = 0;
    //texi Add @var{val} to element @var{i}.
    virtual void accumulate_element(int,double) = 0;
    //texi Return the value of element @var{i}.
    virtual double get_element(int) = 0;
    //texi Sum the result of @var{m} times @var{v} into @code{this}.
    virtual void accumulate_product(SymmSCMatrix* m, SCVector* v) = 0;
    virtual void accumulate_product(SCMatrix* m, SCVector* v) = 0;
    //texi Sum @var{v} into this.
    virtual void accumulate(SCVector*v) = 0;
    //texi Return the dot product.
    virtual double scalar_product(SCVector*) = 0;
    //texi Perform the element operation @var{op} on each element of this.
    virtual void element_op(const RefSCElementOp&) = 0;
    virtual void element_op(const RefSCElementOp2&,
                            SCVector*) = 0;
    virtual void element_op(const RefSCElementOp3&,
                            SCVector*,SCVector*) = 0;
    //texi Print out the vector.
    virtual void print(ostream&);
    virtual void print(const char* title=0,ostream& out=cout, int =10) = 0;
};

//texi The @code{SCMatrix} class is the abstract base class for general
// @code{double} valued n by m matrices.
// For symmetric matrices use @code{SymmSCMatrix} and for
// diagonal matrices use @code{DiagSCMatrix}.
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

    //texi Return the number of rows or columns.
    virtual int nrow();
    virtual int ncol();
    //texi Return the maximum absolute value element.
    virtual double maxabs();
    //texi Assign each element to a random number between -1 and 1
    virtual void randomize();
    //texi Set all elements to @var{val}.
    virtual void assign(double val);
    //texi Assign element @var{i}, @var{j} to
    // @code{@var{m}[@var{i}*nrow()+@var{j}]}.
    virtual void assign(const double* m);
    //texi Assign element @var{i}, @var{j} to @code{@var{m}[@var{i}][@var{j}]}.
    virtual void assign(const double** m);
    //texi Like the @code{assign} members, but these write values
    // to the arguments.
    virtual void convert(double*);
    virtual void convert(double**);
    //texi Make @code{this} have the same elements as @var{m}.
    // The dimensions must match.
    virtual void assign(SCMatrix* m);
    //texi Multiply all elements by @var{val}.
    virtual void scale(double val);
    //texi Scale the diagonal elements by @var{val}.
    virtual void scale_diagonal(double val);
    //texi Shift the diagonal elements by @var{val}.
    virtual void shift_diagonal(double val);
    //texi Make @code{this} equal to the unit matrix.
    virtual void unit();
    //texi Return a matrix with the same dimension and same elements.
    virtual SCMatrix* copy();
    //texi Return a matrix with the same dimension but uninitialized memory.
    virtual SCMatrix* clone();

    // pure virtual functions
    //texi Return the row or column dimension.
    virtual RefSCDimension rowdim() = 0;
    virtual RefSCDimension coldim() = 0;
    //texi Return or modify an element.
    virtual double get_element(int,int) = 0;
    virtual void set_element(int,int,double) = 0;
    virtual void accumulate_element(int,int,double) = 0;
    
    //texi Return a subblock of @code{this}.  The subblock is defined as
    // the rows starting at @code{br} and ending at @code{er}, and the
    // columns beginning at @code{bc} and ending at @code{ec}.
    virtual SCMatrix * get_subblock(int br, int er, int bc, int ec) =0;

    //texi Assign @code{m} to a subblock of @code{this}.
    virtual void assign_subblock(SCMatrix *m, int, int, int, int, int=0, int=0) =0;

    //texi Sum @code{m} into a subblock of @code{this}.
    virtual void accumulate_subblock(SCMatrix *m, int, int, int, int, int=0,int=0) =0;
    
    //texi Return a row or column of @code{this}.
    virtual SCVector * get_row(int i) =0;
    virtual SCVector * get_column(int i) =0;

    //texi Assign @code{v} to a row or column of @code{this}.
    virtual void assign_row(SCVector *v, int i) =0;
    virtual void assign_column(SCVector *v, int i) =0;
    
    //texi Sum @code{v} to a row or column of @code{this}.
    virtual void accumulate_row(SCVector *v, int i) =0;
    virtual void accumulate_column(SCVector *v, int i) =0;
    
    //texi Sum @var{m} into this.
    virtual void accumulate(SCMatrix* m) = 0;
    //texi Sum into @code{this} the products of various vectors or matrices.
    virtual void accumulate_outer_product(SCVector*,SCVector*) = 0;
    virtual void accumulate_product(SCMatrix*,SCMatrix*) = 0;
    virtual void accumulate_product(SCMatrix*,SymmSCMatrix*) = 0;
    virtual void accumulate_product(SCMatrix*,DiagSCMatrix*) = 0;
    virtual void accumulate_product(SymmSCMatrix*,SCMatrix*);
    virtual void accumulate_product(DiagSCMatrix*,SCMatrix*);
    //texi Transpose @code{this}.
    virtual void transpose_this() = 0;
    //texi Return the trace.
    virtual double trace() =0;
    //texi Invert @code{this}.
    virtual double invert_this() = 0;
    //texi Return the determinant of @code{this}.  @code{this} is overwritten.
    virtual double determ_this() = 0;

    //texi Compute the singular value decomposition for @code{this},
    // possibly destroying this.
    virtual void svd_this(SCMatrix *U, DiagSCMatrix *sigma, SCMatrix *V);
    virtual double solve_this(SCVector*) = 0;
    virtual void gen_invert_this() = 0;

    //texi Schmidt orthogonalize @code{this}.  @code{S} is the overlap matrix.
    // @code{n} is the number of columns to orthogonalize.
    virtual void schmidt_orthog(SymmSCMatrix*, int n) =0;
    
    //texi Perform the element operation @var{op} on each element of this.
    virtual void element_op(const RefSCElementOp&) = 0;
    virtual void element_op(const RefSCElementOp2&,
                            SCMatrix*) = 0;
    virtual void element_op(const RefSCElementOp3&,
                            SCMatrix*,SCMatrix*) = 0;
    //texi Print out the matrix.
    virtual void print(ostream&);
    virtual void print(const char* title=0,ostream& out=cout, int =10) = 0;
};

//texi The @code{SymmSCMatrix} class is the abstract base class for symmetric
// @code{double} valued matrices.
class SymmSCMatrix: public SavableState {
#   define CLASSNAME SymmSCMatrix
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SymmSCMatrix();
    SymmSCMatrix(StateIn&);
    void save_data_state(StateOut&);
    //texi Return the maximum absolute value element of this vector.
    virtual double maxabs();
    //texi Assign each element to a random number between -1 and 1
    virtual void randomize();
    //texi Set all elements to @var{val}.
    virtual void assign(double val);
    //texi Assign element @var{i}, @var{j} to
    // @code{@var{m}[@var{i}*(@var{i}+1)/2+@var{j}]}.
    virtual void assign(const double* m);
    //texi Assign element @var{i}, @var{j} to @code{@var{m}[@var{i}][@var{j}]}.
    virtual void assign(const double** m);
    //texi Like the @code{assign} members, but these write values
    // to the arguments.
    virtual void convert(double*);
    virtual void convert(double**);
    //texi Make @code{this} have the same elements as @var{m}.
    // The dimensions must match.
    virtual void assign(SymmSCMatrix* m);
    //texi Multiply all elements by @var{val}.
    virtual void scale(double);
    //texi Scale the diagonal elements by @var{val}.
    virtual void scale_diagonal(double);
    //texi Shift the diagonal elements by @var{val}.
    virtual void shift_diagonal(double);
    //texi Make @code{this} equal to the unit matrix.
    virtual void unit();
    //texi Return the dimension.
    virtual int n();
    //texi Return a matrix with the same dimension and same elements.
    virtual SymmSCMatrix* copy();
    //texi Return a matrix with the same dimension but uninitialized memory.
    virtual SymmSCMatrix* clone();

    // pure virtual functions
    //texi Return the dimension.
    virtual RefSCDimension dim() = 0;
    //texi Return or modify an element.
    virtual double get_element(int,int) = 0;
    virtual void set_element(int,int,double) = 0;
    virtual void accumulate_element(int,int,double) = 0;

    //texi Return a subblock of @code{this}.  The subblock is defined as
    // the rows starting at @code{br} and ending at @code{er}, and the
    // columns beginning at @code{bc} and ending at @code{ec}.
    virtual SCMatrix * get_subblock(int br, int er, int bc, int ec) =0;
    virtual SymmSCMatrix * get_subblock(int br, int er) =0;

    //texi Assign @code{m} to a subblock of @code{this}.
    virtual void assign_subblock(SCMatrix *m, int, int, int, int) =0;
    virtual void assign_subblock(SymmSCMatrix *m, int, int) =0;

    //texi Sum @code{m} into a subblock of @code{this}.
    virtual void accumulate_subblock(SCMatrix *m, int, int, int, int) =0;
    virtual void accumulate_subblock(SymmSCMatrix *m, int, int) =0;

    //texi Return a row of @code{this}.
    virtual SCVector * get_row(int i) =0;

    //texi Assign @code{v} to a row of @code{this}.
    virtual void assign_row(SCVector *v, int i) =0;
    
    //texi Sum @code{v} to a row of @code{this}.
    virtual void accumulate_row(SCVector *v, int i) =0;

    //texi Diagonalize @code{this}, placing the eigenvalues in @var{d}
    // and the eigenvectors in @var{m}.
    virtual void diagonalize(DiagSCMatrix*d,SCMatrix*m) = 0;
    //texi Sum @var{m} into this.
    virtual void accumulate(SymmSCMatrix* m) = 0;
    //texi Sum into @code{this} the products of various vectors or matrices.
    virtual void accumulate_symmetric_product(SCMatrix*) = 0;
    virtual void accumulate_symmetric_sum(SCMatrix*) = 0;
    virtual void accumulate_transform(SCMatrix*,SymmSCMatrix*) = 0;
    virtual void accumulate_transform(SCMatrix*,DiagSCMatrix*) = 0;
    virtual void accumulate_symmetric_outer_product(SCVector*) = 0;
    //texi Return the scalar obtained by multiplying @code{this} on the
    // left and right by @var{v}.
    virtual double scalar_product(SCVector* v) = 0;
    //texi Return the trace.
    virtual double trace() = 0;
    //texi Invert @code{this}.
    virtual double invert_this() = 0;
    //texi Return the determinant of @code{this}.  @code{this} is overwritten.
    virtual double determ_this() = 0;

    virtual double solve_this(SCVector*) = 0;
    virtual void gen_invert_this() = 0;

    //texi Perform the element operation @var{op} on each element of this.
    virtual void element_op(const RefSCElementOp&) = 0;
    virtual void element_op(const RefSCElementOp2&,
                            SymmSCMatrix*) = 0;
    virtual void element_op(const RefSCElementOp3&,
                            SymmSCMatrix*,SymmSCMatrix*) = 0;
    //texi Print out the matrix.
    virtual void print(ostream&);
    virtual void print(const char* title=0,ostream& out=cout, int =10) = 0;
};

//texi The @code{SymmSCMatrix} class is the abstract base class for diagonal
// @code{double} valued matrices.
class DiagSCMatrix: public SavableState {
#   define CLASSNAME DiagSCMatrix
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    DiagSCMatrix();
    DiagSCMatrix(StateIn&);
    void save_data_state(StateOut&);

    //texi Return the maximum absolute value element of this vector.
    virtual double maxabs();
    //texi Assign each element to a random number between -1 and 1
    virtual void randomize();
    //texi Set all elements to @var{val}.
    virtual void assign(double val);
    //texi Assign element @var{i}, @var{i} to
    // @code{@var{m}[@var{i}]}.
    virtual void assign(const double*);
    //texi Like the @code{assign} member, but this write values
    // to the argument.
    virtual void convert(double*);
    //texi Make @code{this} have the same elements as @var{m}.
    // The dimensions must match.
    virtual void assign(DiagSCMatrix*);
    //texi Multiply all elements by @var{val}.
    virtual void scale(double);
    //texi Return the dimension.
    virtual int n();
    //texi Return a matrix with the same dimension and same elements.
    virtual DiagSCMatrix* copy();
    //texi Return a matrix with the same dimension but uninitialized memory.
    virtual DiagSCMatrix* clone();

    // pure virtual functions
    //texi Return the dimension.
    virtual RefSCDimension dim() = 0;
    //texi Return or modify an element.
    virtual double get_element(int) = 0;
    virtual void set_element(int,double) = 0;
    virtual void accumulate_element(int,double) = 0;
    //texi Sum @var{m} into this.
    virtual void accumulate(DiagSCMatrix* m) = 0;
    //texi Return the trace.
    virtual double trace() = 0;
    //texi Return the determinant of @code{this}.  @code{this} is overwritten.
    virtual double determ_this() = 0;
    //texi Invert @code{this}.
    virtual double invert_this() = 0;
    //texi Do a generalized inversion of @code{this}.
    virtual void gen_invert_this() = 0;
    //texi Perform the element operation @var{op} on each element of this.
    virtual void element_op(const RefSCElementOp&) = 0;
    virtual void element_op(const RefSCElementOp2&,
                            DiagSCMatrix*) = 0;
    virtual void element_op(const RefSCElementOp3&,
                            DiagSCMatrix*,DiagSCMatrix*) = 0;
    //texi Print out the matrix.
    virtual void print(ostream&);
    virtual void print(const char* title=0,ostream& out=cout, int =10) = 0;
};

//texi The @code{SCMatrixKit} class produces specialized matrices and
// dimensions.
class SCMatrixKit: public SavableState {
#   define CLASSNAME SCMatrixKit
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SCMatrixKit();
    SCMatrixKit(const RefKeyVal&);
    SCMatrixKit(StateIn&);
    ~SCMatrixKit();
    void save_data_state(StateOut&);

    // this member is default in local.cc
    //texi This returns a @code{LocalSCMatrixKit}.
    static SCMatrixKit* default_matrixkit();

    //texi Return a dimension with length @code{n} and, optionally,
    // name @var{name}.
    virtual SCDimension* dimension(int n, const char* name = 0) = 0;

    //texi Given the dimensions, create matrices or vectors.
    SCMatrix* matrix(const RefSCDimension&,const RefSCDimension&);
    SymmSCMatrix* symmmatrix(const RefSCDimension&);
    DiagSCMatrix* diagmatrix(const RefSCDimension&);
    SCVector* vector(const RefSCDimension&);

    //texi Given the dimensions and a @code{StateIn} object,
    // restore matrices or vectors.
    virtual SCMatrix* restore_matrix(StateIn&,
                                     const RefSCDimension&,
                                     const RefSCDimension&) = 0;
    virtual SymmSCMatrix* restore_symmmatrix(StateIn&,
                                             const RefSCDimension&) = 0;
    virtual DiagSCMatrix* restore_diagmatrix(StateIn&,
                                             const RefSCDimension&) = 0;
    virtual SCVector* restore_vector(StateIn&,
                                     const RefSCDimension&) = 0;
};

#endif
