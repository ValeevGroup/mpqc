
#ifndef _math_scmat_abstract_h
#define _math_scmat_abstract_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/group/message.h>

#include <util/state/state.h>
#include <math/scmat/dim.h>
#include <math/scmat/block.h>
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

DescribedClass_REF_fwddec(SCMatrixKit);

class SCMatrixKit: public DescribedClass {
#   define CLASSNAME SCMatrixKit
#   include <util/class/classda.h>
  protected:
    RefMessageGrp grp_;
    
  public:
    SCMatrixKit();
    SCMatrixKit(const RefKeyVal&);
    ~SCMatrixKit();

    // these members are default in local.cc
    //. This returns a \clsnmref{LocalSCMatrixKit}, unless the
    //. default has been changed with \srccd{set\_default\_matrixkit}.
    static SCMatrixKit* default_matrixkit();
    static void set_default_matrixkit(const RefSCMatrixKit &);

    RefMessageGrp messagegrp() const;

    //. Given the dimensions, create matrices or vectors.
    virtual SCMatrix* matrix(const RefSCDimension&,const RefSCDimension&) = 0;
    virtual SymmSCMatrix* symmmatrix(const RefSCDimension&) = 0;
    virtual DiagSCMatrix* diagmatrix(const RefSCDimension&) = 0;
    virtual SCVector* vector(const RefSCDimension&) = 0;

    //. Given the dimensions and a \clsnmref{StateIn} object,
    //. restore matrices or vectors.
    SCMatrix* restore_matrix(StateIn&,
                             const RefSCDimension&,
                             const RefSCDimension&);
    SymmSCMatrix* restore_symmmatrix(StateIn&,
                                     const RefSCDimension&);
    DiagSCMatrix* restore_diagmatrix(StateIn&,             
                                     const RefSCDimension&);
    SCVector* restore_vector(StateIn&,
                             const RefSCDimension&);
};
DescribedClass_REF_dec(SCMatrixKit);

//. The \clsnmref{SCVector} class is the abstract base class for
//. \srccd{double} valued vectors.
class SCVector: public DescribedClass {
#   define CLASSNAME SCVector
#   include <util/class/classda.h>
  protected:
    RefSCDimension d;
    RefSCMatrixKit kit_;
  public:
    SCVector(const RefSCDimension&, SCMatrixKit *);

    //. Save and restore this in an implementation independent way.
    virtual void save(StateOut&);
    virtual void restore(StateIn&);

    //. Return the \clsnmref{SCMatrixKit} used to create this object.
    RefSCMatrixKit kit() const { return kit_; }

    // concrete functions (some can be overridden)
    //. Return a vector with the same dimension and same elements.
    virtual SCVector* copy();
    //. Return a vector with the same dimension but uninitialized memory.
    virtual SCVector* clone();

    virtual ~SCVector();
    //. Return the length of the vector.
    int n() { return d->n(); }
    //. Return the maximum absolute value element of this vector.
    virtual double maxabs();
    //. Normalize this.
    virtual void normalize();
    //. Assign each element to a random number between -1 and 1
    virtual void randomize();
    //. Assign all elements of this to \vrbl{val}.
    virtual void assign(double val);
    //. Assign element \vrbl{i} to \vrbl{v[i]} for all \vrbl{i}.
    virtual void assign(const double* v);
    //. Assign \vrbl{v[i]} to element \vrbl{i} for all \vrbl{i}.
    virtual void convert(double* v);
    //. Convert an \clsnmref{SCVector} of a different specialization
    //. to this specialization and possibly accumulate the data.
    virtual void convert(SCVector*);
    virtual void convert_accumulate(SCVector*);
    //. Make \srccd{this} have the same elements as \vrbl{v}.
    //. The dimensions must match.
    virtual void assign(SCVector* v);
    //. Multiply each element by \vrbl{val}.
    virtual void scale(double val);

    //. Return the \clsnmref{RefSCDimension} corresponding to this vector.
    RefSCDimension dim() const { return d; }
    //. Set element \vrbl{i} to \vrbl{val}.
    virtual void set_element(int,double) = 0;
    //. Add \vrbl{val} to element \vrbl{i}.
    virtual void accumulate_element(int,double) = 0;
    //. Return the value of element \vrbl{i}.
    virtual double get_element(int) = 0;
    //. Sum the result of \vrbl{m} times \vrbl{v} into \srccd{this}.
    virtual void accumulate_product(SymmSCMatrix* m, SCVector* v);
    virtual void accumulate_product(SCMatrix* m, SCVector* v) = 0;
    //. Sum \vrbl{v} into this.
    virtual void accumulate(SCVector*v) = 0;
    //. Sum \vrbl{m} into this.  One of \vrbl{m}'s dimensions must be 1.
    virtual void accumulate(SCMatrix*m) = 0;
    //. Return the dot product.
    virtual double scalar_product(SCVector*) = 0;
    //. Perform the element operation \vrbl{op} on each element of this.
    virtual void element_op(const RefSCElementOp&) = 0;
    virtual void element_op(const RefSCElementOp2&,
                            SCVector*) = 0;
    virtual void element_op(const RefSCElementOp3&,
                            SCVector*,SCVector*) = 0;
    //. Print out the vector.
    virtual void print(ostream&);
    virtual void print(const char* title=0,ostream& out=cout, int =10) = 0;

    //. Returns the message group used by the matrix kit
    RefMessageGrp messagegrp();
    
    //. Returns iterators for the local (rapidly accessible)
    //. blocks used in this vector.  Only one iterator is allowed
    //. for a matrix is it has \srccd{Accum} or \srccd{Write}
    //. access is allowed.  Multiple \srccd{Read} iterators are permitted.
    virtual RefSCMatrixSubblockIter local_blocks(
        SCMatrixSubblockIter::Access) = 0;
    //. Returns iterators for the all blocks used in this vector.
    virtual RefSCMatrixSubblockIter all_blocks(SCMatrixSubblockIter::Access) = 0;
};

//. The \clsnmref{SCMatrix} class is the abstract base class for general
//. \srccd{double} valued n by m matrices.
//. For symmetric matrices use \clsnmref{SymmSCMatrix} and for
//. diagonal matrices use \clsnmref{DiagSCMatrix}.
class SCMatrix: public DescribedClass {
#   define CLASSNAME SCMatrix
#   include <util/class/classda.h>
  protected:
    RefSCDimension d1,d2;
    RefSCMatrixKit kit_;
  public:
    // concrete functions (some can be overridden)
    SCMatrix(const RefSCDimension&, const RefSCDimension&, SCMatrixKit *);
    virtual ~SCMatrix();

    //. Save and restore this in an implementation independent way.
    virtual void save(StateOut&);
    virtual void restore(StateIn&);

    //. Return the \clsnmref{SCMatrixKit} used to create this object.
    RefSCMatrixKit kit() const { return kit_; }

    //. Return the number of rows or columns.
    int nrow() const { return d1->n(); }
    int ncol() const { return d2->n(); }
    //. Return the maximum absolute value element.
    virtual double maxabs();
    //. Assign each element to a random number between -1 and 1
    virtual void randomize();
    //. Set all elements to \vrbl{val}.
    virtual void assign(double val);
    //. Assign element \vrbl{i}, \vrbl{j} to
    //. \srccd{\vrbl{m}[\vrbl{i}*nrow()+\vrbl{j}]}.
    virtual void assign(const double* m);
    //. Assign element \vrbl{i}, \vrbl{j} to \srccd{\vrbl{m}[\vrbl{i}][\vrbl{j}]}.
    virtual void assign(const double** m);
    //. Like the \srccd{assign} members, but these write values
    //. to the arguments.
    virtual void convert(double*);
    virtual void convert(double**);
    //. Convert an \clsnmref{SCMatrix} of a different specialization
    //. to this specialization and possibly accumulate the data.
    virtual void convert(SCMatrix*);
    virtual void convert_accumulate(SCMatrix*);
    //. Make \srccd{this} have the same elements as \vrbl{m}.
    //. The dimensions must match.
    virtual void assign(SCMatrix* m);
    //. Multiply all elements by \vrbl{val}.
    virtual void scale(double val);
    //. Scale the diagonal elements by \vrbl{val}.
    virtual void scale_diagonal(double val);
    //. Shift the diagonal elements by \vrbl{val}.
    virtual void shift_diagonal(double val);
    //. Make \srccd{this} equal to the unit matrix.
    virtual void unit();
    //. Return a matrix with the same dimension and same elements.
    virtual SCMatrix* copy();
    //. Return a matrix with the same dimension but uninitialized memory.
    virtual SCMatrix* clone();

    // pure virtual functions
    //. Return the row or column dimension.
    RefSCDimension rowdim() const { return d1; }
    RefSCDimension coldim() const { return d2; }
    //. Return or modify an element.
    virtual double get_element(int,int) = 0;
    virtual void set_element(int,int,double) = 0;
    virtual void accumulate_element(int,int,double) = 0;
    
    //. Return a subblock of \srccd{this}.  The subblock is defined as
    //. the rows starting at \srccd{br} and ending at \srccd{er}, and the
    //. columns beginning at \srccd{bc} and ending at \srccd{ec}.
    virtual SCMatrix * get_subblock(int br, int er, int bc, int ec) =0;

    //. Assign \srccd{m} to a subblock of \srccd{this}.
    virtual void assign_subblock(SCMatrix *m, int, int, int, int, int=0, int=0) =0;

    //. Sum \srccd{m} into a subblock of \srccd{this}.
    virtual void accumulate_subblock(SCMatrix *m, int, int, int, int, int=0,int=0) =0;
    
    //. Return a row or column of \srccd{this}.
    virtual SCVector * get_row(int i) =0;
    virtual SCVector * get_column(int i) =0;

    //. Assign \srccd{v} to a row or column of \srccd{this}.
    virtual void assign_row(SCVector *v, int i) =0;
    virtual void assign_column(SCVector *v, int i) =0;
    
    //. Sum \srccd{v} to a row or column of \srccd{this}.
    virtual void accumulate_row(SCVector *v, int i) =0;
    virtual void accumulate_column(SCVector *v, int i) =0;
    
    //. Sum \vrbl{m} into this.
    virtual void accumulate(SCMatrix* m) = 0;
    virtual void accumulate(SymmSCMatrix* m) = 0;
    virtual void accumulate(DiagSCMatrix* m) = 0;
    virtual void accumulate(SCVector*) = 0;
    //. Sum into \srccd{this} the products of various vectors or matrices.
    virtual void accumulate_outer_product(SCVector*,SCVector*) = 0;
    virtual void accumulate_product(SCMatrix*,SCMatrix*) = 0;
    virtual void accumulate_product(SCMatrix*,SymmSCMatrix*);
    virtual void accumulate_product(SCMatrix*,DiagSCMatrix*);
    virtual void accumulate_product(SymmSCMatrix*,SCMatrix*);
    virtual void accumulate_product(DiagSCMatrix*,SCMatrix*);
    //. Transpose \srccd{this}.
    virtual void transpose_this() = 0;
    //. Return the trace.
    virtual double trace() =0;
    //. Invert \srccd{this}.
    virtual double invert_this() = 0;
    //. Return the determinant of \srccd{this}.  \srccd{this} is overwritten.
    virtual double determ_this() = 0;

    //. Compute the singular value decomposition for \srccd{this},
    //. possibly destroying this.
    virtual void svd_this(SCMatrix *U, DiagSCMatrix *sigma, SCMatrix *V);
    virtual double solve_this(SCVector*) = 0;
    virtual void gen_invert_this() = 0;

    //. Schmidt orthogonalize \srccd{this}.  \srccd{S} is the overlap matrix.
    //. \srccd{n} is the number of columns to orthogonalize.
    virtual void schmidt_orthog(SymmSCMatrix*, int n) =0;
    
    //. Perform the element operation \vrbl{op} on each element of this.
    virtual void element_op(const RefSCElementOp&) = 0;
    virtual void element_op(const RefSCElementOp2&,
                            SCMatrix*) = 0;
    virtual void element_op(const RefSCElementOp3&,
                            SCMatrix*,SCMatrix*) = 0;
    //. Print out the matrix.
    virtual void print(ostream&);
    virtual void print(const char* title=0,ostream& out=cout, int =10) = 0;

    //. Returns the message group used by the matrix kit
    RefMessageGrp messagegrp();
    
    //. Returns iterators for the local (rapidly accessible)
    //. blocks used in this matrix.
    virtual RefSCMatrixSubblockIter local_blocks(
        SCMatrixSubblockIter::Access) = 0;
    //. Returns iterators for the all blocks used in this matrix.
    virtual RefSCMatrixSubblockIter all_blocks(
        SCMatrixSubblockIter::Access) = 0;
};

//. The \clsnmref{SymmSCMatrix} class is the abstract base class for symmetric
//. \srccd{double} valued matrices.
class SymmSCMatrix: public DescribedClass {
#   define CLASSNAME SymmSCMatrix
#   include <util/class/classda.h>
  protected:
    RefSCDimension d;
    RefSCMatrixKit kit_;
  public:
    SymmSCMatrix(const RefSCDimension&, SCMatrixKit *);

    //. Return the \clsnmref{SCMatrixKit} object that created this object.
    RefSCMatrixKit kit() const { return kit_; }

    //. Save and restore this in an implementation independent way.
    virtual void save(StateOut&);
    virtual void restore(StateIn&);
    //. Return the maximum absolute value element of this vector.
    virtual double maxabs();
    //. Assign each element to a random number between -1 and 1
    virtual void randomize();
    //. Set all elements to \vrbl{val}.
    virtual void assign(double val);
    //. Assign element \vrbl{i}, \vrbl{j} to
    //. \srccd{\vrbl{m}[\vrbl{i}*(\vrbl{i}+1)/2+\vrbl{j}]}.
    virtual void assign(const double* m);
    //. Assign element \vrbl{i}, \vrbl{j} to \srccd{\vrbl{m}[\vrbl{i}][\vrbl{j}]}.
    virtual void assign(const double** m);
    //. Like the \srccd{assign} members, but these write values
    //. to the arguments.
    virtual void convert(double*);
    virtual void convert(double**);
    //. Convert an \clsnmref{SCSymmSCMatrix} of a different specialization
    //. to this specialization and possibly accumulate the data.
    virtual void convert(SymmSCMatrix*);
    virtual void convert_accumulate(SymmSCMatrix*);
    //. Make \srccd{this} have the same elements as \vrbl{m}.
    //. The dimensions must match.
    virtual void assign(SymmSCMatrix* m);
    //. Multiply all elements by \vrbl{val}.
    virtual void scale(double);
    //. Scale the diagonal elements by \vrbl{val}.
    virtual void scale_diagonal(double);
    //. Shift the diagonal elements by \vrbl{val}.
    virtual void shift_diagonal(double);
    //. Make \srccd{this} equal to the unit matrix.
    virtual void unit();
    //. Return the dimension.
    int n() { return d->n(); }
    //. Return a matrix with the same dimension and same elements.
    virtual SymmSCMatrix* copy();
    //. Return a matrix with the same dimension but uninitialized memory.
    virtual SymmSCMatrix* clone();

    // pure virtual functions
    //. Return the dimension.
    RefSCDimension dim() const { return d; }
    //. Return or modify an element.
    virtual double get_element(int,int) = 0;
    virtual void set_element(int,int,double) = 0;
    virtual void accumulate_element(int,int,double) = 0;

    //. Return a subblock of \srccd{this}.  The subblock is defined as
    //. the rows starting at \srccd{br} and ending at \srccd{er}, and the
    //. columns beginning at \srccd{bc} and ending at \srccd{ec}.
    virtual SCMatrix * get_subblock(int br, int er, int bc, int ec) =0;
    virtual SymmSCMatrix * get_subblock(int br, int er) =0;

    //. Assign \srccd{m} to a subblock of \srccd{this}.
    virtual void assign_subblock(SCMatrix *m, int, int, int, int) =0;
    virtual void assign_subblock(SymmSCMatrix *m, int, int) =0;

    //. Sum \srccd{m} into a subblock of \srccd{this}.
    virtual void accumulate_subblock(SCMatrix *m, int, int, int, int) =0;
    virtual void accumulate_subblock(SymmSCMatrix *m, int, int) =0;

    //. Return a row of \srccd{this}.
    virtual SCVector * get_row(int i) =0;

    //. Assign \srccd{v} to a row of \srccd{this}.
    virtual void assign_row(SCVector *v, int i) =0;
    
    //. Sum \srccd{v} to a row of \srccd{this}.
    virtual void accumulate_row(SCVector *v, int i) =0;

    //. Diagonalize \srccd{this}, placing the eigenvalues in \vrbl{d}
    //. and the eigenvectors in \vrbl{m}.
    virtual void diagonalize(DiagSCMatrix*d,SCMatrix*m) = 0;
    //. Sum \vrbl{m} into this.
    virtual void accumulate(SymmSCMatrix* m) = 0;
    //. Sum into \srccd{this} the products of various vectors or matrices.
    virtual void accumulate_symmetric_sum(SCMatrix*) = 0;
    virtual void accumulate_symmetric_product(SCMatrix*);
    virtual void accumulate_transform(SCMatrix*,SymmSCMatrix*);
    virtual void accumulate_transform(SCMatrix*,DiagSCMatrix*);
    virtual void accumulate_transform(SymmSCMatrix*,SymmSCMatrix*);
    virtual void accumulate_symmetric_outer_product(SCVector*);
    //. Return the scalar obtained by multiplying \srccd{this} on the
    //. left and right by \vrbl{v}.
    virtual double scalar_product(SCVector* v);
    //. Return the trace.
    virtual double trace() = 0;
    //. Invert \srccd{this}.
    virtual double invert_this() = 0;
    //. Return the determinant of \srccd{this}.  \srccd{this} is overwritten.
    virtual double determ_this() = 0;

    virtual double solve_this(SCVector*) = 0;
    virtual void gen_invert_this() = 0;

    //. Perform the element operation \vrbl{op} on each element of this.
    virtual void element_op(const RefSCElementOp&) = 0;
    virtual void element_op(const RefSCElementOp2&,
                            SymmSCMatrix*) = 0;
    virtual void element_op(const RefSCElementOp3&,
                            SymmSCMatrix*,SymmSCMatrix*) = 0;
    //. Print out the matrix.
    virtual void print(ostream&);
    virtual void print(const char* title=0,ostream& out=cout, int =10);

    //. Returns the message group used by the matrix kit
    RefMessageGrp messagegrp();
    
    //. Returns iterators for the local (rapidly accessible)
    //. blocks used in this matrix.
    virtual RefSCMatrixSubblockIter local_blocks(
        SCMatrixSubblockIter::Access) = 0;
    //. Returns iterators for the all blocks used in this matrix.
    virtual RefSCMatrixSubblockIter all_blocks(
        SCMatrixSubblockIter::Access) = 0;
};

//. The \clsnmref{SymmSCMatrix} class is the abstract base class for diagonal
//. \srccd{double} valued matrices.
class DiagSCMatrix: public DescribedClass {
#   define CLASSNAME DiagSCMatrix
#   include <util/class/classda.h>
  protected:
    RefSCDimension d;
    RefSCMatrixKit kit_;
  public:
    DiagSCMatrix(const RefSCDimension&, SCMatrixKit *);

    //. Return the \clsnmref{SCMatrixKit} used to create this object.
    RefSCMatrixKit kit() const { return kit_; }

    //. Save and restore this in an implementation independent way.
    virtual void save(StateOut&);
    virtual void restore(StateIn&);

    //. Return the maximum absolute value element of this vector.
    virtual double maxabs();
    //. Assign each element to a random number between -1 and 1
    virtual void randomize();
    //. Set all elements to \vrbl{val}.
    virtual void assign(double val);
    //. Assign element \vrbl{i}, \vrbl{i} to
    //. \srccd{\vrbl{m}[\vrbl{i}]}.
    virtual void assign(const double*);
    //. Like the \srccd{assign} member, but this write values
    //. to the argument.
    virtual void convert(double*);
    //. Convert an \clsnmref{SCDiagSCMatrix} of a different specialization
    //. to this specialization and possibly accumulate the data.
    virtual void convert(DiagSCMatrix*);
    virtual void convert_accumulate(DiagSCMatrix*);
    //. Make \srccd{this} have the same elements as \vrbl{m}.
    //. The dimensions must match.
    virtual void assign(DiagSCMatrix*);
    //. Multiply all elements by \vrbl{val}.
    virtual void scale(double);
    //. Return the dimension.
    int n() const { return d->n(); }
    //. Return a matrix with the same dimension and same elements.
    virtual DiagSCMatrix* copy();
    //. Return a matrix with the same dimension but uninitialized memory.
    virtual DiagSCMatrix* clone();

    // pure virtual functions
    //. Return the dimension.
    RefSCDimension dim() const { return d; }
    //. Return or modify an element.
    virtual double get_element(int) = 0;
    virtual void set_element(int,double) = 0;
    virtual void accumulate_element(int,double) = 0;
    //. Sum \vrbl{m} into this.
    virtual void accumulate(DiagSCMatrix* m) = 0;
    //. Return the trace.
    virtual double trace() = 0;
    //. Return the determinant of \srccd{this}.  \srccd{this} is overwritten.
    virtual double determ_this() = 0;
    //. Invert \srccd{this}.
    virtual double invert_this() = 0;
    //. Do a generalized inversion of \srccd{this}.
    virtual void gen_invert_this() = 0;
    //. Perform the element operation \vrbl{op} on each element of this.
    virtual void element_op(const RefSCElementOp&) = 0;
    virtual void element_op(const RefSCElementOp2&,
                            DiagSCMatrix*) = 0;
    virtual void element_op(const RefSCElementOp3&,
                            DiagSCMatrix*,DiagSCMatrix*) = 0;
    //. Print out the matrix.
    virtual void print(ostream&);
    virtual void print(const char* title=0,ostream& out=cout, int =10);

    //. Returns the message group used by the matrix kit
    RefMessageGrp messagegrp();
    
    //. Returns iterators for the local (rapidly accessible)
    //. blocks used in this matrix.
    virtual RefSCMatrixSubblockIter local_blocks(
        SCMatrixSubblockIter::Access) = 0;
    //. Returns iterators for the all blocks used in this matrix.
    virtual RefSCMatrixSubblockIter all_blocks(
        SCMatrixSubblockIter::Access) = 0;
};

#endif
