
#ifndef _math_scmat_matrix_h
#define _math_scmat_matrix_h
#ifdef __GNUC__
#pragma interface
#endif

#include <iostream.h>
#include <util/container/array.h>
#include <util/container/set.h>

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

class RefSCMatrix;
class RefSymmSCMatrix;
DescribedClass_named_REF_dec(RefDCSCVector,SCVector);
//. The \clsnm{RefSCVector} class is a smart pointer to an \clsnmref{SCVector}
//. specialization.  Valid indices range from \srccd{0} to \srccd{n-1}.
class RefSCVector: public RefDCSCVector {
    // standard overrides
  public:
    //. Initializes the vector pointer to \srccd{0}.  The reference must
    //. be initialized before it is used.
    RefSCVector();
    //. Make this and \vrbl{v} refer to the same \clsnmref{SCVector}.
    RefSCVector(const RefSCVector& v);
    //. Make this refer to \vrbl{v}.
    RefSCVector(SCVector *v);
    // don't allow automatic conversion from any reference to a
    // described class
    ~RefSCVector();
    //. Make this refer to \vrbl{v}.
    RefSCVector& operator=(SCVector* v);
    //. Make this and \vrbl{v} refer to the same \clsnmref{SCVector}.
    RefSCVector& operator=(const RefSCVector& v);

    // vector specific members
  public:
    //. Create a vector with dimension \vrbl{dim}.  The data values
    //. are undefined.
    RefSCVector(const RefSCDimension& dim,const RefSCMatrixKit&);

    //. Return an l-value that can be used to assign or retrieve an element.
    SCVectordouble operator()(int) const;
    SCVectordouble operator[](int) const;
    //. Add two vectors.
    RefSCVector operator+(const RefSCVector&a) const;
    //. Subtract two vectors.
    RefSCVector operator-(const RefSCVector&a) const;
    //. Scale a vector.
    RefSCVector operator*(double) const;
    //. Return the outer product between this and \vrbl{v}.
    RefSCMatrix outer_product(const RefSCVector& v) const;
    //. The outer product of this with itself is a symmetric matrix.
    RefSymmSCMatrix symmetric_outer_product() const;

    //. These call the \clsnmref{SCMatrix} members of the same name
    //. after checking for references to \srccd{0}.
    void set_element(int i,double val) const;
    void accumulate_element(int i,double val) const;
    double get_element(int) const;
    int n() const;
    RefSCDimension dim() const;
    RefSCMatrixKit kit() const;
    RefSCVector clone() const;
    RefSCVector copy() const;
    double maxabs() const;
    double scalar_product(const RefSCVector&) const;
    double dot(const RefSCVector&) const;
    void normalize() const;
    void randomize() const;
    void assign(const RefSCVector& v) const;
    void assign(double val) const;
    void assign(const double* v) const;
    void convert(double*) const;
    void scale(double val) const;
    void accumulate(const RefSCVector& v) const;
    void accumulate_product(const RefSymmSCMatrix&, const RefSCVector&);
    void accumulate_product(const RefSCMatrix&, const RefSCVector&);
    void element_op(const RefSCElementOp& op) const;
    void element_op(const RefSCElementOp2&,
                    const RefSCVector&) const;
    void element_op(const RefSCElementOp3&,
                    const RefSCVector&,
                    const RefSCVector&) const;
    void print(ostream&out) const;
    void print(const char*title=0, ostream&out=cout, int precision=10) const;
    void save(StateOut&);
    void restore(StateIn&);
};
RefSCVector operator*(double,const RefSCVector&);
ARRAY_dec(RefSCVector);
SET_dec(RefSCVector);

class RefSymmSCMatrix;
class RefDiagSCMatrix;
DescribedClass_named_REF_dec(RefDCSCMatrix,SCMatrix);
//. The \clsnm{RefSCMatrix} class is a smart pointer to an \clsnmref{SCMatrix}
//. specialization.
class RefSCMatrix: public RefDCSCMatrix {
    // standard overrides
  public:
    //. Initializes the matrix pointer to \vrbl{0}.  The reference must
    //. be initialized before it is used.
    RefSCMatrix();
    //. Make this and \srccd{m} refer to the same \srccd{SCMatrix}.
    RefSCMatrix(const RefSCMatrix& m);
    //. Make this refer to \srccd{m}.
     RefSCMatrix(SCMatrix* m);
    ~RefSCMatrix();
    //. Make this refer to \srccd{m}.
    RefSCMatrix& operator=(SCMatrix* m);
    //. Make this and \srccd{m} refer to the same matrix.
    RefSCMatrix& operator=(const RefSCMatrix& m);

    // matrix specific members
  public:
    //. Create a vector with dimension \vrbl{d1} by \vrbl{d2}.
    //. The data values are undefined.
    RefSCMatrix(const RefSCDimension& d1,const RefSCDimension& d2,
                const RefSCMatrixKit&);

    //. Multiply this by a vector and return a vector.
    RefSCVector operator*(const RefSCVector&) const;

    //. Multiply this by a matrix and return a matrix.
    RefSCMatrix operator*(const RefSCMatrix&) const;
    RefSCMatrix operator*(const RefSymmSCMatrix&) const;
    RefSCMatrix operator*(const RefDiagSCMatrix&) const;

    //. Multiply this by a scalar and return the result.
    RefSCMatrix operator*(double) const;

    //. Matrix addition and subtraction.
    RefSCMatrix operator+(const RefSCMatrix&) const;
    RefSCMatrix operator-(const RefSCMatrix&) const;

    //. Return the transpose of this.
    RefSCMatrix t() const;
    //. Return the inverse of this.
    RefSCMatrix i() const;
    //. Return the generalized inverse of this.
    RefSCMatrix gi() const;

    //. These call the \clsnmref{SCMatrix} members of the same name
    //. after checking for references to \srccd{0}.
    RefSCMatrix clone() const;
    RefSCMatrix copy() const;

    RefSCMatrix get_subblock(int br, int er, int bc, int ec);
    void assign_subblock(const RefSCMatrix&, int br, int er, int bc, int ec,
                         int source_br = 0, int source_bc = 0);
    void accumulate_subblock(const RefSCMatrix&, int, int, int, int,
                             int source_br = 0, int source_bc = 0);
    RefSCVector get_row(int);
    RefSCVector get_column(int);
    void assign_row(const RefSCVector&, int);
    void assign_column(const RefSCVector&, int);
    void accumulate_row(const RefSCVector&, int);
    void accumulate_column(const RefSCVector&, int);

    void accumulate_outer_product(const RefSCVector&,const RefSCVector&) const;
    void accumulate_product(const RefSCMatrix&,const RefSCMatrix&) const;
    void assign(const RefSCMatrix&) const;
    void scale(double) const;
    void randomize() const;
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
    RefSCDimension rowdim() const;
    RefSCDimension coldim() const;
    RefSCMatrixKit kit() const;
    void set_element(int,int,double) const;
    void accumulate_element(int,int,double) const;
    double get_element(int,int) const;
    void print(ostream&) const;
    void print(const char*title=0,ostream&out=cout, int =10) const;
    double trace() const;
    void save(StateOut&);
    void restore(StateIn&);

    //. Compute the singular value decomposition.
    //. \srccd{this} = \vrbl{U} \vrbl{sigma} \vrbl{V}\srccd{.t()}.
    //. The dimension of \vrbl{sigma} is the smallest dimension of
    //. \srccd{this}.  \vrbl{U}, \vrbl{V}, and \vrbl{sigma} must already
    //. have the correct dimensions and are overwritten.
    void svd(const RefSCMatrix &U,
             const RefDiagSCMatrix &sigma,
             const RefSCMatrix &V);
    //. Solves \srccd{this} \vrbl{x} = \vrbl{v}.  Overwrites
    //. \vrbl{v} with \vrbl{x}.
    double solve_lin(const RefSCVector& v) const;
    //. Returns the determinant of the referenced matrix.
    double determ() const;
    //. Assign and examine matrix elements.
    SCMatrixdouble operator()(int i,int j) const;
};
//. Allow multiplication with a scalar on the left.
RefSCMatrix operator*(double,const RefSCMatrix&);

DescribedClass_named_REF_dec(RefDCSymmSCMatrix,SymmSCMatrix);
//. The \clsnmref{RefSymmSCMatrix} class is a smart pointer to an
//. \clsnmref{SCSymmSCMatrix} specialization.
class RefSymmSCMatrix: public RefDCSymmSCMatrix {
    // standard overrides
  public:
    //. Initializes the matrix pointer to \vrbl{0}.  The reference must
    //. be initialized before it is used.
    RefSymmSCMatrix();
    //. Make this and \vrbl{m} refer to the same \clsnmref{SCMatrix}.
    RefSymmSCMatrix(const RefSymmSCMatrix& m);
    //. Make this refer to \vrbl{m}.
    RefSymmSCMatrix(SymmSCMatrix *m);
    ~RefSymmSCMatrix();
    //. Make this refer to \vrbl{m}.
    RefSymmSCMatrix& operator=(SymmSCMatrix* m);
    //. Make this and \vrbl{m} refer to the same matrix.
    RefSymmSCMatrix& operator=(const RefSymmSCMatrix& m);

    // matrix specific members
  public:
    //. Create a vector with dimension \vrbl{d} by \vrbl{d}.
    //. The data values are undefined.
    RefSymmSCMatrix(const RefSCDimension& d,const RefSCMatrixKit&);
    //. Multiply this by a matrix and return a matrix.
    RefSCMatrix operator*(const RefSCMatrix&) const;
    //. Multiply this by a vector and return a vector.
    RefSCVector operator*(const RefSCVector&a) const;
    RefSymmSCMatrix operator*(double) const;
    //. Matrix addition and subtraction.
    RefSymmSCMatrix operator+(const RefSymmSCMatrix&) const;
    RefSymmSCMatrix operator-(const RefSymmSCMatrix&) const;
    //. Return the inverse of this.
    RefSymmSCMatrix i() const;
    //. Return the generalized inverse of this.
    RefSymmSCMatrix gi() const;
    //. These call the \clsnmref{SCMatrix} members of the same name
    //. after checking for references to \srccd{0}.
    RefSymmSCMatrix clone() const;
    RefSymmSCMatrix copy() const;
    void set_element(int,int,double) const;
    void accumulate_element(int,int,double) const;
    double get_element(int,int) const;

    RefSCMatrix get_subblock(int br, int er, int bc, int ec);
    RefSymmSCMatrix get_subblock(int br, int er);
    void assign_subblock(const RefSCMatrix&, int br, int er, int bc, int ec);
    void assign_subblock(const RefSymmSCMatrix&, int br, int er);
    void accumulate_subblock(const RefSCMatrix&, int, int, int, int);
    void accumulate_subblock(const RefSymmSCMatrix&, int, int);
    RefSCVector get_row(int);
    void assign_row(const RefSCVector&, int);
    void accumulate_row(const RefSCVector&, int);

    void accumulate_symmetric_outer_product(const RefSCVector&) const;
    double scalar_product(const RefSCVector&) const;
    void accumulate_symmetric_product(const RefSCMatrix&) const;
    void accumulate_symmetric_sum(const RefSCMatrix&) const;
    //. Add a * b * a.t() to this.
    void accumulate_transform(const RefSCMatrix&a,const RefSymmSCMatrix&b) const;
    void accumulate_transform(const RefSCMatrix&,const RefDiagSCMatrix&) const;
    void accumulate_transform(const RefSymmSCMatrix&,const RefSymmSCMatrix&) const;

    void randomize() const;
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
    double trace() const;
    int n() const;
    RefSCDimension dim() const;
    RefSCMatrixKit kit() const;
    void print(ostream&) const;
    void print(const char*title=0,ostream&out=cout, int =10) const;
    void save(StateOut&);
    void restore(StateIn&);

    //. Solves \srccd{this} \vrbl{x} = \vrbl{v}.  Overwrites
    //. \vrbl{v} with \vrbl{x}.
    double solve_lin(const RefSCVector&) const;
    //. Returns the determinant of the referenced matrix.
    double determ() const;
    //. Returns the eigenvalues of the reference matrix.
    RefDiagSCMatrix eigvals() const;
    //. Returns the eigenvectors of the reference matrix.
    RefSCMatrix eigvecs() const;
    //. Sets \vrbl{eigvals} to the eigenvalues and \vrbl{eigvecs}
    //. to the eigenvalues and eigenvectors of the referenced matrix.
    void diagonalize(const RefDiagSCMatrix& eigvals,
                     const RefSCMatrix& eigvecs) const;
    //. Assign and examine matrix elements.
    SymmSCMatrixdouble operator()(int i,int j) const;
};
//. Allow multiplication with a scalar on the left.
RefSymmSCMatrix operator*(double,const RefSymmSCMatrix&);

DescribedClass_named_REF_dec(RefDCDiagSCMatrix,DiagSCMatrix);
//. The \clsnmref{RefDiagSCMatrix} class is a smart pointer to an
//. \clsnmref{DiagSCMatrix} specialization.
class RefDiagSCMatrix: public RefDCDiagSCMatrix {
    // standard overrides
  public:
    //. Initializes the matrix pointer to \vrbl{0}.  The reference must
    //. be initialized before it is used.
    RefDiagSCMatrix();
    //. Make this and \vrbl{m} refer to the same \clsnmref{SCMatrix}.
    RefDiagSCMatrix(const RefDiagSCMatrix& m);
    //. Make this refer to \vrbl{m}.
    RefDiagSCMatrix(DiagSCMatrix *m);
    ~RefDiagSCMatrix();
    //. Make this refer to \vrbl{m}.
    RefDiagSCMatrix& operator=(DiagSCMatrix* m);
    //. Make this and \vrbl{m} refer to the same matrix.
    RefDiagSCMatrix& operator=(const RefDiagSCMatrix & m);

    // matrix specific members
  public:
    //. Create a vector with dimension \vrbl{d} by \vrbl{d}.
    //. The data values are undefined.
    RefDiagSCMatrix(const RefSCDimension&,const RefSCMatrixKit&);
    //. Multiply this by a matrix and return a matrix.
    RefSCMatrix operator*(const RefSCMatrix&) const;
    RefDiagSCMatrix operator*(double) const;
    //. Matrix addition and subtraction.
    RefDiagSCMatrix operator+(const RefDiagSCMatrix&) const;
    RefDiagSCMatrix operator-(const RefDiagSCMatrix&) const;
    //. Return the inverse of this.
    RefDiagSCMatrix i() const;
    //. Return the generalized inverse of this.
    RefDiagSCMatrix gi() const;
    //. These call the \clsnmref{SCMatrix} members of the same name
    //. after checking for references to \srccd{0}.
    RefDiagSCMatrix clone() const;
    RefDiagSCMatrix copy() const;
    void set_element(int,double) const;
    void accumulate_element(int,double) const;
    double get_element(int) const;
    void randomize() const;
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
    RefSCDimension dim() const;
    RefSCMatrixKit kit() const;
    double trace() const;
    void print(ostream&) const;
    void print(const char*title=0,ostream&out=cout, int =10) const;
    void save(StateOut&);
    void restore(StateIn&);
    //. Returns the determinant of the referenced matrix.
    double determ() const;
    //. Assign and examine matrix elements.
    DiagSCMatrixdouble operator()(int i) const;
};
//. Allow multiplication with a scalar on the left.
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
