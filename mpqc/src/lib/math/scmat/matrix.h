
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

DCRef_declare(SCDimension);
SSRef_declare(SCDimension);
//texi
//  The @code{RefSCDimension} class is a smart pointer to an
//   @code{SCDimension} specialization.
class RefSCDimension: public SSRefSCDimension {
    // standard overrides
  public:
    //texi Initializes the dimension pointer to @code{0}.  The
    // reference must be initialized before it is used.
    RefSCDimension();
    //texi Make this and @var{d} refer to the same @code{SCDimension}.
    RefSCDimension(const RefSCDimension& d);
    //texi Make this refer to @var{d}.
    RefSCDimension(SCDimension *d);

    RefSCDimension(const RefDescribedClassBase&);
    ~RefSCDimension();
    //texi Make this refer to @var{d}.
    RefSCDimension& operator=(SCDimension* d);

    RefSCDimension& operator=(const RefDescribedClassBase & c);
    //texi Make this and @var{d} refer to the same @code{SCDimension}.
    RefSCDimension& operator=(const RefSCDimension & d);

    // dimension specific functions
  public:
    //texi Return the dimension.
    operator int();
    int n();
};

class RefSCMatrix;
class RefSymmSCMatrix;
SavableState_named_REF_dec(RefSSSCVector,SCVector);
//texi
//  The @code{RefSCVector} class is a smart pointer to an @code{SCVector}
//  specialization.  Valid indices range from @code{0} to @code{n-1}.
class RefSCVector: public RefSSSCVector {
    // standard overrides
  public:
    //texi Initializes the vector pointer to @code{0}.  The reference must
    // be initialized before it is used.
    RefSCVector();
    //texi Restore the vector's state.
    RefSCVector(StateIn&);
    //texi Make this and @var{v} refer to the same @code{SCVector}.
    RefSCVector(const RefSCVector& v);
    //texi Make this refer to @var{v}.
    RefSCVector(SCVector *v);
    // don't allow automatic conversion from any reference to a
    // described class
    //RefSCVector(RefDescribedClassBase&);
    ~RefSCVector();
    //texi Make this refer to @var{v}.
    RefSCVector& operator=(SCVector* v);
    //RefSCVector& operator=( RefDescribedClassBase & c);
    //texi Make this and @var{v} refer to the same @code{SCVector}.
    RefSCVector& operator=(const RefSCVector& v);

    // vector specific members
  public:
    //texi Create a vector with dimension @var{dim}.  The data values
    // are undefined.
    RefSCVector(const RefSCDimension& dim);

    //texi Return an l-value that can be used to assign or retrieve an element.
    SCVectordouble operator()(int) const;
    SCVectordouble operator[](int) const;
    //texi Add two vectors.
    RefSCVector operator+(const RefSCVector&a) const;
    //texi Subtract two vectors.
    RefSCVector operator-(const RefSCVector&a) const;
    //texi Scale a vector.
    RefSCVector operator*(double) const;
    //texi Return the outer product between this and @var{v}.
    RefSCMatrix outer_product(const RefSCVector& v) const;
    //texi The outer product of this with itself is a symmetric matrix.
    RefSymmSCMatrix symmetric_outer_product() const;

    //texi These call the @code{SCMatrix} members of the same name
    // after checking for references to @code{0}.
    void set_element(int i,double val) const;
    double get_element(int) const;
    int n() const;
    RefSCDimension dim() const;
    RefSCVector clone() const;
    RefSCVector copy() const;
    double maxabs() const;
    double scalar_product(const RefSCVector&) const;
    double dot(const RefSCVector&) const;
    void normalize() const;
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
};
RefSCVector operator*(double,const RefSCVector&);
ARRAY_dec(RefSCVector);
SET_dec(RefSCVector);

class RefSymmSCMatrix;
class RefDiagSCMatrix;
SavableState_named_REF_dec(RefSSSCMatrix,SCMatrix);
//texi
//  The @code{RefSCMatrix} class is a smart pointer to an @code{SCMatrix}
//  specialization.
class RefSCMatrix: public RefSSSCMatrix {
    // standard overrides
  public:
    //texi Initializes the matrix pointer to @var{0}.  The reference must
    // be initialized before it is used.
    RefSCMatrix();
    //texi Restore the matrix's state.
    RefSCMatrix(StateIn&);
    //texi Make this and @var{m} refer to the same @code{SCMatrix}.
    RefSCMatrix(const RefSCMatrix& m);
    //texi Make this refer to @var{m}.
     RefSCMatrix(SCMatrix* m);
    //RefSCMatrix(RefDescribedClassBase&);
    ~RefSCMatrix();
    //texi Make this refer to @var{m}.
    RefSCMatrix& operator=(SCMatrix* m);
    //RefSCMatrix& operator=( RefDescribedClassBase & c);
    //texi Make this and @var{m} refer to the same matrix.
    RefSCMatrix& operator=(const RefSCMatrix& m);

    // matrix specific members
  public:
    //texi Create a vector with dimension @var{d1} by @var{d2}.
    // The data values are undefined.
    RefSCMatrix(const RefSCDimension& d1,const RefSCDimension& d2);
    //texi Multiply this by a vector and return a vector.
    RefSCVector operator*(const RefSCVector&) const;
    //texi Multiply this by a matrix and return a matrix.
    RefSCMatrix operator*(const RefSCMatrix&) const;
    RefSCMatrix operator*(const RefSymmSCMatrix&) const;
    RefSCMatrix operator*(const RefDiagSCMatrix&) const;
    //texi Multiply this by a scalar and return the result.
    RefSCMatrix operator*(double) const;
    //texi Matrix addition and subtraction.
    RefSCMatrix operator+(const RefSCMatrix&) const;
    RefSCMatrix operator-(const RefSCMatrix&) const;
    //texi Return the transpose of this.
    RefSCMatrix t() const;
    //texi Return the inverse of this.
    RefSCMatrix i() const;
    //texi Return the generalized inverse of this.
    RefSCMatrix gi() const;
    //texi These call the @code{SCMatrix} members of the same name
    // after checking for references to @code{0}.
    RefSCMatrix clone() const;
    RefSCMatrix copy() const;
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
    RefSCDimension rowdim() const;
    RefSCDimension coldim() const;
    void set_element(int,int,double) const;
    double get_element(int,int) const;
    void print(ostream&) const;
    void print(const char*title=0,ostream&out=cout, int =10) const;
    double trace() const;

    //texi Solves @code{this} @var{x} = @var{v}.  Overwrites
    // @var{v} with @var{x}.
    double solve_lin(const RefSCVector& v) const;
    //texi Returns the determinant of the referenced matrix.
    double determ() const;
    //texi Assign and examine matrix elements.
    SCMatrixdouble operator()(int i,int j) const;
};
//texi Allow multiplication with a scalar on the left.
RefSCMatrix operator*(double,const RefSCMatrix&);

SavableState_named_REF_dec(RefSSSymmSCMatrix,SymmSCMatrix);
//texi
//  The @code{RefSymmSCMatrix} class is a smart pointer to an
//   @code{SCSymmSCMatrix} specialization.
class RefSymmSCMatrix: public RefSSSymmSCMatrix {
    // standard overrides
  public:
    //texi Initializes the matrix pointer to @var{0}.  The reference must
    // be initialized before it is used.
    RefSymmSCMatrix();
    //texi Restore the matrix's state.
    RefSymmSCMatrix(StateIn&);
    //texi Make this and @var{m} refer to the same @code{SCMatrix}.
    RefSymmSCMatrix(const RefSymmSCMatrix& m);
    //texi Make this refer to @var{m}.
    RefSymmSCMatrix(SymmSCMatrix *m);
    //RefSymmSCMatrix(RefDescribedClassBase&);
    ~RefSymmSCMatrix();
    //texi Make this refer to @var{m}.
    RefSymmSCMatrix& operator=(SymmSCMatrix* m);
    //RefSymmSCMatrix& operator=( RefDescribedClassBase & c);
    //texi Make this and @var{m} refer to the same matrix.
    RefSymmSCMatrix& operator=(const RefSymmSCMatrix& m);

    // matrix specific members
  public:
    //texi Create a vector with dimension @var{d} by @var{d}.
    // The data values are undefined.
    RefSymmSCMatrix(const RefSCDimension& d);
    //texi Multiply this by a matrix and return a matrix.
    RefSCMatrix operator*(const RefSCMatrix&) const;
    //texi Multiply this by a vector and return a vector.
    RefSCVector operator*(const RefSCVector&a) const;
    RefSymmSCMatrix operator*(double) const;
    //texi Matrix addition and subtraction.
    RefSymmSCMatrix operator+(const RefSymmSCMatrix&) const;
    RefSymmSCMatrix operator-(const RefSymmSCMatrix&) const;
    //texi Return the inverse of this.
    RefSymmSCMatrix i() const;
    //texi Return the generalized inverse of this.
    RefSymmSCMatrix gi() const;
    //texi These call the @code{SCMatrix} members of the same name
    // after checking for references to @code{0}.
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
    double trace() const;
    int n() const;
    RefSCDimension dim() const;
    void print(ostream&) const;
    void print(const char*title=0,ostream&out=cout, int =10) const;

    //texi Solves @code{this} @var{x} = @var{v}.  Overwrites
    // @var{v} with @var{x}.
    double solve_lin(const RefSCVector&) const;
    //texi Returns the determinant of the referenced matrix.
    double determ() const;
    //texi Returns the eigenvalues of the reference matrix.
    RefDiagSCMatrix eigvals() const;
    //texi Returns the eigenvectors of the reference matrix.
    RefSCMatrix eigvecs() const;
    //texi Sets @var{eigvals} to the eigenvalues and @var{eigvecs}
    // to the eigenvalues and eigenvectors of the referenced matrix.
    void diagonalize(const RefDiagSCMatrix& eigvals,
                     const RefSCMatrix& eigvecs) const;
    //texi Assign and examine matrix elements.
    SymmSCMatrixdouble operator()(int i,int j) const;
};
//texi Allow multiplication with a scalar on the left.
RefSymmSCMatrix operator*(double,const RefSymmSCMatrix&);

SavableState_named_REF_dec(RefSSDiagSCMatrix,DiagSCMatrix);
//texi
//  The @code{RefDiagSCMatrix} class is a smart pointer to an
//  @code{DiagSCMatrix} specialization.
class RefDiagSCMatrix: public RefSSDiagSCMatrix {
    // standard overrides
  public:
    //texi Initializes the matrix pointer to @var{0}.  The reference must
    // be initialized before it is used.
    RefDiagSCMatrix();
    //texi Restore the matrix's state.
    RefDiagSCMatrix(StateIn&);
    //texi Make this and @var{m} refer to the same @code{SCMatrix}.
    RefDiagSCMatrix(const RefDiagSCMatrix& m);
    //texi Make this refer to @var{m}.
    RefDiagSCMatrix(DiagSCMatrix *m);
    //RefDiagSCMatrix(RefDescribedClassBase&);
    ~RefDiagSCMatrix();
    //texi Make this refer to @var{m}.
    RefDiagSCMatrix& operator=(DiagSCMatrix* m);
    //RefDiagSCMatrix& operator=( RefDescribedClassBase & c);
    //texi Make this and @var{m} refer to the same matrix.
    RefDiagSCMatrix& operator=(const RefDiagSCMatrix & m);

    // matrix specific members
  public:
    //texi Create a vector with dimension @var{d} by @var{d}.
    // The data values are undefined.
    RefDiagSCMatrix(const RefSCDimension&);
    //texi Multiply this by a matrix and return a matrix.
    RefSCMatrix operator*(const RefSCMatrix&) const;
    RefDiagSCMatrix operator*(double) const;
    //texi Matrix addition and subtraction.
    RefDiagSCMatrix operator+(const RefDiagSCMatrix&) const;
    RefDiagSCMatrix operator-(const RefDiagSCMatrix&) const;
    //texi Return the inverse of this.
    RefDiagSCMatrix i() const;
    //texi Return the generalized inverse of this.
    RefDiagSCMatrix gi() const;
    //texi These call the @code{SCMatrix} members of the same name
    // after checking for references to @code{0}.
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
    RefSCDimension dim() const;
    double trace() const;
    void print(ostream&) const;
    void print(const char*title=0,ostream&out=cout, int =10) const;
    //texi Returns the determinant of the referenced matrix.
    double determ() const;
    //texi Assign and examine matrix elements.
    DiagSCMatrixdouble operator()(int i) const;
};
//texi Allow multiplication with a scalar on the left.
RefDiagSCMatrix operator*(double,const RefDiagSCMatrix&);

SavableState_REF_dec(SCMatrixKit);

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
