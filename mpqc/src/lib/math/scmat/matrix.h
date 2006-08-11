//
// matrix.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifndef _math_scmat_matrix_h
#define _math_scmat_matrix_h
#ifdef __GNUC__
#pragma interface
#endif

#include <iostream>

#include <math/scmat/abstract.h>

namespace sc {

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
/** The RefSCVector class is a smart pointer to an SCVector specialization.
 */
class RefSCVector: public Ref<SCVector> {
    // standard overrides
  public:
    /** Initializes the vector pointer to 0.  The reference must
        be initialized before it is used. */
    RefSCVector();
    /// Make this and v refer to the same SCVector.
    RefSCVector(const RefSCVector& v);
    /// Make this refer to v.
    RefSCVector(SCVector *v);
    // don't allow automatic conversion from any reference to a
    // described class
    ~RefSCVector();
    /// Make this refer to v.
    RefSCVector& operator=(SCVector* v);
    /// Make this and v refer to the same SCVector.
    RefSCVector& operator=(const RefSCVector& v);

    // vector specific members
  public:
    /** Create a vector with dimension dim.  The data values
        are undefined. */
    RefSCVector(const RefSCDimension& dim,const Ref<SCMatrixKit>&);

    /// Return an l-value that can be used to assign or retrieve an element.
    SCVectordouble operator()(int) const;
    /// Return an l-value that can be used to assign or retrieve an element.
    SCVectordouble operator[](int) const;
    /// Add two vectors.
    RefSCVector operator+(const RefSCVector&a) const;
    /// Subtract two vectors.
    RefSCVector operator-(const RefSCVector&a) const;
    /// Scale a vector.
    RefSCVector operator*(double) const;
    /// Return the outer product between this and v.
    RefSCMatrix outer_product(const RefSCVector& v) const;
    /// The outer product of this with itself is a symmetric matrix.
    RefSymmSCMatrix symmetric_outer_product() const;

    void set_element(int i,double val) const;
    void accumulate_element(int i,double val) const;
    double get_element(int) const;
    int n() const;
    RefSCDimension dim() const;
    Ref<SCMatrixKit> kit() const;
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
    void element_op(const Ref<SCElementOp>& op) const;
    void element_op(const Ref<SCElementOp2>&,
                    const RefSCVector&) const;
    void element_op(const Ref<SCElementOp3>&,
                    const RefSCVector&,
                    const RefSCVector&) const;
    void print(std::ostream&out) const;
    void print(const char*title=0,
               std::ostream&out=ExEnv::out0(), int precision=10) const;
    void save(StateOut&);
    /// Restores the matrix from StateIn object. The vector must have been initialized already.
    void restore(StateIn&);
};
RefSCVector operator*(double,const RefSCVector&);

class RefSymmSCMatrix;
class RefDiagSCMatrix;
/** The RefSCMatrix class is a smart pointer to an SCMatrix
    specialization.
*/
class RefSCMatrix: public Ref<SCMatrix> {
    // standard overrides
  public:
    /** Initializes the matrix pointer to 0.  The reference must
        be initialized before it is used. */
    RefSCMatrix();
    /// Make this and m refer to the same SCMatrix.
    RefSCMatrix(const RefSCMatrix& m);
    /// Make this refer to m.
     RefSCMatrix(SCMatrix* m);
    ~RefSCMatrix();
    /// Make this refer to m.
    RefSCMatrix& operator=(SCMatrix* m);
    /// Make this and m refer to the same matrix.
    RefSCMatrix& operator=(const RefSCMatrix& m);

    // matrix specific members
  public:
    /** Create a vector with dimension d1 by d2.
        The data values are undefined. */
    RefSCMatrix(const RefSCDimension& d1,const RefSCDimension& d2,
                const Ref<SCMatrixKit>&);

    /// Multiply this by a vector and return a vector.
    RefSCVector operator*(const RefSCVector&) const;

    /// Multiply this by a matrix and return a matrix.
    RefSCMatrix operator*(const RefSCMatrix&) const;
    RefSCMatrix operator*(const RefSymmSCMatrix&) const;
    RefSCMatrix operator*(const RefDiagSCMatrix&) const;

    /// Multiply this by a scalar and return the result.
    RefSCMatrix operator*(double) const;

    /// Matrix addition.
    RefSCMatrix operator+(const RefSCMatrix&) const;
    /// Matrix subtraction.
    RefSCMatrix operator-(const RefSCMatrix&) const;

    /// Return the transpose of this.
    RefSCMatrix t() const;
    /// Return the inverse of this.
    RefSCMatrix i() const;
    /// Return the generalized inverse of this.
    RefSCMatrix gi() const;

    /** These call the SCMatrix members of the same name
        after checking for references to 0. */
    RefSCMatrix clone() const;
    RefSCMatrix copy() const;

    RefSCMatrix get_subblock(int br, int er, int bc, int ec);
    void assign_subblock(const RefSCMatrix&, int br, int er, int bc, int ec,
                         int source_br = 0, int source_bc = 0);
    void accumulate_subblock(const RefSCMatrix&, int, int, int, int,
                             int source_br = 0, int source_bc = 0);
    RefSCVector get_row(int) const;
    RefSCVector get_column(int) const;
    void assign_row(const RefSCVector&, int) const;
    void assign_column(const RefSCVector&, int) const;
    void accumulate_row(const RefSCVector&, int) const;
    void accumulate_column(const RefSCVector&, int) const;

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
    void accumulate(const RefSymmSCMatrix&) const;
    void accumulate(const RefDiagSCMatrix&) const;
    void element_op(const Ref<SCElementOp>&) const;
    void element_op(const Ref<SCElementOp2>&,
                    const RefSCMatrix&) const;
    void element_op(const Ref<SCElementOp3>&,
                    const RefSCMatrix&,
                    const RefSCMatrix&) const;
    int nrow() const;
    int ncol() const;
    RefSCDimension rowdim() const;
    RefSCDimension coldim() const;
    Ref<SCMatrixKit> kit() const;
    void set_element(int,int,double) const;
    void accumulate_element(int,int,double) const;
    double get_element(int,int) const;
    void print(std::ostream&) const;
    void print(const char*title=0,
               std::ostream&out=ExEnv::out0(), int =10) const;
    double trace() const;
    void save(StateOut&);
    /// Restores the matrix from StateIn object. The matrix must have been initialized already.
    void restore(StateIn&);

    /** Compute the singular value decomposition,  this = U sigma V.t().
        The dimension of sigma is the smallest dimension of this.  U, V,
        and sigma must already have the correct dimensions and are
        overwritten. */
    void svd(const RefSCMatrix &U,
             const RefDiagSCMatrix &sigma,
             const RefSCMatrix &V);
    /** Solves this x = v.  Overwrites v with x. */
    double solve_lin(const RefSCVector& v) const;
    /// Returns the determinant of the referenced matrix.
    double determ() const;
    /// Assign and examine matrix elements.
    SCMatrixdouble operator()(int i,int j) const;

    /** If this matrix is blocked return the number of blocks.
     * Otherwise return 1.
     */
    int nblock() const;
    /** If this matrix is blocked return block i.
     * Otherwise return this as block 0.
     */
    RefSCMatrix block(int i) const;
};
/// Allow multiplication with a scalar on the left.
RefSCMatrix operator*(double,const RefSCMatrix&);

/** The RefSymmSCMatrix class is a smart pointer to an SCSymmSCMatrix
    specialization.  */
class RefSymmSCMatrix: public Ref<SymmSCMatrix> {
    // standard overrides
  public:
    /** Initializes the matrix pointer to 0.  The reference must
        be initialized before it is used. */
    RefSymmSCMatrix();
    /// Make this and m refer to the same SCMatrix.
    RefSymmSCMatrix(const RefSymmSCMatrix& m);
    /// Make this refer to m.
    RefSymmSCMatrix(SymmSCMatrix *m);
    ~RefSymmSCMatrix();
    /// Make this refer to m.
    RefSymmSCMatrix& operator=(SymmSCMatrix* m);
    /// Make this and m refer to the same matrix.
    RefSymmSCMatrix& operator=(const RefSymmSCMatrix& m);

    // matrix specific members
  public:
    /** Create a vector with dimension d by d.
        The data values are undefined. */
    RefSymmSCMatrix(const RefSCDimension& d,const Ref<SCMatrixKit>&);
    /// Multiply this by a matrix and return a matrix.
    RefSCMatrix operator*(const RefSCMatrix&) const;
    RefSCMatrix operator*(const RefSymmSCMatrix&) const;
    RefSCMatrix operator*(const RefDiagSCMatrix&) const;
    /// Multiply this by a vector and return a vector.
    RefSCVector operator*(const RefSCVector&a) const;
    RefSymmSCMatrix operator*(double) const;
    /// Matrix addition and subtraction.
    RefSymmSCMatrix operator+(const RefSymmSCMatrix&) const;
    RefSymmSCMatrix operator-(const RefSymmSCMatrix&) const;
    /// Return the inverse of this.
    RefSymmSCMatrix i() const;
    /// Return the generalized inverse of this.
    RefSymmSCMatrix gi() const;
    /** These call the SCMatrix members of the same name after checking for
        references to 0. */
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
    /// Add a * b * a.t() to this.
    void accumulate_transform(const RefSCMatrix&a,const RefSymmSCMatrix&b,
                        SCMatrix::Transform = SCMatrix::NormalTransform) const;
    void accumulate_transform(const RefSCMatrix&a,const RefDiagSCMatrix&b,
                        SCMatrix::Transform = SCMatrix::NormalTransform) const;
    void accumulate_transform(const RefSymmSCMatrix&a,
                              const RefSymmSCMatrix&b) const;

    void randomize() const;
    void assign(const RefSymmSCMatrix&) const;
    void scale(double) const;
    void assign(double) const;
    void assign(const double*) const;
    void assign(const double**) const;
    void convert(double*) const;
    void convert(double**) const;
    void accumulate(const RefSymmSCMatrix&) const;
    void element_op(const Ref<SCElementOp>&) const;
    void element_op(const Ref<SCElementOp2>&,
                    const RefSymmSCMatrix&) const;
    void element_op(const Ref<SCElementOp3>&,
                    const RefSymmSCMatrix&,
                    const RefSymmSCMatrix&) const;
    double trace() const;
    int n() const;
    RefSCDimension dim() const;
    Ref<SCMatrixKit> kit() const;
    void print(std::ostream&) const;
    void print(const char*title=0,
               std::ostream&out=ExEnv::out0(), int =10) const;
    void save(StateOut&);
    /// Restores the matrix from StateIn object. The matrix must have been initialized already.
    void restore(StateIn&);

    /** Solves this x = v.  Overwrites v with x. */
    double solve_lin(const RefSCVector&) const;
    /// Returns the determinant of the referenced matrix.
    double determ() const;
    /// Returns the eigenvalues of the reference matrix.
    RefDiagSCMatrix eigvals() const;
    /// Returns the eigenvectors of the reference matrix.
    RefSCMatrix eigvecs() const;
    /** Sets eigvals to the eigenvalues and eigvecs
        to the eigenvalues and eigenvectors of the referenced matrix.
        The result satisfies eigvecs * eigvals * eigvecs.t() = (*this).  */
    void diagonalize(const RefDiagSCMatrix& eigvals,
                     const RefSCMatrix& eigvecs) const;
    /// Assign and examine matrix elements.
    SymmSCMatrixdouble operator()(int i,int j) const;
    /** If this matrix is blocked return the number of blocks.
     * Otherwise return 1.
     */
    int nblock() const;
    /** If this matrix is blocked return block i.
     * Otherwise return this as block 0.
     */
    RefSymmSCMatrix block(int i) const;
};
/// Allow multiplication with a scalar on the left.
RefSymmSCMatrix operator*(double,const RefSymmSCMatrix&);

/** The RefDiagSCMatrix class is a smart pointer to an DiagSCMatrix
    specialization. */
class RefDiagSCMatrix: public Ref<DiagSCMatrix> {
    // standard overrides
  public:
    /** Initializes the matrix pointer to 0.  The reference must
        be initialized before it is used. */
    RefDiagSCMatrix();
    /// Make this and m refer to the same SCMatrix.
    RefDiagSCMatrix(const RefDiagSCMatrix& m);
    /// Make this refer to m.
    RefDiagSCMatrix(DiagSCMatrix *m);
    ~RefDiagSCMatrix();
    /// Make this refer to m.
    RefDiagSCMatrix& operator=(DiagSCMatrix* m);
    /// Make this and m refer to the same matrix.
    RefDiagSCMatrix& operator=(const RefDiagSCMatrix & m);

    // matrix specific members
  public:
    /** Create a diagonal matrix with dimension d by d.  The data values
        are undefined. */
    RefDiagSCMatrix(const RefSCDimension&,const Ref<SCMatrixKit>&);
    /// Multiply this by a matrix and return a matrix.
    RefSCMatrix operator*(const RefSCMatrix&) const;
    RefSCMatrix operator*(const RefSymmSCMatrix&) const;
    RefDiagSCMatrix operator*(const RefDiagSCMatrix&) const;
    RefDiagSCMatrix operator*(double) const;
    /// Matrix addition and subtraction.
    RefDiagSCMatrix operator+(const RefDiagSCMatrix&) const;
    RefDiagSCMatrix operator-(const RefDiagSCMatrix&) const;
    /// Return the inverse of this.
    RefDiagSCMatrix i() const;
    /// Return the generalized inverse of this.
    RefDiagSCMatrix gi() const;
    /// These call the SCMatrix members of the same name
    /// after checking for references to 0.
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
    void element_op(const Ref<SCElementOp>&) const;
    void element_op(const Ref<SCElementOp2>&,
                    const RefDiagSCMatrix&) const;
    void element_op(const Ref<SCElementOp3>&,
                    const RefDiagSCMatrix&,
                    const RefDiagSCMatrix&) const;
    int n() const;
    RefSCDimension dim() const;
    Ref<SCMatrixKit> kit() const;
    double trace() const;
    void print(std::ostream&) const;
    void print(const char*title=0,
               std::ostream&out=ExEnv::out0(), int =10) const;
    void save(StateOut&);
    /// Restores the matrix from StateIn object. The matrix must have been initialized already.
    void restore(StateIn&);
    /// Returns the determinant of the referenced matrix.
    double determ() const;
    /// Assign and examine matrix elements.
    DiagSCMatrixdouble operator()(int i) const;
    /** If this matrix is blocked return the number of blocks.
     * Otherwise return 1.
     */
    int nblock() const;
    /** If this matrix is blocked return block i.
     * Otherwise return this as block 0.
     */
    RefDiagSCMatrix block(int i) const;
};
/// Allow multiplication with a scalar on the left.
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

}

#ifdef INLINE_FUNCTIONS
#include <math/scmat/matrix_i.h>
#endif

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
