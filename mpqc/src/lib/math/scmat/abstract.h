//
// abstract.h
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

SavableState_REF_fwddec(SCElementOp);
SavableState_REF_fwddec(SCElementOp2);
SavableState_REF_fwddec(SCElementOp3);

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
    /** This returns a LocalSCMatrixKit, unless the
        default has been changed with set_default_matrixkit. */
    static SCMatrixKit* default_matrixkit();
    static void set_default_matrixkit(const RefSCMatrixKit &);

    RefMessageGrp messagegrp() const;

    /// Given the dimensions, create matrices or vectors.
    virtual SCMatrix* matrix(const RefSCDimension&,const RefSCDimension&) = 0;
    virtual SymmSCMatrix* symmmatrix(const RefSCDimension&) = 0;
    virtual DiagSCMatrix* diagmatrix(const RefSCDimension&) = 0;
    virtual SCVector* vector(const RefSCDimension&) = 0;

    /** Given the dimensions and a StateIn object,
        restore matrices or vectors. */
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

/** The SCVector class is the abstract base class for
    double valued vectors. */
class SCVector: public DescribedClass {
#   define CLASSNAME SCVector
#   include <util/class/classda.h>
  protected:
    RefSCDimension d;
    RefSCMatrixKit kit_;
  public:
    SCVector(const RefSCDimension&, SCMatrixKit *);

    /// Save and restore this in an implementation independent way.
    virtual void save(StateOut&);
    virtual void restore(StateIn&);

    /// Return the SCMatrixKit used to create this object.
    RefSCMatrixKit kit() const { return kit_; }

    // concrete functions (some can be overridden)
    /// Return a vector with the same dimension and same elements.
    virtual SCVector* copy();
    /// Return a vector with the same dimension but uninitialized memory.
    virtual SCVector* clone();

    virtual ~SCVector();
    /// Return the length of the vector.
    int n() const { return d->n(); }
    /// Return the maximum absolute value element of this vector.
    virtual double maxabs() const;
    /// Normalize this.
    virtual void normalize();
    /// Assign each element to a random number between -1 and 1
    virtual void randomize();
    /// Assign all elements of this to val.
    void assign(double val) { assign_val(val); }
    /// Assign element i to v[i] for all i.
    void assign(const double* v) { assign_p(v); }
    /** Make this have the same elements as v.  The dimensions must
        match. */
    void assign(SCVector* v) { assign_v(v); }
    /// Overridden to implement the assign functions.
    virtual void assign_val(double val);
    virtual void assign_p(const double* v);
    virtual void assign_v(SCVector *v);
    /// Assign v[i] to element i for all i.
    virtual void convert(double* v) const;
    /** Convert an SCVector of a different specialization
        to this specialization and possibly accumulate the data. */
    virtual void convert(SCVector*);
    virtual void convert_accumulate(SCVector*);
    /// Multiply each element by val.
    virtual void scale(double val);

    /// Return the RefSCDimension corresponding to this vector.
    RefSCDimension dim() const { return d; }
    /// Set element i to val.
    virtual void set_element(int i,double val) = 0;
    /// Add val to element i.
    virtual void accumulate_element(int,double) = 0;
    /// Return the value of element i.
    virtual double get_element(int i) const = 0;
    /// Sum the result of m times v into this.
    void accumulate_product(SymmSCMatrix* m, SCVector* v)
        { accumulate_product_sv(m,v); }
    void accumulate_product(SCMatrix* m, SCVector* v)
        {  accumulate_product_rv(m,v); }
    virtual void accumulate_product_sv(SymmSCMatrix* m, SCVector* v);
    virtual void accumulate_product_rv(SCMatrix* m, SCVector* v) = 0;
    /// Sum v into this.
    virtual void accumulate(const SCVector*v) = 0;
    /// Sum m into this.  One of m's dimensions must be 1.
    virtual void accumulate(const SCMatrix*m) = 0;
    /// Return the dot product.
    virtual double scalar_product(SCVector*) = 0;
    /// Perform the element operation op on each element of this.
    virtual void element_op(const RefSCElementOp&) = 0;
    virtual void element_op(const RefSCElementOp2&,
                            SCVector*) = 0;
    virtual void element_op(const RefSCElementOp3&,
                            SCVector*,SCVector*) = 0;
    /// Print out the vector.
    void print(ostream&o=cout) const;
    void print(const char* title=0,ostream& out=cout, int =10) const;
    virtual void vprint(const char*title=0,ostream&out=cout,int =10) const = 0;

    /// Returns the message group used by the matrix kit
    RefMessageGrp messagegrp() const;
    
    /** Returns iterators for the local (rapidly accessible) blocks used in
        this vector.  Only one iterator is allowed for a matrix is it has
        Accum or Write access is allowed.  Multiple Read iterators are
        permitted. */
    virtual RefSCMatrixSubblockIter local_blocks(
        SCMatrixSubblockIter::Access) = 0;
    /// Returns iterators for the all blocks used in this vector.
    virtual RefSCMatrixSubblockIter all_blocks(SCMatrixSubblockIter::Access) = 0;
};

/** The SCMatrix class is the abstract base class for general double valued
    n by m matrices.  For symmetric matrices use SymmSCMatrix and for
    diagonal matrices use DiagSCMatrix. */
class SCMatrix: public DescribedClass {
#   define CLASSNAME SCMatrix
#   include <util/class/classda.h>
  protected:
    RefSCDimension d1,d2;
    RefSCMatrixKit kit_;
  public:
    // used to control transformations
    enum Transform { NormalTransform = 0, TransposeTransform = 1 };
    
    // concrete functions (some can be overridden)
    SCMatrix(const RefSCDimension&, const RefSCDimension&, SCMatrixKit *);
    virtual ~SCMatrix();

    /// Save and restore this in an implementation independent way.
    virtual void save(StateOut&);
    virtual void restore(StateIn&);

    /// Return the SCMatrixKit used to create this object.
    RefSCMatrixKit kit() const { return kit_; }

    /// Return the number of rows.
    int nrow() const { return d1->n(); }
    /// Return the number of columns.
    int ncol() const { return d2->n(); }
    /// Return the maximum absolute value element.
    virtual double maxabs() const;
    /// Assign each element to a random number between -1 and 1
    virtual void randomize();
    /// Set all elements to val.
    void assign(double val) { assign_val(val); }
    /// Assign element i, j to m[ir*nrow()+j].
    void assign(const double* m) { assign_p(m); }
    /// Assign element i, j to m[i][j].
    void assign(const double** m) { assign_pp(m); }
    /// Make this have the same elements as m. The dimensions must match.
    void assign(SCMatrix* m) { assign_r(m); }
    /// Overridden to implement to assign members.
    virtual void assign_val(double val);
    virtual void assign_p(const double* m);
    virtual void assign_pp(const double** m);
    virtual void assign_r(SCMatrix* m);
    /** Like the assign members, but these write values to the
        arguments. */
    virtual void convert(double*) const;
    virtual void convert(double**) const;
    /** Convert an SCMatrix of a different specialization to this
        specialization and possibly accumulate the data. */
    virtual void convert(SCMatrix*);
    virtual void convert_accumulate(SCMatrix*);
    /// Multiply all elements by val.
    virtual void scale(double val);
    /// Scale the diagonal elements by val.
    virtual void scale_diagonal(double val);
    /// Shift the diagonal elements by val.
    virtual void shift_diagonal(double val);
    /// Make this equal to the unit matrix.
    virtual void unit();
    /// Return a matrix with the same dimension and same elements.
    virtual SCMatrix* copy();
    /// Return a matrix with the same dimension but uninitialized memory.
    virtual SCMatrix* clone();

    // pure virtual functions
    /// Return the row or column dimension.
    RefSCDimension rowdim() const { return d1; }
    RefSCDimension coldim() const { return d2; }
    /// Return or modify an element.
    virtual double get_element(int,int) const = 0;
    virtual void set_element(int,int,double) = 0;
    virtual void accumulate_element(int,int,double) = 0;
    
    /** Return a subblock of this.  The subblock is defined as
        the rows starting at br and ending at er, and the
        columns beginning at bc and ending at ec. */
    virtual SCMatrix * get_subblock(int br, int er, int bc, int ec) =0;

    /// Assign m to a subblock of this.
    virtual void assign_subblock(SCMatrix *m, int, int, int, int, int=0, int=0) =0;

    /// Sum m into a subblock of this.
    virtual void accumulate_subblock(SCMatrix *m, int, int, int, int, int=0,int=0) =0;
    
    /// Return a row or column of this.
    virtual SCVector * get_row(int i) =0;
    virtual SCVector * get_column(int i) =0;

    /// Assign v to a row or column of this.
    virtual void assign_row(SCVector *v, int i) =0;
    virtual void assign_column(SCVector *v, int i) =0;
    
    /// Sum v to a row or column of this.
    virtual void accumulate_row(SCVector *v, int i) =0;
    virtual void accumulate_column(SCVector *v, int i) =0;
    
    /// Sum m into this.
    virtual void accumulate(const SCMatrix* m) = 0;
    virtual void accumulate(const SymmSCMatrix* m) = 0;
    virtual void accumulate(const DiagSCMatrix* m) = 0;
    virtual void accumulate(const SCVector*) = 0;
    /// Sum into this the products of various vectors or matrices.
    virtual void accumulate_outer_product(SCVector*,SCVector*) = 0;
    void accumulate_product(SCMatrix*m1,SCMatrix*m2)
        { accumulate_product_rr(m1,m2); }
    void accumulate_product(SCMatrix*m1,SymmSCMatrix*m2)
        { accumulate_product_rs(m1,m2); }
    void accumulate_product(SCMatrix*m1,DiagSCMatrix*m2)
        { accumulate_product_rd(m1,m2); }
    void accumulate_product(SymmSCMatrix*m1,SCMatrix*m2)
        { accumulate_product_sr(m1,m2); }
    void accumulate_product(DiagSCMatrix*m1,SCMatrix*m2)
        { accumulate_product_dr(m1,m2); }
    void accumulate_product(SymmSCMatrix*m1,SymmSCMatrix*m2)
        { accumulate_product_ss(m1,m2); }
    virtual void accumulate_product_rr(SCMatrix*,SCMatrix*) = 0;
    virtual void accumulate_product_rs(SCMatrix*,SymmSCMatrix*);
    virtual void accumulate_product_rd(SCMatrix*,DiagSCMatrix*);
    virtual void accumulate_product_sr(SymmSCMatrix*,SCMatrix*);
    virtual void accumulate_product_dr(DiagSCMatrix*,SCMatrix*);
    virtual void accumulate_product_ss(SymmSCMatrix*,SymmSCMatrix*);
    /// Transpose this.
    virtual void transpose_this() = 0;
    /// Return the trace.
    virtual double trace() =0;
    /// Invert this.
    virtual double invert_this() = 0;
    /// Return the determinant of this.  this is overwritten.
    virtual double determ_this() = 0;

    /** Compute the singular value decomposition for this,
        possibly destroying this. */
    virtual void svd_this(SCMatrix *U, DiagSCMatrix *sigma, SCMatrix *V);
    virtual double solve_this(SCVector*) = 0;
    virtual void gen_invert_this();

    /** Schmidt orthogonalize this.  S is the overlap matrix.
        n is the number of columns to orthogonalize. */
    virtual void schmidt_orthog(SymmSCMatrix*, int n) =0;
    
    /// Perform the element operation op on each element of this.
    virtual void element_op(const RefSCElementOp&) = 0;
    virtual void element_op(const RefSCElementOp2&,
                            SCMatrix*) = 0;
    virtual void element_op(const RefSCElementOp3&,
                            SCMatrix*,SCMatrix*) = 0;
    /// Print out the matrix.
    void print(ostream&o=cout) const;
    void print(const char* title=0,ostream& out=cout, int =10) const;
    virtual void vprint(const char*title=0,ostream&out=cout,int =10) const = 0;

    /// Returns the message group used by the matrix kit
    RefMessageGrp messagegrp() const;
    
    /** Returns iterators for the local (rapidly accessible)
        blocks used in this matrix. */
    virtual RefSCMatrixSubblockIter local_blocks(
        SCMatrixSubblockIter::Access) = 0;
    /// Returns iterators for the all blocks used in this matrix.
    virtual RefSCMatrixSubblockIter all_blocks(
        SCMatrixSubblockIter::Access) = 0;
};

/** The SymmSCMatrix class is the abstract base class for symmetric
    double valued matrices. */
class SymmSCMatrix: public DescribedClass {
#   define CLASSNAME SymmSCMatrix
#   include <util/class/classda.h>
  protected:
    RefSCDimension d;
    RefSCMatrixKit kit_;
  public:
    SymmSCMatrix(const RefSCDimension&, SCMatrixKit *);

    /// Return the SCMatrixKit object that created this object.
    RefSCMatrixKit kit() const { return kit_; }

    /// Save and restore this in an implementation independent way.
    virtual void save(StateOut&);
    virtual void restore(StateIn&);
    /// Return the maximum absolute value element of this vector.
    virtual double maxabs() const;
    /// Assign each element to a random number between -1 and 1
    virtual void randomize();
    /// Set all elements to val.
    void assign(double val) { assign_val(val); }
    /** Assign element i, j to m[i*(i+1)/2+j]. */
    void assign(const double* m) { assign_p(m); }
    /// Assign element i, j to m[i][j].
    void assign(const double** m) { assign_pp(m); }
    /** Make this have the same elements as m.  The dimensions must
        match. */
    void assign(SymmSCMatrix* m) { assign_s(m); }
    /// Overridden to implement the assign functions
    virtual void assign_val(double val);
    virtual void assign_p(const double* m);
    virtual void assign_pp(const double** m);
    virtual void assign_s(SymmSCMatrix* m);
    /// Like the assign members, but these write values to the arguments.
    virtual void convert(double*) const;
    virtual void convert(double**) const;
    /** Convert an SCSymmSCMatrix of a different specialization
        to this specialization and possibly accumulate the data. */
    virtual void convert(SymmSCMatrix*);
    virtual void convert_accumulate(SymmSCMatrix*);
    /// Multiply all elements by val.
    virtual void scale(double);
    /// Scale the diagonal elements by val.
    virtual void scale_diagonal(double);
    /// Shift the diagonal elements by val.
    virtual void shift_diagonal(double);
    /// Make this equal to the unit matrix.
    virtual void unit();
    /// Return the dimension.
    int n() const { return d->n(); }
    /// Return a matrix with the same dimension and same elements.
    virtual SymmSCMatrix* copy();
    /// Return a matrix with the same dimension but uninitialized memory.
    virtual SymmSCMatrix* clone();

    // pure virtual functions
    /// Return the dimension.
    RefSCDimension dim() const { return d; }
    /// Return or modify an element.
    virtual double get_element(int,int) const = 0;
    virtual void set_element(int,int,double) = 0;
    virtual void accumulate_element(int,int,double) = 0;

    /** Return a subblock of this.  The subblock is defined as the rows
        starting at br and ending at er, and the columns beginning at bc
        and ending at ec. */
    virtual SCMatrix * get_subblock(int br, int er, int bc, int ec) =0;
    virtual SymmSCMatrix * get_subblock(int br, int er) =0;

    /// Assign m to a subblock of this.
    virtual void assign_subblock(SCMatrix *m, int, int, int, int) =0;
    virtual void assign_subblock(SymmSCMatrix *m, int, int) =0;

    /// Sum m into a subblock of this.
    virtual void accumulate_subblock(SCMatrix *m, int, int, int, int) =0;
    virtual void accumulate_subblock(SymmSCMatrix *m, int, int) =0;

    /// Return a row of this.
    virtual SCVector * get_row(int i) =0;

    /// Assign v to a row of this.
    virtual void assign_row(SCVector *v, int i) =0;
    
    /// Sum v to a row of this.
    virtual void accumulate_row(SCVector *v, int i) =0;

    /** Diagonalize this, placing the eigenvalues in d and the eigenvectors
        in m. */
    virtual void diagonalize(DiagSCMatrix*d,SCMatrix*m) = 0;
    /// Sum m into this.
    virtual void accumulate(const SymmSCMatrix* m) = 0;
    /// Sum into this the products of various vectors or matrices.
    virtual void accumulate_symmetric_sum(SCMatrix*) = 0;
    virtual void accumulate_symmetric_product(SCMatrix*);
    virtual void accumulate_transform(SCMatrix*,SymmSCMatrix*,
                            SCMatrix::Transform = SCMatrix::NormalTransform);
    virtual void accumulate_transform(SCMatrix*,DiagSCMatrix*, 
                            SCMatrix::Transform = SCMatrix::NormalTransform);
    virtual void accumulate_transform(SymmSCMatrix*,SymmSCMatrix*);
    virtual void accumulate_symmetric_outer_product(SCVector*);
    /** Return the scalar obtained by multiplying this on the
        left and right by v. */
    virtual double scalar_product(SCVector* v);
    /// Return the trace.
    virtual double trace() = 0;
    /// Invert this.
    virtual double invert_this() = 0;
    /// Return the determinant of this.  this is overwritten.
    virtual double determ_this() = 0;

    virtual double solve_this(SCVector*) = 0;
    virtual void gen_invert_this() = 0;

    /// Perform the element operation op on each element of this.
    virtual void element_op(const RefSCElementOp&) = 0;
    virtual void element_op(const RefSCElementOp2&,
                            SymmSCMatrix*) = 0;
    virtual void element_op(const RefSCElementOp3&,
                            SymmSCMatrix*,SymmSCMatrix*) = 0;
    /// Print out the matrix.
    void print(ostream&o=cout) const;
    void print(const char* title=0,ostream& out=cout, int =10) const;
    virtual void vprint(const char* title=0,ostream& out=cout, int =10) const;

    /// Returns the message group used by the matrix kit
    RefMessageGrp messagegrp() const;
    
    /** Returns iterators for the local (rapidly accessible)
        blocks used in this matrix. */
    virtual RefSCMatrixSubblockIter local_blocks(
        SCMatrixSubblockIter::Access) = 0;
    /// Returns iterators for the all blocks used in this matrix.
    virtual RefSCMatrixSubblockIter all_blocks(
        SCMatrixSubblockIter::Access) = 0;
};

/** The SymmSCMatrix class is the abstract base class for diagonal double
    valued matrices.  */
class DiagSCMatrix: public DescribedClass {
#   define CLASSNAME DiagSCMatrix
#   include <util/class/classda.h>
  protected:
    RefSCDimension d;
    RefSCMatrixKit kit_;
  public:
    DiagSCMatrix(const RefSCDimension&, SCMatrixKit *);

    /// Return the SCMatrixKit used to create this object.
    RefSCMatrixKit kit() const { return kit_; }

    /// Save and restore this in an implementation independent way.
    virtual void save(StateOut&);
    virtual void restore(StateIn&);

    /// Return the maximum absolute value element of this vector.
    virtual double maxabs() const;
    /// Assign each element to a random number between -1 and 1
    virtual void randomize();
    /// Set all elements to val.
    void assign(double val) { assign_val(val); }
    /// Assign element i, i to m[i].
    void assign(const double*p) { assign_p(p); }
    /** Make this have the same elements as m.  The dimensions must
        match. */
    void assign(DiagSCMatrix*d) { assign_d(d); }
    /// Overridden to implement the assign members.
    virtual void assign_val(double val);
    virtual void assign_p(const double*);
    virtual void assign_d(DiagSCMatrix*);
    /// Like the assign member, but this writes values to the argument.
    virtual void convert(double*) const;
    /** Convert an SCDiagSCMatrix of a different specialization
        to this specialization and possibly accumulate the data. */
    virtual void convert(DiagSCMatrix*);
    virtual void convert_accumulate(DiagSCMatrix*);
    /// Multiply all elements by val.
    virtual void scale(double);
    /// Return the dimension.
    int n() const { return d->n(); }
    /// Return a matrix with the same dimension and same elements.
    virtual DiagSCMatrix* copy();
    /// Return a matrix with the same dimension but uninitialized memory.
    virtual DiagSCMatrix* clone();

    // pure virtual functions
    /// Return the dimension.
    RefSCDimension dim() const { return d; }
    /// Return or modify an element.
    virtual double get_element(int) const = 0;
    virtual void set_element(int,double) = 0;
    virtual void accumulate_element(int,double) = 0;
    /// Sum m into this.
    virtual void accumulate(const DiagSCMatrix* m) = 0;
    /// Return the trace.
    virtual double trace() = 0;
    /// Return the determinant of this.  this is overwritten.
    virtual double determ_this() = 0;
    /// Invert this.
    virtual double invert_this() = 0;
    /// Do a generalized inversion of this.
    virtual void gen_invert_this() = 0;
    /// Perform the element operation op on each element of this.
    virtual void element_op(const RefSCElementOp&) = 0;
    virtual void element_op(const RefSCElementOp2&,
                            DiagSCMatrix*) = 0;
    virtual void element_op(const RefSCElementOp3&,
                            DiagSCMatrix*,DiagSCMatrix*) = 0;
    /// Print out the matrix.
    void print(ostream&o=cout) const;
    void print(const char* title=0,ostream& out=cout, int =10) const;
    virtual void vprint(const char* title=0,ostream& out=cout, int =10) const;

    /// Returns the message group used by the matrix kit
    RefMessageGrp messagegrp() const;
    
    /** Returns iterators for the local (rapidly accessible)
        blocks used in this matrix. */
    virtual RefSCMatrixSubblockIter local_blocks(
        SCMatrixSubblockIter::Access) = 0;
    /// Returns iterators for the all blocks used in this matrix.
    virtual RefSCMatrixSubblockIter all_blocks(
        SCMatrixSubblockIter::Access) = 0;
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
