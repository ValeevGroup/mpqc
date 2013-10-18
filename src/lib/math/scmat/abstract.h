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

#include <util/group/message.h>

#include <util/state/state.h>
#include <math/scmat/dim.h>
#include <math/scmat/block.h>
#include <util/misc/xml.h>
#include <iostream>

namespace sc {

class SCMatrix;
class SymmSCMatrix;
class DiagSCMatrix;
class SCVector;

class SCElementOp;
class SCElementOp2;
class SCElementOp3;

class RefSCDimension;

/** The SCMatrixKit abstract class acts as a factory for producing
matrices.  By using one of these, the program makes sure that all of the
matrices are consistent.  */
class SCMatrixKit: public DescribedClass {
  protected:
    Ref<MessageGrp> grp_;

  public:
    SCMatrixKit();
    SCMatrixKit(const Ref<KeyVal>&);
    ~SCMatrixKit();

    // these members are defined in repl.cc
    /** This returns the default matrix kit. It is currently a ReplSCMatrixKit, unless the
        the user has changed it with with set_default_matrixkit. */
    static SCMatrixKit* default_matrixkit();
    /// change the default matrix kit
    static void set_default_matrixkit(const Ref<SCMatrixKit> &);

    Ref<MessageGrp> messagegrp() const;

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


/** The SCVector class is the abstract base class for
    double valued vectors. */
class SCVector: public DescribedClass, public XMLWritable {
  protected:
    RefSCDimension d;
    Ref<SCMatrixKit> kit_;
  public:
    SCVector(const RefSCDimension&, SCMatrixKit *);

    /// Save and restore this in an implementation independent way.
    virtual void save(StateOut&);
    virtual void restore(StateIn&);
    virtual void write_xml(boost::property_tree::ptree& pt);

    /// Return the SCMatrixKit used to create this object.
    Ref<SCMatrixKit> kit() const { return kit_; }

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
    void convert(double* v) const { convert_p(v); }
    void convert(SCVector*v) { convert_v(v); }
    virtual void convert_p(double* v) const;
    virtual void convert_v(SCVector*);
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
    virtual void element_op(const Ref<SCElementOp>&) = 0;
    virtual void element_op(const Ref<SCElementOp2>&,
                            SCVector*) = 0;
    virtual void element_op(const Ref<SCElementOp3>&,
                            SCVector*,SCVector*) = 0;
    /// Print out the vector.
    void print(std::ostream&o=ExEnv::out0()) const;
    void print(const char* title=0,std::ostream&out=ExEnv::out0(),int=10) const;
    virtual void vprint(const char*title=0,std::ostream&out=ExEnv::out0(),
                        int=10) const = 0;

    /// Returns the message group used by the matrix kit
    Ref<MessageGrp> messagegrp() const;

    /** Returns iterators for the local (rapidly accessible) blocks used in
        this vector.  Only one iterator is allowed for a matrix is it has
        Accum or Write access is allowed.  Multiple Read iterators are
        permitted. */
    virtual Ref<SCMatrixSubblockIter> local_blocks(
        SCMatrixSubblockIter::Access) = 0;
    /// Returns iterators for the all blocks used in this vector.
    virtual Ref<SCMatrixSubblockIter> all_blocks(SCMatrixSubblockIter::Access) = 0;
};

/** The SCMatrix class is the abstract base class for general double valued
    n by m matrices.  For symmetric matrices use SymmSCMatrix and for
    diagonal matrices use DiagSCMatrix. */
class SCMatrix: public DescribedClass {
  protected:
    RefSCDimension d1,d2;
    Ref<SCMatrixKit> kit_;
  public:
    /// types of matrix transforms. Only real-valued matrices are assumed.
    enum Transform { NormalTransform = 0, TransposeTransform = 1 };

    // concrete functions (some can be overridden)
    SCMatrix(const RefSCDimension&, const RefSCDimension&, SCMatrixKit *);
    virtual ~SCMatrix();

    /// Save and restore this in an implementation independent way.
    virtual void save(StateOut&);
    virtual void restore(StateIn&);

    /// Return the SCMatrixKit used to create this object.
    Ref<SCMatrixKit> kit() const { return kit_; }

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
    void convert(double*a) const { convert_p(a); }
    void convert(double**a) const { convert_pp(a); }
    void convert(SCMatrix*a) { convert_r(a); }
    virtual void convert_p(double*) const;
    virtual void convert_pp(double**) const;
    virtual void convert_r(SCMatrix*a);
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
    virtual void accumulate_subblock(SCMatrix *m, int br, int er, int bc, int ec,
                                     int source_br = 0, int source_bc = 0) =0;

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
    void accumulate_product(SymmSCMatrix*m1,DiagSCMatrix*m2)
        { accumulate_product_sd(m1,m2); }
    void accumulate_product(DiagSCMatrix*m1,SymmSCMatrix*m2)
        { accumulate_product_ds(m1,m2); }
    virtual void accumulate_product_rr(SCMatrix*,SCMatrix*) = 0;
    virtual void accumulate_product_rs(SCMatrix*,SymmSCMatrix*);
    virtual void accumulate_product_rd(SCMatrix*,DiagSCMatrix*);
    virtual void accumulate_product_sr(SymmSCMatrix*,SCMatrix*);
    virtual void accumulate_product_dr(DiagSCMatrix*,SCMatrix*);
    virtual void accumulate_product_ss(SymmSCMatrix*,SymmSCMatrix*);
    virtual void accumulate_product_sd(SymmSCMatrix*,DiagSCMatrix*);
    virtual void accumulate_product_ds(DiagSCMatrix*,SymmSCMatrix*);
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
    /** Return the generalized inverse of this using SVD decomposition.
    \param condition_number_threshold specifies the extent of truncation of the singular numbers.
    The greater the threshold the fewer singular numbers will be discarded, but numerical errors
    will result for threshold reaching the inverse of the machine precision.
    */
    virtual void gen_invert_this(double condition_number_threshold = 1e8);

    /** Schmidt orthogonalize this.  S is the overlap matrix.
        n is the number of columns to orthogonalize. */
    virtual void schmidt_orthog(SymmSCMatrix*, int n) =0;

    /** Schmidt orthogonalize this.  S is the overlap matrix.  tol is the
        tolerance.  The number of linearly independent vectors is
        returned. */
    virtual int schmidt_orthog_tol(SymmSCMatrix*, double tol, double*res=0)=0;

    /// Perform the element operation op on each element of this.
    virtual void element_op(const Ref<SCElementOp>&) = 0;
    virtual void element_op(const Ref<SCElementOp2>&,
                            SCMatrix*) = 0;
    virtual void element_op(const Ref<SCElementOp3>&,
                            SCMatrix*,SCMatrix*) = 0;
    /// Print out the matrix.
    void print(std::ostream&o=ExEnv::out0()) const;
    void print(const char* title=0,std::ostream& out=ExEnv::out0(),
               int =10) const;
    virtual void vprint(const char*title=0,
                        std::ostream&out=ExEnv::out0(),int =10) const = 0;

    /// Returns the message group used by the matrix kit
    Ref<MessageGrp> messagegrp() const;

    /** Returns iterators for the local (rapidly accessible)
        blocks used in this matrix. */
    virtual Ref<SCMatrixSubblockIter> local_blocks(
        SCMatrixSubblockIter::Access) = 0;
    /// Returns iterators for the all blocks used in this matrix.
    virtual Ref<SCMatrixSubblockIter> all_blocks(
        SCMatrixSubblockIter::Access) = 0;
};

/** The SymmSCMatrix class is the abstract base class for symmetric
    double valued matrices. */
class SymmSCMatrix: public DescribedClass {
  protected:
    RefSCDimension d;
    Ref<SCMatrixKit> kit_;
  public:
    SymmSCMatrix(const RefSCDimension&, SCMatrixKit *);
    ~SymmSCMatrix();

    /// Return the SCMatrixKit object that created this object.
    Ref<SCMatrixKit> kit() const { return kit_; }

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
    void convert(double*a) const { convert_p(a); }
    void convert(double**a) const { convert_pp(a); }
    void convert(SymmSCMatrix*a) { convert_s(a); }
    virtual void convert_p(double*) const;
    virtual void convert_pp(double**) const;
    virtual void convert_s(SymmSCMatrix*);
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
    /** Solve generalized eigensystem for this with metric s, placing the eigenvalues in d and the eigenvectors
        in m.
        \exception AlgorithmException eigensystem could not be solved within available precision. Check condition number of s.
      */
    virtual void eigensystem(SymmSCMatrix*s,DiagSCMatrix*d,SCMatrix*m) = 0;
    /// Sum m into this.
    virtual void accumulate(const SymmSCMatrix* m) = 0;
    /// Sum into a + a.t()
    virtual void accumulate_symmetric_sum(SCMatrix* a) = 0;
    /// Sum into this a * a.t()
    virtual void accumulate_symmetric_product(SCMatrix* a);
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
    /// Return the generalized inverse of this using SVD decomposition. \sa SCMatrix::gen_invert_this()
    virtual void gen_invert_this(double condition_number_threshold = 1e8) = 0;

    //@{ Perform the element operation op on each element of this. Note that the operation is
    ///   only applied to the unique elements of this. For example, to compute the sum of all matrix
    ///   elements you need to do the following:
    ///   \code
    ///   SymmSCMatrix* A;   // presumed initialized elsewhere
    ///   A->element_op(new SCElementScaleDiagonal(0.5));   // scale the diagonal by 1/2
    ///   SCElementSum* sum_op = new SCElementSum;
    ///   A->element_op(sum_op);
    ///   std::cout << "Sum of element of matrix A = " << sum_op->result() * 2.0 << std::endl;
    ///   \endcode
    virtual void element_op(const Ref<SCElementOp>&) = 0;
    virtual void element_op(const Ref<SCElementOp2>&,
                            SymmSCMatrix*) = 0;
    virtual void element_op(const Ref<SCElementOp3>&,
                            SymmSCMatrix*,SymmSCMatrix*) = 0;
    //@}
    /// Print out the matrix.
    void print(std::ostream&o=ExEnv::out0()) const;
    void print(const char* title=0,std::ostream& out=ExEnv::out0(),
               int =10) const;
    virtual void vprint(const char* title=0,
                        std::ostream& out=ExEnv::out0(), int =10) const;

    /// Returns the message group used by the matrix kit
    Ref<MessageGrp> messagegrp() const;

    /** Returns iterators for the local (rapidly accessible)
        blocks used in this matrix. */
    virtual Ref<SCMatrixSubblockIter> local_blocks(
        SCMatrixSubblockIter::Access) = 0;
    /// Returns iterators for the all blocks used in this matrix.
    virtual Ref<SCMatrixSubblockIter> all_blocks(
        SCMatrixSubblockIter::Access) = 0;
};

/** The SymmSCMatrix class is the abstract base class for diagonal double
    valued matrices.  */
class DiagSCMatrix: public DescribedClass {
  protected:
    RefSCDimension d;
    Ref<SCMatrixKit> kit_;
  public:
    DiagSCMatrix(const RefSCDimension&, SCMatrixKit *);
    ~DiagSCMatrix();

    /// Return the SCMatrixKit used to create this object.
    Ref<SCMatrixKit> kit() const { return kit_; }

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
    void assign(DiagSCMatrix*d_a) { assign_d(d_a); }
    /// Overridden to implement the assign members.
    virtual void assign_val(double val);
    virtual void assign_p(const double*);
    virtual void assign_d(DiagSCMatrix*);
    /// Like the assign member, but this writes values to the argument.
    void convert(double*a) const { convert_p(a); }
    void convert(DiagSCMatrix*a) { convert_d(a); }
    virtual void convert_p(double*) const;
    virtual void convert_d(DiagSCMatrix*);
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
    /// Return the generalized inverse of this using SVD decomposition. \sa SCMatrix::gen_invert_this()
    virtual void gen_invert_this(double condition_number_threshold = 1e8) = 0;
    /// Perform the element operation op on each element of this.
    virtual void element_op(const Ref<SCElementOp>&) = 0;
    virtual void element_op(const Ref<SCElementOp2>&,
                            DiagSCMatrix*) = 0;
    virtual void element_op(const Ref<SCElementOp3>&,
                            DiagSCMatrix*,DiagSCMatrix*) = 0;
    /// Print out the matrix.
    void print(std::ostream&o=ExEnv::out0()) const;
    void print(const char* title=0,
               std::ostream& out=ExEnv::out0(), int =10) const;
    virtual void vprint(const char* title=0,
                        std::ostream& out=ExEnv::out0(), int =10) const;

    /// Returns the message group used by the matrix kit
    Ref<MessageGrp> messagegrp() const;

    /** Returns iterators for the local (rapidly accessible)
        blocks used in this matrix. */
    virtual Ref<SCMatrixSubblockIter> local_blocks(
        SCMatrixSubblockIter::Access) = 0;
    /// Returns iterators for the all blocks used in this matrix.
    virtual Ref<SCMatrixSubblockIter> all_blocks(
        SCMatrixSubblockIter::Access) = 0;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
