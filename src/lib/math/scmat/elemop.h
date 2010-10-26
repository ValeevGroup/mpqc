//
// elemop.h
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

#ifndef _math_scmat_elemop_h
#define _math_scmat_elemop_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/state/state.h>
#include <util/group/message.h>
#include <cmath>

namespace sc {

class SCMatrixBlock;
class SCMatrixBlockIter;
class SCMatrixRectBlock;
class SCMatrixLTriBlock;
class SCMatrixDiagBlock;
class SCVectorSimpleBlock;
class SCMatrixRectSubBlock;
class SCMatrixLTriSubBlock;
class SCMatrixDiagSubBlock;
class SCVectorSimpleSubBlock;

class SCMatrix;
class SymmSCMatrix;
class DiagSCMatrix;
class SCVector;

/** Objects of class SCElementOp are used to perform operations on the
    elements of matrices.  When the SCElementOp object is given to the
    element_op member of a matrix, each block the matrix is passed to one
    of the process, process_base, or process_base members. */
class SCElementOp: public SavableState {
  public:
    SCElementOp();
    SCElementOp(StateIn&s): SavableState(s) {}
    virtual ~SCElementOp();
    /** If duplicates of the SCElementOp exist (that is, there is more than
        one node), then if has_collect returns nonzero then collect is
        called with a MessageGrp reference after all of the blocks have
        been processed.  The default return value of has_collect is 0 and
        collect's default action is do nothing.  If defer_collect member is
        called with nonzero, collect will do nothing (this is only used by
        the blocked matrices). */
    virtual int has_collect();
    virtual void defer_collect(int);
    virtual void collect(const Ref<MessageGrp>&);
    /** Multithreaded use of cloneable SCElementOp objects requires that
        data from cloned objects be collected.  The default implementation
        will throw an exception. */
    virtual void collect(const Ref<SCElementOp>&);
    /** By default this returns nonzero.  If the ElementOp specialization
        will change any elements of the matrix, then this must be
        overridden to return nonzero. */
    virtual int has_side_effects();

    /** Returns true if this SCElementOp is threadsafe. The default
     * implementation returns false. */
    virtual bool threadsafe();

    /** Returns true if this SCElementOp supports the cloneable member. The
     * default implmentation returns false. */
    virtual bool cloneable();

    /** Returns a clone of this object.  This is needed for multithreaded
        use of SCElementOp objects that are not thread safe. The default
        implemenation throws an exception. */
    virtual Ref<SCElementOp> clone();

    /** This is the fallback routine to process blocks and is called
        by process_spec members that are not overridden. */
    virtual void process(SCMatrixBlockIter&) = 0;

    /** Lazy matrix implementors can call this member when the
        type of block specialization is unknown.  However, this
        will attempt to dynamic_cast block to a block specialization
        and will thus be less efficient. */
    void process_base(SCMatrixBlock*block);

    /** Matrices should call these members when the type of block is known.
        ElementOp specializations should override these when
        efficiency is important, since these give the most efficient access
        to the elements of the block. */
    virtual void process_spec_rect(SCMatrixRectBlock*);
    virtual void process_spec_ltri(SCMatrixLTriBlock*);
    virtual void process_spec_diag(SCMatrixDiagBlock*);
    virtual void process_spec_vsimp(SCVectorSimpleBlock*);
    virtual void process_spec_rectsub(SCMatrixRectSubBlock*);
    virtual void process_spec_ltrisub(SCMatrixLTriSubBlock*);
    virtual void process_spec_diagsub(SCMatrixDiagSubBlock*);
    virtual void process_spec_vsimpsub(SCVectorSimpleSubBlock*);
};

/** The SCElementOp2 class is very similar to the SCElementOp class except
    that pairs of blocks are treated simultaneously.  The two matrices
    involved must have identical storage layout, which will be the case if
    both matrices are of the same type and dimensions.  */
class SCElementOp2: public SavableState {
  public:
    SCElementOp2();
    SCElementOp2(StateIn&s): SavableState(s) {}
    virtual ~SCElementOp2();
    virtual int has_collect();
    virtual void defer_collect(int);
    virtual int has_side_effects();
    virtual int has_side_effects_in_arg();
    virtual void collect(const Ref<MessageGrp>&);
    virtual void process(SCMatrixBlockIter&,SCMatrixBlockIter&) = 0;
    void process_base(SCMatrixBlock*,SCMatrixBlock*);
    virtual void process_spec_rect(SCMatrixRectBlock*,SCMatrixRectBlock*);
    virtual void process_spec_ltri(SCMatrixLTriBlock*,SCMatrixLTriBlock*);
    virtual void process_spec_diag(SCMatrixDiagBlock*,SCMatrixDiagBlock*);
    virtual void process_spec_vsimp(SCVectorSimpleBlock*,SCVectorSimpleBlock*);
};

/** The SCElementOp3 class is very similar to the SCElementOp class except
    that a triplet of blocks is treated simultaneously.  The three matrices
    involved must have identical storage layout, which will be the case if
    all matrices are of the same type and dimensions.  */
class SCElementOp3: public SavableState {
  public:
    SCElementOp3();
    SCElementOp3(StateIn&s): SavableState(s) {}
    virtual ~SCElementOp3();
    virtual int has_collect();
    virtual void defer_collect(int);
    virtual int has_side_effects();
    virtual int has_side_effects_in_arg1();
    virtual int has_side_effects_in_arg2();
    virtual void collect(const Ref<MessageGrp>&);
    virtual void process(SCMatrixBlockIter&,
                         SCMatrixBlockIter&,
                         SCMatrixBlockIter&) = 0;
    void process_base(SCMatrixBlock*,SCMatrixBlock*,SCMatrixBlock*);
    virtual void process_spec_rect(SCMatrixRectBlock*,
                                   SCMatrixRectBlock*,
                                   SCMatrixRectBlock*);
    virtual void process_spec_ltri(SCMatrixLTriBlock*,
                                   SCMatrixLTriBlock*,
                                   SCMatrixLTriBlock*);
    virtual void process_spec_diag(SCMatrixDiagBlock*,
                                   SCMatrixDiagBlock*,
                                   SCMatrixDiagBlock*);
    virtual void process_spec_vsimp(SCVectorSimpleBlock*,
                                    SCVectorSimpleBlock*,
                                    SCVectorSimpleBlock*);
};

class SCElementScalarProduct: public SCElementOp2 {
  private:
    int deferred_;
    double product;
  public:
    SCElementScalarProduct();
    SCElementScalarProduct(StateIn&);
    ~SCElementScalarProduct();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&,SCMatrixBlockIter&);
    int has_collect();
    void defer_collect(int);
    void collect(const Ref<MessageGrp>&);
    double result();
    void init() { product = 0.0; }
};


class SCDestructiveElementProduct: public SCElementOp2 {
  public:
    SCDestructiveElementProduct();
    SCDestructiveElementProduct(StateIn&);
    ~SCDestructiveElementProduct();
    int has_side_effects();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&,SCMatrixBlockIter&);
};

class SCElementScale: public SCElementOp {
  private:
    double scale;
  public:
    SCElementScale(double a);
    SCElementScale(StateIn&);
    ~SCElementScale();
    int has_side_effects();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&);
};

class SCElementRandomize: public SCElementOp {
  private:
    double assign;
  public:
    SCElementRandomize();
    SCElementRandomize(StateIn&);
    ~SCElementRandomize();
    int has_side_effects();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&);
};

class SCElementAssign: public SCElementOp {
  private:
    double assign;
  public:
    SCElementAssign(double a);
    SCElementAssign(StateIn&);
    ~SCElementAssign();
    int has_side_effects();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&);
};

class SCElementSquareRoot: public SCElementOp {
  public:
    SCElementSquareRoot();
    SCElementSquareRoot(double a);
    SCElementSquareRoot(StateIn&);
    ~SCElementSquareRoot();
    int has_side_effects();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&);
};

class SCElementInvert: public SCElementOp {
  private:
    double threshold_;
    int nbelowthreshold_;
    int deferred_;
  public:
    SCElementInvert(double threshold = 0.0);
    SCElementInvert(StateIn&);
    ~SCElementInvert();
    int has_side_effects();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&);
    int has_collect();
    void defer_collect(int);
    void collect(const Ref<MessageGrp>&);
    void collect(const Ref<SCElementOp>&);
    int result() { return nbelowthreshold_; }
};


class SCElementScaleDiagonal: public SCElementOp {
  private:
    double scale_diagonal;
  public:
    SCElementScaleDiagonal(double a);
    SCElementScaleDiagonal(StateIn&);
    ~SCElementScaleDiagonal();
    int has_side_effects();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&);
};

class SCElementShiftDiagonal: public SCElementOp {
  private:
    double shift_diagonal;
  public:
    SCElementShiftDiagonal(double a);
    SCElementShiftDiagonal(StateIn&);
    ~SCElementShiftDiagonal();
    int has_side_effects();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&);
};

class SCElementMaxAbs: public SCElementOp {
  private:
    int deferred_;
    double r;
  public:
    SCElementMaxAbs();
    SCElementMaxAbs(StateIn&);
    ~SCElementMaxAbs();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&);
    int has_collect();
    void defer_collect(int);
    void collect(const Ref<MessageGrp>&);
    void collect(const Ref<SCElementOp>&);
    double result();
};


class SCElementMinAbs: public SCElementOp {
  private:
    int deferred_;
    double r;
  public:
    // rinit must be greater than the magnitude of the smallest element
    SCElementMinAbs(double rinit);
    SCElementMinAbs(StateIn&);
    ~SCElementMinAbs();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&);
    int has_collect();
    void defer_collect(int);
    void collect(const Ref<MessageGrp>&);
    void collect(const Ref<SCElementOp>&);
    double result();
};


class SCElementSum: public SCElementOp {
  private:
    int deferred_;
    double r;
  public:
    SCElementSum();
    SCElementSum(StateIn&);
    ~SCElementSum();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&);
    int has_collect();
    void defer_collect(int);
    void collect(const Ref<MessageGrp>&);
    void collect(const Ref<SCElementOp>&);
    double result();
    void init() { r = 0.0; }
};

/// Computes k-norm of matrix.
class SCElementKNorm: public SCElementOp {
  private:
    int deferred_;
    double r_;  // result
    unsigned int k_;  // norm parameter

    static double _process(SCMatrixBlockIter& i, int k);
    static double _process1(SCMatrixBlockIter& i);
    static double _process2(SCMatrixBlockIter& i);

  public:
    /// by default compute 2-norm
    SCElementKNorm(unsigned int k = 2);
    SCElementKNorm(StateIn&);
    ~SCElementKNorm();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&);
    int has_collect();
    void defer_collect(int);
    void collect(const Ref<MessageGrp>&);
    void collect(const Ref<SCElementOp>&);
    double result();
    void init() { r_ = 0.0; }
};

class SCElementDot: public SCElementOp {
  private:
    double** avects;
    double** bvects;
    int length;
  public:
    SCElementDot(StateIn&);
    void save_data_state(StateOut&);
    SCElementDot(double**a, double**b, int length);
    void process(SCMatrixBlockIter&);
    int has_side_effects();
};

class SCElementAccumulateSCMatrix: public SCElementOp {
  private:
    SCMatrix *m;
  public:
    SCElementAccumulateSCMatrix(SCMatrix *);
    int has_side_effects();
    void process(SCMatrixBlockIter&);
};

class SCElementAccumulateSymmSCMatrix: public SCElementOp {
  private:
    SymmSCMatrix *m;
  public:
    SCElementAccumulateSymmSCMatrix(SymmSCMatrix *);
    int has_side_effects();
    void process(SCMatrixBlockIter&);
};

class SCElementAccumulateDiagSCMatrix: public SCElementOp {
  private:
    DiagSCMatrix *m;
  public:
    SCElementAccumulateDiagSCMatrix(DiagSCMatrix *);
    int has_side_effects();
    void process(SCMatrixBlockIter&);
};

class SCElementAccumulateSCVector: public SCElementOp {
  private:
    SCVector *m;
  public:
    SCElementAccumulateSCVector(SCVector *);
    int has_side_effects();
    void process(SCMatrixBlockIter&);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
