//
// blocked.h
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

#ifndef _math_scmat_blocked_h
#define _math_scmat_blocked_h

#include <math/scmat/block.h>
#include <math/scmat/elemop.h>
#include <math/scmat/matrix.h>
#include <math/scmat/abstract.h>

namespace sc {

class BlockedSCMatrixKit;
class BlockedSCVector;
class BlockedSCMatrix;
class BlockedSymmSCMatrix;
class BlockedDiagSCMatrix;

/// BlockedSCMatrixKit is a SCMatrixKit that produces blocked matrices
class BlockedSCMatrixKit: public SCMatrixKit {
  private:
    Ref<SCMatrixKit> subkit_;
  public:
    BlockedSCMatrixKit(const Ref<SCMatrixKit>& subkit);
    BlockedSCMatrixKit(const Ref<KeyVal>&);
    ~BlockedSCMatrixKit();
    SCMatrix* matrix(const RefSCDimension&,const RefSCDimension&);
    SymmSCMatrix* symmmatrix(const RefSCDimension&);
    DiagSCMatrix* diagmatrix(const RefSCDimension&);
    SCVector* vector(const RefSCDimension&);

    /// the kit used to implement blocks of the matrices
    Ref<SCMatrixKit> subkit() { return subkit_; }

    /// same as subkit, but recursive
    Ref<SCMatrixKit> subkit_lowest() {
      SCMatrixKit* subkit = subkit_.pointer();
      BlockedSCMatrixKit* bsubkit = dynamic_cast<BlockedSCMatrixKit*>(subkit);
      while (bsubkit) {
        subkit = bsubkit->subkit().pointer();
        bsubkit = dynamic_cast<BlockedSCMatrixKit*>(subkit);
      }
      return subkit;
    }

};


class BlockedSCVector: public SCVector {
    friend class BlockedSCMatrix;
    friend class BlockedSymmSCMatrix;
    friend class BlockedDiagSCMatrix;
  private:
    Ref<SCMatrixKit> subkit;
    RefSCVector *vecs_;

    void resize(SCDimension*);

  public:
    BlockedSCVector(const RefSCDimension&, BlockedSCMatrixKit*);
    ~BlockedSCVector();

    // Save and restore this in an implementation independent way.
    void save(StateOut&);
    void restore(StateIn&);

    void assign_val(double);
    void assign_v(SCVector*);
    void assign_p(const double*);

    double get_element(int) const;
    void set_element(int,double);
    void accumulate_element(int,double);

    void accumulate_product_rv(SCMatrix*,SCVector*);
    void accumulate_product_sv(SymmSCMatrix*,SCVector*);

    void accumulate(const SCVector*);
    void accumulate(const SCMatrix*);
    double scalar_product(SCVector*);

    void element_op(const Ref<SCElementOp>&);
    void element_op(const Ref<SCElementOp2>&,
                    SCVector*);
    void element_op(const Ref<SCElementOp3>&,
                    SCVector*,SCVector*);
    void vprint(const char* title=0,
                std::ostream& out=ExEnv::out0(), int =10) const;

    // BlockedSCVector specific functions
    RefSCDimension dim() const { return d; }
    RefSCDimension dim(int) const;
    int nblocks() const;
    RefSCVector block(int);

    Ref<SCMatrixSubblockIter> local_blocks(SCMatrixSubblockIter::Access);
    Ref<SCMatrixSubblockIter> all_blocks(SCMatrixSubblockIter::Access);
};

/// Blocked SCMatrix
class BlockedSCMatrix: public SCMatrix {
    friend class BlockedSymmSCMatrix;
    friend class BlockedDiagSCMatrix;
    friend class BlockedSCVector;
  private:
    Ref<SCMatrixKit> subkit;
    RefSCMatrix *mats_;
    int nblocks_;

    void resize(SCDimension*, SCDimension*);

  public:
    BlockedSCMatrix(const RefSCDimension&, const RefSCDimension&,
                    BlockedSCMatrixKit*);
    ~BlockedSCMatrix();

    // Save and restore this in an implementation independent way.
    void save(StateOut&);
    void restore(StateIn&);

    void assign_val(double);
    double get_element(int,int) const;
    void set_element(int,int,double);
    void accumulate_element(int,int,double);

    void assign_p(const double*);
    void assign_pp(const double**);
    void convert_p(double*) const;
    void convert_pp(double**) const;

    SCMatrix * get_subblock(int,int,int,int);
    void assign_subblock(SCMatrix*, int,int,int,int,int=0,int=0);
    void accumulate_subblock(SCMatrix*, int,int,int,int,int=0,int=0);

    SCVector * get_row(int i);
    SCVector * get_column(int i);
    void assign_row(SCVector *v, int i);
    void assign_column(SCVector *v, int i);
    void accumulate_row(SCVector *v, int i);
    void accumulate_column(SCVector *v, int i);

    void accumulate_outer_product(SCVector*,SCVector*);
    void accumulate_product_rr(SCMatrix*,SCMatrix*);
    void accumulate_product_rs(SCMatrix*,SymmSCMatrix*);
    void accumulate_product_rd(SCMatrix*,DiagSCMatrix*);
    void accumulate(const SCMatrix*);
    void accumulate(const SymmSCMatrix*);
    void accumulate(const DiagSCMatrix*);
    void accumulate(const SCVector*);

    void transpose_this();
    double invert_this();
    void svd_this(SCMatrix *U, DiagSCMatrix *sigma, SCMatrix *V);
    double solve_this(SCVector*);
    double determ_this();
    double trace();
    /// generalized-invert this. \sa SCMatrix::gen_invert_this()
    void gen_invert_this(double condition_number_threshold = 1e8);
    void schmidt_orthog(SymmSCMatrix*,int);
    int schmidt_orthog_tol(SymmSCMatrix*, double tol, double *res=0);

    void convert_accumulate(SCMatrix*a);

    void element_op(const Ref<SCElementOp>&);
    void element_op(const Ref<SCElementOp2>&,
                    SCMatrix*);
    void element_op(const Ref<SCElementOp3>&,
                    SCMatrix*,SCMatrix*);

    void vprint(const char* title=0,
                std::ostream& out=ExEnv::out0(), int =10) const;

    // BlockedSCMatrix specific functions
    RefSCDimension rowdim() const { return d1; }
    RefSCDimension coldim() const { return d2; }
    RefSCDimension rowdim(int) const;
    RefSCDimension coldim(int) const;
    int nblocks() const;
    RefSCMatrix block(int);

    Ref<SCMatrixSubblockIter> local_blocks(SCMatrixSubblockIter::Access);
    Ref<SCMatrixSubblockIter> all_blocks(SCMatrixSubblockIter::Access);
};

/// Blocked SymmSCMatrix
class BlockedSymmSCMatrix: public SymmSCMatrix {
    friend class BlockedSCMatrix;
    friend class BlockedDiagSCMatrix;
    friend class BlockedSCVector;
  private:
    Ref<SCMatrixKit> subkit;
    RefSymmSCMatrix *mats_;

    void resize(SCDimension*);

  public:
    BlockedSymmSCMatrix(const RefSCDimension&,BlockedSCMatrixKit*);
    ~BlockedSymmSCMatrix();

    // Save and restore this in an implementation independent way.
    void save(StateOut&);
    void restore(StateIn&);

    double get_element(int,int) const;
    void set_element(int,int,double);
    void accumulate_element(int,int,double);
    void scale(double);
    void assign_val(double);
    void assign_s(SymmSCMatrix*m);

    void assign_p(const double*);
    void assign_pp(const double**);
    void convert_p(double*) const;
    void convert_pp(double**) const;

    SCMatrix * get_subblock(int,int,int,int);
    SymmSCMatrix * get_subblock(int,int);
    void assign_subblock(SCMatrix*, int,int,int,int);
    void assign_subblock(SymmSCMatrix*, int,int);
    void accumulate_subblock(SCMatrix*, int,int,int,int);
    void accumulate_subblock(SymmSCMatrix*, int,int);
    SCVector * get_row(int i);
    void assign_row(SCVector *v, int i);
    void accumulate_row(SCVector *v, int i);

    double invert_this();
    double determ_this();
    double trace();
    double solve_this(SCVector*);
    /// generalized-invert this. \sa SymmSCMatrix::gen_invert_this()
    void gen_invert_this(double condition_number_threshold = 1e8);

    double scalar_product(SCVector*);
    void diagonalize(DiagSCMatrix*,SCMatrix*);
    /// like diagonalize(), but with general metric S
    void eigensystem(SymmSCMatrix* S, DiagSCMatrix* e, SCMatrix* V);

    void accumulate(const SymmSCMatrix*);
    void accumulate_symmetric_outer_product(SCVector*);
    void accumulate_symmetric_product(SCMatrix*);
    void accumulate_symmetric_sum(SCMatrix*);
    void accumulate_transform(SCMatrix*,SymmSCMatrix*,
                              SCMatrix::Transform = SCMatrix::NormalTransform);
    void accumulate_transform(SCMatrix*,DiagSCMatrix*,
                              SCMatrix::Transform = SCMatrix::NormalTransform);
    void accumulate_transform(SymmSCMatrix*,SymmSCMatrix*);

    void convert_accumulate(SymmSCMatrix*a);

    void element_op(const Ref<SCElementOp>&);
    void element_op(const Ref<SCElementOp2>&,
                    SymmSCMatrix*);
    void element_op(const Ref<SCElementOp3>&,
                    SymmSCMatrix*,SymmSCMatrix*);

    void vprint(const char* title=0,
                std::ostream& out=ExEnv::out0(), int =10) const;

    // BlockedSymmSCMatrix specific functions
    RefSCDimension dim() const { return d; }
    RefSCDimension dim(int) const;
    int nblocks() const;
    RefSymmSCMatrix block(int);

    Ref<SCMatrixSubblockIter> local_blocks(SCMatrixSubblockIter::Access);
    Ref<SCMatrixSubblockIter> all_blocks(SCMatrixSubblockIter::Access);
};

/// Blocked DiagSCMatrix
class BlockedDiagSCMatrix: public DiagSCMatrix {
    friend class BlockedSCMatrix;
    friend class BlockedSymmSCMatrix;
    friend class BlockedSCVector;
  private:
    Ref<SCMatrixKit> subkit;
    RefDiagSCMatrix *mats_;

    void resize(SCDimension*);

  public:
    BlockedDiagSCMatrix(const RefSCDimension&,BlockedSCMatrixKit*);
    ~BlockedDiagSCMatrix();

    // Save and restore this in an implementation independent way.
    void save(StateOut&);
    void restore(StateIn&);

    double get_element(int) const;
    void set_element(int,double);
    void accumulate_element(int,double);
    void accumulate(const DiagSCMatrix*);

    double invert_this();
    double determ_this();
    double trace();
    /// generalized-invert this. \sa DiagSCMatrix::gen_invert_this()
    void gen_invert_this(double condition_number_threshold = 1e8);

    void convert_accumulate(DiagSCMatrix*a);

    void element_op(const Ref<SCElementOp>&);
    void element_op(const Ref<SCElementOp2>&,
                    DiagSCMatrix*);
    void element_op(const Ref<SCElementOp3>&,
                    DiagSCMatrix*,DiagSCMatrix*);
    void vprint(const char* title=0,
                std::ostream& out=ExEnv::out0(), int =10) const;

    // BlockedDiagSCMatrix specific functions
    RefSCDimension dim() const { return d; }
    RefSCDimension dim(int) const;
    int nblocks() const;
    RefDiagSCMatrix block(int);

    Ref<SCMatrixSubblockIter> local_blocks(SCMatrixSubblockIter::Access);
    Ref<SCMatrixSubblockIter> all_blocks(SCMatrixSubblockIter::Access);
};

class BlockedSCElementOp : public SCElementOp {
  private:
    int current_block_;

  public:
    BlockedSCElementOp();
    void working_on(int);
    int current_block() const;
};

class BlockedSCElementOp2 : public SCElementOp2 {
  private:
    int current_block_;

  public:
    BlockedSCElementOp2();
    void working_on(int);
    int current_block() const;
};

class BlockedSCElementOp3 : public SCElementOp3 {
  private:
    int current_block_;

  public:
    BlockedSCElementOp3();
    void working_on(int);
    int current_block() const;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
