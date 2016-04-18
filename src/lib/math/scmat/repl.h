//
// repl.h
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

#ifndef _math_scmat_repl_h
#define _math_scmat_repl_h

#include <util/group/message.h>

#include <math/scmat/block.h>
#include <math/scmat/matrix.h>
#include <math/scmat/abstract.h>

namespace sc {

/** The ReplSCMatrixKit produces matrices that work in a many processor
environment.  A copy of the entire matrix is stored on each node.
*/
class ReplSCMatrixKit: public SCMatrixKit {
  public:
    ReplSCMatrixKit();
    ReplSCMatrixKit(const Ref<KeyVal>&);
    ~ReplSCMatrixKit();
    SCMatrix* matrix(const RefSCDimension&,const RefSCDimension&);
    SymmSCMatrix* symmmatrix(const RefSCDimension&);
    DiagSCMatrix* diagmatrix(const RefSCDimension&);
    SCVector* vector(const RefSCDimension&);
};


class ReplSCMatrixListSubblockIter: public SCMatrixListSubblockIter {
  protected:
    Ref<MessageGrp> grp_;
    double *data_;
    int ndata_;
  public:
    ReplSCMatrixListSubblockIter(Access,
                             const Ref<SCMatrixBlockList> &list,
                             const Ref<MessageGrp> &grp,
                             double *data, int ndata);
    ~ReplSCMatrixListSubblockIter();
};

class ReplSCVector: public SCVector {
    friend class ReplSCMatrix;
    friend class ReplSymmSCMatrix;
    friend class ReplDiagSCMatrix;
  protected:
    Ref<SCMatrixBlockList> blocklist;
    double* vector;
    void init_blocklist();
    void before_elemop();
    void after_elemop();
  public:
    ReplSCVector(const RefSCDimension&,ReplSCMatrixKit*);
    ~ReplSCVector();
    void assign_val(double);
    void assign_v(SCVector*);
    void assign_p(const double*);

    void set_element(int,double);
    void accumulate_element(int,double);
    double get_element(int) const;
    void accumulate_product_sv(SymmSCMatrix*,SCVector*);
    void accumulate_product_rv(SCMatrix*,SCVector*);
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

    // return a pointer to the data for fast access
    double *get_data() { return vector; }

    Ref<SCMatrixSubblockIter> local_blocks(SCMatrixSubblockIter::Access);
    Ref<SCMatrixSubblockIter> all_blocks(SCMatrixSubblockIter::Access);

    Ref<ReplSCMatrixKit> skit();
};

class ReplSCMatrix: public SCMatrix {
    friend class ReplSymmSCMatrix;
    friend class ReplDiagSCMatrix;
    friend class ReplSCVector;
  protected:
    Ref<SCMatrixBlockList> blocklist;
    double* matrix;
    double** rows;
  protected:
    // utility functions
    size_t compute_offset(int,int) const;
    void init_blocklist();

    void before_elemop();
    void after_elemop();
  public:
    ReplSCMatrix(const RefSCDimension&,const RefSCDimension&,
                 ReplSCMatrixKit*);
    ~ReplSCMatrix();

    // implementations and overrides of virtual functions
    void assign_val(double);
    double get_element(int,int) const;
    void set_element(int,int,double);
    void accumulate_element(int,int,double);
    SCMatrix * get_subblock(int,int,int,int);
    void assign_subblock(SCMatrix*, int,int,int,int,int=0,int=0);
    void accumulate_subblock(SCMatrix*, int,int,int,int,int=0,int=0);
    SCVector * get_row(int i);
    SCVector * get_column(int i);
    void assign_row(SCVector *v, int i);
    void assign_column(SCVector *v, int i);
    void accumulate_row(SCVector *v, int i);
    void accumulate_column(SCVector *v, int i);
    void assign_p(const double*);
    void assign_pp(const double**);
    void convert_p(double*) const;
    void convert_pp(double**) const;

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
    void schmidt_orthog(SymmSCMatrix*,int);
    int schmidt_orthog_tol(SymmSCMatrix*, double tol, double *res=0);
    void element_op(const Ref<SCElementOp>&);
    void element_op(const Ref<SCElementOp2>&,
                    SCMatrix*);
    void element_op(const Ref<SCElementOp3>&,
                    SCMatrix*,SCMatrix*);
    void vprint(const char* title=0,
                std::ostream& out=ExEnv::out0(), int =10) const;

    // return a pointer to the data for fast access
    double *get_data() { return matrix; }
    double **get_rows() { return rows; }

    Ref<SCMatrixSubblockIter> local_blocks(SCMatrixSubblockIter::Access);
    Ref<SCMatrixSubblockIter> all_blocks(SCMatrixSubblockIter::Access);

    Ref<ReplSCMatrixKit> skit();
};

/// Replicated SymmSCMatrix
class ReplSymmSCMatrix: public SymmSCMatrix {
    friend class ReplSCMatrix;
    friend class ReplDiagSCMatrix;
    friend class ReplSCVector;
  protected:
    Ref<SCMatrixBlockList> blocklist;
    double* matrix;
    double** rows;
  protected:
    // utility functions
    size_t compute_offset(int,int) const;
    void init_blocklist();

    void before_elemop();
    void after_elemop();
  public:
    ReplSymmSCMatrix(const RefSCDimension&, ReplSCMatrixKit*);
    ~ReplSymmSCMatrix();

    // implementations and overrides of virtual functions
    void assign_val(double);
    void assign_s(SymmSCMatrix*);
    void assign_p(const double*);
    void assign_pp(const double**);
    void convert_p(double*) const;
    void convert_pp(double**) const;
    double get_element(int,int) const;
    void set_element(int,int,double);
    void accumulate_element(int,int,double);
    void scale(double);

    SCMatrix * get_subblock(int,int,int,int);
    SymmSCMatrix * get_subblock(int,int);
    void assign_subblock(SCMatrix*, int,int,int,int);
    void assign_subblock(SymmSCMatrix*, int,int);
    void accumulate_subblock(SCMatrix*, int,int,int,int);
    void accumulate_subblock(SymmSCMatrix*, int,int);
    SCVector * get_row(int i);
    void assign_row(SCVector *v, int i);
    void accumulate_row(SCVector *v, int i);

    void accumulate_product_rr(SCMatrix*,SCMatrix*);
    void accumulate(const SymmSCMatrix*);
    double invert_this();
    double solve_this(SCVector*);
    double trace();
    double determ_this();
    void gen_invert_this(double condition_number_threshold = 1e8);

    double scalar_product(SCVector*);
    void diagonalize(DiagSCMatrix*,SCMatrix*);
    void eigensystem(SymmSCMatrix*,DiagSCMatrix*,SCMatrix*);
    void accumulate_symmetric_outer_product(SCVector*);
    void accumulate_symmetric_product(SCMatrix*);
    void accumulate_symmetric_sum(SCMatrix*);
    void accumulate_transform(SCMatrix*,SymmSCMatrix*,
                              SCMatrix::Transform = SCMatrix::NormalTransform);
    void accumulate_transform(SCMatrix*,DiagSCMatrix*,
                              SCMatrix::Transform = SCMatrix::NormalTransform);
    void accumulate_transform(SymmSCMatrix*,SymmSCMatrix*);
    void element_op(const Ref<SCElementOp>&);
    void element_op(const Ref<SCElementOp2>&,
                    SymmSCMatrix*);
    void element_op(const Ref<SCElementOp3>&,
                    SymmSCMatrix*,SymmSCMatrix*);
    void vprint(const char* title=0,
                std::ostream& out=ExEnv::out0(), int =10) const;

    // return a pointer to the data for fast access
    double *get_data() { return matrix; }
    double **get_rows() { return rows; }

    Ref<SCMatrixSubblockIter> local_blocks(SCMatrixSubblockIter::Access);
    Ref<SCMatrixSubblockIter> all_blocks(SCMatrixSubblockIter::Access);

    Ref<ReplSCMatrixKit> skit();
};

/// Replicated DiagSCMatrix
class ReplDiagSCMatrix: public DiagSCMatrix {
    friend class ReplSCMatrix;
    friend class ReplSymmSCMatrix;
    friend class ReplSCVector;
  protected:
    Ref<SCMatrixBlockList> blocklist;
    void init_blocklist();
    double* matrix;

    void before_elemop();
    void after_elemop();
  public:
    ReplDiagSCMatrix(const RefSCDimension&, ReplSCMatrixKit*);
    ~ReplDiagSCMatrix();

    // implementations and overrides of virtual functions
    void assign_val(double);
    double get_element(int) const;
    void set_element(int,double);
    void accumulate_element(int,double);
    void accumulate(const DiagSCMatrix*);
    double invert_this();
    double determ_this();
    double trace();
    void gen_invert_this(double condition_number_threshold = 1e8);

    void element_op(const Ref<SCElementOp>&);
    void element_op(const Ref<SCElementOp2>&,
                    DiagSCMatrix*);
    void element_op(const Ref<SCElementOp3>&,
                    DiagSCMatrix*,DiagSCMatrix*);
    void vprint(const char* title=0,
                std::ostream& out=ExEnv::out0(), int =10) const;

    // return a pointer to the data for fast access
    double *get_data() { return matrix; }

    Ref<SCMatrixSubblockIter> local_blocks(SCMatrixSubblockIter::Access);
    Ref<SCMatrixSubblockIter> all_blocks(SCMatrixSubblockIter::Access);

    Ref<ReplSCMatrixKit> skit();
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
