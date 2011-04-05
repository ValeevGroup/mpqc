//
// local.h
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

#ifndef _math_scmat_local_h
#define _math_scmat_local_h

#include <math/scmat/block.h>
#include <math/scmat/matrix.h>
#include <math/scmat/abstract.h>

namespace sc {

class LocalSCMatrixKit;
class LocalSCVector;
class LocalSCMatrix;
class LocalSymmSCMatrix;
class LocalDiagSCMatrix;

/** The LocalSCMatrixKit produces matrices that work in a single processor
environment.  */
class LocalSCMatrixKit: public SCMatrixKit {
  public:
    LocalSCMatrixKit();
    LocalSCMatrixKit(const Ref<KeyVal>&);
    ~LocalSCMatrixKit();
    SCMatrix* matrix(const RefSCDimension&,const RefSCDimension&);
    SymmSCMatrix* symmmatrix(const RefSCDimension&);
    DiagSCMatrix* diagmatrix(const RefSCDimension&);
    SCVector* vector(const RefSCDimension&);
};

class LocalSCVector: public SCVector {
    friend class LocalSCMatrix;
    friend class LocalSymmSCMatrix;
    friend class LocalDiagSCMatrix;
  private:
    Ref<SCVectorSimpleBlock> block;

    void resize(int);
  public:
    LocalSCVector();
    LocalSCVector(const RefSCDimension&,LocalSCMatrixKit*);
    ~LocalSCVector();
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
    double *get_data();

    Ref<SCMatrixSubblockIter> local_blocks(SCMatrixSubblockIter::Access);
    Ref<SCMatrixSubblockIter> all_blocks(SCMatrixSubblockIter::Access);
};

class LocalSCMatrix: public SCMatrix {
    friend class LocalSymmSCMatrix;
    friend class LocalDiagSCMatrix;
    friend class LocalSCVector;
  private:
    Ref<SCMatrixRectBlock> block;
    double** rows;
  private:
    // utility functions
    int compute_offset(int,int) const;
    void resize(int,int);
  public:
    LocalSCMatrix(const RefSCDimension&,const RefSCDimension&,
                  LocalSCMatrixKit*);
    ~LocalSCMatrix();

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
    double *get_data();
    double **get_rows();

    Ref<SCMatrixSubblockIter> local_blocks(SCMatrixSubblockIter::Access);
    Ref<SCMatrixSubblockIter> all_blocks(SCMatrixSubblockIter::Access);
};

/// Local SymmSCMatrix
class LocalSymmSCMatrix: public SymmSCMatrix {
    friend class LocalSCMatrix;
    friend class LocalDiagSCMatrix;
    friend class LocalSCVector;
  private:
    Ref<SCMatrixLTriBlock> block;
    double** rows;
  private:
    // utility functions
    int compute_offset(int,int) const;
    void resize(int n);
  public:
    LocalSymmSCMatrix(const RefSCDimension&, LocalSCMatrixKit*);
    ~LocalSymmSCMatrix();

    // implementations and overrides of virtual functions
    double get_element(int,int) const;
    void set_element(int,int,double);
    void accumulate_element(int,int,double);

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
    double *get_data();
    double **get_rows();

    Ref<SCMatrixSubblockIter> local_blocks(SCMatrixSubblockIter::Access);
    Ref<SCMatrixSubblockIter> all_blocks(SCMatrixSubblockIter::Access);
};

/// Local DiagSCMatrix
class LocalDiagSCMatrix: public DiagSCMatrix {
    friend class LocalSCMatrix;
    friend class LocalSymmSCMatrix;
    friend class LocalSCVector;
  private:
    Ref<SCMatrixDiagBlock> block;
    void resize(int n);
  public:
    LocalDiagSCMatrix(const RefSCDimension&, LocalSCMatrixKit*);
    ~LocalDiagSCMatrix();

    // implementations and overrides of virtual functions
    void save_data_state(StateOut&);
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
    double *get_data();

    Ref<SCMatrixSubblockIter> local_blocks(SCMatrixSubblockIter::Access);
    Ref<SCMatrixSubblockIter> all_blocks(SCMatrixSubblockIter::Access);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
