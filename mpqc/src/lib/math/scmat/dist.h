//
// dist.h
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

#ifdef __GNUC__
#pragma interface
#endif

#ifndef _math_scmat_dist_h
#define _math_scmat_dist_h

#include <util/group/message.h>
#include <util/group/mstate.h>

#include <math/scmat/block.h>
#include <math/scmat/matrix.h>
#include <math/scmat/abstract.h>

class DistSCMatrixKit: public SCMatrixKit {
#   define CLASSNAME DistSCMatrixKit
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  public:
    DistSCMatrixKit(const RefMessageGrp &grp = 0);
    DistSCMatrixKit(const RefKeyVal&);
    ~DistSCMatrixKit();
    SCMatrix* matrix(const RefSCDimension&,const RefSCDimension&);
    SymmSCMatrix* symmmatrix(const RefSCDimension&);
    DiagSCMatrix* diagmatrix(const RefSCDimension&);
    SCVector* vector(const RefSCDimension&);
};
DescribedClass_REF_dec(DistSCMatrixKit);

class DistSCVector: public SCVector {
    friend class DistSCMatrix;
    friend class DistSymmSCMatrix;
    friend class DistDiagSCMatrix;
#   define CLASSNAME DistSCVector
#   include <util/class/classd.h>
  protected:
    RefSCMatrixBlockList blocklist;

    void init_blocklist();
    double *find_element(int i) const;
    int element_to_node(int i) const;
    int block_to_node(int) const;
    RefSCMatrixBlock block_to_block(int) const;
    void error(const char *);
  public:
    DistSCVector(const RefSCDimension&, DistSCMatrixKit*);
    ~DistSCVector();
    void assign_p(const double*);
    void assign_v(SCVector*a);
    void convert(double* v) const;
    void convert(SCVector *);

    void set_element(int,double);
    void accumulate_element(int,double);
    double get_element(int) const;
    void accumulate(const SCVector*);
    void accumulate(const SCMatrix*m);
    double scalar_product(SCVector*);
    void accumulate_product_rv(SCMatrix *, SCVector *);
    void element_op(const RefSCElementOp&);
    void element_op(const RefSCElementOp2&,
                    SCVector*);
    void element_op(const RefSCElementOp3&,
                    SCVector*,SCVector*);
    void vprint(const char* title=0,ostream& out=cout, int =10) const;

    RefSCMatrixSubblockIter local_blocks(SCMatrixSubblockIter::Access);
    RefSCMatrixSubblockIter all_blocks(SCMatrixSubblockIter::Access);

    RefDistSCMatrixKit skit();
};

class DistSCMatrix: public SCMatrix {
    friend class DistSymmSCMatrix;
    friend class DistDiagSCMatrix;
    friend DistSCVector;
#   define CLASSNAME DistSCMatrix
#   include <util/class/classd.h>
  protected:
    RefSCMatrixBlockList blocklist;

    int vecoff;
    int nvec;
    double **vec;
  protected:
    // utility functions
    void init_blocklist();
    void error(const char *);
    double *find_element(int i, int j) const;
    int element_to_node(int i, int j) const;
    int block_to_node(int,int) const;
    RefSCMatrixBlock block_to_block(int, int) const;
    RefSCBlockInfo rowblocks() const { return d1->blocks(); }
    RefSCBlockInfo colblocks() const { return d2->blocks(); }
    
    enum VecOp {CopyFromVec, CopyToVec, AccumFromVec, AccumToVec};
    enum Form { Row, Col } form;
    void create_vecform(Form, int nvec = -1);
    void delete_vecform();
    void vecform_op(VecOp op, int *ivec = 0);
    void vecform_zero();
  public:
    DistSCMatrix(const RefSCDimension&, const RefSCDimension&,
                 DistSCMatrixKit*);
    ~DistSCMatrix();

    // implementations and overrides of virtual functions
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
    void accumulate(const SCMatrix*);
    void accumulate(const SymmSCMatrix*);
    void accumulate(const DiagSCMatrix*);
    void accumulate(const SCVector*);
    void transpose_this();
    double invert_this();
    double solve_this(SCVector*);
    double determ_this();
    double trace();
    void gen_invert_this();
    void schmidt_orthog(SymmSCMatrix*,int);
    void element_op(const RefSCElementOp&);
    void element_op(const RefSCElementOp2&,
                    SCMatrix*);
    void element_op(const RefSCElementOp3&,
                    SCMatrix*,SCMatrix*);
    void vprint(const char* title=0,ostream& out=cout, int =10);
    void vprint(const char* title=0,ostream& out=cout, int =10) const;

    RefSCMatrixSubblockIter local_blocks(SCMatrixSubblockIter::Access);
    RefSCMatrixSubblockIter all_blocks(SCMatrixSubblockIter::Access);

    RefDistSCMatrixKit skit();
};

class DistSymmSCMatrix: public SymmSCMatrix {
    friend class DistSCMatrix;
    friend class DistDiagSCMatrix;
    friend DistSCVector;
#   define CLASSNAME DistSymmSCMatrix
#   include <util/class/classd.h>
  protected:
    RefSCMatrixBlockList blocklist;
  protected:
    // utility functions
    void init_blocklist();
    double *find_element(int i, int j) const;
    int element_to_node(int i, int j) const;
    int block_to_node(int,int) const;
    RefSCMatrixBlock block_to_block(int, int) const;

    void error(const char *msg);
  public:
    DistSymmSCMatrix(const RefSCDimension&, DistSCMatrixKit*);
    ~DistSymmSCMatrix();

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
    void gen_invert_this();

    void diagonalize(DiagSCMatrix*,SCMatrix*);
    void accumulate_symmetric_sum(SCMatrix*);
    void element_op(const RefSCElementOp&);
    void element_op(const RefSCElementOp2&,
                    SymmSCMatrix*);
    void element_op(const RefSCElementOp3&,
                    SymmSCMatrix*,SymmSCMatrix*);

    virtual void convert_accumulate(SymmSCMatrix*);

    RefSCMatrixSubblockIter local_blocks(SCMatrixSubblockIter::Access);
    RefSCMatrixSubblockIter all_blocks(SCMatrixSubblockIter::Access);

    RefDistSCMatrixKit skit();
};

class DistDiagSCMatrix: public DiagSCMatrix {
    friend DistSCMatrix;
    friend DistSymmSCMatrix;
    friend DistSCVector;
#   define CLASSNAME DistDiagSCMatrix
#   include <util/class/classd.h>
  protected:
    RefSCMatrixBlockList blocklist;

    void init_blocklist();
    double *find_element(int i) const;
    int element_to_node(int i) const;
    int block_to_node(int) const;
    RefSCMatrixBlock block_to_block(int) const;
    void error(const char *msg);
  public:
    DistDiagSCMatrix(const RefSCDimension&, DistSCMatrixKit*);
    ~DistDiagSCMatrix();

    // implementations and overrides of virtual functions
    double get_element(int) const;
    void set_element(int,double);
    void accumulate_element(int,double);
    void accumulate(const DiagSCMatrix*);
    double invert_this();
    double determ_this();
    double trace();
    void gen_invert_this();

    void element_op(const RefSCElementOp&);
    void element_op(const RefSCElementOp2&,
                    DiagSCMatrix*);
    void element_op(const RefSCElementOp3&,
                    DiagSCMatrix*,DiagSCMatrix*);

    RefSCMatrixSubblockIter local_blocks(SCMatrixSubblockIter::Access);
    RefSCMatrixSubblockIter all_blocks(SCMatrixSubblockIter::Access);

    RefDistSCMatrixKit skit();
};

class DistSCMatrixListSubblockIter: public SCMatrixListSubblockIter {
  protected:
    RefMessageGrp grp_;
    StateSend out_;
    StateRecv in_;
    int step_;
    RefSCMatrixBlockList locallist_;

    void maybe_advance_list();
    void advance_list();
  public:
    DistSCMatrixListSubblockIter(Access,
                                 const RefSCMatrixBlockList &locallist,
                                 const RefMessageGrp &grp);
    void begin();
    void next();
    ~DistSCMatrixListSubblockIter();
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
