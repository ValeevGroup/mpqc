//
// util.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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

#include <typeinfo>
#include <vector>

#include <math/scmat/util.h>
#include <util/group/thread.h>
#include <math/scmat/matrix.h>
#include <math/scmat/blocked.h>

namespace sc {

class BlockOpThread: public Thread {
    int me_;
    int n_;
    Ref<SCElementOp> op_;
    Ref<SCMatrixBlockList> blocklist_;
  public:
    BlockOpThread(int me,
                  int n,
                  const Ref<SCElementOp>& op,
                  const Ref<SCMatrixBlockList> &blocklist);
    void run();
};

}

using namespace sc;

BlockOpThread::BlockOpThread(int me,
                             int n,
                             const Ref<SCElementOp>& op,
                             const Ref<SCMatrixBlockList> &blocklist)
{
  me_ = me;
  n_ = n;
  op_ = op;
  blocklist_ = blocklist;
}

void
BlockOpThread::run()
{
  unsigned long count = 0;
  SCMatrixBlockListIter i;
  for (i = blocklist_->begin(); i != blocklist_->end(); i++,count++) {
      if (count%n_ == me_) {
          op_->process_base(i.block());
        }
    }
}

void
sc::scmat_perform_op_on_blocks(const Ref<SCElementOp>& op,
                               const Ref<SCMatrixBlockList> &blocklist)
{
  Ref<ThreadGrp> thr = ThreadGrp::get_default_threadgrp();

  for (int i=0; i<thr->nthread(); i++) {
      thr->add_thread(i,0);
    }

  Ref<SCElementOp> *ops = new Ref<SCElementOp>[thr->nthread()];

  int nthread;
  if (op->threadsafe()) {
      nthread = thr->nthread();
      for (int i=0; i<nthread; i++) ops[i] = op;
    }
  else if (op->cloneable()) {
      nthread = thr->nthread();
      ops[0] = op;
      for (int i=1; i<nthread; i++) ops[i] = op->clone();
    }
  else {
      ops[0] = op;
      nthread = 1;
    }

  for (int i=0; i<nthread; i++) {
      thr->add_thread(i, new BlockOpThread(i,nthread,ops[i],blocklist));
    }

  thr->start_threads();
  thr->wait_threads();
  thr->delete_threads();

  if (!op->threadsafe() && op->cloneable() && op->has_collect()) {
      for (int i = 1; i < nthread; i++) {
          ops[0]->collect(ops[i]);
        }
    }

  delete[] ops;
}

void
sc::canonicalize_column_phases(RefSCMatrix& A) {
  Ref<BlockedSCMatrix> A_blkd = dynamic_cast<BlockedSCMatrix*>(A.pointer());
  if (A_blkd.null()) { // if matrix is nonblocked, use SCElementOp
    RefDiagSCMatrix U = A.kit()->diagmatrix(A.coldim());
    const int ncol = A.ncol();

    // find max abs element in each column
    typedef SCElementFindExtremum<sc::abs_less<double>, SCMatrixIterationRanges::Columns> ElOp;
    Ref<ElOp> find_maxabs_op = new ElOp;
    A.element_op(find_maxabs_op);
    std::vector<SCElement> maxabs_elems = find_maxabs_op->result();

    // scale each column, if needed
    for(int c=0; c<ncol; ++c) {
      const double phase_correction = (maxabs_elems[c].value < 0.0) ? -1.0 : 1.0;  // max magnitude, not most positive!!!
      //ExEnv::out0() << indent << "col = " << c << ": max_value = " << maxabs_elems[c].value << " phase = " << phase_correction << std::endl;
      U.set_element(c, phase_correction);
    }
    A = A * U;
  }
  else { // with blocks don't be fancy, do serial but safe
    const double soft_zero = 1e-12;

    const int nrow = A.nrow();
    const int ncol = A.ncol();
    // in each column
    for(int c=0; c<ncol; ++c) {
      // find the largest-magnitude element
      int rmax = -1;
      int valmax = 0.0;
      for(int r=0; r<nrow; ++r) {
        const double value = fabs(A.get_element(r,c));
        if (value - valmax > soft_zero) {
          valmax = value;
          rmax = r;
        }
      }
      // if the largest-magnitude element is negative, scale the column by -1
      if (A.get_element(rmax,c) < -soft_zero) {
        for(int r=0; r<nrow; ++r) {
          A.set_element(r,c, A.get_element(r,c) * (-1.0));
        }
      }
    }

  }
}
