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

#include <math/scmat/util.h>
#include <util/group/thread.h>

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

