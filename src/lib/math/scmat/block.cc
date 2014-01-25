//
// block.cc
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

#include <iostream>
#include <string.h>
#include <util/state/stateio.h>
#include <util/misc/consumableresources.h>
#include <math/scmat/block.h>
#include <math/scmat/blkiter.h>
#include <math/scmat/elemop.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////////////
// SCMatrixBlock member functions

static ClassDesc SCMatrixBlock_cd(
  typeid(SCMatrixBlock),"SCMatrixBlock",1,"public SavableState",
  0, 0, 0);

SCMatrixBlock::SCMatrixBlock()
{
  blocki = blockj = -1;
}

SCMatrixBlock::SCMatrixBlock(StateIn&s):
  SavableState(s)
{
  s.get(blocki);
  s.get(blockj);
}

void
SCMatrixBlock::save_data_state(StateOut&s)
{
  s.put(blocki);
  s.put(blockj);
}

SCMatrixBlock::~SCMatrixBlock()
{
}

SCMatrixBlock *
SCMatrixBlock::deepcopy() const
{
  ExEnv::errn() << "SCMatrixBlock of type " << class_name()
       << " cannot be deep copied" << endl;
  abort();
  return 0;
}

double *
SCMatrixBlock::dat()
{
  ExEnv::errn() << "SCMatrixBlock of type " << class_name()
       << " cannot provide internal data" << endl;
  abort();
  return 0;
}

int
SCMatrixBlock::ndat() const
{
  ExEnv::errn() << "SCMatrixBlock of type " << class_name()
       << " cannot provide size of internal data" << endl;
  abort();
  return 0;
}

/////////////////////////////////////////////////////////////////////////////
// SCMatrixBlockListLink member functions

SCMatrixBlockListLink::SCMatrixBlockListLink(SCMatrixBlock* b,
                                             SCMatrixBlockListLink* l)
{
  block(b);
  next(l);
}

SCMatrixBlockListLink::~SCMatrixBlockListLink()
{
  if (_block) _block->dereference();
  if (_block->nreference() == 0) delete _block;

  for (SCMatrixBlockListLink *nexti, *i=_next; i; i = nexti) {
      nexti = i->_next;
      i->_next = 0;
      delete i;
    }
}

void
SCMatrixBlockListLink::block(SCMatrixBlock* b)
{
  _block = b;
  if (_block) _block->reference();
}

/////////////////////////////////////////////////////////////////////////////
// SCMatrixBlockList member functions

static ClassDesc SCMatrixBlockList_cd(
  typeid(SCMatrixBlockList),"SCMatrixBlockList",1,"public SavableState",
  0, 0, create<SCMatrixBlockList>);

SCMatrixBlockList::SCMatrixBlockList()
{
  _begin = 0;
}

SCMatrixBlockList::SCMatrixBlockList(StateIn& s):
  SavableState(s)
{
  int i, count;
  Ref<SCMatrixBlock> b;
  s.get(count);
  _begin = 0;
  for (i=0; i<count; i++) {
      b << SavableState::restore_state(s);
      append(b);
    }
}

SCMatrixBlockList::~SCMatrixBlockList()
{
  if (_begin) delete _begin;
}

void
SCMatrixBlockList::save_data_state(StateOut&s)
{
  int count = 0;
  SCMatrixBlockListIter i;
  for (i = begin(); i != end(); i++) count++;
  s.put(count);
  for (i = begin(); i != end(); i++) {
      i.block()->save_state(s);
    }
}

void
SCMatrixBlockList::insert(SCMatrixBlock* b)
{
  _begin = new SCMatrixBlockListLink(b, _begin);
}

void
SCMatrixBlockList::append(SCMatrixBlock* b)
{
  if (_begin == 0) {
      _begin = new SCMatrixBlockListLink(b);
    }
  else {
      SCMatrixBlockListLink* i;
      for (i = _begin; i->next() != 0; i = i->next());
      i->next(new SCMatrixBlockListLink(b));
    }
}

SCMatrixBlockList *
SCMatrixBlockList::deepcopy()
{
  SCMatrixBlockListIter i;
  SCMatrixBlockList *ret = new SCMatrixBlockList();
  for (i=begin(); i!=end(); i++) {
      ret->append(i.block()->deepcopy());
    }
  return ret;
}

/////////////////////////////////////////////////////////////////////////////
// SCMatrixRectBlock member functions

static ClassDesc SCMatrixRectBlock_cd(
  typeid(SCMatrixRectBlock),"SCMatrixRectBlock",1,"public SCMatrixBlock",
  0, 0, create<SCMatrixRectBlock>);

SCMatrixRectBlock::SCMatrixRectBlock(int is, int ie, int js, int je):
  istart(is),
  jstart(js),
  iend(ie),
  jend(je)
{
  data = allocate<double>((ie-is)*(je-js));
}

SCMatrixRectBlock::SCMatrixRectBlock(StateIn&s):
  SCMatrixBlock(s)
{
  s.get(istart);
  s.get(jstart);
  s.get(iend);
  s.get(jend);
  s.get(data);
}

void
SCMatrixRectBlock::save_data_state(StateOut&s)
{
  SCMatrixBlock::save_data_state(s);
  s.put(istart);
  s.put(jstart);
  s.put(iend);
  s.put(jend);
  s.put(data,(iend-istart)*(jend-jstart));
}

SCMatrixBlock *
SCMatrixRectBlock::deepcopy() const
{
  SCMatrixRectBlock *ret = new SCMatrixRectBlock(istart,iend,jstart,jend);
  ret->blocki = blocki;
  ret->blockj = blockj;
  memcpy(ret->data, data, sizeof(double)*ndat());
  return ret;
}

double *
SCMatrixRectBlock::dat()
{
  return data;
}

int
SCMatrixRectBlock::ndat() const
{
  return (iend-istart)*(jend-jstart);
}

SCMatrixRectBlock::~SCMatrixRectBlock()
{
  deallocate(data);
}

void
SCMatrixRectBlock::process(SCElementOp*op)
{
  SCMatrixRectBlockIter i(this);
  op->process(i);
}

void
SCMatrixRectBlock::process(SCElementOp2*op,
                           SCMatrixBlock* b)
{
  SCMatrixRectBlockIter i(this);
  SCMatrixRectBlockIter j((SCMatrixRectBlock*)b);
  op->process(i,j);
}

void
SCMatrixRectBlock::process(SCElementOp3*op,
                           SCMatrixBlock* b1, SCMatrixBlock* b2)
{
  SCMatrixRectBlockIter i(this);
  SCMatrixRectBlockIter j((SCMatrixRectBlock*)b1);
  SCMatrixRectBlockIter k((SCMatrixRectBlock*)b2);
  op->process(i,j,k);
}

/////////////////////////////////////////////////////////////////////////////
// SCMatrixRectSubBlock member functions

static ClassDesc SCMatrixRectSubBlock_cd(
  typeid(SCMatrixRectSubBlock),"SCMatrixRectSubBlock",1,"public SCMatrixBlock",
  0, 0, create<SCMatrixRectSubBlock>);

SCMatrixRectSubBlock::SCMatrixRectSubBlock(int is, int ie, int istr,
                                           int js, int je, double* d):
  istart(is),
  jstart(js),
  iend(ie),
  jend(je),
  istride(istr),
  data(d)
{
}

SCMatrixRectSubBlock::SCMatrixRectSubBlock(StateIn&s):
  SCMatrixBlock(s)
{
  s.get(istart);
  s.get(istride);
  s.get(jstart);
  s.get(iend);
  s.get(jend);
  data = 0;
}

void
SCMatrixRectSubBlock::save_data_state(StateOut&s)
{
  SCMatrixBlock::save_data_state(s);
  s.put(istart);
  s.put(istride);
  s.put(jstart);
  s.put(iend);
  s.put(jend);
}

SCMatrixRectSubBlock::~SCMatrixRectSubBlock()
{
}

void
SCMatrixRectSubBlock::process(SCElementOp*op)
{
  SCMatrixRectSubBlockIter i(this);
  op->process(i);
}

void
SCMatrixRectSubBlock::process(SCElementOp2*op,
                              SCMatrixBlock* b)
{
  SCMatrixRectSubBlockIter i(this);
  SCMatrixRectSubBlockIter j((SCMatrixRectSubBlock*)b);
  op->process(i,j);
}

void
SCMatrixRectSubBlock::process(SCElementOp3*op,
                              SCMatrixBlock* b1,
                              SCMatrixBlock* b2)
{
  SCMatrixRectSubBlockIter i(this);
  SCMatrixRectSubBlockIter j((SCMatrixRectSubBlock*)b1);
  SCMatrixRectSubBlockIter k((SCMatrixRectSubBlock*)b2);
  op->process(i,j,k);
}

/////////////////////////////////////////////////////////////////////////////
// SCMatrixLTriBlock member functions

static ClassDesc SCMatrixLTriBlock_cd(
  typeid(SCMatrixLTriBlock),"SCMatrixLTriBlock",1,"public SCMatrixBlock",
  0, 0, create<SCMatrixLTriBlock>);

SCMatrixLTriBlock::SCMatrixLTriBlock(int s,int e):
  start(s),
  end(e)
{
  data = allocate<double>(((e-s)*(e-s+1))/2);
}

SCMatrixLTriBlock::SCMatrixLTriBlock(StateIn&s):
  SCMatrixBlock(s)
{
  s.get(start);
  s.get(end);
  s.get(data);
}

void
SCMatrixLTriBlock::save_data_state(StateOut&s)
{
  SCMatrixBlock::save_data_state(s);
  s.put(start);
  s.put(end);
  s.put(data,((end-start)*(end-start+1))/2);
}

SCMatrixBlock *
SCMatrixLTriBlock::deepcopy() const
{
  SCMatrixLTriBlock *ret = new SCMatrixLTriBlock(start,end);
  ret->blocki = blocki;
  ret->blockj = blockj;
  memcpy(ret->data, data, sizeof(double)*ndat());
  return ret;
}

double *
SCMatrixLTriBlock::dat()
{
  return data;
}

int
SCMatrixLTriBlock::ndat() const
{
  return ((end-start)*(end-start+1))/2;
}

SCMatrixLTriBlock::~SCMatrixLTriBlock()
{
  deallocate(data);
}

void
SCMatrixLTriBlock::process(SCElementOp*op)
{
  SCMatrixLTriBlockIter i(this);
  op->process(i);
}

void
SCMatrixLTriBlock::process(SCElementOp2*op,
                           SCMatrixBlock* b)
{
  SCMatrixLTriBlockIter i(this);
  SCMatrixLTriBlockIter j((SCMatrixLTriBlock*)b);
  op->process(i,j);
}

void
SCMatrixLTriBlock::process(SCElementOp3*op,
                           SCMatrixBlock* b1, SCMatrixBlock* b2)
{
  SCMatrixLTriBlockIter i(this);
  SCMatrixLTriBlockIter j((SCMatrixLTriBlock*)b1);
  SCMatrixLTriBlockIter k((SCMatrixLTriBlock*)b2);
  op->process(i,j,k);
}

/////////////////////////////////////////////////////////////////////////////
// SCMatrixLTriSubBlock member functions

static ClassDesc SCMatrixLTriSubBlock_cd(
  typeid(SCMatrixLTriSubBlock),"SCMatrixLTriSubBlock",1,"public SCMatrixBlock",
  0, 0, create<SCMatrixLTriSubBlock>);

SCMatrixLTriSubBlock::SCMatrixLTriSubBlock(int is, int ie,
                                           int js, int je,
                                           double*d):
  istart(is),
  iend(ie),
  jstart(js),
  jend(je),
  data(d)
{
}

SCMatrixLTriSubBlock::SCMatrixLTriSubBlock(StateIn&s):
  SCMatrixBlock(s)
{
  s.get(istart);
  s.get(iend);
  s.get(jstart);
  s.get(jend);
  data = 0;
}

void
SCMatrixLTriSubBlock::save_data_state(StateOut&s)
{
  SCMatrixBlock::save_data_state(s);
  s.put(istart);
  s.put(iend);
  s.put(jstart);
  s.put(jend);
}

SCMatrixLTriSubBlock::~SCMatrixLTriSubBlock()
{
}

void
SCMatrixLTriSubBlock::process(SCElementOp*op)
{
  SCMatrixLTriSubBlockIter i(this);
  op->process(i);
}

void
SCMatrixLTriSubBlock::process(SCElementOp2*op,
                              SCMatrixBlock* b)
{
  SCMatrixLTriSubBlockIter i(this);
  SCMatrixLTriSubBlockIter j((SCMatrixLTriSubBlock*)b);
  op->process(i,j);
}

void
SCMatrixLTriSubBlock::process(SCElementOp3*op,
                              SCMatrixBlock* b1,
                              SCMatrixBlock* b2)
{
  SCMatrixLTriSubBlockIter i(this);
  SCMatrixLTriSubBlockIter j((SCMatrixLTriSubBlock*)b1);
  SCMatrixLTriSubBlockIter k((SCMatrixLTriSubBlock*)b2);
  op->process(i,j,k);
}

/////////////////////////////////////////////////////////////////////////////
// SCMatrixDiagBlock member functions

static ClassDesc SCMatrixDiagBlock_cd(
  typeid(SCMatrixDiagBlock),"SCMatrixDiagBlock",1,"public SCMatrixBlock",
  0, 0, create<SCMatrixDiagBlock>);

SCMatrixDiagBlock::SCMatrixDiagBlock(int s, int e):
  istart(s),
  jstart(s),
  iend(e)
{
  data = allocate<double>(e-s);
}

SCMatrixDiagBlock::SCMatrixDiagBlock(int is, int ie,int js):
  istart(is),
  jstart(js),
  iend(ie)
{
  data = allocate<double>(ie-is);
}

SCMatrixDiagBlock::SCMatrixDiagBlock(StateIn&s):
  SCMatrixBlock(s)
{
  s.get(istart);
  s.get(jstart);
  s.get(iend);
  s.get(data);
}

void
SCMatrixDiagBlock::save_data_state(StateOut&s)
{
  SCMatrixBlock::save_data_state(s);
  s.put(istart);
  s.put(jstart);
  s.put(iend);
  s.put(data,iend-istart);
}

SCMatrixBlock *
SCMatrixDiagBlock::deepcopy() const
{
  SCMatrixDiagBlock *ret = new SCMatrixDiagBlock(istart,iend,jstart);
  ret->blocki = blocki;
  ret->blockj = blockj;
  memcpy(ret->data, data, sizeof(double)*ndat());
  return ret;
}

double *
SCMatrixDiagBlock::dat()
{
  return data;
}

int
SCMatrixDiagBlock::ndat() const
{
  return iend-istart;
}

SCMatrixDiagBlock::~SCMatrixDiagBlock()
{
  deallocate(data);
}

void
SCMatrixDiagBlock::process(SCElementOp*op)
{
  SCMatrixDiagBlockIter i(this);
  op->process(i);
}

void
SCMatrixDiagBlock::process(SCElementOp2*op,
                           SCMatrixBlock* b)
{
  SCMatrixDiagBlockIter i(this);
  SCMatrixDiagBlockIter j((SCMatrixDiagBlock*)b);
  op->process(i,j);
}

void
SCMatrixDiagBlock::process(SCElementOp3*op,
                           SCMatrixBlock* b1, SCMatrixBlock* b2)
{
  SCMatrixDiagBlockIter i(this);
  SCMatrixDiagBlockIter j((SCMatrixDiagBlock*)b1);
  SCMatrixDiagBlockIter k((SCMatrixDiagBlock*)b2);
  op->process(i,j,k);
}

/////////////////////////////////////////////////////////////////////////////
// SCMatrixDiagSubBlock member functions

static ClassDesc SCMatrixDiagSubBlock_cd(
  typeid(SCMatrixDiagSubBlock),"SCMatrixDiagSubBlock",1,"public SCMatrixBlock",
  0, 0, create<SCMatrixDiagSubBlock>);

SCMatrixDiagSubBlock::SCMatrixDiagSubBlock(int s, int e, int o,
                                           double* d):
  istart(s),
  jstart(s),
  iend(e),
  offset(o),
  data(d)
{
}

SCMatrixDiagSubBlock::SCMatrixDiagSubBlock(int is, int ie,
                                           int js, int o, double* d):
  istart(is),
  jstart(js),
  iend(ie),
  offset(o),
  data(d)
{
}

SCMatrixDiagSubBlock::SCMatrixDiagSubBlock(StateIn&s):
  SCMatrixBlock(s)
{
  s.get(istart);
  s.get(jstart);
  s.get(iend);
  s.get(offset);
  data = 0;
}

void
SCMatrixDiagSubBlock::save_data_state(StateOut&s)
{
  SCMatrixBlock::save_data_state(s);
  s.put(istart);
  s.put(jstart);
  s.put(iend);
  s.put(offset);
}

SCMatrixDiagSubBlock::~SCMatrixDiagSubBlock()
{
}

void
SCMatrixDiagSubBlock::process(SCElementOp*op)
{
  SCMatrixDiagSubBlockIter i(this);
  op->process(i);
}

void
SCMatrixDiagSubBlock::process(SCElementOp2*op,
                              SCMatrixBlock* b)
{
  SCMatrixDiagSubBlockIter i(this);
  SCMatrixDiagSubBlockIter j((SCMatrixDiagSubBlock*)b);
  op->process(i,j);
}

void
SCMatrixDiagSubBlock::process(SCElementOp3*op,
                              SCMatrixBlock* b1,
                              SCMatrixBlock* b2)
{
  SCMatrixDiagSubBlockIter i(this);
  SCMatrixDiagSubBlockIter j((SCMatrixDiagSubBlock*)b1);
  SCMatrixDiagSubBlockIter k((SCMatrixDiagSubBlock*)b2);
  op->process(i,j,k);
}

/////////////////////////////////////////////////////////////////////////////
// SCVectorSimpleBlock member functions

static ClassDesc SCVectorSimpleBlock_cd(
  typeid(SCVectorSimpleBlock),"SCVectorSimpleBlock",1,"public SCMatrixBlock",
  0, 0, create<SCVectorSimpleBlock>);

SCVectorSimpleBlock::SCVectorSimpleBlock(int s, int e):
  istart(s),
  iend(e)
{
  data = allocate<double>(e-s);
}

SCVectorSimpleBlock::SCVectorSimpleBlock(StateIn&s):
  SCMatrixBlock(s)
{
  s.get(istart);
  s.get(iend);
  s.get(data);
}

void
SCVectorSimpleBlock::save_data_state(StateOut&s)
{
  SCMatrixBlock::save_data_state(s);
  s.put(istart);
  s.put(iend);
  s.put(data,iend-istart);
}

SCMatrixBlock *
SCVectorSimpleBlock::deepcopy() const
{
  SCVectorSimpleBlock *ret = new SCVectorSimpleBlock(istart,iend);
  ret->blocki = blocki;
  ret->blockj = blockj;
  memcpy(ret->data, data, sizeof(double)*ndat());
  return ret;
}

double *
SCVectorSimpleBlock::dat()
{
  return data;
}

int
SCVectorSimpleBlock::ndat() const
{
  return iend-istart;
}

SCVectorSimpleBlock::~SCVectorSimpleBlock()
{
  deallocate(data);
}

void
SCVectorSimpleBlock::process(SCElementOp*op)
{
  SCVectorSimpleBlockIter i(this);
  op->process(i);
}

void
SCVectorSimpleBlock::process(SCElementOp2*op,
                             SCMatrixBlock* b)
{
  SCVectorSimpleBlockIter i(this);
  SCVectorSimpleBlockIter j((SCVectorSimpleBlock*)b);
  op->process(i,j);
}

void
SCVectorSimpleBlock::process(SCElementOp3*op,
                             SCMatrixBlock* b1,
                             SCMatrixBlock* b2)
{
  SCVectorSimpleBlockIter i(this);
  SCVectorSimpleBlockIter j((SCVectorSimpleBlock*)b1);
  SCVectorSimpleBlockIter k((SCVectorSimpleBlock*)b2);
  op->process(i,j,k);
}

/////////////////////////////////////////////////////////////////////////////
// SCVectorSimpleSubBlock member functions

static ClassDesc SCVectorSimpleSubBlock_cd(
  typeid(SCVectorSimpleSubBlock),"SCVectorSimpleSubBlock",1,"public SCMatrixBlock",
  0, 0, create<SCVectorSimpleSubBlock>);

SCVectorSimpleSubBlock::SCVectorSimpleSubBlock(int s, int e, int o,
                                               double* d):
  istart(s),
  iend(e),
  offset(o),
  data(d)
{
}

SCVectorSimpleSubBlock::SCVectorSimpleSubBlock(StateIn&s):
  SCMatrixBlock(s)
{
  s.get(istart);
  s.get(iend);
  s.get(offset);
  data = 0;
}

void
SCVectorSimpleSubBlock::save_data_state(StateOut&s)
{
  SCMatrixBlock::save_data_state(s);
  s.put(istart);
  s.put(iend);
  s.put(offset);
}

SCVectorSimpleSubBlock::~SCVectorSimpleSubBlock()
{
}

void
SCVectorSimpleSubBlock::process(SCElementOp*op)
{
  SCVectorSimpleSubBlockIter i(this);
  op->process(i);
}

void
SCVectorSimpleSubBlock::process(SCElementOp2*op,
                                SCMatrixBlock* b)
{
  SCVectorSimpleSubBlockIter i(this);
  SCVectorSimpleSubBlockIter j((SCVectorSimpleSubBlock*)b);
  op->process(i,j);
}

void
SCVectorSimpleSubBlock::process(SCElementOp3*op,
                                SCMatrixBlock* b1,
                                SCMatrixBlock* b2)
{
  SCVectorSimpleSubBlockIter i(this);
  SCVectorSimpleSubBlockIter j((SCVectorSimpleSubBlock*)b1);
  SCVectorSimpleSubBlockIter k((SCVectorSimpleSubBlock*)b2);
  op->process(i,j,k);
}

///////////////////////////////////////////////////////////////////////
// SCMatrixSubblockIter

SCMatrixSubblockIter::~SCMatrixSubblockIter()
{
}

///////////////////////////////////////////////////////////////////////
// SCMatrixSimpleSubblockIter

SCMatrixSimpleSubblockIter::SCMatrixSimpleSubblockIter(
    Access access_,
    const Ref<SCMatrixBlock> &b):
  SCMatrixSubblockIter(access_)
{
  block_ = b;
}

void
SCMatrixSimpleSubblockIter::begin()
{
  if (block_) ready_ = 1;
  else ready_ = 0;
}

int
SCMatrixSimpleSubblockIter::ready()
{
  return ready_;
}

void
SCMatrixSimpleSubblockIter::next()
{
  ready_ = 0;
}

SCMatrixBlock *
SCMatrixSimpleSubblockIter::block()
{
  return block_.pointer();
}

///////////////////////////////////////////////////////////////////////
// SCMatrixListSubblockIter

SCMatrixListSubblockIter::SCMatrixListSubblockIter(
    Access access,
    const Ref<SCMatrixBlockList> &list
    ):
  SCMatrixSubblockIter(access),
  list_(list)
{
}

void
SCMatrixListSubblockIter::begin()
{
  iter_ = list_->begin();
}

int
SCMatrixListSubblockIter::ready()
{
  return iter_ != list_->end();
}

void
SCMatrixListSubblockIter::next()
{
  iter_++;
}

SCMatrixBlock *
SCMatrixListSubblockIter::block()
{
  return iter_.block();
}

///////////////////////////////////////////////////////////////////////
// SCMatrixNullSubblockIter

SCMatrixNullSubblockIter::SCMatrixNullSubblockIter():
  SCMatrixSubblockIter(None)
{
}

SCMatrixNullSubblockIter::SCMatrixNullSubblockIter(Access access):
  SCMatrixSubblockIter(access)
{
}

void
SCMatrixNullSubblockIter::begin()
{
}

int
SCMatrixNullSubblockIter::ready()
{
  return 0;
}

void
SCMatrixNullSubblockIter::next()
{
}

SCMatrixBlock *
SCMatrixNullSubblockIter::block()
{
  return 0;
}

///////////////////////////////////////////////////////////////////////
// SCMatrixCompositeSubblockIter

SCMatrixCompositeSubblockIter::SCMatrixCompositeSubblockIter(
    Ref<SCMatrixSubblockIter>& i1,
    Ref<SCMatrixSubblockIter>& i2):
  SCMatrixSubblockIter(None)
{
  niters_ = 0;
  if (i1) { niters_++; }
  if (i2) { niters_++; }
  iters_ = new Ref<SCMatrixSubblockIter>[niters_];
  iiter_ = 0;
  if (i1) { iters_[iiter_] = i1; iiter_++; }
  if (i2) { iters_[iiter_] = i2; iiter_++; }

  if (niters_) access_ = iters_[0]->access();
  for (int i=0; i<niters_; i++) {
      if (iters_[i]->access() != access_) {
          ExEnv::errn() << "SCMatrixCompositeSubblockIter: access not compatible"
               << endl;
          abort();
        }
    }
}

SCMatrixCompositeSubblockIter::SCMatrixCompositeSubblockIter(
    Access access_,
    int niters):
  SCMatrixSubblockIter(access_)
{
  niters_ = niters;
  iters_ = new Ref<SCMatrixSubblockIter>[niters_];
}

SCMatrixCompositeSubblockIter::~SCMatrixCompositeSubblockIter()
{
  delete[] iters_;
}

void
SCMatrixCompositeSubblockIter::set_iter(int i,
                                        const Ref<SCMatrixSubblockIter>& iter)
{
  iters_[i] = iter;
  if (iters_[i]->access() != access_) {
      ExEnv::errn() << "SCMatrixCompositeSubblockIter: access not compatible"
           << endl;
      abort();
    }
}

void
SCMatrixCompositeSubblockIter::begin()
{
  if (niters_ == 0) return;
  iiter_ = 0;
  iters_[iiter_]->begin();
  while (!iters_[iiter_]->ready()) {
      if (iiter_ < niters_-1) {
          iiter_++;
          iters_[iiter_]->begin();
        }
      else break;
    }
}

int
SCMatrixCompositeSubblockIter::ready()
{
  return iters_[iiter_]->ready();
}

void
SCMatrixCompositeSubblockIter::next()
{
  iters_[iiter_]->next();
  while (!iters_[iiter_]->ready()) {
      if (iiter_ < niters_-1) {
          iiter_++;
          iters_[iiter_]->begin();
        }
      else break;
    }
}

SCMatrixBlock *
SCMatrixCompositeSubblockIter::block()
{
  return iters_[iiter_]->block();
}

///////////////////////////////////////////////////////////////////////
// SCMatrixJointSubblockIter

SCMatrixJointSubblockIter::SCMatrixJointSubblockIter(
    const Ref<SCMatrixSubblockIter>& i1,
    const Ref<SCMatrixSubblockIter>& i2,
    const Ref<SCMatrixSubblockIter>& i3,
    const Ref<SCMatrixSubblockIter>& i4,
    const Ref<SCMatrixSubblockIter>& i5):
  SCMatrixSubblockIter(None)
{
  niters_ = 0;
  if (i1) { niters_++; }
  if (i2) { niters_++; }
  if (i3) { niters_++; }
  if (i4) { niters_++; }
  if (i5) { niters_++; }
  iters_ = new Ref<SCMatrixSubblockIter>[niters_];
  int i = 0;
  if (i1) { iters_[i] = i1; i++; }
  if (i2) { iters_[i] = i2; i++; }
  if (i3) { iters_[i] = i3; i++; }
  if (i4) { iters_[i] = i4; i++; }
  if (i5) { iters_[i] = i5; i++; }
}

SCMatrixJointSubblockIter::~SCMatrixJointSubblockIter()
{
  delete[] iters_;
}

void
SCMatrixJointSubblockIter::begin()
{
  for (int i=0; i<niters_; i++) {
      iters_[i]->begin();
    }
}

int
SCMatrixJointSubblockIter::ready()
{
  int nready = 0;
  for (int i=0; i<niters_; i++) {
      nready += (iters_[i]->ready()?1:0);
    }

  if (nready == niters_)
    return 1;
  else if (!nready)
    return 0;

  ExEnv::errn() << "SCMatrixJointSubblockIter: incompatible iterators" << endl;
  abort();
  return 0;
}

void
SCMatrixJointSubblockIter::next()
{
  for (int i=0; i<niters_; i++) {
      iters_[i]->next();
    }
}

SCMatrixBlock *
SCMatrixJointSubblockIter::block()
{
  return block(0);
}

SCMatrixBlock *
SCMatrixJointSubblockIter::block(int b)
{
  return iters_[b]->block();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
