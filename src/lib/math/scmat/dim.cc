//
// dim.cc
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

#include <cmath>
#include <numeric>
#include <assert.h>

#include <util/keyval/keyval.h>
#include <math/scmat/dim.h>
#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <iomanip>

using namespace std;
using namespace sc;

static void
fail(const char *s)
{
  ExEnv::errn() << indent << "math/scmat/dim.cc: " << s << endl;
  abort();
}

/////////////////////////////////////////////////////////////////////////////
// SCBlockInfo member functions

static ClassDesc SCBlockInfo_cd(
  typeid(SCBlockInfo),"SCBlockInfo",1,"public SavableState",
  0, create<SCBlockInfo>, create<SCBlockInfo>);

SCBlockInfo::SCBlockInfo(int n, int nblocks, const int *blocksizes):
  subdims_(0)
{
  n_ = n;
  nblocks_ = nblocks;
  size_ = 0;

  // validate input
  if (blocksizes) {
    const int nn = std::accumulate(blocksizes, blocksizes+nblocks, 0);
    MPQC_ASSERT(n == nn);
  }

  if (n_ != 0 && nblocks_ == 0) {
      nblocks_ = 1;
      size_ = new int[1];
      size_[0] = n;
    }
  else if (nblocks_ == 0) {
      size_ = 0;
    }
  else {
      int i;
      size_ = new int[nblocks_];
      if (blocksizes) {
          for (i=0; i<nblocks_; i++) {
              size_[i] = blocksizes[i];
            }
        }
      else {
          int nper = n/nblocks;
          int nleft = n%nblocks;
          for (i=0; i<nblocks_; i++) {
              size_[i] = nper;
              if (i<nleft) size_[i]++;
            }
        }
    }

  init_start();
}

SCBlockInfo::SCBlockInfo(const Ref<KeyVal>&keyval):
  subdims_(0)
{
  nblocks_ = keyval->count("sizes");
  n_ = 0;
  size_ = new int[nblocks_];
  for (int i=0; i<nblocks_; i++) {
      size_[i] = keyval->intvalue("sizes",i);
      n_ += size_[i];
    }
  int nsubdims = keyval->count("subdims");
  if (nsubdims) {
      if (nblocks_ != 0 && nsubdims != nblocks_) {
          fail("SCBlockInfo(const Ref<KeyVal>&): nsubdims != nblocks");
        }
      subdims_ = new RefSCDimension[nsubdims];
      for (int i=0; i<nsubdims; i++) {
          subdims_[i] << keyval->describedclassvalue("subdims",i);
        }
      if (nblocks_ == 0) {
          delete[] size_;
          size_ = new int[nsubdims];
          for (int i=0; i<nsubdims; i++) {
              size_[i] = subdims_[i].n();
              n_ += size_[i];
            }
          nblocks_ = nsubdims;
        }
    }
  init_start();
}

SCBlockInfo::SCBlockInfo(StateIn&s):
  SavableState(s)
{
  s.get(n_);
  s.get(nblocks_);
  s.get(size_);
  s.get(start_);
  int have_subdims;
  s.get(have_subdims);
  if (have_subdims) {
      subdims_ = new RefSCDimension[nblocks_];
      for (int i=0; i<nblocks_; i++) {
          subdims_[i] << SavableState::restore_state(s);
        }
    }
  else {
      subdims_ = 0;
    }
}

void
SCBlockInfo::save_data_state(StateOut&s)
{
  s.put(n_);
  s.put(nblocks_);
  s.put(size_,nblocks_);
  s.put(start_,nblocks_);
  if (subdims_) {
      s.put(1);
      for (int i=0; i<nblocks_; i++) {
          SavableState::save_state(subdims_[i].pointer(),s);
        }
    }
  else {
      s.put(0);
    }
}

void
SCBlockInfo::init_start()
{
  start_ = 0;
  if (nblocks_) {
      start_ = new int[nblocks_];
      start_[0] = 0;
      for (int i=1; i<nblocks_; i++) {
          start_[i] = start_[i-1] + size_[i-1];
        }
    }
}

int
SCBlockInfo::equiv(SCBlockInfo *bi)
{
  if (bi == 0) return 0;
  if (n_ != bi->n_) return 0;
  if (nblocks_ != bi->nblocks_) return 0;
  int i;
  for (i=0; i<nblocks_; i++) {
      if (start_[i] != bi->start_[i]) return 0;
    }
  if (subdims_) {
      if (bi->subdims_) {
          for (i=0; i<nblocks_; i++) {
              if (subdims_[i].nonnull()) {
                  if (bi->subdims_[i].nonnull()) {
                      if (!subdims_[i]->equiv(bi->subdims_[i])) return 0;
                    }
                  else return 0;
                }
              else if (bi->subdims_[i].nonnull()) return 0;
            }
        }
      else {
          return 0;
        }
    }
  else if (bi->subdims_) {
      return 0;
    }
  return 1;
}

void
SCBlockInfo::elem_to_block(int elem, int &block, int &offset)
{
  for (int i=nblocks_-1; i>=0; i--) {
      if (start_[i] <= elem) {
          block = i;
          offset = elem-start_[i];
          return;
        }
    }
  fail("SCBlockInfo::elem_to_block: couldn't find block");
}

SCBlockInfo::~SCBlockInfo()
{
  delete[] start_;
  delete[] size_;
  delete[] subdims_;
}

RefSCDimension
SCBlockInfo::subdim(int i)
{
  if (i >= nblocks_ || i < 0) {
      fail("SCBlockInfo::subdim: bad block index");
    }
  if (!subdims_) return 0;
  return subdims_[i];
}

void
SCBlockInfo::set_subdim(int i, const RefSCDimension &dim)
{
  if (!subdims_) {
      subdims_ = new RefSCDimension[nblocks_];
    }
  if (i >= nblocks_ || i < 0) {
      fail("SCBlockInfo::set_subdim: bad block index");
    }
  if (size(i) != dim.n()) {
      fail("SCBlockInfo::set_subdim: size mismatch");
    }
  subdims_[i] = dim;
}

void
SCBlockInfo::print(ostream&o) const
{
  indent(o); o << "nblocks = " << nblocks_ << endl;
  indent(o); o << "sizes =";
  for (int i=0; i<nblocks_; i++) {
      o << " " << size_[i];
    }
  o << endl;
  if (subdims_) {
      for (int i=0; i<nblocks_; i++) {
          indent(o); o << "subdimension " << i << ":" << endl;
          incindent(o);
          subdims_[i].print(o);
          decindent(o);
        }
    }
}

/////////////////////////////////////////////////////////////////////////
// SCDimension members

static ClassDesc SCDimension_cd(
  typeid(SCDimension),"SCDimension",1,"public SavableState",
  0, create<SCDimension>, create<SCDimension>);

SCDimension::SCDimension(const char* name):
  n_(0)
{
  if (name) name_ = name;
}

SCDimension::SCDimension(int n, const char* name):
  n_(n)
{
  if (name) name_ = name;
  blocks_ = new SCBlockInfo(n, 1);
}

SCDimension::SCDimension(const Ref<SCBlockInfo>& b, const char* name):
  n_(b->nelem()), blocks_(b)
{
  if (name) name_ = name;
}

SCDimension::SCDimension(int n,
                         int nblocks, const int *blocksizes,
                         const char* name):
  n_(n)
{
  if (name) name_ = name;
  blocks_ = new SCBlockInfo(n, nblocks, blocksizes);
}

SCDimension::SCDimension(const Ref<KeyVal>& keyval)
{
  blocks_ << keyval->describedclassvalue("blocks");
  n_ = keyval->intvalue("n");
  if (blocks_.null()) {
      if (keyval->error() != KeyVal::OK) {
          fail("SCDimension(const Ref<KeyVal>&): missing input");
        }
      blocks_ = new SCBlockInfo(n_);
    }
  else {
      if (n_ != 0 && n_ != blocks_->nelem()) {
          fail("SCDimension(const Ref<KeyVal>&): inconsistent sizes");
        }
      n_ = blocks_->nelem();
    }
  name_ = keyval->stringvalue("name");
}

SCDimension::SCDimension(StateIn&s):
  SavableState(s)
{
  s.get(name_);
  s.get(n_);
  blocks_ << SavableState::restore_state(s);
}

void
SCDimension::save_data_state(StateOut&s)
{
  s.put(name_);
  s.put(n_);
  SavableState::save_state(blocks_.pointer(), s);
}

SCDimension::~SCDimension()
{
}

int
SCDimension::equiv(const SCDimension *a) const
{
  if (n_ != a->n_) return 0;

  if (!blocks_->equiv(a->blocks_)) return 0;

  return 1;
}

void
SCDimension::print(ostream&o) const
{
  indent(o); o << "n = " << n_;
  if (!name_.empty()) {
      o << ", name = " << name_;
    }
  o << endl;
  incindent(o);
  if (blocks_.nonnull()) blocks_->print(o);
  decindent(o);
}

/////////////////////////////////////////////////////////////////////////////
// RefSCDimension member functions

void
RefSCDimension::print(ostream&o) const
{
  if (null()) { indent(o); o << "n = 0" << endl; }
  else pointer()->print(o);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
