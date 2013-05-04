//
// memmsg.cc
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

#ifndef _util_group_memmsg_cc
#define _util_group_memmsg_cc

#include <util/group/memmsg.h>

using namespace std;
using namespace sc;

static ClassDesc MsgMemoryGrp_cd(
  typeid(MsgMemoryGrp),"MsgMemoryGrp",1,"public MemoryGrp",
  0, 0, 0);

MsgMemoryGrp::MsgMemoryGrp(const Ref<MessageGrp> &msg)
{
  msg_ = msg;
  n_ = msg->n();
  me_ = msg->me();
}

MsgMemoryGrp::MsgMemoryGrp(const Ref<KeyVal> &keyval):
  MemoryGrp(keyval)
{
  Ref<MessageGrp> msg; msg << keyval->describedclassvalue("message");
  if (msg.null()) {
      msg = MessageGrp::get_default_messagegrp();
    }
  if (msg.null()) {
      ExEnv::errn() << "MsgMemoryGrp(const Ref<KeyVal>&): couldn't find MessageGrp"
           << endl;
      abort();
    }

  msg_ = msg;
  n_ = msg->n();
  me_ = msg->me();
  offsets_ = new distsize_t[n_ + 1];
  std::fill(offsets_, offsets_ + n_ + 1, distsize_t(0));
}

MsgMemoryGrp::~MsgMemoryGrp()
{
}

void
MsgMemoryGrp::set_localsize(size_t localsize)
{
  if (offsets_ != 0) delete[] offsets_;

  offsets_ = new distsize_t[n_ + 1];
  size_t *sizes = new size_t[n_];

  for (size_t i=0; i<n_; i++) sizes[i] = 0;
  sizes[me_] = localsize;

  msg_->sum(sizes, n_);

  offsets_[0] = 0;
  for (int i=1; i<=n_; i++) {
      offsets_[i] = sizes[i-1] + offsets_[i-1];
      if (offsets_[i] < offsets_[i-1]) {
          ExEnv::errn() << "MsgMemoryGrp::set_localsize: distsize_t cannot handle biggest size" << endl;
          abort();
        }
    }

  delete[] sizes;
}

void
MsgMemoryGrp::sync()
{
  msg_->sync();
#if 0
  // debug code
  for (int i=0; i<n(); i++) {
      if (i == me()) {
          cout << "Data on node " << i << ":" << endl;
          for (size_t j=0; j<localsize()/sizeof(double); j++) {
              double dat = ((double*)localdata())[j];
              if (fabs(dat) > 1.e-12) {
                  cout << " " << i << " " << offset(me()) + j << " "
                       << setw(6)
                       << dat << endl;
                }
            }
        }
      cout.flush();
      msg_->sync();
    }
  // end debug code
#endif
}

#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
