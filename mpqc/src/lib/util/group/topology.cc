//
// topology.cc
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

#ifdef __GNUG__
#pragma implementation
#endif

#include <util/group/message.h>
#include <util/group/topology.h>

using namespace sc;

static ClassDesc GlobalMsgIter_cd(
  typeid(GlobalMsgIter),"GlobalMsgIter",1,"public DescribedClass",
  0, 0, 0);

GlobalMsgIter::GlobalMsgIter(int nproc, int me, int root)
{
  nproc_ = nproc;
  me_ = me;
  root_ = root;
  forwards();
}

GlobalMsgIter::~GlobalMsgIter()
{
}

static ClassDesc MachineTopology_cd(
  typeid(MachineTopology),"MachineTopology",1,"public DescribedClass",
  0, 0, 0);

MachineTopology::MachineTopology()
{
}

MachineTopology::MachineTopology(const Ref<KeyVal>&)
{
}

MachineTopology::~MachineTopology()
{
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
