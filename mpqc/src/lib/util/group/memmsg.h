//
// memmsg.h
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

#ifndef _util_group_memmsg_h
#define _util_group_memmsg_h

#include <util/group/message.h>
#include <util/group/memory.h>

namespace sc {

/** A MsgMemoryGrp that initializes its data using a messagegrp. */
class MsgMemoryGrp: public MemoryGrp {

  protected:
    Ref<MessageGrp> msg_;
  public:
    MsgMemoryGrp(const Ref<MessageGrp>& msg);
    MsgMemoryGrp(const Ref<KeyVal>& keyval);
    ~MsgMemoryGrp();
    void set_localsize(size_t localsize);

    void sync();
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
