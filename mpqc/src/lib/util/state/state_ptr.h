//
// state_ptr.h
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

#ifndef _util_state_state_ptr_h
#define _util_state_state_ptr_h

#ifdef __GNUC__
#pragma interface
#endif

class StateData {
  public:
    void* ptr;
    int num;
    int size;
    int type;
    int offset;
    int next_reference;
    int can_refer;
  private:
    void init();
  public:
    StateData(int n);
    StateData(void *p);
    StateData(int n, void *p);
};

class StateDataNum: public StateData {
  public:
    StateDataNum(int n):StateData(n) {}
    StateDataNum(int n,void*p):StateData(n,p) {}
    int operator==(const StateDataNum&n) const { return num == n.num; }
    int compare(const StateDataNum&n) const;
  };

inline int StateDataNum::compare(const StateDataNum&p) const
{
  return (p.num == num)?0:((p.num<num)?1:-1);
}

class StateDataPtr: public StateData {
  public:
    StateDataPtr(void*p):StateData(p) {}
    int operator==(const StateDataPtr&p) const { return ptr == p.ptr; }
    int compare(const StateDataPtr&p) const;
  };

inline int StateDataPtr::compare(const StateDataPtr&p) const
{
  return (p.ptr == ptr)?0:((p.ptr<ptr)?1:-1);
}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
