//
// state_ptr.cc
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <util/class/class.h>
#include <util/state/state.h>
#include <util/state/stateptrImplSet.h>
#include <util/state/statenumImplSet.h>

// Returns the number of the new object, if this object is new.
// If the object is found then zero is returned.
int StateIn::getpointer(void**p)
{
  int refnum;
  get(refnum);
  if (refnum == 0) {
    *p = 0;
    //cout << "StateOut::getpointer: pointer is 0" << endl;
    return 0;
    }
  StateDataNum num(refnum);
  Pix ind = (ps_?ps_->seek(num):0);
  //cout << "StateOut::getpointer: looking for " << refnum << " and got "
  //     << (int)ind << endl;
  if (ind == 0) {
    *p = 0;
    return refnum;
    }
  else {
    *p = ((*this->ps_)(ind)).ptr();
    //cout << "StateOut::getpointer: pointer is made 0x"
    //     << setbase(16) << *p << endl;
    return 0;
    }
  }

void StateIn::nextobject(int objnum)
{
  _nextobject = objnum;
}

void StateIn::havepointer(void*p)
{
  if (_nextobject) {
      havepointer(_nextobject,p);
      _nextobject = 0;
    }
}

void StateIn::havepointer(int objnum,void*p)
{
  StateDataNum num(objnum,p);
  if (ps_) ps_->add(num);
}

// Returns 0 if the object has already been written.
// Returns 1 if the object must yet be written.
int StateOut::putpointer(void*p)
{
  if (p == 0) {
    put(0);
    return 0;
    }
  StateDataPtr dp(p);
  Pix ind = (ps_?ps_->seek(dp):0);
  //cout << "StateOut::putpointer: ind = " << (int)ind << " for 0x"
  //     << setbase(16) << p << endl;
  if (ind == 0) {
      if (ps_) {
          dp.assign_num(next_pointer_number++);
          ps_->add(dp);
        }
      put(dp.num());
      return 1;
    }
  else {
      put((*this->ps_)(ind).num());
      return 0;
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
