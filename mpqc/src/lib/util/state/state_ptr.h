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

#ifndef _libqc_state_ptr_h
#define _libqc_state_ptr_h

#ifdef __GNUC__
#pragma interface
#endif

class StateDataNum {
  private:
    int num_;
    void* ptr_;
  public:
    StateDataNum(int num):num_(num),ptr_(0) {}
    StateDataNum(int num,void*p):num_(num),ptr_(p) {}
    int num() const { return num_; }
    void* ptr() const { return ptr_; }
    void assign_ptr(void*p) { ptr_ = p; }
    int operator==(const StateDataNum&num) const { return num_ == num.num_; }
    int compare(const StateDataNum&num) const;
  };

inline int StateDataNum::compare(const StateDataNum&p) const
{
  return (p.num_ == num_)?0:((p.num_<num_)?1:-1);
}

class StateDataPtr {
  private:
    int num_;
    void* ptr_;
  public:
    StateDataPtr(void*p):ptr_(p) {}
    void* ptr() const { return ptr_; }
    int num() const { return num_; }
    void assign_num(int n) { num_ = n; }
    int operator==(const StateDataPtr&p) const { return ptr_ == p.ptr_; }
    int compare(const StateDataPtr&p) const;
  };

inline int StateDataPtr::compare(const StateDataPtr&p) const
{
  return (p.ptr_ == ptr_)?0:((p.ptr_<ptr_)?1:-1);
}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
