//
// ptrset.h
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

#ifndef _ptrset_h
#define _ptrset_h

#include <util/container/voidptrSet.h>
#include <util/container/voidptrAVLSet.h>

template <class Type>
class PtrSet
{
 private:

  VoidPtrSet* setimpl;

 public:
                       PtrSet():setimpl(new VoidPtrAVLSet) {}
  virtual              ~PtrSet() { delete setimpl; }

  int                   length() { return setimpl->length(); }
  int                   empty() { return setimpl->empty(); }

  virtual Pix           add(Type* item) { return setimpl->add((void*)item); }
  virtual void          del(Type* item) { setimpl->del((void*)item); }
  virtual int           contains(Type*  item)
                            { return setimpl->contains((void*)item); }

  virtual void          clear() { setimpl->clear(); }

  virtual Pix           first() { return setimpl->first(); }
  virtual void          next(Pix& i) { return setimpl->next(i); }
  virtual Type*         operator () (Pix i) { return (*setimpl)(i); }

  virtual int           owns(Pix i) { return setimpl->owns(i); }
  virtual Pix           seek(Type*  item) { return setimpl->seek((void*)item); }

  void                  operator |= (PtrSet<Type>& b)
                            { setimpl->operator|=(*b.setimpl); }
  void                  operator -= (PtrSet<Type>& b)
                            { setimpl->operator-=(*b.setimpl); }
  void                  operator &= (PtrSet<Type>& b)
                            { setimpl->operator&=(*b.setimpl); }

  int                   operator == (PtrSet<Type>& b)
                            { return setimpl->operator==(*b.setimpl); }
  int                   operator != (PtrSet<Type>& b)
                            { return setimpl->operator!=(*b.setimpl); }
  int                   operator <= (PtrSet<Type>& b)
                            { return setimpl->operator<=(*b.setimpl); }

};

#endif
