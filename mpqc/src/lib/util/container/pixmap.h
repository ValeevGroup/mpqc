//
// pixmap.h
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

#ifndef _pixmap_h
#define _pixmap_h

#include <util/container/pixvpMap.h>
#include <util/container/pixvpRAVLMap.h>

template <class Type>
class PixMap
{
 private:
  // Since PixVoidPtrMap does use const members, trickery must be used.
  inline PixVoidPtrMap& cimpl() const { return *((PixVoidPtrMap*)&impl); }

  PixVoidPtrRAVLMap impl;
  Pix guess;
 public:
  inline PixMap():guess((Pix)10000),impl((VoidPtr)0) {}
  virtual ~PixMap() {
      for (Pix i=first(); i!=0; next(i)) { delete (Type*) impl.contents(i); }
    }
  inline int length() const { return cimpl().length(); }
  inline int empty() const { return cimpl().empty(); }
  inline int contains(Pix key) const { return cimpl().contains(key); }
  inline void clear() { impl.clear(); }
  inline Type& operator[](Pix key) { 
      if (impl.seek(key) == 0) {
	  impl[key] = (void*) new Type;
	}
      Type* tmp = (Type*)(impl[key].getptr());
      return *tmp;
    }
  inline const Type& operator[](Pix key) const { 
      if (cimpl().seek(key) == 0) {
	  cimpl()[key] = (void*) new Type;
	}
      Type* tmp = (Type*)(cimpl()[key].getptr());
      return *tmp;
    }
  inline void del(Pix key) { impl.del(key); }
  inline Pix first() const { return cimpl().first(); }
  inline void next(Pix&i) const { cimpl().next(i); }
  inline Type& contents(Pix i) { return *((Type*)impl.contents(i).getptr()); }
  inline int owns(Pix i) const { return cimpl().owns(i); }
  inline Pix seek(Pix key) const { return cimpl().seek(key); }
  inline Pix key(Pix pix) const { return cimpl().key(pix); }

  // this finds a key that isn't already in use
  inline Pix newkey() {
    while (impl.seek(guess)) ((int)guess)++;
    return guess;
    }
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
