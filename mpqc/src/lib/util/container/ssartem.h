//
// ssartem.h
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
#pragma interface
#endif

template <class Type>
class SSBArray: public Array<Type> {
  public:
    SSBArray() {}
    SSBArray(const Array<Type>&a): Array<Type>(a) {}
    SSBArray(Type* data,int size): Array<Type>(data,size) {}
    SSBArray(int size): Array<Type>(size) {}
    SSBArray(StateIn&s) {
      s.get(_length);
      if (_length) s.get(_array);
      _managed=1;
    }
    void save_object_state(StateOut&s) {
        s.put(_length);
        if (_length) s.put(_array,_length);
      }
};

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
