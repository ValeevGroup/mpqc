//
// ssar2tem.h
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

template <class Type>
class SSBArray2: public Array2<Type> {
  public:
    SSBArray2() {}
    SSBArray2(const Array2<Type> &a): Array2<Type>(a) {}
    SSBArray2(Type* data,int size0,int size1): Array2<Type>(data,size0,size1){}
    SSBArray2(int size0,int size1): Array2<Type>(size0,size1) {}
    SSBArray2(StateIn&s) {
        s.get(_length0);
        s.get(_length1);
        if (_length0&&_length1) s.get(_array);
        _managed=1;
      }
    void save_object_state(StateOut&s) {
        s.put(_length0);
        s.put(_length1);
        if (_length0&&_length1) s.put(_array,_length0*_length1);
      }
};

// ///////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
