//
// asettmpl.h --- template for the Arrayset class
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

template <class Type>
class Arrayset : public Set<Type>
{
 public:
    Arrayset() {}
    Arrayset(const Arrayset<Type>&a): Set<Type>(a) {}
    ~Arrayset() {}
    Type& operator[](int i)
    {
      if (i<0 || i>=nelement) {
          cerr << "Arrayset::operator[] out of range: "
               << i << " (nelement = " << nelement << ")" << endl;
          abort();
        };
      return element[i];
    }
    const Type& operator[](int i) const
    {
      if (i<0 || i>=nelement) {
          cerr << "Arrayset::operator[] out of range: "
               << i << " (nelement = " << nelement << ")" << endl;
          abort();
        };
      return element[i];
    }
    int iseek(const Type &item)
    {
      for (int i=0; i<nelement; i++) {
          if (item == element[i]) return i;
        }
      return -1;
    }
};

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
