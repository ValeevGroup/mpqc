//
// avlaset.h
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

#include <util/container/pixintRAVLMap.h>
#include <util/container/intpixRAVLMap.h>

template <class Type>
class Arrayset : public AVLSet<Type>
{
  private:
    int nelement;
    PixintRAVLMap pixtoint;
    intPixRAVLMap inttopix;
  public:
    Arrayset(): pixtoint(-1), inttopix(0), nelement(0) {}
    Arrayset(const Arrayset<Type>&a):
      AVLSet<Type>(a),pixtoint(-1),inttopix(0),nelement(0) {}
    ~Arrayset() {}
    Type& operator[](int i)
    {
      return operator()(inttopix[i]);
    }
    const Type& operator[](int i) const
    {
      // must cast away the constness
      Arrayset<Type> *ncthis = (Arrayset<Type>*) this;
      return ncthis->operator()(ncthis->inttopix[i]);
    }
    int iseek(const Type &item)
    {
      return pixtoint[seek(item)];
    }
    Pix add(Type& item) {
        Pix r = AVLSet<Type>::add(item);
        pixtoint[r] = nelement;
        inttopix[nelement] = r;
        nelement++;
      }
    void del(Type& item) {
        Pix p = seek(item);
        int i = pixtoint[p];
        AVLSet<Type>::del(item);
        pixtoint.del(p);
        inttopix.del(i);
      }
};

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
