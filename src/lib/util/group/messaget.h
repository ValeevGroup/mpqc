//
// messaget.h
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

#ifndef _util_group_messaget_h
#define _util_group_messaget_h

#ifdef __GNUC__
#pragma implementation
#endif

#include <util/group/message.h>

namespace sc {

template <class T>
void
GrpSumReduce<T>::reduce(T*target, T*data, int nelement)
{
  for (int i=0; i<nelement; i++) {
      target[i] += data[i];
    }
}

template <class T>
void
GrpMinReduce<T>::reduce(T*target, T*data, int nelement)
{
  for (int i=0; i<nelement; i++) {
      if (target[i] > data[i]) target[i] = data[i];
    }
}

template <class T>
void
GrpMaxReduce<T>::reduce(T*target, T*data, int nelement)
{
  for (int i=0; i<nelement; i++) {
      if (target[i] < data[i]) target[i] = data[i];
    }
}

template <class T>
void
GrpArithmeticAndReduce<T>::reduce(T*target, T*data, int nelement)
{
  for (int i=0; i<nelement; i++) {
      target[i] = target[i] & data[i];
    }
}

template <class T>
void
GrpArithmeticOrReduce<T>::reduce(T*target, T*data, int nelement)
{
  for (int i=0; i<nelement; i++) {
      target[i] = target[i] | data[i];
    }
}

template <class T>
void
GrpArithmeticXOrReduce<T>::reduce(T*target, T*data, int nelement)
{
  for (int i=0; i<nelement; i++) {
      target[i] = target[i] ^ data[i];
    }
}

template <class T>
void
GrpProductReduce<T>::reduce(T*target, T*data, int nelement)
{
  for (int i=0; i<nelement; i++) {
      target[i] *= data[i];
    }
}

template <class T>
void
GrpFunctionReduce<T>::reduce(T*target, T*data, int nelement)
{
  (*func_)(target,data,nelement);
}

}

#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
