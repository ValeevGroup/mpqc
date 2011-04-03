//
// tbgrad.h --- definition of the abstract two-electron gradient builder
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
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

#ifndef _chemistry_qc_scf_tbgrad_h
#define _chemistry_qc_scf_tbgrad_h

#include <util/group/thread.h>

namespace sc {

template<class T>
class TBGrad : public Thread {
  protected:
    T& contribution;
    double exchange_fraction;

  public:
    TBGrad(T&t, double ex = 1.0) : contribution(t), exchange_fraction(ex) {}
    virtual ~TBGrad() {}

    inline void set_scale(double& coulombscale, double& exchangescale,
                          int i, int j, int k, int l) const
    {
      double scale = 1.0;

      if ((i!=k)||(j!=l))
        scale *= 2.0;

      if (i!=j)
        scale *= 2.0;

      coulombscale = 0.5*scale;
      exchangescale = -0.25*scale * exchange_fraction;

      if (k!=l)
        coulombscale *= 2.0;

      if ((k!=l)&&(i==j))
        exchangescale *= 2.0;
    }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
