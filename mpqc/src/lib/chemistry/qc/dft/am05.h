//
// am05.h -- derived from:
// functional.h --- definition of the dft functional
//
// Copyright (C) 1997 Limit Point Systems, Inc.
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

#ifndef _chemistry_qc_dft_am05_h
#define _chemistry_qc_dft_am05_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/dft/functional.h>

namespace sc {

/** Implements the Perdew-Burke-Ernzerhof (PBE) correlation functional.

    John P. Perdew, Kieron Burke, and Yue Wang, Phys. Rev. B, 54(23),
    pp. 16533-16539, 1996.

    John P. Perdew, Kieron Burke, and Matthias Ernzerhof, Phys. Rev. Lett.,
    77(18), pp. 3865-3868, 1996.
*/
class AM05Functional: public DenFunctional {
  protected:
    void am05xc(double rho,double gam,
                double &fxc,double &dfxcdrho,double &dfxcdgamma);
    double am05_lambertw(double z);
    void xs(double rho,double &ex, double &vx);
    void pw(double rhoa,double rhob,
            double &ec,double &vca,double &vcb);
    void cpbe_lsd(double rhoa,double rhob,
                  double &eps,double &vca,double &vcb);
    void pbe_gcor(double a,double a1,
                  double b1,double b2,double b3,double b4,
                  double rtrs,
                  double &gg,double &ggrs);
  public:
    AM05Functional();
    AM05Functional(const Ref<KeyVal> &);
    AM05Functional(StateIn &);
    ~AM05Functional();
    void save_data_state(StateOut &);
    int need_density_gradient();
    int need_density_hessian();
    void point(const PointInputData&, PointOutputData&);
    void set_spin_polarized(int);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
