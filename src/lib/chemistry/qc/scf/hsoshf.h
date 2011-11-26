//
// hsoshf.h --- definition of the high-spin open shell Hartree-Fock SCF class
//
// Copyright (C) 1997 Limit Point Systems, Inc.
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

#ifndef _chemistry_qc_scf_hsoshf_h
#define _chemistry_qc_scf_hsoshf_h

#include <chemistry/qc/scf/hsosscf.h>

namespace sc {

// //////////////////////////////////////////////////////////////////////////

/// HSOSHF is a Hartree-Fock specialization of HSOSSCF.
class HSOSHF: public HSOSSCF {
  public:
    HSOSHF(StateIn&);
    HSOSHF(const Ref<KeyVal>&);
    ~HSOSHF();

    void save_data_state(StateOut&);

    void print(std::ostream&o=ExEnv::out0()) const;

    void two_body_energy(double &ec, double &ex);

    int value_implemented() const;
  protected:
    void ao_fock(double accuracy);
    void two_body_deriv(double*);
  private:
    bool analytic_gradient_implemented() const;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
