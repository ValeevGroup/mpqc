//
// uhf.h --- definition of the unrestricted Hartree-Fock class
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

#ifndef _chemistry_qc_scf_uhf_h
#define _chemistry_qc_scf_uhf_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/scf/uscf.h>

////////////////////////////////////////////////////////////////////////////

class UHF: public UnrestrictedSCF {
#   define CLASSNAME UHF
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    UHF(StateIn&);
    UHF(const RefKeyVal&);
    ~UHF();

    void save_data_state(StateOut&);

    void print(ostream&o=cout);

    void two_body_energy(double &ec, double &ex);

    int value_implemented();
    int gradient_implemented();

  protected:
    void ao_fock();
    void two_body_deriv(double*);
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
