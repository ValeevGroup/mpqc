//
// clks.h --- definition of the closed shell Kohn-Sham SCF class
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

#ifndef _chemistry_qc_dft_clks_h
#define _chemistry_qc_dft_clks_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/scf/clscf.h>
#include <chemistry/qc/dft/integrator.h>
#include <chemistry/qc/dft/functional.h>

namespace sc {

// //////////////////////////////////////////////////////////////////////////

/**
   This provides a Kohn-Sham implementation for closed-shell systems.
 */
class CLKS: public CLSCF {
  protected:
    Ref<DenIntegrator> integrator_;
    Ref<DenFunctional> functional_;
    RefSymmSCMatrix vxc_;
    
  public:
    CLKS(StateIn&);
    /**
       This KeyVal constructor reads the following keywords:
        <dl>

       <dt><tt>integrator</tt><dd>Specifies the DenIntegrator that will be
       used to integrate the density functional.  The default is
       RadialAngularIntegrator.

       <dt><tt>functional</tt><dd>Specifies the DenFunctional that will be
       used to compute the exchange/correlation contribution.  This is no
       default.

       </dl>
    */
    CLKS(const Ref<KeyVal>&);
    ~CLKS();

    void save_data_state(StateOut&);

    void print(std::ostream&o=ExEnv::out0()) const;

    void two_body_energy(double &ec, double &ex);

    int value_implemented() const;
    int gradient_implemented() const;

    RefSymmSCMatrix density();
  protected:
    void ao_fock(double accuracy);
    double exc_;
    double scf_energy();
    Ref<SCExtrapData> initial_extrap_data();
    Ref<SCExtrapData> extrap_data();
    RefSymmSCMatrix effective_fock();

    void init_vector();
    void done_vector();
    
    void two_body_deriv(double*);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
