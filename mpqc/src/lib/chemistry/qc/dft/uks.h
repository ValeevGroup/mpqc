//
// uks.h --- definition of the unrestricted Kohn-Sham class
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

#ifndef _chemistry_qc_scf_uks_h
#define _chemistry_qc_scf_uks_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/scf/uscf.h>
#include <chemistry/qc/dft/integrator.h>
#include <chemistry/qc/dft/functional.h>

namespace sc {

// //////////////////////////////////////////////////////////////////////////

/**
   This provides a Kohn-Sham implementation for unrestricted-orbital
   open-shell systems.
 */
class UKS: public UnrestrictedSCF {
  protected:
    Ref<DenIntegrator> integrator_;
    Ref<DenFunctional> functional_;
    RefSymmSCMatrix vaxc_;
    RefSymmSCMatrix vbxc_;

  public:
    UKS(StateIn&);
    /** This KeyVal constructor is used to construct UKS
        objects from the input.

        The keywords used by this constructor are listed below.  The KeyVal
        constructor for the parent class, UnrestrictedSCF, will also be
        called, so consult the documentation for
        UnrestrictedSCF(const Ref<KeyVal>&) for additional keywords
        that will be read.

        <dl>

       <dt><tt>integrator</tt><dd>Specifies the DenIntegrator that will be
       used to integrate the density functional.  The default is
       RadialAngularIntegrator.

       <dt><tt>functional</tt><dd>Specifies the DenFunctional that will be
       used to compute the exchange/correlation contribution.  This is no
       default.

       </dl>
    */
    UKS(const Ref<KeyVal>&);
    ~UKS();

    void save_data_state(StateOut&);

    void print(std::ostream&o=ExEnv::out0()) const;

    void two_body_energy(double &ec, double &ex);

    int value_implemented() const;
    int gradient_implemented() const;

  protected:
    double exc_;
    
    void ao_fock(double accuracy);
    double scf_energy();
    Ref<SCExtrapData> extrap_data();
    void two_body_deriv(double*);

    void init_vector();
    void done_vector();
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
