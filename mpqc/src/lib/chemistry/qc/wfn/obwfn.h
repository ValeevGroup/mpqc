//
// obwfn.h
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

#ifndef _chemistry_qc_wfn_obwfn_h
#define _chemistry_qc_wfn_obwfn_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/wfn/wfn.h>
#include <chemistry/qc/wfn/hcore.h>

SavableState_REF_fwddec(OneBodyWavefunction);
class OneBodyWavefunction: public Wavefunction {
#   define CLASSNAME OneBodyWavefunction
#   include <util/state/stated.h>
#   include <util/class/classda.h>
 protected:
    ResultRefSymmSCMatrix density_;
    AccResultRefSCMatrix eigenvectors_;
    int nirrep_;
    int *nvecperirrep_;

    void init_sym_info();
 public:
    OneBodyWavefunction(StateIn&);
    OneBodyWavefunction(const RefKeyVal&);
    ~OneBodyWavefunction();

    void save_data_state(StateOut&);

    virtual RefSCMatrix eigenvectors() = 0;
    virtual double occupation(int irrep, int vectornum) = 0;
    double occupation(int vectornum);

    virtual RefSCMatrix projected_eigenvectors(const RefOneBodyWavefunction&);
    virtual RefSCMatrix hcore_guess();

    double orbital(const SCVector3& r, int iorb);
    double orbital_density(const SCVector3& r, int iorb, double* orbval = 0);

    double density(const SCVector3&);
    RefSymmSCMatrix density();

    void print(ostream&o=cout);
};
SavableState_REF_dec(OneBodyWavefunction);

// This is useful as an initial guess for other one body wavefunctions
class HCoreWfn: public OneBodyWavefunction {
#   define CLASSNAME HCoreWfn
#   define HAVE_STATEIN_CTOR
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  private:
    RefAccumHCore accumh;

    int nirrep_;
    int *docc;
    int *socc;
    
    void compute();

  public:
    HCoreWfn(StateIn&);
    HCoreWfn(const RefKeyVal&);
    ~HCoreWfn();

    void save_data_state(StateOut&);

    double occupation(int irrep, int vectornum);

    RefSCMatrix eigenvectors();
};
SavableState_REF_dec(HCoreWfn);

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")
// End:
