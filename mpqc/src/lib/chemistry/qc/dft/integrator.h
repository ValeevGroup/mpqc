//
// integrator.h --- definition of the electron density integrator
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

#ifndef _chemistry_qc_dft_integrator_h
#define _chemistry_qc_dft_integrator_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/state/state.h>
#include <chemistry/qc/dft/functional.h>

class DenIntegrator: virtual_base public SavableState {
#   define CLASSNAME DenIntegrator
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    RefWavefunction wfn_;

    double value_;

    int spin_polarized_;

    double *bs_values_;
    double *bsg_values_;
    double *alpha_dmat_;
    double *beta_dmat_;
    double *alpha_vmat_; // lower triangle of xi_i(r) v(r) xi_j(r) integrals
    double *beta_vmat_; // lower triangle of xi_i(r) v(r) xi_j(r) integrals
    int need_density_; // specialization must set to 1 if it needs density_
    double density_;
    int nbasis_;
    int compute_potential_integrals_; // 1 if potential integrals are needed

    int need_gradient_;

    void get_density(double *dmat, double &den, double grad[3]);
    void init_integration(const RefDenFunctional &func,
                          const RefSymmSCMatrix& densa,
                          const RefSymmSCMatrix& densb);
    double do_point(const SCVector3 &r, const RefDenFunctional &,
                    double weight);
  public:
    DenIntegrator();
    DenIntegrator(const RefKeyVal &);
    DenIntegrator(StateIn &);
    ~DenIntegrator();
    void save_data_state(StateOut &);

    virtual void set_wavefunction(const RefWavefunction &);
    RefWavefunction wavefunction() const { return wfn_; }
    double value() const { return value_; }

    // Call with non zero if the potential integrals are to be computed.
    // They can be returned with the vmat() member.
    void set_compute_potential_integrals(int);
    const double *alpha_vmat() const { return alpha_vmat_; }
    const double *beta_vmat() const { return beta_vmat_; }

    virtual void integrate(const RefDenFunctional &,
                           const RefSymmSCMatrix& densa =0,
                           const RefSymmSCMatrix& densb =0) = 0;
};
SavableState_REF_dec(DenIntegrator);

// Based on C.W. Murray, et al. Mol. Phys. 78, No. 4, 997-1014, 1993
class Murray93Integrator: public DenIntegrator {
#   define CLASSNAME Murray93Integrator
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    int nr_;
    int ntheta_;
    int nphi_;
    int Ktheta_;
    
  public:
    Murray93Integrator();
    Murray93Integrator(const RefKeyVal &);
    Murray93Integrator(StateIn &);
    ~Murray93Integrator();
    void save_data_state(StateOut &);

    void integrate(const RefDenFunctional &,
                   const RefSymmSCMatrix& densa =0,
                   const RefSymmSCMatrix& densb =0);

    void print(ostream & =cout) const;
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
