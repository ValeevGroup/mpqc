//
// clscf.h --- definition of the closed shell SCF class
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

#ifndef _chemistry_qc_scf_clscf_h
#define _chemistry_qc_scf_clscf_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/scf/scf.h>

namespace sc {

// //////////////////////////////////////////////////////////////////////////

/** The CLSCF class is a base for classes implementing a self-consistent
procedure for closed-shell molecules. */
class CLSCF: public SCF {
  protected:
    Ref<PointGroup> most_recent_pg_;
    int user_occupations_;
    int tndocc_;
    int nirrep_;
    int *initial_ndocc_;
    int *ndocc_;

    ResultRefSymmSCMatrix cl_fock_;

  public:
    CLSCF(StateIn&);
    /** The KeyVal constructor:
        <dl>

        <dt><tt>total_charge</tt><dd> This floating point number
        gives the total charge
        of the molecule.  The default is 0.

        <dt><tt>docc</tt><dd> This vector of integers gives the total
        number of doubly occupied orbitals of each irreducible
        representation.  By default, this will be chosen to make the
        molecule uncharged and the electrons will be distributed among the
        irreducible representations according to the orbital energies.

        </dl> */
    CLSCF(const Ref<KeyVal>&);
    ~CLSCF();

    void save_data_state(StateOut&);

    void print(std::ostream&o=ExEnv::out0()) const;

    double occupation(int irrep, int vectornum);

    int n_fock_matrices() const;
    /** Returns closed-shell Fock matrix in AO basis (excluding XC contribution in KS DFT).
	Use effective_fock() if you want the full KS Fock matrix.
	Argument i must be 0.
    */
    RefSymmSCMatrix fock(int i);
    /** Returns closed-shell Fock matrix in MO basis (including XC contribution for KS DFT). */
    RefSymmSCMatrix effective_fock();

    RefSymmSCMatrix density();

    void symmetry_changed();
    
    // returns 0
    int spin_polarized();

  protected:
    // these are temporary data, so they should not be checkpointed
    RefSymmSCMatrix cl_dens_;
    RefSymmSCMatrix cl_dens_diff_;
    RefSymmSCMatrix cl_gmat_;

    void set_occupations(const RefDiagSCMatrix& evals);

    // scf things
    void init_vector();
    void done_vector();
    void reset_density();
    double new_density();
    double scf_energy();

    Ref<SCExtrapData> initial_extrap_data();
    Ref<SCExtrapData> extrap_data();
    
    // gradient things
    void init_gradient();
    void done_gradient();

    RefSymmSCMatrix lagrangian();
    RefSymmSCMatrix gradient_density();

    // hessian things
    void init_hessian();
    void done_hessian();

    // The Hartree-Fock derivatives
    void two_body_deriv_hf(double*grad,double exchange_fraction);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
