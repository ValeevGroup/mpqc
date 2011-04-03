//
// ossscf.h --- definition of the open shell singlet SCF class
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

#ifndef _chemistry_qc_scf_ossscf_h
#define _chemistry_qc_scf_ossscf_h

#include <chemistry/qc/scf/scf.h>

namespace sc {

// //////////////////////////////////////////////////////////////////////////

class OSSSCF: public SCF {
 protected:
    int user_occupations_;
    int tndocc_;
    int nirrep_;
    int *ndocc_;
    int osa_;
    int osb_;

    ResultRefSymmSCMatrix cl_fock_;
    ResultRefSymmSCMatrix op_focka_;
    ResultRefSymmSCMatrix op_fockb_;

  public:
    OSSSCF(StateIn&);
    OSSSCF(const Ref<KeyVal>&);
    ~OSSSCF();

    void save_data_state(StateOut&);

    void print(std::ostream&o=ExEnv::out0()) const;

    double occupation(int ir, int vectornum);
    double alpha_occupation(int irrep, int vectornum);
    double beta_occupation(int irrep, int vectornum);

    int n_fock_matrices() const;
    RefSymmSCMatrix fock(int);
    RefSymmSCMatrix effective_fock();
    RefSymmSCMatrix density();
    RefSymmSCMatrix alpha_density();
    RefSymmSCMatrix beta_density();

    void symmetry_changed();
    
    int spin_polarized();

  protected:
    // these are temporary data, so they should not be checkpointed
    RefSymmSCMatrix cl_dens_;
    RefSymmSCMatrix cl_dens_diff_;
    RefSymmSCMatrix cl_gmat_;
    RefSymmSCMatrix op_densa_;
    RefSymmSCMatrix op_densa_diff_;
    RefSymmSCMatrix op_gmata_;
    RefSymmSCMatrix op_densb_;
    RefSymmSCMatrix op_densb_diff_;
    RefSymmSCMatrix op_gmatb_;

    RefSymmSCMatrix cl_hcore_;
    
    void set_occupations(const RefDiagSCMatrix& evals);

    // scf things
    void init_vector();
    void done_vector();
    void reset_density();
    double new_density();
    double scf_energy();

    Ref<SCExtrapData> extrap_data();
    
    // gradient things
    void init_gradient();
    void done_gradient();

    RefSymmSCMatrix lagrangian();
    RefSymmSCMatrix gradient_density();

    // hessian things
    void init_hessian();
    void done_hessian();
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
