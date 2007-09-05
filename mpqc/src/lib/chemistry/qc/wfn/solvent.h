//
// solvent.h
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

#ifndef _chemistry_qc_wfn_solvent_h
#define _chemistry_qc_wfn_solvent_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/solvent/bem.h>
#include <chemistry/qc/wfn/wfn.h>
#include <chemistry/qc/wfn/accum.h>

namespace sc {

/** <em>WARNING: The BEMSolventH class is not thoroughly tested.  It should
    only be used by developers wishing to test or fix it.</em>

    This specialization of AccumH computes the
    contribution to the energy one body Hamiltonian from
    a solvent using a polarizable continuum model.
*/
class BEMSolventH: public AccumH {
  private:
    double gamma_;
    int onebody_;
    int normalize_q_;
    int separate_surf_charges_;
    int y_equals_j_;
    int integrate_nelectron_;

    Ref<Wavefunction> wfn_;
    Ref<BEMSolvent> solvent_;

    double **charge_positions_;
    double **normals_;
    double *efield_dot_normals_;
    double *charges_;
    double *charges_n_;
    double enucsurf_;
    double eelecsurf_;
    double esurfsurf_;
    double escalar_;
    double ecavitation_;
    double edisprep_;

  public:
    BEMSolventH(StateIn&);
    BEMSolventH(const Ref<KeyVal>&);
    virtual ~BEMSolventH();

    void save_data_state(StateOut&);

    void init(const Ref<Wavefunction>&);
    void accum(const RefSymmSCMatrix& h);
    void done();
    void print_summary();

    double e();
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
