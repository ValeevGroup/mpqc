//
// wfn.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
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

#ifndef _chemistry_qc_wfn_wfn_h
#define _chemistry_qc_wfn_wfn_h

#ifdef __GNUC__
#pragma interface
#endif

#include <iostream.h>

#include <util/misc/compute.h>
#include <math/scmat/matrix.h>
#include <math/scmat/vector3.h>
#include <chemistry/molecule/energy.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/integral.h>

class Wavefunction: public MolecularEnergy {
#   define CLASSNAME Wavefunction
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  private:
    RefSCDimension basisdim_;
    RefSCMatrixKit basiskit_;
    
    ResultRefSymmSCMatrix overlap_;
    ResultRefSymmSCMatrix hcore_;
    ResultRefSCMatrix natural_orbitals_;
    ResultRefDiagSCMatrix natural_density_;

    double * bs_values;
    double * bsg_values;

    RefGaussianBasisSet gbs_;
    RefIntegral integral_;

    int print_nao_;
    int print_npa_;
    
  public:
    Wavefunction(StateIn&);
    Wavefunction(const RefKeyVal&);
    virtual ~Wavefunction();

    void save_data_state(StateOut&);

    double density(const SCVector3&);
    double density_gradient(const SCVector3&,double*);
    double natural_orbital(const SCVector3& r, int iorb);
    double natural_orbital_density(const SCVector3& r,
                                   int orb, double* orbval = 0);
    double orbital(const SCVector3& r, int iorb, const RefSCMatrix& orbs);

    double orbital_density(const SCVector3& r,
                           int iorb,
                           const RefSCMatrix& orbs,
                           double* orbval = 0);

    virtual RefSymmSCMatrix density() = 0;
    virtual RefSCMatrix natural_orbitals();
    virtual RefDiagSCMatrix natural_density();

    // returns the ao to nao transformation matrix
    virtual RefSCMatrix nao();

    virtual RefSymmSCMatrix overlap();
    virtual RefSymmSCMatrix core_hamiltonian();

    RefSCDimension basis_dimension();
    RefSCMatrixKit basis_matrixkit();
    RefGaussianBasisSet basis();
    RefIntegral integral();

    // returns a matrix which transforms AO's to orthogonal AO's
    // can be overridden, but defaults to S^-1/2
    virtual RefSymmSCMatrix ao_to_orthog_ao();

    void print(ostream& = cout);
};
SavableState_REF_dec(Wavefunction);

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
