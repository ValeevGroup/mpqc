//
// densval.cc
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

#include <util/misc/formio.h>
#include <util/keyval/keyval.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/petite.h>

#include <chemistry/qc/wfn/wfn.h>

using namespace sc;

// Function for returning electron charge density at a point
double Wavefunction::density(const SCVector3&r)
{
  int nbasis = basis()->nbasis();
  if (!bs_values) bs_values=new double[nbasis];

  // compute the basis set values
  GaussianBasisSet::ValueData *valdat
      = new GaussianBasisSet::ValueData(basis(), integral_);
  basis()->values(r,valdat,bs_values);
  delete valdat;

  //for (int i=0; i<nbasis; i++)
  //     ExEnv::out0() << indent
  //          << scprintf("bs_values[%d] = %12.8f\n",i,bs_values[i]);

  // Assuming this will be called many times for the same wavefunction,
  // it is more efficient to force the computation of the natural
  // orbitals now and use them.
  // get the natural orbitals and density
  RefSCMatrix nos
      = integral()->petite_list()->evecs_to_AO_basis(natural_orbitals());
  RefDiagSCMatrix nd = natural_density();
  
  // loop over natural orbitals adding contributions to elec_density
  double elec_density=0.0;
  for (int i=0; i<nbasis; i++) {
      double tmp = 0.0;
      for (int j=0; j<nbasis; j++) {
          tmp += nos.get_element(j,i)*bs_values[j];
        }
      elec_density += nd.get_element(i)*tmp*tmp;
    }

  return elec_density;
}     

// Function for returning electron charge density at a point.
// The grad at that point is also computed and put into double grad[3].
double Wavefunction::density_gradient(const SCVector3&r,double*grad)
{
  int nbasis = basis()->nbasis();
  if (!bs_values) bs_values=new double[nbasis];
  if (!bsg_values) bsg_values=new double[nbasis*3];

  // compute the grad values and get the basis set values at the
  // same time
  GaussianBasisSet::ValueData *valdat
      = new GaussianBasisSet::ValueData(basis(), integral_);
  basis()->grad_values(r,valdat,bsg_values,bs_values);
  delete valdat;

  //for (int i=0; i<nbasis; i++)
  //     ExEnv::out0() << indent
  //          << scprintf("bs_values[%d] = % 12.8f\n",i,bs_values[i]);

  // get the natural orbitals and density
  RefSCMatrix nos
      = integral()->petite_list()->evecs_to_AO_basis(natural_orbitals());
  RefDiagSCMatrix nd = natural_density();
    
  // loop over natural orbitals adding contributions to elec_density
  double elec_density=0.0;
  grad[0] = grad[1] = grad[2] = 0.0;
  for (int i=0; i<nbasis; i++) {
      double tmp = 0.0;
      int j;
      for (j=0; j<nbasis; j++) {
          tmp += nos.get_element(j,i)*bs_values[j];
        }
      elec_density += nd.get_element(i)*tmp*tmp;
      double tmpg[3];
      tmpg[0] = tmpg[1] = tmpg[2] = 0.0;
      for (j=0; j<nbasis; j++) {
          tmpg[0] += nos.get_element(j,i)*bsg_values[j*3+0];
          tmpg[1] += nos.get_element(j,i)*bsg_values[j*3+1];
          tmpg[2] += nos.get_element(j,i)*bsg_values[j*3+2];
        }
      grad[0] += nd.get_element(i)*tmpg[0]*tmp;
      grad[1] += nd.get_element(i)*tmpg[1]*tmp;
      grad[2] += nd.get_element(i)*tmpg[2]*tmp;
    }

  grad[0] *= 2.0;
  grad[1] *= 2.0;
  grad[2] *= 2.0;

  return elec_density;
}     

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
