//
// hcorewfn.cc
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
#include <util/state/stateio.h>

#include <chemistry/qc/wfn/obwfn.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/petite.h>

/////////////////////////////////////////////////////////////////////////

SavableState_REF_def(HCoreWfn);

#define CLASSNAME HCoreWfn
#define PARENTS public OneBodyWavefunction
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>

void *
HCoreWfn::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = OneBodyWavefunction::_castdown(cd);
  return do_castdowns(casts,cd);
}

HCoreWfn::HCoreWfn(StateIn& s) :
  maybe_SavableState(s)
  OneBodyWavefunction(s)
{
  s.get(nirrep_);
  s.get(docc);
  s.get(socc);
}

HCoreWfn::HCoreWfn(const RefKeyVal&keyval):
  OneBodyWavefunction(keyval)
{
  CharacterTable ct = molecule()->point_group()->char_table();

  nirrep_ = ct.ncomp();
  docc = new int[nirrep_];
  socc = new int[nirrep_];

  for (int i=0; i < nirrep_; i++) {
    docc[i]=0;
    socc[i]=0;

    if (keyval->exists("docc",i))
      docc[i] = keyval->intvalue("docc",i);
    if (keyval->exists("socc",i))
      socc[i] = keyval->intvalue("socc",i);
  }
}

HCoreWfn::~HCoreWfn()
{
  if (docc) {
    delete[] docc;
    docc=0;
  }
  if (socc) {
    delete[] socc;
    socc=0;
  }
}

void
HCoreWfn::save_data_state(StateOut&s)
{
  OneBodyWavefunction::save_data_state(s);
  s.put(nirrep_);
  s.put(docc,nirrep_);
  s.put(socc,nirrep_);
}

RefSCMatrix
HCoreWfn::eigenvectors()
{
  if (!eigenvectors_.computed()) {
    eigenvectors_=hcore_guess();
    eigenvectors_.computed() = 1;
  }
  
  return eigenvectors_.result_noupdate();
}

RefDiagSCMatrix
HCoreWfn::eigenvalues()
{
  if (!eigenvalues_.computed()) {
    eigenvalues_=core_hamiltonian().eigvals();
    eigenvalues_.computed() = 1;
  }
  
  return eigenvalues_.result_noupdate();
}

RefSymmSCMatrix
HCoreWfn::density()
{
  if (!density_.computed()) {
    RefDiagSCMatrix mo_density(oso_dimension(), basis_matrixkit());
    BlockedDiagSCMatrix *modens
      = BlockedDiagSCMatrix::castdown(mo_density.pointer());
    if (!modens) {
      cerr << node0 << indent << "HCoreWfn::density: basis weirdness" << endl;
      abort();
    }

    modens->assign(0.0);
    for (int iblock=0; iblock < modens->nblocks(); iblock++) {
      RefDiagSCMatrix modens_ib = modens->block(iblock);
      int i;
      for (i=0; i < docc[iblock]; i++)
        modens_ib->set_element(i, 2.0);
      for ( ; i < docc[iblock]+socc[iblock]; i++)
        modens_ib->set_element(i, 1.0);
    }

    RefSymmSCMatrix dens(so_dimension(), basis_matrixkit());
    dens->assign(0.0);
    dens->accumulate_transform(eigenvectors(), mo_density);

    density_ = dens;
    density_.computed() = 1;
  }
      
  return density_.result_noupdate();   
}

double
HCoreWfn::occupation(int ir, int i)
{
  if (i < docc[ir])
    return 2.0;
  else if (i < docc[ir]+socc[ir])
    return 1.0;
  else
    return 0.0;
}

int
HCoreWfn::spin_polarized()
{
  return 0;
}

int
HCoreWfn::spin_unrestricted()
{
  return 0;
}

void
HCoreWfn::compute()
{
  eigenvectors();
  set_energy(0.0);
  set_actual_value_accuracy(desired_value_accuracy());
  return;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
