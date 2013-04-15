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

#include <math/scmat/blocked.h>

#include <chemistry/qc/wfn/femo.h>
#include <chemistry/qc/wfn/obwfn.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/petite.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////////

static ClassDesc HCoreWfn_cd(
  typeid(HCoreWfn),"HCoreWfn",1,"public OneBodyWavefunction",
  0, create<HCoreWfn>, create<HCoreWfn>);

HCoreWfn::HCoreWfn(StateIn& s) :
  SavableState(s),
  OneBodyWavefunction(s)
{
  s.get(nirrep_);
  s.get(docc_);
  s.get(socc_);
  s.get(user_occ_);
  s.get(total_charge_);
}

HCoreWfn::HCoreWfn(const Ref<KeyVal>&keyval):
  OneBodyWavefunction(keyval)
{
  CharacterTable ct = molecule()->point_group()->char_table();

  nirrep_ = ct.ncomp();
  docc_ = new int[nirrep_];
  socc_ = new int[nirrep_];

  user_occ_ = 0;
  total_charge_ = keyval->intvalue("total_charge");

  int nuclear_charge = int(molecule()->nuclear_charge());
  int computed_charge = nuclear_charge;

  for (int i=0; i < nirrep_; i++) {
    docc_[i]=0;
    socc_[i]=0;

    if (keyval->exists("docc",i)) {
      docc_[i] = keyval->intvalue("docc",i);
      computed_charge -= 2;
      user_occ_ = 1;
      }
    if (keyval->exists("socc",i)) {
      socc_[i] = keyval->intvalue("socc",i);
      computed_charge -= 1;
      user_occ_ = 1;
      }
  }
  if (!keyval->exists("total_charge")) {
    if (user_occ_) total_charge_ = computed_charge;
    else total_charge_ = 0;
    }
  else if (total_charge_ != computed_charge && user_occ_) {
    ExEnv::err0() << indent
                 << "ERROR: HCoreWfn: total_charge != computed_charge"
                 << endl;
    abort();
    }
  if (total_charge_ > nuclear_charge) {
    ExEnv::err0() << indent
                 << "ERROR: HCoreWfn: total_charge > nuclear_charge"
                 << endl;
    abort();
  }
}

HCoreWfn::~HCoreWfn()
{
  delete[] docc_;
  delete[] socc_;
}

void
HCoreWfn::save_data_state(StateOut&s)
{
  OneBodyWavefunction::save_data_state(s);
  s.put(nirrep_);
  s.put(docc_,nirrep_);
  s.put(socc_,nirrep_);
  s.put(user_occ_);
  s.put(total_charge_);
}

RefSCMatrix
HCoreWfn::oso_eigenvectors()
{
  if (!oso_eigenvectors_.computed() || !eigenvalues_.computed()) {
    RefSymmSCMatrix hcore_oso(oso_dimension(), basis_matrixkit());
    hcore_oso->assign(0.0);
    hcore_oso->accumulate_transform(so_to_orthog_so(), core_hamiltonian());

    if (debug_ > 1) {
      core_hamiltonian().print("hcore in SO basis");
    }

    if (debug_ > 1) {
      hcore_oso.print("hcore in ortho SO basis");
    }

    RefSCMatrix vec(oso_dimension(), oso_dimension(), basis_matrixkit());
    RefDiagSCMatrix val(oso_dimension(), basis_matrixkit());

    hcore_oso.diagonalize(val,vec);

    if (debug_ > 1) {
      val.print("hcore eigenvalues in ortho SO basis");
      vec.print("hcore eigenvectors in ortho SO basis");
    }
    oso_eigenvectors_=vec;
    oso_eigenvectors_.computed() = 1;

    eigenvalues_ = val;
    eigenvalues_.computed() = 1;

    if (!user_occ_) {
      int nelectron = int(molecule()->nuclear_charge()) - total_charge_;
      fill_occ(val, nelectron, docc_, socc_);

      ExEnv::out0() << indent << "docc = [";
      for (int i=0; i<nirrep_; i++) ExEnv::out0() << " " << docc_[i];
      ExEnv::out0() << "]" << endl;

      ExEnv::out0() << indent << "socc = [";
      for (int i=0; i<nirrep_; i++) ExEnv::out0() << " " << socc_[i];
      ExEnv::out0() << "]" << endl;
      }
  }
  
  return oso_eigenvectors_.result_noupdate();
}

RefDiagSCMatrix
HCoreWfn::eigenvalues()
{
  if (!eigenvalues_.computed()) {
    oso_eigenvectors();
  }
  
  return eigenvalues_.result_noupdate();
}

RefSymmSCMatrix
HCoreWfn::density()
{
  if (!density_.computed()) {
    // Make sure the eigenvectors are computed before
    // the MO density is computed, otherwise, docc and
    // socc might not have been initialized.
    oso_eigenvectors();

    RefDiagSCMatrix mo_density(oso_dimension(), basis_matrixkit());
    BlockedDiagSCMatrix *modens
      = dynamic_cast<BlockedDiagSCMatrix*>(mo_density.pointer());
    if (!modens) {
      ExEnv::err0() << indent
                   << "HCoreWfn::density: wrong MO matrix kit" << endl;
      abort();
    }

    modens->assign(0.0);
    for (int iblock=0; iblock < modens->nblocks(); iblock++) {
      RefDiagSCMatrix modens_ib = modens->block(iblock);
      int i;
      for (i=0; i < docc_[iblock]; i++)
        modens_ib->set_element(i, 2.0);
      for ( ; i < docc_[iblock]+socc_[iblock]; i++)
        modens_ib->set_element(i, 1.0);
    }

    RefSymmSCMatrix dens(so_dimension(), basis_matrixkit());
    dens->assign(0.0);
    dens->accumulate_transform(so_to_orthog_so().t() * mo_to_orthog_so(),
                               mo_density);

    if (debug_ > 1) {
      mo_density.print("MO Density");
      dens.print("SO Density");
      ExEnv::out0()
                   << indent << "Nelectron(MO) = " << mo_density.trace()
                   << endl
                   << indent << "Nelectron(SO) = "
                   << (overlap()*dens).trace()
                   << endl;
    }

    density_ = dens;
    density_.computed() = 1;
  }
      
  return density_.result_noupdate();   
}

double
HCoreWfn::occupation(int ir, int i)
{
  if (i < docc_[ir])
    return 2.0;
  else if (i < docc_[ir]+socc_[ir])
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
  double e = (density()*core_hamiltonian()).trace();
  set_energy(e);
  set_actual_value_accuracy(desired_value_accuracy());
  return;
}

int
HCoreWfn::value_implemented() const
{
  return 1;
}

void
HCoreWfn::fill_occ(const RefDiagSCMatrix &evals,int nelectron,int *docc, int *socc)
{
  HundsFEMOSeeker femoseeker(nelectron, HundsFEMOSeeker::tolerance, true, evals);
  Ref<FEMO> femo = femoseeker.result();
  for(int g=0; g<nirrep_; ++g) {
    const int na = femo->nalpha(g);
    const int nb = femo->nbeta(g);
    docc[g] = nb;
    socc[g] = na - nb;
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
