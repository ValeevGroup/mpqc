//
// eht.cc
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

#include <chemistry/qc/wfn/eht.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/petite.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////////

static ClassDesc ExtendedHuckelWfn_cd(
  typeid(ExtendedHuckelWfn),"ExtendedHuckelWfn",1,"public OneBodyWavefunction",
  0, create<ExtendedHuckelWfn>, create<ExtendedHuckelWfn>);

ExtendedHuckelWfn::ExtendedHuckelWfn(StateIn& s) :
  SavableState(s),
  OneBodyWavefunction(s)
{
  s.get(nirrep_);
  s.get(docc_);
  s.get(socc_);
  s.get(user_occ_);
  s.get(total_charge_);
}

ExtendedHuckelWfn::ExtendedHuckelWfn(const Ref<KeyVal>&keyval):
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
      computed_charge -= 2*docc_[i];
      user_occ_ = 1;
      }
    if (keyval->exists("socc",i)) {
      socc_[i] = keyval->intvalue("socc",i);
      computed_charge -= 1*socc_[i];
      user_occ_ = 1;
      }
  }
  if (!keyval->exists("total_charge")) {
    if (user_occ_) total_charge_ = computed_charge;
    else total_charge_ = 0;
    }
  else if (total_charge_ != computed_charge && user_occ_) {
    ExEnv::err0() << indent
                 << "ERROR: ExtendedHuckelWfn: total_charge != computed_charge"
                 << endl;
    abort();
    }
  if (total_charge_ > nuclear_charge) {
    ExEnv::err0() << indent
                 << "ERROR: ExtendedHuckelWfn: total_charge > nuclear_charge"
                 << endl;
    abort();
  }
}

ExtendedHuckelWfn::~ExtendedHuckelWfn()
{
  delete[] docc_;
  delete[] socc_;
}

void
ExtendedHuckelWfn::save_data_state(StateOut&s)
{
  OneBodyWavefunction::save_data_state(s);
  s.put(nirrep_);
  s.put(docc_,nirrep_);
  s.put(socc_,nirrep_);
  s.put(user_occ_);
  s.put(total_charge_);
}

RefSymmSCMatrix
ExtendedHuckelWfn::h_eht_oso()
{
  Ref<PetiteList> pl = integral()->petite_list();

  // Compute H in the AO basis
  double K = 1.75;
  Ref<AtomInfo> atominfo = molecule()->atominfo();
  RefSymmSCMatrix h_ao = pl->to_AO_basis(overlap());
  int natom = basis()->ncenter();
  int funcoff1 = 0;
  for (int atom1=0; atom1<natom; atom1++) {
      int nbasis1 = basis()->nbasis_on_center(atom1);
      double I1 = atominfo->ip(molecule()->Z(atom1));
      int funcoff2 = 0;
      for (int atom2=0; atom2<=atom1; atom2++) {
          int nbasis2 = basis()->nbasis_on_center(atom2);
          double I2 = atominfo->ip(molecule()->Z(atom2));
          for (int func1=0; func1<nbasis1; func1++) {
              if (atom1 == atom2) nbasis2 = func1 + 1;
              for (int func2=0; func2<nbasis2; func2++) {
                int i1 = funcoff1+func1;
                int i2 = funcoff2+func2;
                  double val = h_ao(i1,i2);
                if (atom1 == atom2 && func1 == func2) {
                  // The overlap integral is not a part of the diagonal
                  // element in standard EHT formulae.  It is here though,
                  // since basis functions are not always normalized (some
                  // d shell components for example).
                  val *= -I1;
                }
                else {
                  val *= -0.5*K*(I1+I2);
                }
                h_ao(i1,i2) = val;
              }
          }
          funcoff2 += nbasis2;
      }
      funcoff1 += nbasis1;
  }

  if (debug_ > 1) h_ao.print("h in the AO basis");

  // Compute H in the SO basis
  RefSymmSCMatrix h_so = pl->to_SO_basis(h_ao);

  if (debug_ > 1) {
    pl->to_AO_basis(overlap()).print("S in AO basis");
    overlap().print("S in SO basis");
    pl->aotoso().print("AO to SO transform");
    h_so.print("h in the SO basis");
  }

  // Compute H in the OSO basis
  RefSymmSCMatrix h_oso(oso_dimension(), basis_matrixkit());
  h_oso->assign(0.0);
  h_oso->accumulate_transform(so_to_orthog_so(),h_so);

  return h_oso;
}

RefSCMatrix
ExtendedHuckelWfn::oso_eigenvectors()
{
  if (!oso_eigenvectors_.computed() || !eigenvalues_.computed()) {
    RefSymmSCMatrix h_oso = h_eht_oso();

    if (debug_ > 1) {
      h_oso.print("h in ortho SO basis");
    }

    RefSCMatrix vec(oso_dimension(), oso_dimension(), basis_matrixkit());
    RefDiagSCMatrix val(oso_dimension(), basis_matrixkit());

    h_oso.diagonalize(val,vec);

    if (debug_ > 1) {
      val.print("h eigenvalues in ortho SO basis");
      vec.print("h eigenvectors in ortho SO basis");
    }
    oso_eigenvectors_=vec;
    oso_eigenvectors_.computed() = 1;

    eigenvalues_ = val;
    eigenvalues_.computed() = 1;

    if (!user_occ_) {
      int nelectron = int(molecule()->nuclear_charge()) - total_charge_;
      int ndocc = nelectron/2;
      int nsocc = nelectron%2;
      fill_occ(val, ndocc, docc_, nsocc, socc_);

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
ExtendedHuckelWfn::eigenvalues()
{
  if (!eigenvalues_.computed()) {
    oso_eigenvectors();
  }
  
  return eigenvalues_.result_noupdate();
}

RefSymmSCMatrix
ExtendedHuckelWfn::density()
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
                   << "ExtendedHuckelWfn::density: wrong MO matrix kit" << endl;
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

    if (debug_ > 1) mo_density.print("MO Density");

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
ExtendedHuckelWfn::occupation(int ir, int i)
{
  if (!eigenvalues_.computed()) {
    oso_eigenvectors();
  }
  
  if (i < docc_[ir])
    return 2.0;
  else if (i < docc_[ir]+socc_[ir])
    return 1.0;
  else
    return 0.0;
}

int
ExtendedHuckelWfn::spin_polarized()
{
  return 0;
}

int
ExtendedHuckelWfn::spin_unrestricted()
{
  return 0;
}

void
ExtendedHuckelWfn::compute()
{
  double e = (density()*core_hamiltonian()).trace();
  set_energy(e);
  set_actual_value_accuracy(desired_value_accuracy());
  return;
}

int
ExtendedHuckelWfn::value_implemented() const
{
  return 1;
}

void
ExtendedHuckelWfn::fill_occ(const RefDiagSCMatrix &evals,int ndocc,int *docc,
                   int nsocc, int *socc)
{
  BlockedDiagSCMatrix *bval
    = require_dynamic_cast<BlockedDiagSCMatrix*>(evals.pointer(),
                                           "ExtendedHuckelWfn: getting occupations");
  int nblock = bval->nblocks();
  if (nblock != nirrep_) {
    ExEnv::errn() << "ERROR: ExtendedHuckelWfn: fill_occ: nblock != nirrep" << endl
                 << "  nblock = " << nblock << endl
                 << "  nirrep = " << nirrep_ << endl;
    abort();
  }
  memset(docc,0,sizeof(docc[0])*nblock);
  memset(socc,0,sizeof(socc[0])*nblock);
  for (int i=0; i<ndocc; i++) {
    double lowest;
    int lowest_j = -1;
    for (int j=0; j<nblock; j++) {
      RefDiagSCMatrix block = bval->block(j);
      if (block.null()) continue;
      double current = block->get_element(docc[j]);
      if (lowest_j < 0 || lowest > current) {
        lowest = current;
        lowest_j = j;
      }
    }
    docc[lowest_j]++;
  }
  for (int i=0; i<nsocc; i++) {
    double lowest;
    int lowest_j = -1;
    for (int j=0; j<nblock; j++) {
      double current = bval->block(j)->get_element(docc[j]+socc[j]);
      if (lowest_j < 0 || lowest > current) {
        lowest = current;
        lowest_j = j;
      }
    }
    socc[lowest_j]++;
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
