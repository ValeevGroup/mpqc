//
// wfn.cc
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <stdlib.h>
#include <math.h>

#include <iostream>

#include <util/keyval/keyval.h>
#include <util/misc/timer.h>
#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/symmint.h>
#include <chemistry/qc/intv3/intv3.h>

#include <chemistry/qc/wfn/wfn.h>

using namespace std;
using namespace sc;

#define CHECK_SYMMETRIZED_INTEGRALS 0

static ClassDesc Wavefunction_cd(
  typeid(Wavefunction),"Wavefunction",5,"public MolecularEnergy",
  0, 0, 0);

Wavefunction::Wavefunction(const Ref<KeyVal>&keyval):
  // this will let molecule be retrieved from basis
  // MolecularEnergy(new AggregateKeyVal(keyval,
  //                                     new PrefixKeyVal("basis", keyval))),
  MolecularEnergy(keyval),
  overlap_(this),
  hcore_(this),
  natural_orbitals_(this),
  natural_density_(this),
  bs_values(0),
  bsg_values(0)
{
  overlap_.compute() = 0;
  hcore_.compute() = 0;
  natural_orbitals_.compute() = 0;
  natural_density_.compute() = 0;

  overlap_.computed() = 0;
  hcore_.computed() = 0;
  natural_orbitals_.computed() = 0;
  natural_density_.computed() = 0;

  print_nao_ = keyval->booleanvalue("print_nao");
  print_npa_ = keyval->booleanvalue("print_npa");
  KeyValValuedouble lindep_tol_def(1.e-8);
  lindep_tol_ = keyval->doublevalue("lindep_tol", lindep_tol_def);
  if (keyval->exists("symm_orthog")) {
      ExEnv::out0() << indent
                   << "WARNING: using obsolete \"symm_orthog\" keyword"
                   << endl;
      if (keyval->booleanvalue("symm_orthog")) {
        orthog_method_ = OverlapOrthog::Symmetric;
      }
      else {
        orthog_method_ = OverlapOrthog::Canonical;
      }
  }
  else {
    char *orthog_name = keyval->pcharvalue("orthog_method");
    if (!orthog_name) {
      orthog_method_ = OverlapOrthog::Symmetric;
    }
    else if (::strcmp(orthog_name, "canonical") == 0) {
      orthog_method_ = OverlapOrthog::Canonical;
    }
    else if (::strcmp(orthog_name, "symmetric") == 0) {
      orthog_method_ = OverlapOrthog::Symmetric;
    }
    else if (::strcmp(orthog_name, "gramschmidt") == 0) {
      orthog_method_ = OverlapOrthog::GramSchmidt;
    }
    else {
      ExEnv::errn() << "ERROR: bad orthog_method: \""
                   << orthog_name << "\"" << endl;
      abort();
    }
    delete[] orthog_name;
  }

  debug_ = keyval->intvalue("debug");

  gbs_ = require_dynamic_cast<GaussianBasisSet*>(
    keyval->describedclassvalue("basis").pointer(),
    "Wavefunction::Wavefunction\n"
    );

  integral_ << keyval->describedclassvalue("integrals");
  if (integral_.null()) {
    Integral* default_intf = Integral::get_default_integral();
    integral_ = default_intf->clone();
  }
  
  integral_->set_basis(gbs_);
  Ref<PetiteList> pl = integral_->petite_list();

  sodim_ = pl->SO_basisdim();
  aodim_ = pl->AO_basisdim();
  basiskit_ = gbs_->so_matrixkit();
}

Wavefunction::Wavefunction(StateIn&s):
  SavableState(s),
  MolecularEnergy(s),
  overlap_(this),
  hcore_(this),
  natural_orbitals_(this),
  natural_density_(this),
  bs_values(0),
  bsg_values(0)
{
  debug_ = 0;

  overlap_.compute() = 0;
  hcore_.compute() = 0;
  natural_orbitals_.compute() = 0;
  natural_density_.compute() = 0;

  overlap_.computed() = 0;
  hcore_.computed() = 0;
  natural_orbitals_.computed() = 0;
  natural_density_.computed() = 0;

  if (s.version(::class_desc<Wavefunction>()) >= 2) {
    s.get(print_nao_);
    s.get(print_npa_);
  }
  else {
    print_nao_ = 0;
    print_npa_ = 0;
  }

  if (s.version(::class_desc<Wavefunction>()) >= 5) {
    int orthog_enum;
    s.get(orthog_enum);
    orthog_method_ = (OverlapOrthog::OrthogMethod) orthog_enum;
  }
  else if (s.version(::class_desc<Wavefunction>()) >= 3) {
    int symm_orthog;
    s.get(symm_orthog);
    if (symm_orthog) orthog_method_ = OverlapOrthog::Symmetric;
    else orthog_method_ = OverlapOrthog::Canonical;
  }
  else {
    orthog_method_ = OverlapOrthog::Symmetric;
  }

  if (s.version(::class_desc<Wavefunction>()) >= 4) {
    s.get(lindep_tol_);
  }
  else {
    lindep_tol_ = 1.e-6;
  }

  gbs_ << SavableState::restore_state(s);
  integral_ << SavableState::restore_state(s);

  integral_->set_basis(gbs_);
  Ref<PetiteList> pl = integral_->petite_list();

  sodim_ = pl->SO_basisdim();
  aodim_ = pl->AO_basisdim();
  basiskit_ = gbs_->so_matrixkit();
}

void
Wavefunction::symmetry_changed()
{
  MolecularEnergy::symmetry_changed();

  Ref<PetiteList> pl = integral_->petite_list();
  sodim_ = pl->SO_basisdim();
  aodim_ = pl->AO_basisdim();
  overlap_.result_noupdate() = 0;
  basiskit_ = gbs_->so_matrixkit();

  orthog_ = 0;
}

Wavefunction::~Wavefunction()
{
  if (bs_values) {
    delete[] bs_values;
    bs_values=0;
  }
  if (bsg_values) {
    delete[] bsg_values;
    bsg_values=0;
  }
}

void
Wavefunction::save_data_state(StateOut&s)
{
  MolecularEnergy::save_data_state(s);

  // overlap and hcore integrals are cheap so don't store them.
  // same goes for natural orbitals

  s.put(print_nao_);
  s.put(print_npa_);
  s.put((int)orthog_method_);
  s.put(lindep_tol_);

  SavableState::save_state(gbs_.pointer(),s);
  SavableState::save_state(integral_.pointer(),s);
}

double
Wavefunction::charge()
{
  return molecule()->nuclear_charge() - nelectron();
}

RefSymmSCMatrix
Wavefunction::ao_density()
{
  return integral()->petite_list()->to_AO_basis(density());
}

RefSCMatrix
Wavefunction::natural_orbitals()
{
  if (!natural_orbitals_.computed()) {
      RefSymmSCMatrix dens = density();

      // transform the density into an orthogonal basis
      RefSCMatrix ortho = so_to_orthog_so();

      RefSymmSCMatrix densortho(oso_dimension(), basis_matrixkit());
      densortho.assign(0.0);
      densortho.accumulate_transform(so_to_orthog_so(),dens);

      RefSCMatrix natorb(oso_dimension(), oso_dimension(),
                         basis_matrixkit());
      RefDiagSCMatrix natden(oso_dimension(), basis_matrixkit());

      densortho.diagonalize(natden, natorb);

      // natorb is the ortho SO to NO basis transform
      // make natural_orbitals_ the SO to the NO basis transform
      natural_orbitals_ = so_to_orthog_so().t() * natorb;
      natural_density_ = natden;

      natural_orbitals_.computed() = 1;
      natural_density_.computed() = 1;
    }

  return natural_orbitals_.result_noupdate();
}

RefDiagSCMatrix
Wavefunction::natural_density()
{
  if (!natural_density_.computed()) {
      natural_orbitals();
    }

  return natural_density_.result_noupdate();
}

RefSymmSCMatrix
Wavefunction::overlap()
{
  if (!overlap_.computed()) {
    integral()->set_basis(gbs_);
    Ref<PetiteList> pl = integral()->petite_list();

#if ! CHECK_SYMMETRIZED_INTEGRALS
    // first form skeleton s matrix
    RefSymmSCMatrix s(basis()->basisdim(), basis()->matrixkit());
    Ref<SCElementOp> ov =
      new OneBodyIntOp(new SymmOneBodyIntIter(integral()->overlap(), pl));

    s.assign(0.0);
    s.element_op(ov);
    ov=0;

    if (debug_ > 1) s.print("AO skeleton overlap");

    // then symmetrize it
    RefSymmSCMatrix sb(so_dimension(), basis_matrixkit());
    pl->symmetrize(s,sb);

    overlap_ = sb;
#else
    ExEnv::out0() << "Checking symmetrized overlap" << endl;

    RefSymmSCMatrix s(basis()->basisdim(), basis()->matrixkit());
    Ref<SCElementOp> ov =
      new OneBodyIntOp(new OneBodyIntIter(integral()->overlap()));

    s.assign(0.0);
    s.element_op(ov);
    ov=0;

    overlap_ = pl->to_SO_basis(s);

    //// use petite list to form saopl

    // form skeleton Hcore in AO basis
    RefSymmSCMatrix saopl(basis()->basisdim(), basis()->matrixkit());
    saopl.assign(0.0);

    ov = new OneBodyIntOp(new SymmOneBodyIntIter(integral_->overlap(), pl));
    saopl.element_op(ov);
    ov=0;

    // also symmetrize using pl->symmetrize():
    RefSymmSCMatrix spl(so_dimension(), basis_matrixkit());
    pl->symmetrize(saopl,spl);

    //// compare the answers

    int n = overlap_.result_noupdate().dim().n();
    int me = MessageGrp::get_default_messagegrp()->me();
    for (int i=0; i<n; i++) {
      for (int j=0; j<=i; j++) {
        double val1 = overlap_.result_noupdate().get_element(i,j);
        double val2 = spl.get_element(i,j);
        if (me == 0) {
          if (fabs(val1-val2) > 1.0e-6) {
            ExEnv::out0() << "bad overlap vals for " << i << " " << j
                         << ": " << val1 << " " << val2 << endl;
          }
        }
      }
    }
#endif

    overlap_.computed() = 1;
  }

  return overlap_.result_noupdate();
}

RefSymmSCMatrix
Wavefunction::core_hamiltonian()
{
  if (!hcore_.computed()) {
    integral()->set_basis(gbs_);
    Ref<PetiteList> pl = integral()->petite_list();

#if ! CHECK_SYMMETRIZED_INTEGRALS
    // form skeleton Hcore in AO basis
    RefSymmSCMatrix hao(basis()->basisdim(), basis()->matrixkit());
    hao.assign(0.0);

    Ref<SCElementOp> hc =
      new OneBodyIntOp(new SymmOneBodyIntIter(integral_->kinetic(), pl));
    hao.element_op(hc);
    hc=0;

    Ref<OneBodyInt> nuc = integral_->nuclear();
    nuc->reinitialize();
    hc = new OneBodyIntOp(new SymmOneBodyIntIter(nuc, pl));
    hao.element_op(hc);
    hc=0;

    // now symmetrize Hso
    RefSymmSCMatrix h(so_dimension(), basis_matrixkit());
    pl->symmetrize(hao,h);

    hcore_ = h;
#else
    ExEnv::out0() << "Checking symmetrized hcore" << endl;

    RefSymmSCMatrix hao(basis()->basisdim(), basis()->matrixkit());
    hao.assign(0.0);

    Ref<SCElementOp> hc =
      new OneBodyIntOp(new OneBodyIntIter(integral_->kinetic()));
    hao.element_op(hc);
    hc=0;

    Ref<OneBodyInt> nuc = integral_->nuclear();
    nuc->reinitialize();
    hc = new OneBodyIntOp(new OneBodyIntIter(nuc));
    hao.element_op(hc);
    hc=0;

    hcore_ = pl->to_SO_basis(hao);

    //// use petite list to form haopl

    // form skeleton Hcore in AO basis
    RefSymmSCMatrix haopl(basis()->basisdim(), basis()->matrixkit());
    haopl.assign(0.0);

    hc = new OneBodyIntOp(new SymmOneBodyIntIter(integral_->kinetic(), pl));
    haopl.element_op(hc);
    hc=0;

    nuc = integral_->nuclear();
    nuc->reinitialize();
    hc = new OneBodyIntOp(new SymmOneBodyIntIter(nuc, pl));
    haopl.element_op(hc);
    hc=0;

    // also symmetrize using pl->symmetrize():
    RefSymmSCMatrix h(so_dimension(), basis_matrixkit());
    pl->symmetrize(haopl,h);

    //// compare the answers

    int n = hcore_.result_noupdate().dim().n();
    int me = MessageGrp::get_default_messagegrp()->me();
    for (int i=0; i<n; i++) {
      for (int j=0; j<=i; j++) {
        double val1 = hcore_.result_noupdate().get_element(i,j);
        double val2 = h.get_element(i,j);
        if (me == 0) {
          if (fabs(val1-val2) > 1.0e-6) {
            ExEnv::outn() << "bad hcore vals for " << i << " " << j
                         << ": " << val1 << " " << val2 << endl;
          }
        }
      }
    }
#endif

    hcore_.computed() = 1;
  }

  return hcore_.result_noupdate();
}

// returns the orthogonalization matrix
RefSCMatrix
Wavefunction::so_to_orthog_so()
{
  if (orthog_.null()) init_orthog();
  return orthog_->basis_to_orthog_basis();
}

RefSCMatrix
Wavefunction::so_to_orthog_so_inverse()
{
  if (orthog_.null()) init_orthog();
  return orthog_->basis_to_orthog_basis_inverse();
}

Ref<GaussianBasisSet>
Wavefunction::basis() const
{
  return gbs_;
}

Ref<Integral>
Wavefunction::integral()
{
  return integral_;
}

RefSCDimension
Wavefunction::so_dimension()
{
  return sodim_;
}

RefSCDimension
Wavefunction::ao_dimension()
{
  return aodim_;
}

RefSCDimension
Wavefunction::oso_dimension()
{
  if (orthog_.null()) init_orthog();
  return orthog_->orthog_dim();
}

Ref<SCMatrixKit>
Wavefunction::basis_matrixkit()
{
  return basiskit_;
}

void
Wavefunction::print(ostream&o) const
{
  MolecularEnergy::print(o);
  basis()->print_brief(o);
  // the other stuff is a wee bit too big to print
  if (print_nao_ || print_npa_) {
    tim_enter("NAO");
    RefSCMatrix naos = ((Wavefunction*)this)->nao();
    tim_exit("NAO");
    if (print_nao_) naos.print("NAO");
  }
}
    
RefSymmSCMatrix
Wavefunction::alpha_density()
{
  if (!spin_polarized()) {
    RefSymmSCMatrix result = density().copy();
    result.scale(0.5);
    return result;
  }
  ExEnv::errn() << class_name() << "::alpha_density not implemented" << endl;
  abort();
  return 0;
}

RefSymmSCMatrix
Wavefunction::beta_density()
{
  if (!spin_polarized()) {
    RefSymmSCMatrix result = density().copy();
    result.scale(0.5);
    return result;
  }
  ExEnv::errn() << class_name() << "::beta_density not implemented" << endl;
  abort();
  return 0;
}

RefSymmSCMatrix
Wavefunction::alpha_ao_density()
{
  return integral()->petite_list()->to_AO_basis(alpha_density());
}

RefSymmSCMatrix
Wavefunction::beta_ao_density()
{
  return integral()->petite_list()->to_AO_basis(beta_density());
}

void
Wavefunction::obsolete()
{
  orthog_ = 0;

  MolecularEnergy::obsolete();
}

void
Wavefunction::copy_orthog_info(const Ref<Wavefunction>&wfn)
{
  if (orthog_.nonnull()) {
    ExEnv::errn() << "WARNING: Wavefunction: orthogonalization info changing"
                 << endl;
  }
  if (wfn->orthog_.null())
    wfn->init_orthog();
  orthog_ = wfn->orthog_->copy();
}

void
Wavefunction::init_orthog()
{
  orthog_ = new OverlapOrthog(
    orthog_method_,
    overlap(),
    basis_matrixkit(),
    lindep_tol_,
    debug_
    );
}

double
Wavefunction::min_orthog_res()
{
  return orthog_->min_orthog_res();
}

double
Wavefunction::max_orthog_res()
{
  if (orthog_.null()) init_orthog();
  return orthog_->max_orthog_res();
}

OverlapOrthog::OrthogMethod
Wavefunction::orthog_method() const
{
  return orthog_method_;
}

double
Wavefunction::lindep_tol() const
{
  return lindep_tol_;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
