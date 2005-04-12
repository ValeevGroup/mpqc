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

#include <stdexcept>
#include <iostream>
#include <iterator>

#include <util/keyval/keyval.h>
#include <util/misc/timer.h>
#include <util/misc/formio.h>
#include <util/misc/autovec.h>
#include <util/state/stateio.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/symmint.h>
#include <chemistry/qc/intv3/intv3.h>

#include <chemistry/qc/wfn/wfn.h>

using namespace std;
using namespace sc;

#define CHECK_SYMMETRIZED_INTEGRALS 0

/////////////////////////////////////////////////////////////////////////

// This maps a TwoBodyThreeCenterInt into a OneBodyInt
namespace sc {
class ChargeDistInt: public OneBodyInt {
    Ref<TwoBodyThreeCenterInt> tbint_;
    Ref<Molecule> mol_;
    Ref<GaussianBasisSet> ebasis0_;
    Ref<GaussianBasisSet> ebasis1_;
    Ref<GaussianBasisSet> atom_basis_;
    std::vector<int> i_cd_;
    const double *coef_;
  public:
    // The electronic basis are bs1 and bs2 in tbint and the
    // nuclear basis is bs3.
    ChargeDistInt(const Ref<TwoBodyThreeCenterInt>& tbint,
                  const double *coef);

    void compute_shell(int,int);

    bool cloneable();
};

ChargeDistInt::ChargeDistInt(const Ref<TwoBodyThreeCenterInt>& tbint,
                             const double *coef):
  OneBodyInt(tbint->integral(),tbint->basis1(),tbint->basis2()),
  tbint_(tbint),
  coef_(coef)
{
  ebasis0_ = tbint->basis1();
  ebasis1_ = tbint->basis2();
  atom_basis_ = tbint->basis3();
  mol_ = atom_basis_->molecule();
  buffer_ = new double[tbint->basis1()->max_nfunction_in_shell()
                       *tbint->basis2()->max_nfunction_in_shell()];

  for (int i=0; i<atom_basis_->ncenter(); i++) {
    if (atom_basis_->nshell_on_center(i) > 0) i_cd_.push_back(i);
  }
}

void
ChargeDistInt::compute_shell(int ish,int jsh)
{
//   std::cout << "starting " << ish << " " << jsh << std::endl;
  int nijbf
    = ebasis0_->shell(ish).nfunction()
    * ebasis1_->shell(jsh).nfunction();
  memset(buffer_,0,nijbf*sizeof(buffer_[0]));
  const double *tbbuffer = tbint_->buffer();
  int ksh = 0;
  int coef_offset = 0;
  for (int ii=0; ii<i_cd_.size(); ii++) {
    int i = i_cd_[ii];
    int nshell = atom_basis_->nshell_on_center(i);
    for (int j=0; j<nshell; j++, ksh++) {
      std::cout << "computing " << ish << " " << jsh << " " << ksh << std::endl;
      tbint_->compute_shell(ish,jsh,ksh);
      int nbfk = atom_basis_->shell(i).nfunction();
      for (int ijbf=0; ijbf<nijbf; ijbf++) {
        for (int kbf=0; kbf<nbfk; kbf++) {
          buffer_[ijbf]
            -= coef_[coef_offset+kbf] * tbbuffer[ijbf*nbfk + kbf];
//           std::cout << "  adding "
//                     << coef_[coef_offset+kbf] * tbbuffer[ijbf*nbfk + kbf]
//                     << " = "
//                     << coef_[coef_offset+kbf]
//                     << " * "
//                     << tbbuffer[ijbf*nbfk + kbf]
//                     << " at location "
//                     << ijbf
//                     << std::endl;
        }
      }
    }
    coef_offset += nshell;
  }
//   std::cout << "done with " << ish << " " << jsh << std::endl;
}

bool
ChargeDistInt::cloneable()
{
  // not cloneable because tbint is not cloneable
  return false;
}

}

/////////////////////////////////////////////////////////////////////////

static ClassDesc Wavefunction_cd(
  typeid(Wavefunction),"Wavefunction",7,"public MolecularEnergy",
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

  atom_basis_ << keyval->describedclassvalue("atom_basis");
  if (atom_basis_.nonnull()) {
    atom_basis_coef_ = new double[atom_basis_->nbasis()];
    for (int i=0; i<atom_basis_->nbasis(); i++) {
      atom_basis_coef_[i] = keyval->doublevalue("atom_basis_coef",i);
    }
    scale_atom_basis_coef();

  }
  else {
    atom_basis_coef_ = 0;
  }

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

  if (s.version(::class_desc<Wavefunction>()) >= 6) {
    Ref<KeyVal> original_override = s.override();
    Ref<AssignedKeyVal> matrixkit_override = new AssignedKeyVal;
    matrixkit_override->assign("matrixkit", gbs_->so_matrixkit().pointer());
    Ref<KeyVal> new_override
      = new AggregateKeyVal(matrixkit_override.pointer(),
                            original_override.pointer());
    s.set_override(new_override);
    orthog_ << SavableState::restore_state(s);
    s.set_override(original_override);
  }

  if (s.version(::class_desc<Wavefunction>()) >= 7) {
    atom_basis_ << SavableState::restore_state(s);
    if (atom_basis_.nonnull()) {
      s.get_array_double(atom_basis_coef_, atom_basis_->nbasis());
    }
  }
  else {
    atom_basis_coef_ = 0;
  }

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
  SavableState::save_state(orthog_.pointer(),s);
  SavableState::save_state(atom_basis_.pointer(), s);
  if (atom_basis_.nonnull()) {
    s.put_array_double(atom_basis_coef_,atom_basis_->nbasis());
  }
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
      densortho.accumulate_transform(so_to_orthog_so_inverse().t(),dens);

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

    if (atom_basis_.null()) {
      Ref<OneBodyInt> nuc = integral_->nuclear();
      nuc->reinitialize();
      hc = new OneBodyIntOp(new SymmOneBodyIntIter(nuc, pl));
      hao.element_op(hc);
      hc=0;
    }
    else {
      // we have an atom_basis, so some nuclear charges will be computed
      // with the two electron integral code and some will be computed
      // with the point charge code

      // find out which atoms are point charges and which ones are charge
      // distributions
      std::vector<int> i_pc; // point charge centers
      for (int i=0; i<atom_basis_->ncenter(); i++) {
        if (atom_basis_->nshell_on_center(i) == 0) i_pc.push_back(i);
      }
      int n_pc = i_pc.size();

      // initialize the point charge data
      auto_vec<double> pc_charges(new double[n_pc]);
      auto_vec<double*> pc_xyz(new double*[n_pc]);
      auto_vec<double> pc_xyz0(new double[n_pc*3]);
      pc_xyz[0] = pc_xyz0.get();
      for (int i=1; i<n_pc; i++) pc_xyz[i] = pc_xyz[i-1] + 3;
      for (int i=0; i<n_pc; i++) {
        pc_charges[i] = molecule()->charge(i_pc[i]);
        for (int j=0; j<3; j++) pc_xyz[i][j] = molecule()->r(i_pc[i],j);
      }
      Ref<PointChargeData> pc_data
        = new PointChargeData(n_pc, pc_xyz.get(), pc_charges.get());

      // compute the point charge contributions
      Ref<OneBodyInt> pc_int = integral_->point_charge(pc_data);
      hc = new OneBodyIntOp(new SymmOneBodyIntIter(pc_int,pl));
      hao.element_op(hc);
      hc=0;
      pc_int=0;

      // compute the charge distribution contributions
      // H_ij += sum_A -q_A sum_k N_A_k (ij|A_k)
      integral()->set_basis(gbs_,gbs_,atom_basis_);
      Ref<TwoBodyThreeCenterInt> cd_tbint
        = integral_->electron_repulsion3();
      Ref<OneBodyInt> cd_int = new ChargeDistInt(cd_tbint, atom_basis_coef_);
      hc = new OneBodyIntOp(new SymmOneBodyIntIter(cd_int,pl));
      hao.element_op(hc);
      integral()->set_basis(gbs_);
      hc=0;
      cd_int=0;
      cd_tbint=0;
    }

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

//     std::cout << "null atom_basis" << std::endl;
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

// Compute lists of centers that are point charges and lists that
// are charge distributions.
void
Wavefunction::set_up_charge_types(
  std::vector<int> &q_pc,
  std::vector<int> &q_cd,
  std::vector<int> &n_pc,
  std::vector<int> &n_cd)
{
  q_pc.clear();
  q_cd.clear();
  n_pc.clear();
  n_cd.clear();

  for (int i=0; i<atom_basis_->ncenter(); i++) {
    bool is_Q = (molecule()->atom_symbol(i) == "Q");
    if (atom_basis_->nshell_on_center(i) > 0) {
      if (is_Q) q_cd.push_back(i);
      else      n_cd.push_back(i);
    }
    else {
      if (is_Q) q_pc.push_back(i);
      else      n_pc.push_back(i);
    }
  }
}

double
Wavefunction::nuclear_repulsion_energy()
{
  if (atom_basis_.null()) return molecule()->nuclear_repulsion_energy();

  double nucrep = 0.0;

  std::vector<int> q_pc, q_cd, n_pc, n_cd;
  set_up_charge_types(q_pc,q_cd,n_pc,n_cd);

  // compute point charge - point charge contribution
  nucrep += nuc_rep_pc_pc(n_pc, n_pc, true /* i > j  */);
  nucrep += nuc_rep_pc_pc(q_pc, n_pc, false /* all i j */);
  if (molecule()->include_qq()) {
    nucrep += nuc_rep_pc_pc(q_pc, q_pc, true /* i > j */);
  }

  // compute point charge - charge distribution contribution
  nucrep += nuc_rep_pc_cd(n_pc, n_cd);
  nucrep += nuc_rep_pc_cd(q_pc, n_cd);
  nucrep += nuc_rep_pc_cd(n_pc, q_cd);
  if (molecule()->include_qq()) {
    nucrep += nuc_rep_pc_cd(q_pc, q_cd);
  }

  // compute the charge distribution - charge distribution contribution
  nucrep += nuc_rep_cd_cd(n_cd, n_cd, true /* i > j  */);
  nucrep += nuc_rep_cd_cd(q_cd, n_cd, false /* all i j  */);
  if (molecule()->include_qq()) {
    nucrep += nuc_rep_cd_cd(q_cd, q_cd, true /* i > j */);
  }

  return nucrep;
}

double
Wavefunction::nuc_rep_pc_pc(const std::vector<int>&c1,
                            const std::vector<int>&c2,
                            bool uniq)
{
  double e = 0.0;

  if (c1.size() == 0 || c2.size() == 0) return e;

  for (int ii=0; ii<c1.size(); ii++) {
    int i = c1[ii];
    SCVector3 ai(molecule()->r(i));
    double Zi = molecule()->charge(i);
    int jfence = (uniq?ii:c2.size());
    for (int jj=0; jj<jfence; jj++) {
      int j = c2[jj];
      SCVector3 aj(molecule()->r(j));
      e += Zi * molecule()->charge(j) / ai.dist(aj);
    }
  }

  return e;
}

double
Wavefunction::nuc_rep_pc_cd(const std::vector<int>&pc,
                            const std::vector<int>&cd)
{
  double e = 0.0;

  if (pc.size() == 0 || cd.size() == 0) return e;

  integral()->set_basis(atom_basis());

  sc::auto_vec<double> charges(new double[pc.size()]);
  sc::auto_vec<double*> positions(new double*[pc.size()]);
  sc::auto_vec<double> xyz(new double[pc.size()*3]);
  for (int i=0; i<pc.size(); i++) {
    positions.get()[i] = &xyz.get()[i*3];
  }

  for (int ii=0; ii<pc.size(); ii++) {
    int i=pc[ii];
    charges.get()[ii] = molecule()->charge(i);
    for (int j=0; j<3; j++) positions.get()[ii][j] = molecule()->r(i,j);
  }

  Ref<PointChargeData> pcdata = new PointChargeData(pc.size(),
                                                    positions.get(),
                                                    charges.get());
  Ref<OneBodyOneCenterInt> ob = integral()->point_charge1(pcdata);

  const double *buffer = ob->buffer();

  for (int ii=0,icoef=0; ii<cd.size(); ii++) {
    int icenter = cd[ii];
    int joff = atom_basis()->shell_on_center(icenter,0);
    for (int j=0; j<atom_basis()->nshell_on_center(icenter); j++) {
      int jsh = j + joff;
      ob->compute_shell(jsh);
      int nfunc = atom_basis()->shell(jsh).nfunction();
      for (int k=0; k<nfunc; k++,icoef++) {
        e -= atom_basis_coef_[icoef] * buffer[k];
      }
    }
  }

  integral()->set_basis(basis());

  return e;
}

double
Wavefunction::nuc_rep_cd_cd(const std::vector<int>&c1,
                            const std::vector<int>&c2,
                            bool uniq)
{
  double e = 0.0;

  if (c1.size() == 0 || c2.size() == 0) return e;

  integral()->set_basis(atom_basis());

  Ref<TwoBodyTwoCenterInt> tb = integral()->electron_repulsion2();

  const double *buffer = tb->buffer();

  for (int ii=0; ii<c1.size(); ii++) {
    int icenter = c1[ii];
    int inshell = atom_basis()->nshell_on_center(icenter);
    int ish = atom_basis()->shell_on_center(icenter,0);
    for (int iish=0; iish<inshell; iish++,ish++) {
      int infunc = atom_basis()->shell(ish).nfunction();
      int ifuncoff = atom_basis()->shell_to_function(ish);
      int jjfence = (uniq?ii:c2.size());
      for (int jj=0; jj<jjfence; jj++) {
        int jcenter = c2[jj];
        int jnshell = atom_basis()->nshell_on_center(jcenter);
        int jsh = atom_basis()->shell_on_center(jcenter,0);
        for (int jjsh=0; jjsh<jnshell; jjsh++,jsh++) {
          int jnfunc = atom_basis()->shell(jsh).nfunction();
          int jfuncoff = atom_basis()->shell_to_function(jsh);
          tb->compute_shell(ish,jsh);
          for (int ifunc=0, ijfunc=0; ifunc<infunc; ifunc++) {
            for (int jfunc=0; jfunc<jnfunc; jfunc++, ijfunc++) {
              e += atom_basis_coef_[ifuncoff+ifunc]
                * atom_basis_coef_[jfuncoff+jfunc]
                * buffer[ijfunc];
            }
          }
        }
      }
    }
  }

  integral()->set_basis(basis());

  return e;
}

void
Wavefunction::nuclear_repulsion_energy_gradient(double *g)
{
  int natom = molecule()->natom();
  sc::auto_vec<double*> gp(new double*[natom]);
  for (int i=0; i<natom; i++) {
    gp.get()[i] = &g[i*3];
  }
  nuclear_repulsion_energy_gradient(gp.get());
}

void
Wavefunction::nuclear_repulsion_energy_gradient(double **g)
{
  if (atom_basis_.null()) {
    int natom = molecule()->natom();
    for (int i=0; i<natom; i++) {
      molecule()->nuclear_repulsion_1der(i,g[i]);
    }
    return;
  }

  // zero the gradient
  for (int i=0; i<molecule()->natom(); i++) {
    for (int j=0; j<3; j++) g[i][j] = 0.0;
  }

  // compute charge types
  std::vector<int> q_pc, q_cd, n_pc, n_cd;
  set_up_charge_types(q_pc,q_cd,n_pc,n_cd);

  // compute point charge - point charge contribution
  nuc_rep_grad_pc_pc(g, n_pc, n_pc, true /* i > j  */);
  nuc_rep_grad_pc_pc(g, q_pc, n_pc, false /* all i j */);
  if (molecule()->include_qq()) {
    nuc_rep_grad_pc_pc(g, q_pc, q_pc, true /* i > j */);
  }

  // compute point charge - charge distribution contribution
  nuc_rep_grad_pc_cd(g, n_pc, n_cd);
  nuc_rep_grad_pc_cd(g, q_pc, n_cd);
  nuc_rep_grad_pc_cd(g, n_pc, q_cd);
  if (molecule()->include_qq()) {
    nuc_rep_grad_pc_cd(g, q_pc, q_cd);
  }

  // compute the charge distribution - charge distribution contribution
  nuc_rep_grad_cd_cd(g, n_cd, n_cd, true /* i > j  */);
  nuc_rep_grad_cd_cd(g, q_cd, n_cd, false /* all i j  */);
  if (molecule()->include_qq()) {
    nuc_rep_grad_cd_cd(g, q_cd, q_cd, true /* i > j */);
  }

  // note: the electronic terms still need to be done in
  // a new hcore_deriv implemented in Wavefunction.
  throw std::runtime_error("Wavefunction::nuclear_repulsion_energy_gradient: not done");

}

void
Wavefunction::nuc_rep_grad_pc_pc(double **grad,
                                 const std::vector<int>&c1,
                                 const std::vector<int>&c2,
                                 bool uniq)
{
  if (c1.size() == 0 || c2.size() == 0) return;

  for (int ii=0; ii<c1.size(); ii++) {
    int i = c1[ii];
    SCVector3 ai(molecule()->r(i));
    double Zi = molecule()->charge(i);
    int jfence = (uniq?ii:c2.size());
    for (int jj=0; jj<jfence; jj++) {
      int j = c2[jj];
      SCVector3 aj(molecule()->r(j));
      double Zj = molecule()->charge(j);
      SCVector3 rdiff = ai - aj;
      double r2 = rdiff.dot(rdiff);
      double factor = - Zi * Zj / (r2*sqrt(r2));
      for (int k=0; k<3; k++) {
        grad[i][k] += factor * rdiff[k];
        grad[j][k] -= factor * rdiff[k];
      }
    }
  }
}

void
Wavefunction::nuc_rep_grad_pc_cd(double **grad,
                                 const std::vector<int>&pc,
                                 const std::vector<int>&cd)
{
  if (pc.size() == 0 || cd.size() == 0) return;

  throw std::runtime_error("Wavefunction::nuclear_repulsion_energy_gradient: not done");
}

void
Wavefunction::nuc_rep_grad_cd_cd(double **grad,
                                 const std::vector<int>&c1,
                                 const std::vector<int>&c2,
                                 bool uniq)
{
  if (c1.size() == 0 || c2.size() == 0) return;

  throw std::runtime_error("Wavefunction::nuclear_repulsion_energy_gradient: not done");
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

Ref<Molecule>
Wavefunction::molecule() const
{
  return basis()->molecule();
}

Ref<GaussianBasisSet>
Wavefunction::atom_basis() const
{
  return atom_basis_;
}

const double *
Wavefunction::atom_basis_coef() const
{
  return atom_basis_coef_;
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
  ExEnv::out0() << indent << "Electronic basis:" << std::endl;
  ExEnv::out0() << incindent;
  basis()->print_brief(o);
  ExEnv::out0() << decindent;
  if (atom_basis_.nonnull()) {
    ExEnv::out0() << indent << "Nuclear basis:" << std::endl;
    ExEnv::out0() << incindent;
    atom_basis_->print_brief(o);
    ExEnv::out0() << decindent;
  }
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

void
Wavefunction::set_orthog_method(const OverlapOrthog::OrthogMethod& omethod)
{
  if (orthog_method_ != omethod) {
    orthog_method_ = omethod;
    init_orthog();
    obsolete();
  }
}

double
Wavefunction::lindep_tol() const
{
  return lindep_tol_;
}

void
Wavefunction::set_lindep_tol(double lindep_tol)
{
  if (lindep_tol_ != lindep_tol) {
    lindep_tol_ = lindep_tol;
    init_orthog();
    obsolete();
  }
}

void
Wavefunction::scale_atom_basis_coef()
{
  std::vector<int> i_cd;
  for (int i=0; i<atom_basis_->ncenter(); i++) {
    if (atom_basis_->nshell_on_center(i) > 0) i_cd.push_back(i);
  }

  if (atom_basis_->max_angular_momentum() > 0) {
    // Only s functions will preserve the full symmetry.
    // Can only normalize functions that don't integrate to zero.
    throw std::runtime_error("ChargeDistInt: max am larger than 0");
  }

  int coef_offset = 0;
  int icoef = 0;
  for (int ii=0; ii<i_cd.size(); ii++) {
    int i = i_cd[ii];
    int nshell = atom_basis_->nshell_on_center(i);
    int nfunc = 0;
    if (nshell > 0) {
      double raw_charge = 0.0;
      for (int jj=0, icoef=coef_offset; jj<nshell; jj++) {
        int j = atom_basis_->shell_on_center(i,jj);
        const GaussianShell &shell = atom_basis_->shell(j);
        // loop over the contractions
        // the number of functions in each contraction is 1
        // since only s functions are allowed
        for (int k=0; k<shell.ncontraction(); k++, icoef++) {
          for (int l=0; l<shell.nprimitive(); l++) {
            double alpha = shell.exponent(l);
            double con_coef = shell.coefficient_unnorm(k,l);
            // The exponent is halved because the normalization
            // formula is for the s function squared.
            double integral = ::pow(M_PI/alpha,1.5);
            raw_charge += atom_basis_coef_[icoef] * con_coef * integral;
//             std::cout << "con_coef = " << con_coef
//                       << " integral = " << integral
//                       << std::endl;
          }
        }
        nfunc += shell.ncontraction();
      }
      double charge = atom_basis_->molecule()->charge(i);
      double correction = charge/raw_charge;
      for (int icoef=coef_offset; icoef<coef_offset+nfunc; icoef++) {
        atom_basis_coef_[icoef] = correction * atom_basis_coef_[icoef];
      }
    }
    coef_offset += nshell;
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
