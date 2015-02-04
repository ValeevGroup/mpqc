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

#include <stdlib.h>
#include <math.h>

#include <stdexcept>
#include <iostream>
#include <iterator>

#include <util/keyval/keyval.h>
#include <util/misc/regtime.h>
#include <util/misc/formio.h>
#include <util/misc/autovec.h>
#include <util/misc/xmlwriter.h>
#include <util/state/stateio.h>
#include <util/misc/scexception.h>
#include <chemistry/qc/basis/uncontract.h>
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

    bool cloneable() const;
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
ChargeDistInt::cloneable() const
{
  // not cloneable because tbint is not cloneable
  return false;
}

}

/////////////////////////////////////////////////////////////////////////

static ClassDesc Wavefunction_cd(
  typeid(Wavefunction),"Wavefunction",9,"public MolecularEnergy",
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

  KeyValValuedouble lindep_tol_def(OverlapOrthog::default_lindep_tol());
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
    std::string orthog_name = keyval->stringvalue("orthog_method");
    if (orthog_name.empty()) {
      orthog_method_ = OverlapOrthog::default_orthog_method();
    }
    else if (orthog_name == "canonical") {
      orthog_method_ = OverlapOrthog::Canonical;
    }
    else if (orthog_name == "symmetric") {
      orthog_method_ = OverlapOrthog::Symmetric;
    }
    else if (orthog_name == "gramschmidt") {
      orthog_method_ = OverlapOrthog::GramSchmidt;
    }
    else {
      ExEnv::errn() << "ERROR: bad orthog_method: \""
                   << orthog_name << "\"" << endl;
      abort();
    }
  }

  debug_ = keyval->intvalue("debug");

  gbs_ = require_dynamic_cast<GaussianBasisSet*>(
    keyval->describedclassvalue("basis").pointer(),
    "Wavefunction::Wavefunction\n"
    );
  if (gbs_.null())
    throw InputError("Wavefunction::Wavefunction -- basis is missing");

  atom_basis_ << keyval->describedclassvalue("atom_basis");
  if (atom_basis_) {
    atom_basis_coef_ = new double[atom_basis_->nbasis()];
    for (int i=0; i<atom_basis_->nbasis(); i++) {
      atom_basis_coef_[i] = keyval->doublevalue("atom_basis_coef",i);
    }
    scale_atom_basis_coef();

  }
  else {
    atom_basis_coef_ = 0;
  }

  momentum_basis_ << keyval->describedclassvalue("momentum_basis");
  dk_ = keyval->intvalue("dk");
  if (dk_ > 0 && momentum_basis_.null()) {
    momentum_basis_ = new UncontractedBasisSet(gbs_);
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

  magnetic_moment_ = aodim_.n() + 1;

  // post-construction validation
  {
    if (electric_field()) {
      if (fabs(electric_field().dot(electric_field())) > 1e-13 &&
          not Wavefunction::nonzero_efield_supported())
        throw FeatureNotImplemented("External electric field is specified, \
but this Wavefunction cannot be computed in its presence",
              __FILE__, __LINE__, this->class_desc());
    }
  }
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
    if (atom_basis_) {
      s.get_array_double(atom_basis_coef_, atom_basis_->nbasis());
    }
  }
  else {
    atom_basis_coef_ = 0;
  }

  if (s.version(::class_desc<Wavefunction>()) >= 8) {
    s.get(dk_);
    momentum_basis_ << SavableState::restore_state(s);
  }
  else {
    dk_ = 0;
  }

  if (s.version(::class_desc<Wavefunction>()) >= 9) {
    s.get(magnetic_moment_);
  }
  else {
    magnetic_moment_ = aodim_.n() + 1;
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
  magnetic_moment_ = aodim_.n() + 1;
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
  if (atom_basis_) {
    s.put_array_double(atom_basis_coef_,atom_basis_->nbasis());
  }

  s.put(dk_);
  SavableState::save_state(momentum_basis_.pointer(), s);

  s.put(magnetic_moment_);
}

ptree&
Wavefunction::write_xml(
    ptree& parent,
    const XMLWriter& writer
)
{
  ptree& child = get_my_ptree(parent);

  if(natural_orbitals_.computed()){
    writer.insert_child(
        child, natural_orbitals_.result_noupdate(), "natural_orbitals"
    );
  }
  if(natural_density_.computed()){
    writer.insert_child(
        child, natural_density_.result_noupdate(), "natural_density"
    );
  }
  child.put("total_charge", total_charge());
  child.put("spin_polarized", (bool)spin_polarized());
  writer.insert_child(child, basis(), "basis");
  if(atom_basis().nonnull() and not atom_basis()->equiv(basis())){
    ptree& atom_basis_tree = writer.insert_child(
        child, atom_basis(), "atom_basis"
    );
    atom_basis_tree.put("comment",
        "basis set describing the nuclear charge distributions"
    );
  }
  return MolecularEnergy::write_xml(parent, writer);
}

double
Wavefunction::total_charge() const
{
  return molecule()->total_charge() - const_cast<Wavefunction*>(this)->nelectron();
}

double
Wavefunction::magnetic_moment() const
{
  Wavefunction* this_ptr_nonconst = const_cast<Wavefunction*>(this);
  if (magnetic_moment_ > aodim_.n()) // magnetic moment greater than the number of states means it has not been computed yet.
    magnetic_moment_ = (this_ptr_nonconst->alpha_density() * this_ptr_nonconst->overlap()).trace() -
                       (this_ptr_nonconst->beta_density() * this_ptr_nonconst->overlap()).trace();
  return magnetic_moment_;
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

void
Wavefunction::core_hamiltonian_dk2_contrib(const RefSymmSCMatrix &h_pbas,
                                           const RefDiagSCMatrix &E,
                                           const RefDiagSCMatrix &K,
                                           const RefDiagSCMatrix &p2,
                                           const RefDiagSCMatrix &p2K2,
                                           const RefDiagSCMatrix &p2K2_inv,
                                           const RefSymmSCMatrix &AVA_pbas,
                                           const RefSymmSCMatrix &BpVpB_pbas)
{
  RefSCDimension pdim = AVA_pbas.dim();

  RefSymmSCMatrix AVA_prime = AVA_pbas.clone();
  RefSymmSCMatrix BpVpB_prime = BpVpB_pbas.clone();
  int npbas = pdim.n();
  for (int i=0; i<npbas; i++) {
    double Ei = E(i);
    for (int j=0; j<=i; j++) {
      double Ej = E(j);
      AVA_prime(i,j) = AVA_pbas(i,j)/(Ei+Ej);
      BpVpB_prime(i,j) = BpVpB_pbas(i,j)/(Ei+Ej);
    }
  }

  RefSCMatrix h_contrib;
  h_contrib
    =
    - 1.0 * BpVpB_prime * E * AVA_prime
    - 0.5 * BpVpB_prime * AVA_prime * E
    - 0.5 * AVA_prime * BpVpB_prime * E
    + 0.5 * AVA_prime * p2K2 * AVA_prime * E
    + 0.5 * BpVpB_prime * p2K2_inv * BpVpB_prime * E
    + 0.5 * AVA_prime * (p2K2 * E) * AVA_prime
    + 0.5 * BpVpB_prime * (p2K2_inv * E) * BpVpB_prime
    ;

  h_pbas.accumulate_symmetric_sum(h_contrib);
}

void
Wavefunction::core_hamiltonian_dk3_contrib(const RefSymmSCMatrix &h_pbas,
                                           const RefDiagSCMatrix &E,
                                           const RefDiagSCMatrix &B,
                                           const RefDiagSCMatrix &p2K2_inv,
                                           const RefSCMatrix &so_to_p,
                                           const RefSymmSCMatrix &pxVp)
{
  RefSCDimension p_oso_dim = so_to_p.rowdim();
  Ref<SCMatrixKit> p_so_kit = so_to_p.kit();
  int noso = p_oso_dim.n();
  RefSymmSCMatrix pxVp_pbas(p_oso_dim, p_so_kit);
  pxVp_pbas.assign(0.0);
  pxVp_pbas.accumulate_transform(so_to_p, pxVp);

  RefSymmSCMatrix BpxVpB_prime(p_oso_dim, p_so_kit);
  for (int i=0; i<noso; i++) {
    double Ei = E(i);
    for (int j=0; j<=i; j++) {
      double Ej = E(j);
      BpxVpB_prime(i,j) = pxVp_pbas(i,j)*B(i)*B(j)/(Ei+Ej);
    }
  }

  RefSCDimension pdim = E.dim();

  RefSCMatrix h_contrib;
  h_contrib
    =
    - 0.5 * BpxVpB_prime * E * p2K2_inv * BpxVpB_prime
    - 0.5 * BpxVpB_prime * p2K2_inv * BpxVpB_prime * E
    ;

  h_pbas.accumulate_symmetric_sum(h_contrib);
}

RefSymmSCMatrix
Wavefunction::core_hamiltonian_dk(int dk,
                                  const Ref<GaussianBasisSet> &bas,
                                  const Ref<GaussianBasisSet> &p_bas)
{
#define DK_DEBUG 0

  ExEnv::out0() << indent << "Including DK" << dk
                << (dk==1?" (free particle projection)":"")
                << (dk==2?" (Douglas-Kroll-Hess)":"")
                << (dk==3?" (complete spin-free Douglas-Kroll)":"")
                << " terms in the one body Hamiltonian."
                << std::endl;

  if (atom_basis_) {
    throw FeatureNotImplemented("atom_basis given and dk > 0",
                                __FILE__, __LINE__, class_desc());
  }
  if (dk > 2) {
    throw FeatureNotImplemented("dk must be 0, 1, or 2",
                                __FILE__, __LINE__, class_desc());
  }
  if (dk > 0 && gradient_needed()) {
    throw FeatureNotImplemented("gradients not available for dk",
                                __FILE__, __LINE__, class_desc());
  }

  // The one electron integrals will be computed in the momentum basis.
  integral()->set_basis(p_bas);

  Ref<PetiteList> p_pl = integral()->petite_list();

  RefSCDimension p_so_dim = p_pl->SO_basisdim();
  RefSCDimension p_ao_dim = p_pl->AO_basisdim();
  Ref<SCMatrixKit> p_kit = p_bas->matrixkit();
  Ref<SCMatrixKit> p_so_kit = p_bas->so_matrixkit();

  // Compute the overlap in the momentum basis.
  RefSymmSCMatrix S_skel(p_ao_dim, p_kit);
  S_skel.assign(0.0);
  Ref<SCElementOp> hc =
    new OneBodyIntOp(new SymmOneBodyIntIter(integral()->overlap(), p_pl));
  S_skel.element_op(hc);
  hc=0;
  RefSymmSCMatrix S(p_so_dim, p_so_kit);
  p_pl->symmetrize(S_skel,S);

  ExEnv::out0() << indent
                << "The momentum basis is:"
                << std::endl;
  ExEnv::out0() << incindent;
  p_bas->print_brief(ExEnv::out0());
  ExEnv::out0() << decindent;

  ExEnv::out0() << indent
                << "Orthogonalizing the momentum basis"
                << std::endl;
  Ref<OverlapOrthog> p_orthog
    = new OverlapOrthog(orthog_method_, S, p_so_kit,
                        lindep_tol_, debug_);

  RefSCDimension p_oso_dim = p_orthog->orthog_dim();

  // form skeleton Hcore in the momentum basis
  RefSymmSCMatrix T_skel(p_ao_dim, p_kit);
  T_skel.assign(0.0);

  hc =
    new OneBodyIntOp(new SymmOneBodyIntIter(integral()->kinetic(), p_pl));
  T_skel.element_op(hc);
  hc=0;

  // finish constructing the kinetic energy integrals,
  // for which the skeleton is in hao
  RefSymmSCMatrix T(p_so_dim, p_so_kit);
  p_pl->symmetrize(T_skel,T);
  T_skel = 0;

  // Transform T into an orthogonal basis
  RefSymmSCMatrix T_oso(p_oso_dim, p_so_kit);
  T_oso.assign(0.0);
  T_oso.accumulate_transform(p_orthog->basis_to_orthog_basis(),T);

  // diagonalize the T integrals to get a momentum basis
  RefDiagSCMatrix Tval(p_oso_dim, p_so_kit);
  RefSCMatrix Tvec(p_oso_dim, p_oso_dim, p_so_kit);
  // Tvec * Tval * Tvec.t() = T_oso
  T_oso.diagonalize(Tval,Tvec);

  T_oso = 0;

#if DK_DEBUG
  T.print("T");
  Tval.print("Tval");
#endif

  // Compute the kinematic factors
  RefDiagSCMatrix A(p_oso_dim, p_so_kit);
  RefDiagSCMatrix B(p_oso_dim, p_so_kit);
  RefDiagSCMatrix E(p_oso_dim, p_so_kit);
  RefDiagSCMatrix K(p_oso_dim, p_so_kit);
  RefDiagSCMatrix p2(p_oso_dim, p_so_kit);
  RefDiagSCMatrix Emc2(p_oso_dim, p_so_kit);
  const double c = 137.0359895; // speed of light in a vacuum in a.u.
  int noso = p_oso_dim.n();
  for (int i=0; i<noso; i++) {
    double T_val = Tval(i);
    // momentum basis sets with near linear dependencies may
    // have T_val approximately equal to zero.  These can be
    // negative, which will cause a SIGFPE in sqrt.
    if (T_val < DBL_EPSILON) T_val = 0.0;
    double p = sqrt(2.0*T_val);
    double E_val = c * sqrt(p*p+c*c);
    double A_val = sqrt((E_val+c*c)/(2.0*E_val));
    double K_val = c/(E_val+c*c);
    double B_val = A_val * K_val;
    double Emc2_val = c*c*p*p/(E_val + c*c); // = E - mc^2
    A(i) = A_val;
    B(i) = B_val;
    E(i) = E_val;
    K(i) = K_val;
    p2(i) = p*p;
    Emc2(i) = Emc2_val;
  }

#if DK_DEBUG
  A.print("A");
  B.print("B");
  E.print("E");
  K.print("K");
  Emc2.print("Emc2");
#endif

  // Construct the transform from the coordinate to the momentum
  // representation in the momentum basis
  RefSCMatrix so_to_p
    = Tvec.t() * p_orthog->basis_to_orthog_basis();

#if DK_DEBUG
  so_to_p.print("so_to_p");
#endif

  // compute the V integrals
  Ref<OneBodyInt> V_obi = integral_->nuclear();
  V_obi->reinitialize();
  hc = new OneBodyIntOp(new SymmOneBodyIntIter(V_obi, p_pl));
  RefSymmSCMatrix V_skel(p_ao_dim, p_kit);
  V_skel.assign(0.0);
  V_skel.element_op(hc);
  V_obi = 0;
  hc = 0;
  RefSymmSCMatrix V(p_so_dim, p_so_kit);
  p_pl->symmetrize(V_skel,V);
  V_skel = 0;

  // include contributions from external electric field, if needed
  if (electric_field()) {
    RefSymmSCMatrix mu_so;
    RefSymmSCMatrix mu(p_ao_dim, p_kit);
    mu.assign(0.0);

    double E[3];  for(int xyz=0; xyz<3; ++xyz) E[xyz] = electric_field().get_element(xyz);
    Ref<GaussianBasisSet> bs = p_bas;
    const int nshell = bs->nshell();
    integral()->set_basis(bs,bs);
    Ref<OneBodyInt> m1_ints = integral()->dipole(0);
    for(int sh1=0; sh1<nshell; sh1++) {
      int bf1_offset = bs->shell_to_function(sh1);
      int nbf1 = bs->shell(sh1).nfunction();

      int sh2max = sh1;
      for(int sh2=0; sh2<=sh2max; sh2++) {
        int bf2_offset = bs->shell_to_function(sh2);
        int nbf2 = bs->shell(sh2).nfunction();

        m1_ints->compute_shell(sh1,sh2);
        const double *m1intsptr = m1_ints->buffer();

        int bf1_index = bf1_offset;
        for(int bf1=0; bf1<nbf1; bf1++, bf1_index++, m1intsptr+=3*nbf2) {
          int bf2_index = bf2_offset;
          const double *ptr1 = m1intsptr;
          int bf2max;
          if (sh1 == sh2)
            bf2max = bf1;
          else
            bf2max = nbf2-1;
          for(int bf2=0; bf2<=bf2max; bf2++, bf2_index++) {

            // V = - mu \dot E; additional -1 accounts for the negative charge of the electron
            const double V = (ptr1[0] * E[0] + ptr1[1] * E[1] + ptr1[2] * E[2]);
            mu.set_element(bf1_index, bf2_index, V);
            ptr1 += 3;

          }
        }
      }
    }
    m1_ints = 0;

    const int nbasis = bs->nbasis();
    for(int bf1=0; bf1<nbasis; bf1++)
      for(int bf2=0; bf2<=bf1; bf2++) {
        mu(bf2,bf1) = mu(bf1,bf2);
      }

    mu_so = p_pl->to_SO_basis(mu);
    mu = 0;
    V.accumulate(mu_so);
  }

#if DK_DEBUG
  V.print("V");
#endif

  // transform V to the momentum basis
  RefSymmSCMatrix V_pbas(p_oso_dim, p_so_kit);
  V_pbas.assign(0.0);
  V_pbas.accumulate_transform(so_to_p, V);

  // compute the p.Vp integrals
  Ref<OneBodyInt> pVp_obi = integral()->p_dot_nuclear_p();
  hc = new OneBodyIntOp(new SymmOneBodyIntIter(pVp_obi, p_pl));
  RefSymmSCMatrix pVp_skel(p_ao_dim, p_kit);
  pVp_skel.assign(0.0);
  pVp_skel.element_op(hc);
#if DK_DEBUG
  const double *buf = pVp_obi->buffer();
  for (int I=0,Ii=0; I<p_bas->nshell(); I++) {
    for (int i=0; i<p_bas->shell(I).nfunction(); i++,Ii++) {
      for (int J=0,Jj=0; J<p_bas->nshell(); J++) {
        pVp_obi->compute_shell(I,J);
        int ij = i*p_bas->shell(J).nfunction();
        for (int j=0; j<p_bas->shell(J).nfunction(); j++,ij++,Jj++) {
          std::cout << "pVp["<<Ii<<"]["<<Jj<<"][0]= " << buf[ij]
                    << std::endl;
        }
      }
    }
  }
#endif
  pVp_obi = 0;
  hc = 0;
  RefSymmSCMatrix pVp(p_so_dim, p_so_kit);
  p_pl->symmetrize(pVp_skel,pVp);
  pVp_skel = 0;

#if DK_DEBUG
  pVp.print("pVp");
  (-2.0*T).print("-2*T");
#endif

  // transform p.Vp to the momentum basis
  RefSymmSCMatrix pVp_pbas(p_oso_dim, p_so_kit);
  pVp_pbas.assign(0.0);
  pVp_pbas.accumulate_transform(so_to_p, pVp);

  RefSymmSCMatrix AVA_pbas(p_oso_dim, p_so_kit);
  RefSymmSCMatrix BpVpB_pbas(p_oso_dim, p_so_kit);
  for (int i=0; i<noso; i++) {
    for (int j=0; j<=i; j++) {
      AVA_pbas(i,j) = V_pbas(i,j)*A(i)*A(j);
      BpVpB_pbas(i,j) = pVp_pbas(i,j)*B(i)*B(j);
    }
  }

  V_pbas = 0;
  pVp_pbas = 0;

  // form the momentum basis hamiltonian
  RefSymmSCMatrix h_pbas(p_oso_dim, p_so_kit);
  h_pbas = AVA_pbas + BpVpB_pbas;

  // Add the kinetic energy
  for (int i=0; i<noso; i++) {
    h_pbas(i,i) = h_pbas(i,i) + Emc2(i);
  }

  if (dk_ > 1) {
    RefDiagSCMatrix p2K2 = p2*K*K;
    RefDiagSCMatrix p2K2_inv = p2K2->clone();

    for (int i=0; i<noso; i++) {
      double p2K2_val = p2K2(i);
      if (fabs(p2K2_val) > DBL_EPSILON) p2K2_inv(i) = 1.0/p2K2_val;
      else p2K2_inv(i) = 0.0;
    }

    core_hamiltonian_dk2_contrib(h_pbas, E, K, p2, p2K2, p2K2_inv,
                                 AVA_pbas, BpVpB_pbas);

    if (dk_ > 2) {
      Ref<OneBodyInt> pxVp_obi = integral()->p_cross_nuclear_p();
      Ref<SCElementOp3> hc3;
      hc3 = new OneBody3IntOp(new SymmOneBodyIntIter(pxVp_obi, p_pl));
      RefSymmSCMatrix pxVp_x_skel(p_ao_dim, p_kit);
      RefSymmSCMatrix pxVp_y_skel(p_ao_dim, p_kit);
      RefSymmSCMatrix pxVp_z_skel(p_ao_dim, p_kit);
      pxVp_x_skel.assign(0.0);
      pxVp_y_skel.assign(0.0);
      pxVp_z_skel.assign(0.0);
      pxVp_x_skel.element_op(hc3,pxVp_y_skel,pxVp_z_skel);
      RefSymmSCMatrix pxVp_x(p_so_dim, p_so_kit);
      RefSymmSCMatrix pxVp_y(p_so_dim, p_so_kit);
      RefSymmSCMatrix pxVp_z(p_so_dim, p_so_kit);
      p_pl->symmetrize(pxVp_x_skel,pxVp_x);
      p_pl->symmetrize(pxVp_y_skel,pxVp_y);
      p_pl->symmetrize(pxVp_z_skel,pxVp_z);
      pxVp_x_skel = 0;
      pxVp_y_skel = 0;
      pxVp_z_skel = 0;

      core_hamiltonian_dk3_contrib(h_pbas,
                                   E, B,
                                   p2K2_inv,
                                   so_to_p,
                                   pxVp_x);
      core_hamiltonian_dk3_contrib(h_pbas,
                                   E, B,
                                   p2K2_inv,
                                   so_to_p,
                                   pxVp_y);
      core_hamiltonian_dk3_contrib(h_pbas,
                                   E, B,
                                   p2K2_inv,
                                   so_to_p,
                                   pxVp_z);
    }

  }

#if DK_DEBUG
  h_pbas.print("h_pbas");
#endif

  AVA_pbas = 0;
  BpVpB_pbas = 0;
  A = 0;
  B = 0;
  E = 0;
  K = 0;
  Emc2 = 0;

  // Construct the transform from the momentum representation to the
  // coordinate representation in the momentum basis
  RefSCMatrix p_to_so
    = p_orthog->basis_to_orthog_basis_inverse() * Tvec;

  // Construct the transform from the momentum basis to the
  // coordinate basis.
  integral()->set_basis(bas,p_bas);
  Ref<PetiteList> pl = integral()->petite_list();
  RefSCMatrix S_ao_p(pl->AO_basisdim(), p_ao_dim, p_kit);
  S_ao_p.assign(0.0);
  hc = new OneBodyIntOp(integral()->overlap());
  S_ao_p.element_op(hc);
  hc=0;
  // convert s_ao_p into the so ao and so p basis
  RefSCMatrix blocked_S_ao_p(pl->AO_basisdim(), p_pl->AO_basisdim(), p_so_kit);
  blocked_S_ao_p->convert(S_ao_p);
  RefSCMatrix S_ao_p_so_l = pl->sotoao() * blocked_S_ao_p;
  RefSCMatrix S_ao_p_so = S_ao_p_so_l * p_pl->aotoso();
  S_ao_p_so_l = 0;

  // transform h_pbas back to the so basis
  RefSymmSCMatrix h_dk_so(pl->SO_basisdim(),bas->so_matrixkit());
  h_dk_so.assign(0.0);
  h_dk_so.accumulate_transform(S_ao_p_so
                               *p_orthog->overlap_inverse()
                               *p_to_so, h_pbas);

  // Compute the overlap in bas
  integral()->set_basis(bas);
  RefSymmSCMatrix S_bas;
  {
    Ref<SCMatrixKit> kit = bas->matrixkit();
    Ref<SCMatrixKit> so_kit = bas->so_matrixkit();
    RefSCDimension so_dim = pl->SO_basisdim();
    RefSCDimension ao_dim = pl->AO_basisdim();
    RefSymmSCMatrix S_skel(ao_dim, kit);
    S_skel.assign(0.0);
    hc = new OneBodyIntOp(new SymmOneBodyIntIter(integral()->overlap(), pl));
    S_skel.element_op(hc);
    hc=0;
    S_bas = so_kit->symmmatrix(so_dim);
    pl->symmetrize(S_skel, S_bas);
  }

  integral()->set_basis(basis());

#if DK_DEBUG
  {
    RefSCMatrix tmp = S_ao_p_so * p_orthog->overlap_inverse() * S_ao_p_so.t();
    tmp.print("S(OBS,pbasis) * S(pbasis)^-1 * S(pbasis,OBS)");
    ExEnv::out0() << indent << " trace = " << tmp.trace() << endl;
    S_bas.print("S(OBS)");

    ExEnv::out0() << indent << " trace = " << S_bas.trace() << endl;

    ExEnv::out0() << indent << "nso = " << pl->SO_basisdim()->n() << endl;
  }
#endif

  // Check to see if the momentum basis spans the coordinate basis.  The
  // following approach seems reasonable, but a more careful mathematical
  // analysis would be desirable.
  const double S_ao_projected_trace
    = (S_ao_p_so * p_orthog->overlap_inverse() * S_ao_p_so.t()).trace()
    / pl->SO_basisdim()->n();
  const double S_ao_trace = S_bas.trace() / pl->SO_basisdim()->n();
  const double completeness_diagnostic = S_ao_projected_trace / S_ao_trace;
  ExEnv::out0() << indent
                << "Tr(basis overlap projected into momentum basis)/Tr(basis overlap) = "
                << completeness_diagnostic
                << std::endl;
  if (fabs(1.0 - completeness_diagnostic)>lindep_tol_) {
    ExEnv::out0() << indent
                  << "WARNING: the momentum basis does not span the orbital basis"
                  << std::endl;
  }

#if DK_DEBUG
  S_ao_p_so.print("S_ao_p_so");
  p_to_so.print("p_to_so");
  //(p_to_so*so_to_p).print("p_to_so*so_to_p");
  (S_ao_p_so*S.gi()*p_to_so).print("S_ao_p_so*S.gi()*p_to_so");
#endif

#if DK_DEBUG
  (T+V).print("T+V");
  h_dk_so.print("h_dk_so");
#endif

  return h_dk_so;
}

RefSymmSCMatrix
Wavefunction::core_hamiltonian_for_basis(
  const Ref<GaussianBasisSet> &basis,
  const Ref<GaussianBasisSet> &p_basis)
{
  RefSymmSCMatrix hcore;

  if (dk_ > 0) {
    hcore = core_hamiltonian_dk(dk_, basis, p_basis);
  }
  else {
    hcore = core_hamiltonian_nr(basis);
  }

  return hcore;
}

RefSymmSCMatrix
Wavefunction::core_hamiltonian()
{
  if (!hcore_.computed()) {
    integral()->set_basis(gbs_);

    if (dk_ > 0) {
      hcore_ = core_hamiltonian_dk(dk_,gbs_,momentum_basis_);
    }
    else {
      hcore_ = core_hamiltonian_nr(gbs_);
    }

    hcore_.computed() = 1;
  }

  return hcore_.result_noupdate();
}

RefSymmSCMatrix
Wavefunction::core_hamiltonian_nr(const Ref<GaussianBasisSet> &bas)
{
    RefSymmSCMatrix hcore;

    integral()->set_basis(bas);

    Ref<PetiteList> pl = integral()->petite_list();

    // form skeleton Hcore in AO basis
    RefSymmSCMatrix hao(bas->basisdim(), bas->matrixkit());
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
      integral()->set_basis(bas,bas,atom_basis_);
      Ref<TwoBodyThreeCenterInt> cd_tbint
        = integral_->electron_repulsion3();
      Ref<OneBodyInt> cd_int = new ChargeDistInt(cd_tbint, atom_basis_coef_);
      hc = new OneBodyIntOp(new SymmOneBodyIntIter(cd_int,pl));
      hao.element_op(hc);
      hc=0;
      cd_int=0;
      cd_tbint=0;
    }

    // include contributions from external electric field, if needed
    RefSymmSCMatrix mu_so;
    if (electric_field()) {
      {
        RefSymmSCMatrix mu = hao.clone();
        mu.assign(0.0);

        double E[3];  for(int xyz=0; xyz<3; ++xyz) E[xyz] = electric_field().get_element(xyz);
        const int nshell = bas->nshell();
        integral()->set_basis(bas,bas);
        Ref<OneBodyInt> m1_ints = integral()->dipole(0);
        for(int sh1=0; sh1<nshell; sh1++) {
          int bf1_offset = bas->shell_to_function(sh1);
          int nbf1 = bas->shell(sh1).nfunction();

          int sh2max = sh1;
          for(int sh2=0; sh2<=sh2max; sh2++) {
            int bf2_offset = bas->shell_to_function(sh2);
            int nbf2 = bas->shell(sh2).nfunction();

            m1_ints->compute_shell(sh1,sh2);
            const double *m1intsptr = m1_ints->buffer();

            int bf1_index = bf1_offset;
            for(int bf1=0; bf1<nbf1; bf1++, bf1_index++, m1intsptr+=3*nbf2) {
              int bf2_index = bf2_offset;
              const double *ptr1 = m1intsptr;
              int bf2max;
              if (sh1 == sh2)
                bf2max = bf1;
              else
                bf2max = nbf2-1;
              for(int bf2=0; bf2<=bf2max; bf2++, bf2_index++) {

                // V = - mu \dot E; additional -1 accounts for the negative charge of the electron
                const double V = (ptr1[0] * E[0] + ptr1[1] * E[1] + ptr1[2] * E[2]);
                mu.set_element(bf1_index, bf2_index, V);
                ptr1 += 3;

              }
            }
          }
        }
        m1_ints = 0;

        const int nbasis = bas->nbasis();
        for(int bf1=0; bf1<nbasis; bf1++)
          for(int bf2=0; bf2<=bf1; bf2++) {
            mu(bf2,bf1) = mu(bf1,bf2);
          }

        mu_so = pl->to_SO_basis(mu);

      }
    }
    // now symmetrize Hso
    RefSymmSCMatrix h(pl->SO_basisdim(), bas->so_matrixkit());
    pl->symmetrize(hao,h);

    if (electric_field()) {
      h.accumulate(mu_so);
    }

    hcore = h;
    integral()->set_basis(basis());

    return hcore;
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
  // include the energy of nuclei in the presence of electric field, if needed
  double ext_efield_contribution = 0.0;
  if (electric_field()) {
    const int natoms = molecule()->natom();
    for(int a=0; a<natoms; ++a) {
      ext_efield_contribution -= electric_field()->get_element(0) * molecule()->charge(a) * molecule()->r(a, 0);
      ext_efield_contribution -= electric_field()->get_element(1) * molecule()->charge(a) * molecule()->r(a, 1);
      ext_efield_contribution -= electric_field()->get_element(2) * molecule()->charge(a) * molecule()->r(a, 2);
    }
  }

  if (atom_basis_.null()) return ext_efield_contribution + molecule()->nuclear_repulsion_energy();

  double nucrep = ext_efield_contribution;

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

Ref<GaussianBasisSet>
Wavefunction::momentum_basis() const
{
  return momentum_basis_;
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
  o << indent << "Electronic basis:" << std::endl;
  o << incindent;
  basis()->print_brief(o);
  o << decindent;
  if (atom_basis_) {
    o << indent << "Nuclear basis:" << std::endl;
    o << incindent;
    atom_basis_->print_brief(o);
    o << decindent;
  }
  if (momentum_basis_) {
    o << indent << "Momentum basis:" << std::endl;
    o << incindent;
    momentum_basis_->print_brief(o);
    o << decindent;
  }
  o << sc::indent << "Integral factory = " << const_cast<Wavefunction*>(this)->integral()->class_name() << std::endl;
  o << indent << "magnetic moment = " << magnetic_moment() << std::endl;
  // the other stuff is a wee bit too big to print
  if (print_nao_ || print_npa_) {
    Timer tim("NAO");
    RefSCMatrix naos = ((Wavefunction*)this)->nao();
    tim.exit("NAO");
    if (print_nao_) naos.print("NAO", o);
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
  magnetic_moment_ = aodim_.n() + 1;

  MolecularEnergy::obsolete();
}

void
Wavefunction::copy_orthog_info(const Ref<Wavefunction>&wfn)
{
  if (orthog_) {
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
    purge();
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

bool
Wavefunction::nonzero_efield_supported() const {
  // support efields in C1 symmetry only
  if (molecule()->point_group()->char_table().order() == 1)
    return true;
  return false;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
