//
// energy.cc
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

#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <math/scmat/local.h>
#include <util/keyval/keyval.h>
#include <chemistry/molecule/energy.h>

/////////////////////////////////////////////////////////////////
// MolecularEnergy

SavableState_REF_def(MolecularEnergy);

#define CLASSNAME MolecularEnergy
#define PARENTS public Function
#define VERSION 4
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
MolecularEnergy::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = Function::_castdown(cd);
  return do_castdowns(casts,cd);
}

MolecularEnergy::MolecularEnergy(const MolecularEnergy& mole):
  Function(mole)
{
  print_molecule_when_changed_ = mole.print_molecule_when_changed_;
  mc_ = mole.mc_;
  moldim_ = mole.moldim_;
  mol_ = mole.mol_;
  initial_pg_ = new PointGroup(mol_->point_group());
}

MolecularEnergy::MolecularEnergy(const RefKeyVal&keyval):
  Function(keyval,1.0e-6,1.0e-6,1.0e-4)
{
  print_molecule_when_changed_
      = keyval->booleanvalue("print_molecule_when_changed");
  if (keyval->error() != KeyVal::OK) print_molecule_when_changed_ = 1;

  mol_ = keyval->describedclassvalue("molecule");
  if (mol_.null()) {
      cerr << indent << "MolecularEnergy(Keyval): no molecule found"
           << endl;
      abort();
    }

  initial_pg_ = new PointGroup(mol_->point_group());

  hess_ = keyval->describedclassvalue("hessian");

  guesshess_ = keyval->describedclassvalue("guess_hessian");

  moldim_ = new SCDimension(3 * mol_->natom(), "3Natom");

  // the molecule coordinate object needs moldim_
  // so constract a keyval that has it
  RefAssignedKeyVal assignedkeyval = new AssignedKeyVal;
  RefDescribedClass dc = moldim_;
  assignedkeyval->assign("natom3", dc);
  dc = matrixkit();
  assignedkeyval->assign("matrixkit", dc);
  RefKeyVal asskeyval(assignedkeyval.pointer());
  RefKeyVal aggkeyval = new AggregateKeyVal(asskeyval, keyval);

  // Don't bother with internal coordinates if there is only 1 atom
  if (mol_->natom() > 1) {
      mc_  = aggkeyval->describedclassvalue("coor");
    }

  RefSCDimension dim;
  if (mc_.null()) {
      dim = moldim_;
    }
  else {
      dim = mc_->dim();
    }
  set_dimension(dim);

  do_value(1);
  do_gradient(0);
  do_hessian(0);

  molecule_to_x();
}

MolecularEnergy::~MolecularEnergy()
{
}

MolecularEnergy::MolecularEnergy(StateIn&s):
  SavableState(s),
  Function(s)
{
  mc_.restore_state(s);
  moldim_.restore_state(s);
  mol_.restore_state(s);
  if (s.version(static_class_desc()) >= 2) s.get(print_molecule_when_changed_);
  else print_molecule_when_changed_ = 1;
  if (s.version(static_class_desc()) >= 3) {
      hess_.restore_state(s);
      guesshess_.restore_state(s);
    }
  if (s.version(static_class_desc()) >= 4) initial_pg_.restore_state(s);
  else initial_pg_ = new PointGroup(mol_->point_group());
}

MolecularEnergy&
MolecularEnergy::operator=(const MolecularEnergy& mole)
{
  Function::operator=(mole);
  mc_ = mole.mc_;
  moldim_ = mole.moldim_;
  mol_ = mole.mol_;
  print_molecule_when_changed_ = mole.print_molecule_when_changed_;
  initial_pg_ = new PointGroup(mole.initial_pg_);
  return *this;
}

void
MolecularEnergy::save_data_state(StateOut&s)
{
  Function::save_data_state(s);
  mc_.save_state(s);
  moldim_.save_state(s);
  mol_.save_state(s);
  s.put(print_molecule_when_changed_);
  hess_.save_state(s);
  guesshess_.save_state(s);
  initial_pg_.save_state(s);
}

void
MolecularEnergy::failure(const char * msg)
{
  cerr << node0 << indent << "MolecularEnergy::failure: " << msg << endl;
  abort();
}

void
MolecularEnergy::set_energy(double e)
{
  set_value(e);
}

double
MolecularEnergy::energy()
{
  return value();
}

void
MolecularEnergy::set_gradient(RefSCVector&g)
{
  cartesian_gradient_ = g.copy();
  if (mc_.null()) {
    Function::set_gradient(g);
  } else {
    RefSCVector grad(dimension(), matrixkit());
    mc_->to_internal(grad,g);
    Function::set_gradient(grad);
  }
}

void
MolecularEnergy::set_hessian(RefSymmSCMatrix&h)
{
  cartesian_hessian_ = h.copy();
  if (mc_.null()) {
    Function::set_hessian(h);
  } else {
    RefSymmSCMatrix hess(dimension(), matrixkit());
    mc_->to_internal(hess,h);
    Function::set_hessian(hess);
  }
}

void
MolecularEnergy::x_to_molecule()
{
  RefSCVector x = get_x_no_copy();

  if (mc_.null()) {
    int c = 0;
    
    for (int i=0; i<mol_->natom(); i++) {
      mol_->r(i,0) = x(c); c++;
      mol_->r(i,1) = x(c); c++;
      mol_->r(i,2) = x(c); c++;
    }
  } else {
    mc_->to_cartesian(get_x_no_copy());
  }
}

void
MolecularEnergy::molecule_to_x()
{
  if (mc_.null()) {
    RefSCVector cartesian(moldim(),matrixkit());
    int c = 0;
    for (int i=0; i < mol_->natom(); i++) {
      cartesian(c) = mol_->r(i,0); c++;
      cartesian(c) = mol_->r(i,1); c++;
      cartesian(c) = mol_->r(i,2); c++;
    }
    Function::set_x(cartesian);
  } else {
    mc_->to_internal(get_x_reference());
  }
}

void
MolecularEnergy::set_x(const RefSCVector&v)
{
  Function::set_x(v);
  x_to_molecule();
  if (print_molecule_when_changed_) {
      cout << node0 << endl << indent << class_name()
           << ": changing atomic coordinates:" << endl;
      molecule()->print();
    }
}

RefSCVector
MolecularEnergy::get_cartesian_x()
{
  RefSCVector cartesian(moldim(),matrixkit());
  int c = 0;
  for (int i=0; i < mol_->natom(); i++) {
      cartesian(c) = mol_->r(i,0); c++;
      cartesian(c) = mol_->r(i,1); c++;
      cartesian(c) = mol_->r(i,2); c++;
    }
  return cartesian;
}

RefSCVector
MolecularEnergy::get_cartesian_gradient()
{
  gradient();
  if (cartesian_gradient_.null()) {
      cerr << "MolecularEnergy::get_cartesian_gradient(): "
           << "cartesian gradient not available"
           << endl;
      abort();
    }
  return cartesian_gradient_;
}

RefSymmSCMatrix
MolecularEnergy::get_cartesian_hessian()
{
  hessian();
  if (cartesian_hessian_.null()) {
      cerr << "MolecularEnergy::get_cartesian_hessian(): "
           << "cartesian hessian not available"
           << endl;
      abort();
    }
  return cartesian_hessian_;
}

RefSCDimension
MolecularEnergy::moldim() const
{
  return moldim_;
}

RefMolecule
MolecularEnergy::molecule() const
{
  return mol_;
}

void
MolecularEnergy::guess_hessian(RefSymmSCMatrix&hessian)
{
  if (guesshess_.nonnull()) {
      int nullmole = (guesshess_->energy() == 0);
      this->reference();
      if (nullmole) guesshess_->set_energy(this);
      RefSymmSCMatrix xhess = guesshess_->cartesian_hessian();
      if (nullmole) guesshess_->set_energy(0);
      this->dereference();
      if (mc_.nonnull()) {
          mc_->to_internal(hessian, xhess);
        }
      else {
          hessian.assign(xhess);
        }
    }
  else if (mc_.nonnull()) {
      mc_->guess_hessian(hessian);
    }
  else {
      Function::guess_hessian(hessian);
    }
}

RefSymmSCMatrix
MolecularEnergy::inverse_hessian(RefSymmSCMatrix&hessian)
{
  if (mc_.nonnull()) {
      return mc_->inverse_hessian(hessian);
    }
  else {
      return Function::inverse_hessian(hessian);
    }
}

RefSymmSCMatrix
MolecularEnergy::hessian()
{
  if (hess_.null()) return hessian_.result();

  if (hessian_.computed()) return hessian_.result();

  int nullmole = (hess_->energy() == 0);
  this->reference();
  if (nullmole) hess_->set_energy(this);
  RefSymmSCMatrix xhess = hess_->cartesian_hessian();
  if (nullmole) hess_->set_energy(0);
  this->dereference();
  set_hessian(xhess);
  return hessian_.result();
}

int
MolecularEnergy::hessian_implemented() const
{
  return hess_.nonnull();
}

void
MolecularEnergy::symmetry_changed()
{
  obsolete();
}

RefNonlinearTransform
MolecularEnergy::change_coordinates()
{
  if (!mc_) return 0;
  RefNonlinearTransform t = mc_->change_coordinates();
  do_change_coordinates(t);
  return t;
}

void
MolecularEnergy::print_natom_3(const RefSCVector &v,
                               const char *title, ostream&o) const
{
  int precision = 10;
  int lwidth = precision + 4;
  int n = v.n()/3;
  if (title) {
    o << node0 << indent << title << endl;
    o << node0 << incindent;
  }
  for (int i=0,ii=0; i<n; i++) {
    o << node0 << indent
      << scprintf("%4d %3s",
                  i+1,AtomInfo::symbol(molecule()->Z(i)));
    for (int j=0; j<3; j++,ii++) {
      o << node0 << scprintf(" % *.*f", lwidth,precision,double(v(ii)));
    }
    o << node0 << endl;
  }
  if (title) {
    o << node0 << decindent;
  }
  o.flush();
}

void
MolecularEnergy::print_natom_3(double **vn3,
                               const char *title, ostream&o) const
{
  int precision = 10;
  int lwidth = precision + 4;
  int n = molecule()->natom();
  if (title) {
    o << node0 << indent << title << endl;
    o << node0 << incindent;
  }
  for (int i=0; i<n; i++) {
    o << node0 << indent
      << scprintf("%4d %3s",
                  i+1,AtomInfo::symbol(molecule()->Z(i)));
    for (int j=0; j<3; j++) {
      o << node0 << scprintf(" % *.*f", lwidth,precision,double(vn3[i][j]));
    }
    o << node0 << endl;
  }
  if (title) {
    o << node0 << decindent;
  }
  o.flush();
}

void
MolecularEnergy::print(ostream&o) const
{
  Function::print(o);
  if (mc_.nonnull()) {
      o << node0 << indent << "Molecular Coordinates:\n" << incindent;
      mc_->print(o);
      o << node0 << decindent;
    }
  else {
      o << node0 << indent << "Molecule:\n" << incindent;
      mol_->print(o);
      o << node0 << decindent << endl;
    }
}

/////////////////////////////////////////////////////////////////
// SumMolecularEnergy

SavableState_REF_def(SumMolecularEnergy);
#define CLASSNAME SumMolecularEnergy
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public MolecularEnergy
#include <util/state/statei.h>
#include <util/class/classi.h>

void *
SumMolecularEnergy::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = MolecularEnergy::_castdown(cd);
  return do_castdowns(casts,cd);
}

SumMolecularEnergy::SumMolecularEnergy(const RefKeyVal &keyval):
  MolecularEnergy(keyval)
{
  n_ = keyval->count("mole");
  mole_ = new RefMolecularEnergy[n_];
  coef_ = new double[n_];
  for (int i=0; i<n_; i++) {
      mole_[i] = keyval->describedclassvalue("mole",i);
      coef_[i] = keyval->intvalue("coef",i);
      if (mole_[i].null()
          || mole_[i]->molecule()->natom() != molecule()->natom()) {
          cerr << "SumMolecularEnergy: a mole is null or has a molecule"
               << " with the wrong number of atoms" << endl;
          abort();
        }
    }
}

SumMolecularEnergy::SumMolecularEnergy(StateIn&s):
  MolecularEnergy(s)
{
  s.get(n_);
  coef_ = new double[n_];
  mole_ = new RefMolecularEnergy[n_];
  s.get_array_double(coef_,n_);
  for (int i=0; i<n_; i++) {
      mole_[i].restore_state(s);
    }
}

void
SumMolecularEnergy::save_data_state(StateOut&s)
{
  MolecularEnergy::save_data_state(s);
  s.put(n_);
  s.put_array_double(coef_,n_);
  for (int i=0; i<n_; i++) {
      mole_[i].save_state(s);
    }
}

SumMolecularEnergy::~SumMolecularEnergy()
{
  delete[] mole_;
  delete[] coef_;
}

int
SumMolecularEnergy::value_implemented() const
{
  for (int i=0; i<n_; i++) {
      if (!mole_[i]->value_implemented()) return 0;
    }
  return 1;
}

int
SumMolecularEnergy::gradient_implemented() const
{
  for (int i=0; i<n_; i++) {
      if (!mole_[i]->gradient_implemented()) return 0;
    }
  return 1;
}

int
SumMolecularEnergy::hessian_implemented() const
{
  for (int i=0; i<n_; i++) {
      if (!mole_[i]->hessian_implemented()) return 0;
    }
  return 1;
}

void
SumMolecularEnergy::set_x(const RefSCVector&v)
{
  MolecularEnergy::set_x(v);
  for (int i=0; i<n_; i++) {
      mole_[i]->set_x(v);
    }
}

void
SumMolecularEnergy::compute()
{
  int i;

  int *old_do_value = new int[n_];
  int *old_do_gradient = new int[n_];
  int *old_do_hessian = new int[n_];

  for (i=0; i<n_; i++)
      old_do_value[i] = mole_[i]->do_value(value_.compute());
  for (i=0; i<n_; i++)
      old_do_gradient[i]=mole_[i]->do_gradient(gradient_.compute());
  for (i=0; i<n_; i++)
      old_do_hessian[i] = mole_[i]->do_hessian(hessian_.compute());

  cout << node0 << indent
       << "SumMolecularEnergy: compute" << endl;

  cout << incindent;

  if (value_needed()) {
      double val = 0.0;
      for (i=0; i<n_; i++) {
          val += coef_[i] * mole_[i]->value();
        }
      cout << node0 << endl << indent
           << "SumMolecularEnergy =" << endl;
      for (i=0; i<n_; i++) {
          cout << node0 << indent
               << scprintf("  %c % 16.12f * % 16.12f",
                           (i==0?' ':'+'),
                           coef_[i], mole_[i]->value())
               << endl;
        }
      cout << node0 << indent
           << scprintf("  = % 16.12f", val) << endl;
      set_energy(val);
    }
  if (gradient_needed()) {
      RefSCVector gradientvec = matrixkit()->vector(moldim());
      gradientvec->assign(0.0);
      for (i=0; i<n_; i++)
          gradientvec.accumulate(coef_[i] * mole_[i]->gradient());
      set_gradient(gradientvec);
    }
  if (hessian_needed()) {
      RefSymmSCMatrix hessianmat = matrixkit()->symmmatrix(moldim());
      hessianmat->assign(0.0);
      for (i=0; i<n_; i++)
          hessianmat.accumulate(coef_[i] * mole_[i]->hessian());
      set_hessian(hessianmat);
    }

  cout << decindent;

  for (i=0; i<n_; i++) mole_[i]->do_value(old_do_value[i]);
  for (i=0; i<n_; i++) mole_[i]->do_gradient(old_do_gradient[i]);
  for (i=0; i<n_; i++) mole_[i]->do_hessian(old_do_hessian[i]);

  delete[] old_do_value;
  delete[] old_do_gradient;
  delete[] old_do_hessian;
}

/////////////////////////////////////////////////////////////////
// MolEnergyConvergence

SavableState_REF_def(MolEnergyConvergence);
#define CLASSNAME MolEnergyConvergence
#define VERSION 3
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public Convergence
#include <util/state/statei.h>
#include <util/class/classi.h>

void *
MolEnergyConvergence::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = Convergence::_castdown(cd);
  return do_castdowns(casts,cd);
}

MolEnergyConvergence::MolEnergyConvergence()
{
  set_defaults();
}

MolEnergyConvergence::MolEnergyConvergence(StateIn&s):
  SavableState(s),
  Convergence(s)
{
  if (s.version(static_class_desc()) >= 2) s.get(cartesian_);
  if (s.version(static_class_desc()) >= 3) mole_.restore_state(s);
}

MolEnergyConvergence::MolEnergyConvergence(const RefKeyVal&keyval)
{
  mole_ = keyval->describedclassvalue("energy");
  if (mole_.null()) {
      cerr << "MolEnergyConvergence(const RefKeyVal&keyval): "
           << "require an energy keyword of type MolecularEnergy"
           << endl;
      abort();
    }

  cartesian_ = keyval->booleanvalue("cartesian");
  if (keyval->error() != KeyVal::OK) cartesian_ = 1;

  use_max_disp_ = keyval->exists("max_disp");
  use_max_grad_ = keyval->exists("max_grad");
  use_rms_disp_ = keyval->exists("rms_disp");
  use_rms_grad_ = keyval->exists("rms_grad");
  use_graddisp_ = keyval->exists("graddisp");
  if (use_max_disp_) max_disp_ = keyval->doublevalue("max_disp");
  if (use_max_grad_) max_grad_ = keyval->doublevalue("max_grad");
  if (use_rms_disp_) rms_disp_ = keyval->doublevalue("rms_disp");
  if (use_rms_grad_) rms_grad_ = keyval->doublevalue("rms_grad");
  if (use_graddisp_) graddisp_ = keyval->doublevalue("graddisp");

  if (!use_max_disp_ && !use_max_grad_
      && !use_rms_disp_ && !use_rms_grad_
      && !use_graddisp_) {
      set_defaults();
    }
}

MolEnergyConvergence::~MolEnergyConvergence()
{
}

void
MolEnergyConvergence::save_data_state(StateOut&s)
{
  Convergence::save_data_state(s);
  s.put(cartesian_);
  mole_.save_state(s);
}

void
MolEnergyConvergence::set_defaults()
{
  use_max_disp_ = 1;
  use_max_grad_ = 1;
  use_rms_disp_ = 0;
  use_rms_grad_ = 0;
  use_graddisp_ = 1;
  max_disp_ = 1.0e-4;
  max_grad_ = 1.0e-4;
  graddisp_ = 1.0e-4;
}

void
MolEnergyConvergence::get_x(const RefFunction &f)
{
  RefMolecularEnergy m(f);
  if (cartesian_ && m.nonnull() && m->molecularcoor().nonnull()) {
      x_ = m->get_cartesian_x();
    }
  else {
      x_ = f->get_x();
    }
}


void
MolEnergyConvergence::set_nextx(const RefSCVector& x)
{
  if (cartesian_ && mole_.nonnull() && mole_->molecularcoor().nonnull()) {
      RefMolecule mol = new Molecule(*(mole_->molecule().pointer()));
      mole_->molecularcoor()->to_cartesian(mol, x);
      nextx_ = mole_->matrixkit()->vector(mole_->moldim());
      int c = 0;
      for (int i=0; i < mol->natom(); i++) {
          nextx_(c) = mol->r(i,0); c++;
          nextx_(c) = mol->r(i,1); c++;
          nextx_(c) = mol->r(i,2); c++;
        }
    }
  else if (mole_.null()) {
      // this only happens after restoring state from old versions
      // of MolEnergyConvergence
      nextx_ = 0;
    }
  else {
      nextx_ = x.copy();
    }
}

void
MolEnergyConvergence::get_grad(const RefFunction &f)
{
  RefMolecularEnergy m(f);
  if (cartesian_ && m.nonnull() && m->molecularcoor().nonnull()) {
      RefSCVector cartesian_grad = m->get_cartesian_gradient()->copy();
      if (m->molecularcoor()->nconstrained()) {
          // convert the gradient to internal coordinates and back
          // this will project out the fixed coordinates
          RefSCVector internal_grad(m->dimension(), m->matrixkit());
          m->molecularcoor()->to_internal(internal_grad,cartesian_grad);
          m->molecularcoor()->to_cartesian(cartesian_grad,internal_grad);
        }
      grad_ = cartesian_grad;
    }
  else {
      grad_ = f->gradient();
    }
}

int
MolEnergyConvergence::converged()
{
  return Convergence::converged();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
