
#ifdef __GNUC__
#pragma implementation
#endif

#include <stdlib.h>
#include <math.h>

#include <util/misc/formio.h>
#include <math/scmat/local.h>
#include <util/keyval/keyval.h>
#include <chemistry/molecule/energy.h>

/////////////////////////////////////////////////////////////////
// MolecularEnergy

SavableState_REF_def(MolecularEnergy);

#define CLASSNAME MolecularEnergy
#define PARENTS public Function
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
  mc_ = mole.mc_;
  moldim_ = mole.moldim_;
  mol_ = mole.mol_;
}

MolecularEnergy::MolecularEnergy(const RefKeyVal&keyval):
  Function(keyval)
{
  if (!keyval->exists("value_accuracy")) {
      value_.set_desired_accuracy(1.0e-6);
    }
  if (!keyval->exists("gradient_accuracy")) {
      gradient_.set_desired_accuracy(1.0e-6);
    }
  if (!keyval->exists("hessian_accuracy")) {
      hessian_.set_desired_accuracy(1.0e-4);
    }

  mol_ = keyval->describedclassvalue("molecule");

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
  mc_  = aggkeyval->describedclassvalue("coor");

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
  Function(s)
  maybe_SavableState(s)
{
  mc_.restore_state(s);
  moldim_.restore_state(s);
  mol_.restore_state(s);
}

MolecularEnergy&
MolecularEnergy::operator=(const MolecularEnergy& mole)
{
  Function::operator=(mole);
  mc_ = mole.mc_;
  moldim_ = mole.moldim_;
  mol_ = mole.mol_;
  return *this;
}

void
MolecularEnergy::save_data_state(StateOut&s)
{
  Function::save_data_state(s);
  mc_.save_state(s);
  moldim_.save_state(s);
  mol_.save_state(s);
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
      mol_->operator[](i)[0] = x(c); c++;
      mol_->operator[](i)[1] = x(c); c++;
      mol_->operator[](i)[2] = x(c); c++;
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
      cartesian(c) = mol_->operator[](i)[0]; c++;
      cartesian(c) = mol_->operator[](i)[1]; c++;
      cartesian(c) = mol_->operator[](i)[2]; c++;
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
}

RefSCVector
MolecularEnergy::get_cartesian_x()
{
  RefSCVector cartesian(moldim(),matrixkit());
  int c = 0;
  for (int i=0; i < mol_->natom(); i++) {
      cartesian(c) = mol_->operator[](i)[0]; c++;
      cartesian(c) = mol_->operator[](i)[1]; c++;
      cartesian(c) = mol_->operator[](i)[2]; c++;
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

RefSCDimension
MolecularEnergy::moldim()
{
  return moldim_;
}

RefMolecule
MolecularEnergy::molecule()
{
  return mol_;
}

void
MolecularEnergy::guess_hessian(RefSymmSCMatrix&hessian)
{
  if (mc_.nonnull()) {
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

void
MolecularEnergy::print(ostream&o)
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
// MolEnergyConvergence

SavableState_REF_def(MolEnergyConvergence);
#define CLASSNAME MolEnergyConvergence
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
  Convergence(s)
{
}

MolEnergyConvergence::MolEnergyConvergence(const RefKeyVal&keyval)
{
  cartesian_ = keyval->booleanvalue("cartesian");

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
  if (m.nonnull() && cartesian_) {
      x_ = m->get_cartesian_x();
    }
  else {
      x_ = f->get_x();
    }
}

void
MolEnergyConvergence::get_nextx(const RefFunction &f)
{
  RefMolecularEnergy m(f);
  if (m.nonnull() && cartesian_) {
      nextx_ = m->get_cartesian_x();
    }
  else {
      nextx_ = f->get_x();
    }
}

void
MolEnergyConvergence::get_grad(const RefFunction &f)
{
  RefMolecularEnergy m(f);
  if (m.nonnull() && cartesian_) {
      grad_ = m->get_cartesian_gradient();
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
// eval: (c-set-style "CLJ")
