
#ifdef __GNUC__
#pragma implementation
#endif

#include <stdlib.h>
#include <math.h>

#include <util/misc/formio.h>
#include <math/scmat/local.h>
#include <util/keyval/keyval.h>
#include <chemistry/molecule/energy.h>

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

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
