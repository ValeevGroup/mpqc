
#ifdef __GNUC__
#pragma implementation
#endif

extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
}

#include <math/scmat/local.h>
#include <util/keyval/keyval.h>
#include <chemistry/molecule/energy.h>

SavableState_REF_def(MolecularEnergy);

#define CLASSNAME MolecularEnergy
#define PARENTS public NLP2
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
MolecularEnergy::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = NLP2::_castdown(cd);
  return do_castdowns(casts,cd);
}

MolecularEnergy::MolecularEnergy(const MolecularEnergy& mole):
  NLP2(mole),
  _energy(_value)
{
  _mc = mole._mc;
  _moldim = mole._moldim;
  _mol = mole._mol;
}

MolecularEnergy::MolecularEnergy(const RefKeyVal&keyval):
  NLP2(keyval),
  _energy(_value)
{
  _mol = keyval->describedclassvalue("molecule");

  _moldim = matrixkit()->dimension(3 * _mol->natom(), "3Natom");

  // the molecule coordinate object needs _moldim
  // so constract a keyval that has it
  RefAssignedKeyVal assignedkeyval = new AssignedKeyVal;
  RefDescribedClass dc = _moldim;
  assignedkeyval->assign("natom3", dc);
  dc = matrixkit();
  assignedkeyval->assign("matrixkit", dc);
  RefKeyVal asskeyval(assignedkeyval.pointer());
  RefKeyVal aggkeyval = new AggregateKeyVal(asskeyval, keyval);
  _mc  = aggkeyval->describedclassvalue("coor");

  RefSCDimension dim;
  if (_mc.null()) {
      dim = _moldim;
    }
  else {
      dim = _mc->dim();
    }
  set_dimension(dim);

  _energy.compute() = 1;
  _gradient.compute() = 0;
  _hessian.compute() = 0;

  molecule_to_x();
}

MolecularEnergy::~MolecularEnergy()
{
}

MolecularEnergy::MolecularEnergy(StateIn&s):
  NLP2(s),
  _energy(_value)
  maybe_SavableState(s)
{
  _mc.restore_state(s);
  _moldim.restore_state(s);
  _mol.restore_state(s);
}

MolecularEnergy&
MolecularEnergy::operator=(const MolecularEnergy& mole)
{
  NLP2::operator=(mole);
  _mc = mole._mc;
  _moldim = mole._moldim;
  _mol = mole._mol;
  return *this;
}

void
MolecularEnergy::save_data_state(StateOut&s)
{
  NLP2::save_data_state(s);
  _mc.save_state(s);
  _moldim.save_state(s);
  _mol.save_state(s);
}

void
MolecularEnergy::failure(const char * msg)
{
  fprintf(stderr,"MolecularEnergy::failure: \"%s\"\n",msg);
  abort();
}

void
MolecularEnergy::set_energy(double e)
{
  _energy.result_noupdate() = e;
  _energy.computed() = 1;
}

double
MolecularEnergy::energy()
{
  return _energy;
}

void
MolecularEnergy::set_gradient(RefSCVector&g)
{
  if (_mc.null()) {
      _gradient.result_noupdate() = g;
    }
  else {
      _mc->to_internal(_gradient.result_noupdate(),g);
    }
  _gradient.computed() = 1;
}

void
MolecularEnergy::set_hessian(RefSymmSCMatrix&h)
{
  if (_mc.null()) {
      _hessian.result_noupdate() = h;
    }
  else {
      _mc->to_internal(_hessian.result_noupdate(),h);
    }
  _hessian.computed() = 1;
}

void
MolecularEnergy::x_to_molecule()
{

  if (_mc.null()) {
      int c = 0;
      for (int i=0; i<_mol->natom(); i++) {
          _mol->operator[](i)[0] = _x(c); c++;
          _mol->operator[](i)[1] = _x(c); c++;
          _mol->operator[](i)[2] = _x(c); c++;
        }
    }
  else {
      _mc->to_cartesian(_x);
    }
}

void
MolecularEnergy::molecule_to_x()
{
  if (_mc.null()) {
      RefSCVector cartesian(_moldim);
      int c = 0;
      for (int i=0; i<_mol->natom(); i++) {
          cartesian(c) = _mol->operator[](i)[0]; c++;
          cartesian(c) = _mol->operator[](i)[1]; c++;
          cartesian(c) = _mol->operator[](i)[2]; c++;
        }
      _x = cartesian;
    }
  else {
      _mc->to_internal(_x);
    }

  obsolete();
}

void
MolecularEnergy::set_x(const RefSCVector&v)
{
  NLP0::set_x(v);
  x_to_molecule();
}

RefMolecule
MolecularEnergy::molecule()
{
  return _mol;
}

void
MolecularEnergy::guess_hessian(RefSymmSCMatrix&hessian)
{
  if (_mc.nonnull()) {
      _mc->guess_hessian(hessian);
    }
  else {
      NLP2::guess_hessian(hessian);
    }
}

RefSymmSCMatrix
MolecularEnergy::inverse_hessian(RefSymmSCMatrix&hessian)
{
  if (_mc.nonnull()) {
      return _mc->inverse_hessian(hessian);
    }
  else {
      return NLP2::inverse_hessian(hessian);
    }
}

void
MolecularEnergy::print(SCostream&o)
{
  NLP2::print(o);
  if (_mc.nonnull()) {
      o.indent(); o << "Molecular Coordinates:\n";
      o++;
      _mc->print(o);
      o--;
    }
  else {
      o.indent(); o << "Molecule:\n";
      o++;
      _mol->print(o);
      o--;
    }
}
