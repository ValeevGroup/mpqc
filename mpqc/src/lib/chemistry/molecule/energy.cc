
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
#define PARENTS virtual public NLP2
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
MolecularEnergy::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = NLP2::_castdown(cd);
  return do_castdowns(casts,cd);
}

static RefKeyVal ugly_CTOR_hack_keyval(0);
static KeyVal&
ugly_CTOR_hack_get_keyval(KeyVal&keyval)
{
  if (ugly_CTOR_hack_keyval.nonnull()) {
      fprintf(stderr,"MolecularEnergy KeyVal CTOR recursively called:"
              " this is not yet supported--aborting.\n");
      abort();
    }

  if (!keyval.exists("dimension")) {
      // put the correct dimension in the input
      RefSCDimension dim;
      if (!keyval.exists("coor")) {
          RefMolecule mol = keyval.describedclassvalue("molecule");
          dim = new LocalSCDimension(mol->natom()*3);
        }
      else {
          RefMolecularCoor coor = keyval.describedclassvalue("coor");
          dim = coor->dim();
        }
      AssignedKeyVal * assignedkeyval = new AssignedKeyVal;
      RefDescribedClass dc = dim;
      assignedkeyval->assign("dimension",dc);
      ugly_CTOR_hack_keyval = new AggregateKeyVal(*assignedkeyval,keyval);

      return *ugly_CTOR_hack_keyval;
    }
  else {
      return keyval;
    }
}
MolecularEnergy::MolecularEnergy(KeyVal&keyval):
  NLP2(ugly_CTOR_hack_get_keyval(keyval)),
  _energy(_value)
{
  ugly_CTOR_hack_keyval = 0;

  _mc  = keyval.describedclassvalue("coor");

  _mol = keyval.describedclassvalue("molecule");

  _moldim = _mol->dim_natom3();
  
  _energy.compute() = 1;
  _gradient.compute() = 0;
  _hessian.compute() = 0;

  molecule_to_x();
}

MolecularEnergy::MolecularEnergy(RefMolecule&mol):
  NLP2(mol->dim_natom3()),
  _mol(mol),
  _mc(0),
  _energy(_value)
{
  _moldim = mol->dim_natom3();
  
  _energy.compute() = 1;
  _gradient.compute() = 0;
  _hessian.compute() = 0;

  molecule_to_x();
}

MolecularEnergy::MolecularEnergy(RefMolecule&mol,RefMolecularCoor&mc):
  NLP2(mc->dim()),
  _mol(mol),
  _mc(mc),
  _energy(_value)
{
  _moldim = new LocalSCDimension(mol->natom()*3);
  
  _energy.compute() = 1;
  _gradient.compute() = 0;
  _hessian.compute() = 0;

  molecule_to_x();
}

MolecularEnergy::~MolecularEnergy()
{
}

MolecularEnergy::MolecularEnergy(StateIn&s):
  SavableState(s,class_desc_),
  NLP2(s),
  _energy(_value)
{
  _mc->restore_state(s);
}

void
MolecularEnergy::save_data_state(StateOut&s)
{
  _mc->save_state(s);
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
  if (_mc == 0) {
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
  if (_mc == 0) {
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

  if (_mc == 0) {
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
  RefSCVector cartesian(_moldim);

  if (_mc == 0) {
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
MolecularEnergy::set_x(RefSCVector&v)
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
