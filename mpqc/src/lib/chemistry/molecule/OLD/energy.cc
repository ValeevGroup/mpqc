
extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
}

#include <math/newmat7/newmat.h>
#include <util/keyval/keyval.h>
#include "molecule.h"
#include "energy.h"
#include "coor.h"

MolecularEnergy::MolecularEnergy(Molecule&mol):
  NLP2(mol.natom()*3),
  _mol(mol),
  _mc(0),
  _energy(fvalue),
  _gradient(grad),
  _hessian(Hessian),
  _do_energy(1),
  _do_gradient(0),
  _do_hessian(0)
{
  molecule_to_X();
}

MolecularEnergy::MolecularEnergy(Molecule&mol,MolecularCoor&mc):
  NLP2(mc.dim()),
  _mol(mol),
  _mc(&mc),
  _energy(fvalue),
  _gradient(grad),
  _hessian(Hessian),
  _do_gradient(0),
  _do_hessian(0)
{
  molecule_to_X();
}

MolecularEnergy::~MolecularEnergy()
{
}

void MolecularEnergy::failure(const char * msg)
{
  fprintf(stderr,"MolecularEnergy::failure: \"%s\"\n",msg);
  abort();
}

void MolecularEnergy::Eval()
{
  hessian(); gradient(); energy();
}

double MolecularEnergy::EvalF()
{
  return energy();
}

ColumnVector MolecularEnergy::EvalG()
{
  ColumnVector result = gradient();
  return result;
}

SymmetricMatrix MolecularEnergy::EvalH()
{
  SymmetricMatrix result = hessian();
  return result;
}

void MolecularEnergy::set_energy(double e)
{
  _energy = e;
  _have_energy = 1;
}

double MolecularEnergy::energy()
{
  if (!_have_energy) {
      int old = do_energy(1);
      compute();
      do_energy(old);
    }
  if (!_have_energy) {
      failure("could not compute energy");
    }
  return _energy;
}

void MolecularEnergy::set_gradient(ColumnVector&g)
{
  if (_mc == 0) {
      _gradient = g;
    }
  else {
      _mc->to_internal(_gradient,g);
    }
  _have_gradient = 1;
}

const ColumnVector& MolecularEnergy::gradient()
{
  if (!_have_gradient) {
      int old = do_gradient(1);
      compute();
      do_gradient(old);
    }
  if (!_have_gradient) {
      failure("could not compute gradient");
    }
  return _gradient;
}

void MolecularEnergy::set_hessian(SymmetricMatrix&h)
{
  if (_mc == 0) {
      _hessian = h;
    }
  else {
      _mc->to_internal(_hessian,h);
    }
  _have_hessian = 1;
}

const SymmetricMatrix& MolecularEnergy::hessian()
{
  if (!_have_hessian) {
      int old = do_hessian(1);
      compute();
      do_hessian(old);
    }
  if (!_have_hessian) {
      failure("could not compute hessian");
    }
  return _hessian;
}

int MolecularEnergy::do_energy()
{
  return _do_energy;
}

int MolecularEnergy::do_gradient()
{
  return _do_gradient;
}

int MolecularEnergy::do_hessian()
{
  return _do_hessian;
}

int MolecularEnergy::do_energy(int f)
{
  int old = _do_energy;
  _do_energy = f;
  return old;
}

int MolecularEnergy::do_gradient(int f)
{
  int old = _do_gradient;
  _do_gradient = f;
  return old;
}

int MolecularEnergy::do_hessian(int f)
{
  int old = _do_hessian;
  _do_hessian = f;
  return old;
}

void MolecularEnergy::x_changed()
{
  _have_energy = 0;
  _have_gradient = 0;
  _have_hessian = 0;
}

void MolecularEnergy::X_to_molecule()
{
  ColumnVector cartesian(_mol.natom()*3);

  if (_mc == 0) {
      cartesian = xc;
    }
  else {
      _mc->to_cartesian(cartesian,xc);
    }

  //printf("xc:\n");
  //Print(xc);

  //printf("cartesian:\n");
  //Print(cartesian);

  int c = 1;
  for (int i=0; i<_mol.natom(); i++) {
      _mol[i][0] = cartesian(c); c++;
      _mol[i][1] = cartesian(c); c++;
      _mol[i][2] = cartesian(c); c++;
    }
}

void MolecularEnergy::molecule_to_X()
{
  ColumnVector cartesian(_mol.natom()*3);

  int c = 1;
  for (int i=0; i<_mol.natom(); i++) {
      cartesian(c) = _mol[i][0]; c++;
      cartesian(c) = _mol[i][1]; c++;
      cartesian(c) = _mol[i][2]; c++;
    }

  if (_mc == 0) {
      xc = cartesian;
    }
  else {
      _mc->to_internal(xc,cartesian);
    }

  //printf("xc:\n");
  //Print(xc);

  //printf("cartesian:\n");
  //Print(cartesian);
}

