
#ifdef __GNUC__
#pragma implementation
#endif

extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
}

#include <iostream.h>
#include <util/keyval/keyval.h>
#include <math/scmat/local.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/integral/integralv2.h>
#include <chemistry/qc/wfn/wfn.h>

SavableState_REF_def(Wavefunction);

#define CLASSNAME Wavefunction
#define PARENTS public MolecularEnergy
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
Wavefunction::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = MolecularEnergy::_castdown(cd);
  return do_castdowns(casts,cd);
}

Wavefunction::Wavefunction(const Wavefunction& wfn) :
  MolecularEnergy(wfn),
  _overlap(wfn._overlap,this),
  _natural_orbitals(wfn._natural_orbitals,this),
  _natural_density(wfn._natural_density,this),
  bs_values(0),
  bsg_values(0)
{
  _basisdim = wfn._basisdim;
  _gbs = wfn._gbs;
}

Wavefunction::Wavefunction(const RefKeyVal&keyval):
  MolecularEnergy(keyval),
  _overlap(this),
  _natural_orbitals(this),
  _natural_density(this),
  bs_values(0),
  bsg_values(0)
{
  _overlap.compute() = 0;
  _natural_orbitals.compute() = 0;
  _natural_density.compute() = 0;

  _gbs
    = GaussianBasisSet::require_castdown(keyval->describedclassvalue("basis").pointer(),
                                         "Wavefunction::Wavefunction\n");

  _basisdim = new LocalSCDimension(_gbs->nbasis(), "AO-basis");
}

Wavefunction::~Wavefunction()
{
  if (bs_values) delete[] bs_values;
  if (bsg_values) delete[] bsg_values;
}

Wavefunction::Wavefunction(StateIn&s):
  MolecularEnergy(s),
  _overlap(s,this),
  _natural_orbitals(s,this),
  _natural_density(s,this),
  bs_values(0),
  bsg_values(0)
  maybe_SavableState(s)
{
  _gbs.restore_state(s);
  _basisdim.restore_state(s);
}

Wavefunction&
Wavefunction::operator=(const Wavefunction& wfn)
{
  MolecularEnergy::operator=(wfn);
  _basisdim = wfn._basisdim;
  _gbs = wfn._gbs;
  _overlap = wfn._overlap;
  _natural_orbitals = wfn._natural_orbitals;
  _natural_density = wfn._natural_density;
  return *this;
}

void
Wavefunction::save_data_state(StateOut&s)
{
  MolecularEnergy::save_data_state(s);
  _overlap.save_data_state(s);
  _natural_orbitals.save_data_state(s);
  _natural_density.save_data_state(s);
  _gbs.save_state(s);
  _basisdim.save_state(s);
}

RefSCMatrix
Wavefunction::natural_orbitals()
{
  if (!_natural_orbitals.computed()) {
      RefSymmSCMatrix dens = density();

      // transform the density into an orthogonal basis
      RefSCMatrix ortho;
      RefSCMatrix orthoi;
      basis()->ortho(ortho,orthoi);

      RefSymmSCMatrix densortho(_basisdim);
      densortho.assign(0.0);
      densortho.accumulate_transform(orthoi.t(),dens);

      RefSCMatrix natorb(_basisdim,_basisdim);
      RefDiagSCMatrix natden(_basisdim);
      _natural_orbitals = natorb;
      _natural_density = natden;

      densortho.diagonalize(_natural_density,_natural_orbitals);

      // _natural_orbitals is the ortho to NO basis transform
      // make _natural_orbitals the AO to the NO basis transform
      _natural_orbitals = ortho * _natural_orbitals;

      _natural_orbitals.computed() = 1;
      _natural_density.computed() = 1;
    }

  return _natural_orbitals;
}

RefDiagSCMatrix
Wavefunction::natural_density()
{
  if (!_natural_density.computed()) {
      RefSymmSCMatrix dens = density();

      RefSCMatrix natorb(_basisdim,_basisdim);
      RefDiagSCMatrix natden(_basisdim);
      _natural_orbitals = natorb;
      _natural_density = natden;

      dens.diagonalize(_natural_density,_natural_orbitals);

      _natural_orbitals.computed() = 1;
      _natural_density.computed() = 1;
    }

  return _natural_density;
}

RefSymmSCMatrix
Wavefunction::overlap()
{
  if (!_overlap.computed()) {
    RefSymmSCMatrix s(basis_dimension());
    s.assign(0.0);

    RefSCElementOp op = new GaussianOverlapIntv2(basis());
    s.element_op(op);
    op=0;

    _overlap = s;
    _overlap.computed() = 1;
  }

  return _overlap;
}

RefGaussianBasisSet
Wavefunction::basis()
{
  return _gbs;
}

void
Wavefunction::print(SCostream&o)
{
  MolecularEnergy::print(o);
  // the other stuff is a wee bit too big to print
}
