
extern "C" {
#include <stdio.h>
#include <stdlib.h>
}

#include <math.h>
#include <math/newmat7/newmat.h>
#include <math/newmat7/newmatap.h>
#include <util/keyval/keyval.h>
#include <chemistry/molecule/molecule.h>
#include "wfn.h"
#include <ostream.h>

void Wavefunction::init(KeyVal&keyval)
{
  x_changed();
}

Wavefunction::Wavefunction(KeyVal&keyval,
                           Molecule&mol,
                           GaussianBasisSet&gbs):
MolecularEnergy(mol),_gbs(gbs)
{
  init(keyval);
}

Wavefunction::Wavefunction(KeyVal&keyval,
                           Molecule&mol,
                           GaussianBasisSet&gbs,
                           MolecularCoor&mc):
MolecularEnergy(mol,mc),_gbs(gbs)
{
  init(keyval);
}

Wavefunction::~Wavefunction()
{
  if (bs_values) delete[] bs_values;
  if (bsg_values) delete[] bsg_values;
}

void Wavefunction::x_changed()
{
  MolecularEnergy::x_changed();
  _have_natural_orbitals = 0;
  _have_natural_density = 0;
}

const Matrix& Wavefunction::natural_orbitals()
{
  if (!_have_natural_orbitals) {
      //cout.width(8);
      //cout.precision(4);
      //cout.setf(ios::fixed,ios::floatfield);
      const SymmetricMatrix& dens = density();

      // transform the density into an orthogonal basis
      const GaussianBasisSet& gbs = basis();
      int nbasis = gbs.nbasis();
      Matrix ortho;
      Matrix orthoi;
      gbs.ortho(ortho,orthoi);

      SymmetricMatrix densortho(nbasis);
      densortho << orthoi.t()*dens*orthoi;

      _natural_orbitals.ReDimension(nbasis,nbasis);
      _natural_density.ReDimension(nbasis);

      EigenValues(densortho,_natural_density,_natural_orbitals);

      // _natural_orbitals is the ortho to NO basis transform
      // make _natural_orbitals the AO to the NO basis transform
      _natural_orbitals = ortho * _natural_orbitals;

      _have_natural_orbitals = 1;
      _have_natural_density = 1;
    }

  return _natural_orbitals;
}

const DiagonalMatrix& Wavefunction::natural_density()
{
  if (!_have_natural_density) {
      const SymmetricMatrix& dens = density();
      int nbasis = basis().nbasis();
      _natural_orbitals.ReDimension(nbasis,nbasis);
      _natural_density.ReDimension(nbasis);
      EigenValues(dens,_natural_density,_natural_orbitals);
      _have_natural_orbitals = 1;
      _have_natural_density = 1;
    }

  return _natural_density;
}

const GaussianBasisSet& Wavefunction::basis()
{
  return _gbs;
}
