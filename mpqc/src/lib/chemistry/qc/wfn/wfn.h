
#ifndef _chemistry_qc_wfn_wfn_h
#define _chemistry_qc_wfn_wfn_h

#include <stdio.h>
#include <chemistry/molecule/energy.h>
#include <math/newmat7/newmat.h>
#include <math/topology/point.h>
#include <chemistry/qc/basis/basis.h>
class Molecule;
class MolecularCoor;
class ColumnVector;
class SymmetricMatrix;
class Matrix;

class Wavefunction: public MolecularEnergy
{
 private:
  int _have_natural_orbitals;
  int _have_natural_density;
  Matrix _natural_orbitals;
  DiagonalMatrix _natural_density;

  double* bs_values;
  double* bsg_values;
  GaussianBasisSet& _gbs;
  void init(KeyVal&);
 public:
  Wavefunction(KeyVal&,Molecule&,GaussianBasisSet&);
  Wavefunction(KeyVal&,Molecule&,GaussianBasisSet&,MolecularCoor&);
  virtual ~Wavefunction();

  void x_changed();

  void print(FILE* =stdout);
  double density(cart_point&);
  double density_gradient(cart_point&,double*);
  double natural_orbital(cart_point& r, int iorb);
  double natural_orbital_density(cart_point& r, int iorb, double* orbval = 0);
  double orbital(cart_point& r, int iorb, const Matrix& orbs);
  double orbital_density(cart_point& r,
                         int iorb,
                         const Matrix& orbs,
                         double* orbval = 0);

  virtual const SymmetricMatrix& density() = 0;
  const Matrix& natural_orbitals();
  const DiagonalMatrix& natural_density();
  const GaussianBasisSet& basis();
};

#endif
