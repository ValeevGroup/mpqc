
#include "obwfn.h"


void OneBodyWavefunction::init(KeyVal&keyval)
{
  x_changed();
}

OneBodyWavefunction::OneBodyWavefunction(KeyVal&keyval,
                           Molecule&mol,
                           GaussianBasisSet&gbs):
  Wavefunction(keyval,mol,gbs)
{
  init(keyval);
}

OneBodyWavefunction::OneBodyWavefunction(KeyVal&keyval,
                           Molecule&mol,
                           GaussianBasisSet&gbs,
                           MolecularCoor&mc):
  Wavefunction(keyval,mol,gbs,mc)
{
  init(keyval);
}

OneBodyWavefunction::~OneBodyWavefunction()
{
}

void OneBodyWavefunction::x_changed()
{
  Wavefunction::x_changed();
  _have_density = 0;
}


const SymmetricMatrix& OneBodyWavefunction::density()
{

  if (!_have_density) {
      const Matrix& vec = eigenvectors();
      Matrix ortho;
      Matrix orthoi;
      basis().ortho(ortho,orthoi);
      int nbasis = basis().nbasis();

      _density.ReDimension(nbasis);
      _density = 0.0;
      for (int k=0; k<nbasis; k++) {
          double occ = occupation(k);
          if (occ == 0.0) continue;
          for (int i=0; i<nbasis; i++) {
              for (int j=0; j<=i; j++) {
                  _density.element(i,j)
                    += occ*vec.element(i,k)*vec.element(j,k);
                }
            }
        }

      _have_density = 1;
    }

  return _density;
}

// Function for returning an orbital value at a point
double OneBodyWavefunction::orbital(cart_point& r, int iorb)
{
  return Wavefunction::orbital(r,iorb,eigenvectors());
}

// Function for returning an orbital value at a point
double OneBodyWavefunction::orbital_density(cart_point& r,
                                            int iorb,
                                            double* orbval)
{
  return Wavefunction::orbital_density(r,iorb,eigenvectors(),orbval);
}
