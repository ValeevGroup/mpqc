
#ifdef __GNUC__
#pragma implementation
#endif

#include "obwfn.h"

#define CLASSNAME OneBodyWavefunction
#define PARENTS public Wavefunction
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
OneBodyWavefunction::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = Wavefunction::_castdown(cd);
  return do_castdowns(casts,cd);
}

OneBodyWavefunction::OneBodyWavefunction(KeyVal&keyval):
  Wavefunction(keyval),
  _density(this)
{
}

OneBodyWavefunction::~OneBodyWavefunction()
{
}

OneBodyWavefunction::OneBodyWavefunction(StateIn&s):
  SavableState(s,OneBodyWavefunction::class_desc_),
  Wavefunction(s),
  _density(this)
{
  abort();
}

void
OneBodyWavefunction::save_data_state(StateOut&s)
{
  Wavefunction::save_data_state(s);
  abort();
}

double OneBodyWavefunction::density(cart_point&c)
{
  return Wavefunction::density(c);
}

RefSymmSCMatrix
OneBodyWavefunction::density()
{

  if (!_density.computed()) {
      RefSCMatrix vec = eigenvectors();
      RefSCMatrix ortho(dimension(),dimension());
      RefSCMatrix orthoi(dimension(),dimension());
      basis()->ortho(ortho,orthoi);
      int nbasis = basis()->nbasis();

      RefSymmSCMatrix newdensity(dimension());
      _density = newdensity;
      newdensity.assign(0.0);
      for (int k=0; k<nbasis; k++) {
          double occ = occupation(k);
          if (occ == 0.0) continue;
          for (int i=0; i<nbasis; i++) {
              for (int j=0; j<=i; j++) {
                  newdensity.set_element(i,j,
                             newdensity.get_element(i,j)
                             + occ*vec.get_element(i,k)*vec.get_element(j,k));
                }
            }
        }

      _density.computed() = 1;
    }

  return _density;
}

// Function for returning an orbital value at a point
double
OneBodyWavefunction::orbital(cart_point& r, int iorb)
{
  return Wavefunction::orbital(r,iorb,eigenvectors());
}

// Function for returning an orbital value at a point
double
OneBodyWavefunction::orbital_density(cart_point& r,
                                            int iorb,
                                            double* orbval)
{
  return Wavefunction::orbital_density(r,iorb,eigenvectors(),orbval);
}

void
OneBodyWavefunction::print(SCostream&o)
{
  Wavefunction::print(o);
}
