
extern "C" {
#include <stdio.h>
}

#include <chemistry/molecule/molecule.h>
#include "density.h"

ElectronDensity::ElectronDensity(Wavefunction&wfn):
  Volume(3),
  _wfn(wfn)
{
}

ElectronDensity::~ElectronDensity()
{
}

void
ElectronDensity::compute()
{
  const ColumnVector& cv = NLP0::GetX();
  Point3 r;
  r[0] = cv.element(0);
  r[1] = cv.element(1);
  r[2] = cv.element(2);
  // do_gradient will automatically cause the value to be computed
  if (do_gradient()) {
      DVector d(3);
      set_value(_wfn.density_gradient(r,d.pointer()));
      set_gradient(d);
    }
  else if (do_value()) {
      set_value(_wfn.density(r));
    }
  if (do_hessian()) {
      fprintf(stderr,"ElectronDensity::compute(): "
              " isn't yet implemented\n");
      abort();
    }
}

// make a wild guess about the bounding box
void
ElectronDensity::boundingbox(double valuemin,
                             double valuemax,
                             Point& p1, Point& p2)
{
  Molecule& mol = _wfn.molecule();

  if (mol.natom() == 0) {
      for (int i=0; i<3; i++) p1[i] = p2[i] = 0.0;
    }

  int i;
  for (i=0; i<3; i++) p1[i] = p2[i] = mol[0][i];
  for (i=1; i<mol.natom(); i++) {
      for (int j=0; j<3; j++) {
          if (mol[i][j] < p1[i]) p1[i] = mol[i][j];
          if (mol[i][j] > p2[i]) p2[i] = mol[i][j];
        }
    }
  for (i=0; i<3; i++) {
      p1[i] -= 3.0;
      p2[i] += 3.0;
    }
}
