
#ifdef __GNUC__
#pragma implementation
#endif

extern "C" {
#include <stdio.h>
}

#include <math/scmat/local.h>
#include <math/scmat/vector3.h>
#include <chemistry/molecule/molecule.h>
#include "density.h"

ElectronDensity::ElectronDensity(Wavefunction&wfn):
  Volume(new LocalSCDimension(3)),
  _wfn(wfn)
{
}

ElectronDensity::~ElectronDensity()
{
}

void
ElectronDensity::compute()
{
  RefSCVector cv = get_x();
  RefSCVector r(dimension());
  r[0] = cv.get_element(0);
  r[1] = cv.get_element(1);
  r[2] = cv.get_element(2);
  cart_point rc;
  rc[0] = r[0]; rc[1] = r[1]; rc[2] = r[2];
  // do_gradient will automatically cause the value to be computed
  if (do_gradient()) {
      RefSCVector d(dimension());
      double v[3];
      set_value(_wfn.density_gradient(rc,v));
      d.assign(v);
      set_gradient(d);
    }
  else if (do_value()) {
      set_value(_wfn.density(rc));
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
                             RefSCVector& p1, RefSCVector& p2)
{
  Molecule& mol = *_wfn.molecule();

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
      p1[i] = p1[i] - 3.0;
      p2[i] = p2[i] + 3.0;
    }
}
