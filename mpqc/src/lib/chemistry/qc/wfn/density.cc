
#ifdef __GNUC__
#pragma implementation
#endif

extern "C" {
#include <stdio.h>
}

#include <math/scmat/local.h>
#include <math/scmat/vector3.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/wfn/density.h>

ElectronDensity::ElectronDensity(Wavefunction&wfn):
  Volume(),
  wfn_(wfn)
{
}

ElectronDensity::~ElectronDensity()
{
}

void
ElectronDensity::compute()
{
  SCVector3 r;
  get_x(r);
  // do_gradient will automatically cause the value to be computed
  if (do_gradient()) {
      double v[3];
      set_value(wfn_.density_gradient(r,v));
      SCVector3 d(v);
      set_gradient(d);
    }
  else if (do_value()) {
      set_value(wfn_.density(r));
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
  Molecule& mol = *wfn_.molecule();

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
