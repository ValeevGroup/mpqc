
#include <math/scmat/offset.h>

#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/tbint.h>

#include <chemistry/qc/scf/scf.h>

//////////////////////////////////////////////////////////////////////////////

static void
set_scale(double& coulombscale, double& exchangescale,
          int i, int j, int k, int l)
{
  double scale = 1.0;

  if ((i!=k)||(j!=l))
    scale *= 2.0;

  if (i!=j)
    scale *= 2.0;

  coulombscale = 0.5*scale;
  exchangescale = -0.25*scale;

  if (k!=l)
    coulombscale *= 2.0;

  if ((k!=l)&&(i==j))
    exchangescale *= 2.0;
}

void
SCF::compute_gradient(const RefSCVector& gradient)
{
  init_gradient();
  gradient.assign(0.0);

  // handy things
  GaussianBasisSet& gbs = *basis().pointer();

  // form one electron contribution
  RefSymmSCMatrix lag = lagrangian();
  RefSymmSCMatrix dens = gradient_density();

  RefOneBodyDerivInt overlap_der = integral()->overlap_deriv();
  RefOneBodyDerivInt hcore_der = integral()->hcore_deriv();
  for (int ish=0; ish < gbs.nshell(); ish++) {
    GaussianShell& gsi = gbs(ish);
      
    int istart = gbs.shell_to_function(ish);
    int iend = istart + gsi.nfunction();

    for (int jsh=0; jsh <= ish; jsh++) {
      GaussianShell& gsj = gbs(jsh);

      int jstart = gbs.shell_to_function(jsh);
      int jend = jstart + gsj.nfunction();

      DerivCenters cent;
      int x;

      overlap_der->compute_shell(ish,jsh,cent);
      const double *buf = overlap_der->buffer();
      int index=0;
      for (x=0; x<cent.n(); x++) {
        double dx, dy, dz;
        dx = dy = dz = 0.0;
        for (int i=istart; i < iend; i++) {
          for (int j=jstart; j < jend; j++) {
            double ewij = lag.get_element(i,j);
            dx += buf[index++] * ewij;
            dy += buf[index++] * ewij;
            dz += buf[index++] * ewij;
          }
        }
        if (ish != jsh) {
          dx *= 2.0;
          dy *= 2.0;
          dz *= 2.0;
        }
        gradient.accumulate_element(x*3+0, dx);
        gradient.accumulate_element(x*3+1, dy);
        gradient.accumulate_element(x*3+2, dz);
        if (cent.has_omitted_center()) {
          int o = cent.omitted_center();
          gradient.accumulate_element(o*3+0, -dx);
          gradient.accumulate_element(o*3+1, -dy);
          gradient.accumulate_element(o*3+2, -dz);
        }
      }

      hcore_der->compute_shell(ish,jsh,cent);
      buf = hcore_der->buffer();
      index=0;
      for (x=0; x<cent.n(); x++) {
        double dx, dy, dz;
        dx = dy = dz = 0.0;
        for (int i=istart; i < iend; i++) {
          for (int j=jstart; j < jend; j++) {
            double dij = dens.get_element(i,j);
            dx += buf[index++] * dij;
            dy += buf[index++] * dij;
            dz += buf[index++] * dij;
          }
        }
        if (ish != jsh) {
          dx *= 2.0;
          dy *= 2.0;
          dz *= 2.0;
        }
        gradient.accumulate_element(x*3+0, dx);
        gradient.accumulate_element(x*3+1, dy);
        gradient.accumulate_element(x*3+2, dz);
        if (cent.has_omitted_center()) {
          int o = cent.omitted_center();
          gradient.accumulate_element(o*3+0, -dx);
          gradient.accumulate_element(o*3+1, -dy);
          gradient.accumulate_element(o*3+2, -dz);
        }
      }
    }
  }

  lag=0;
  dens=0;
  overlap_der=0;
  hcore_der=0;
  
  gradient.print("one electron contribution");

  RefTwoBodyDerivInt tbints = integral()->electron_repulsion_deriv();
  
  // now the two electron part
  for (int i=0; i < gbs.nshell(); i++) {
    for (int j=0; j <= i; j++) {
      for (int k=0; k <= i; k++) {
        for (int l=0; l <= ((i==k)?j:k); l++) {
          DerivCenters cent;
          tbints->compute_shell(i,j,k,l,cent);
          const double * buf = tbints->buffer();
          
          double coulombscale, exchangescale;

          set_scale(coulombscale, exchangescale, i, j, k, l);
          make_gradient_contribution();
        }
      }
    }
  }
  
  // clean up
  tbints=0;
  
  done_gradient();
}

void
SCF::compute_hessian(const RefSymmSCMatrix& hessian)
{
  init_hessian();

  done_hessian();
}
