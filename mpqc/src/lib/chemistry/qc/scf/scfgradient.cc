
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

  RefOneBodyDerivInt obints = integral()->deriv();
  for (int x=0; x < molecule()->natom(); x++) {
    for (int ish=0; ish < gbs.nshell(); ish++) {
      GaussianShell& gsi = gbs(ish);
      
      int istart = gbs.shell_to_function(ish);
      int iend = istart + gsi.nfunction();

      for (int jsh=0; jsh <= ish; jsh++) {
        GaussianShell& gsj = gbs(jsh);

        int jstart = gbs.shell_to_function(jsh);
        int jend = jstart + gsj.nfunction();

        obints->compute_overlap_shell(x,ish,jsh);

        const double *buf = obints->buffer();
        
        int index=0;
        double dx=0, dy=0, dz=0;
        for (int i=istart; i < iend; i++) {
          for (int j=jstart; j < jend; j++) {
            double ewij = lag.get_element(i,j);
            
            dx += buf[index++] * ewij;
            dy += buf[index++] * ewij;
            dz += buf[index++] * ewij;
          }
        }

        obints->compute_hcore_shell(x,ish,jsh);
        buf = obints->buffer();
        
        index=0;
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
      }
    }
  }

  lag=0;
  dens=0;
  obints=0;
  
  gradient.print("one electron contribution");

  RefTwoBodyDerivInt tbints = integral()->electron_repulsion_deriv();
  
  // now the two electron part
  for (int i=0; i < gbs.nshell(); i++) {
    for (int j=0; j <= i; j++) {
      for (int k=0; k <= i; k++) {
        for (int l=0; l <= ((i==k)?j:k); l++) {
          tbints->compute_shell(i,j,k,l);
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
