
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

//////////////////////////////////////////////////////////////////////////////

static void
ob_gradient(const RefOneBodyDerivInt& derint, const RefSCVector& gradient,
            const RefSymmSCMatrix& density, const RefGaussianBasisSet& gbs_)
{
  GaussianBasisSet& gbs = *gbs_.pointer();
  
  for (int ish=0; ish < gbs.nshell(); ish++) {
    GaussianShell& gsi = gbs(ish);
      
    int istart = gbs.shell_to_function(ish);
    int iend = istart + gsi.nfunction();

    for (int jsh=0; jsh <= ish; jsh++) {
      GaussianShell& gsj = gbs(jsh);

      int jstart = gbs.shell_to_function(jsh);
      int jend = jstart + gsj.nfunction();

      DerivCenters cent;

      derint->compute_shell(ish,jsh,cent);
      const double *buf = derint->buffer();

      int index=0;
      for (int x=0; x < cent.n(); x++) {
        double dx=0, dy=0, dz=0;
        for (int i=istart; i < iend; i++) {
          for (int j=jstart; j < jend; j++) {
            double ewij = density.get_element(i,j);
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
        gradient.accumulate_element(cent.atom(x)*3+0, dx);
        gradient.accumulate_element(cent.atom(x)*3+1, dy);
        gradient.accumulate_element(cent.atom(x)*3+2, dz);
        if (cent.has_omitted_center()) {
          int o = cent.omitted_atom();
          if (cent.atom(x) != o) {
            gradient.accumulate_element(o*3+0, -dx);
            gradient.accumulate_element(o*3+1, -dy);
            gradient.accumulate_element(o*3+2, -dz);
          }
          x++;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void
SCF::compute_gradient(const RefSCVector& gradient)
{
  init_gradient();

  // handy things
  GaussianBasisSet& gbs = *basis().pointer();
  Molecule& mol = *molecule().pointer();

  // do the nuclear contribution
  gradient.assign(0.0);
  double xyz[3];
  for (int x=0; x < mol.natom(); x++) {
    mol.nuclear_repulsion_1der(x, xyz);
    for (int x1=0; x1 < 3; x1++)
      gradient.accumulate_element(x*3+x1,xyz[x1]);
  }
  gradient.print("nuclear-nuclear term");

  // form overlap contribution
  gradient.assign(0.0);
  RefSymmSCMatrix dens = lagrangian();
  dens.print("lag");
  RefOneBodyDerivInt derint = integral()->overlap_deriv();
  ob_gradient(derint, gradient, dens, basis());
  gradient.print("overlap term");
  
  // other one electron contributions
  gradient.assign(0.0);
  dens = gradient_density();
  dens.print("density");
  derint = integral()->hcore_deriv();
  ob_gradient(derint, gradient, dens, basis());

  dens=0;
  derint=0;
  
  gradient.print("one electron contribution");

#if 0
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
#endif
}

void
SCF::compute_hessian(const RefSymmSCMatrix& hessian)
{
}
