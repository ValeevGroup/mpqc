
#include <math/scmat/offset.h>

#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/tbint.h>

#include <chemistry/qc/scf/scf.h>

//////////////////////////////////////////////////////////////////////////////

static void
ob_gradient(const RefOneBodyDerivInt& derint, const RefSCVector& gradient,
            const RefSymmSCMatrix& density, const RefGaussianBasisSet& gbs_)
{
  GaussianBasisSet& gbs = *gbs_.pointer();
  Molecule& mol = *gbs_->molecule().pointer();
  
  for (int x=0; x < mol.natom(); x++) {
    for (int ish=0; ish < gbs.nshell(); ish++) {
      GaussianShell& gsi = gbs(ish);
      
      int istart = gbs.shell_to_function(ish);
      int iend = istart + gsi.nfunction();

      for (int jsh=0; jsh <= ish; jsh++) {
        GaussianShell& gsj = gbs(jsh);

        int jstart = gbs.shell_to_function(jsh);
        int jend = jstart + gsj.nfunction();

        derint->compute_shell(ish,jsh,x);
        const double *buf = derint->buffer();

        int index=0;
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

        gradient.accumulate_element(x*3+0, dx);
        gradient.accumulate_element(x*3+1, dy);
        gradient.accumulate_element(x*3+2, dz);
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
  //gradient.print("nuclear-nuclear term");

  // form overlap contribution
  RefSymmSCMatrix dens = lagrangian();
  RefOneBodyDerivInt derint = integral()->overlap_deriv();
  RefSCVector obint = gradient.clone();
  obint.assign(0.0);
  ob_gradient(derint, gradient, dens, basis());
  //dens.print("lag");
  //obint.print("overlap term");
  gradient.accumulate(obint);
  
  // other one electron contributions
  obint.assign(0.0);
  dens = gradient_density();
  derint = integral()->hcore_deriv();
  ob_gradient(derint, gradient, dens, basis());
  //dens.print("density");
  //obint.print("one electron term");
  gradient.accumulate(obint);

  //gradient.print("one electron contribution");

  dens=0;
  derint=0;
  
  // now calculate two electron contribution
  obint.assign(0.0);
  two_body_deriv(obint);
  //obint.print("two electron contribution");
  
  gradient.accumulate(obint);
  
  // clean up
  obint=0;
  
  done_gradient();
}

//////////////////////////////////////////////////////////////////////////////

void
SCF::compute_hessian(const RefSymmSCMatrix& hessian)
{
}
