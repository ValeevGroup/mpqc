
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/scf/scf.h>

//////////////////////////////////////////////////////////////////////////////

static void
nuc_repulsion(const RefSCVector& g, const RefMolecule& m)
{
  // handy things
  Molecule& mol = *m.pointer();

  for (int x=0; x < mol.natom(); x++) {
    double xyz[3];
    mol.nuclear_repulsion_1der(x, xyz);
    for (int x1=0, x3=x*3; x1 < 3; x1++,x3++)
      g.accumulate_element(x3, xyz[x1]);
  }
}

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

  // do the nuclear contribution
  gradient.assign(0.0);
  nuc_repulsion(gradient, molecule());
  //gradient.print("nuclear contribution");

  // form overlap contribution
  RefSymmSCMatrix dens = lagrangian();
  RefOneBodyDerivInt derint = integral()->overlap_deriv();
  RefSCVector obint = gradient.clone();
  obint.assign(0.0);
  ob_gradient(derint, obint, dens, basis());
  gradient.accumulate(obint);
  //dens.print("lagrangian");
  //obint.print("overlap contribution");
  
  // other one electron contributions
  obint.assign(0.0);
  dens = gradient_density();
  derint = integral()->hcore_deriv();
  ob_gradient(derint, obint, dens, basis());
  gradient.accumulate(obint);
  //dens.print("dens");
  //obint.print("one electron contribution");

  dens=0;
  derint=0;
  
  // now calculate two electron contribution
  obint.assign(0.0);
  two_body_deriv(obint);
  //obint.print("two electron contribution");
  
  gradient.accumulate(obint);
  
  obint=0;
  
  done_gradient();
}

//////////////////////////////////////////////////////////////////////////////

void
SCF::compute_hessian(const RefSymmSCMatrix& hessian)
{
}
