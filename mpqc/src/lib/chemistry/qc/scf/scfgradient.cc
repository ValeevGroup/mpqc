
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

#include <chemistry/qc/basis/petite.h>

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

  obint=0;
  //dens=0;
  derint=0;
  
  RefTwoBodyDerivInt tbints = integral()->electron_repulsion_deriv();

  RefSCVector tbint = gradient.clone();
  tbint.assign(0.0);

  RefPetiteList rpl = integral()->petite_list();
  PetiteList& pl = *rpl.pointer();

  RefSCElementMaxAbs m = new SCElementMaxAbs();
  dens.element_op(m);
  double pmax = m->result();
  m=0;

  int Pmax = (int) (log(pmax)/log(2.0));
  int PPmax = (int) (log(6.0*pmax*pmax)/log(2.0));
  int threshold = (int) (log(desired_gradient_accuracy())/log(2.0));
  
  // now the two electron part
  for (int i=0; i < gbs.nshell(); i++) {
    if (!pl.in_p1(i))
      continue;
    
    int ni=gbs(i).nfunction();
    int fi=gbs.shell_to_function(i);
    
    for (int j=0; j <= i; j++) {
      int ij=i_offset(i)+j;
      if (!pl.in_p2(ij))
        continue;
      
      if (tbints->log2_shell_bound(i,j,-1,-1)+PPmax < threshold)
        continue;
      
      int nj=gbs(j).nfunction();
      int fj=gbs.shell_to_function(j);
    
      for (int k=0; k <= i; k++) {
        int nk=gbs(k).nfunction();
        int fk=gbs.shell_to_function(k);
    
        for (int l=0; l <= ((i==k)?j:k); l++) {
          if (tbints->log2_shell_bound(i,j,k,l)+PPmax < threshold)
            continue;
          
          int kl=i_offset(k)+l;
          int qijkl;
          if (!(qijkl=pl.in_p4(ij,kl,i,j,k,l)))
            continue;
          
          int nl=gbs(l).nfunction();
          int fl=gbs.shell_to_function(l);

          DerivCenters cent;
          tbints->compute_shell(i,j,k,l,cent);

          const double * buf = tbints->buffer();
          
          double coulombscale, exchangescale;

          set_scale(coulombscale, exchangescale, i, j, k, l);

          int indijkl=0;
          int ii,jj,kk,ll;
          int ip,jp,kp,lp;
          int nx=cent.n();
          if (cent.has_omitted_center()) nx--;
          for (int x=0; x < nx; x++) {
            for (int ixyz=0; ixyz < 3; ixyz++) {
              for (ip=0, ii=fi; ip < ni; ip++, ii++) {
                for (jp=0, jj=fj; jp < nj; jp++, jj++) {
                  for (kp=0, kk=fk; kp < nk; kp++, kk++) {
                    for (lp=0, ll=fl; lp < nl; lp++, ll++, indijkl++) {
                      double contrib;

                      contrib = coulombscale*buf[indijkl]*qijkl*
                                             dens.get_element(ii,jj)*
                                             dens.get_element(kk,ll);

                      tbint.accumulate_element(ixyz+cent.atom(x)*3,contrib);
                      tbint.accumulate_element(ixyz+cent.omitted_atom()*3,
                                               -contrib);

                      contrib = exchangescale*buf[indijkl]*qijkl*
                                              dens.get_element(ii,kk)*
                                              dens.get_element(jj,ll);

                      tbint.accumulate_element(ixyz+cent.atom(x)*3,contrib);
                      tbint.accumulate_element(ixyz+cent.omitted_atom()*3,
                                               -contrib);

                      if (i!=j && k!=l) {
                        contrib = exchangescale*buf[indijkl]*qijkl*
                                              dens.get_element(ii,ll)*
                                              dens.get_element(jj,kk);

                        tbint.accumulate_element(ixyz+cent.atom(x)*3,contrib);
                        tbint.accumulate_element(ixyz+cent.omitted_atom()*3,
                                                 -contrib);
                      }
                    }
                  }
                }
              }
            }

            if (cent.atom(x) == cent.omitted_atom())
              x++;
          }
        }
      }
    }
  }
  
  //tbint.print("two electron contribution");

  RefSCVector sym2ei = tbint.copy();
  CharacterTable ct = molecule()->point_group().char_table();
  SymmetryOperation so;
  
  for (int alpha=0; alpha < molecule()->natom(); alpha++) {
    for (int g=1; g < ct.order(); g++) {
      so = ct.symm_operation(g);
      int ap = pl.atom_map(alpha,g);

      sym2ei.accumulate_element(alpha*3+0,
                               tbint.get_element(ap*3+0)*so(0,0) +
                               tbint.get_element(ap*3+1)*so(1,0) +
                               tbint.get_element(ap*3+2)*so(2,0));
      sym2ei.accumulate_element(alpha*3+1,
                               tbint.get_element(ap*3+0)*so(0,1) +
                               tbint.get_element(ap*3+1)*so(1,1) +
                               tbint.get_element(ap*3+2)*so(2,1));
      sym2ei.accumulate_element(alpha*3+2,
                               tbint.get_element(ap*3+0)*so(0,2) +
                               tbint.get_element(ap*3+1)*so(1,2) +
                               tbint.get_element(ap*3+2)*so(2,2));
    }
  }
    
  sym2ei.scale(1.0/ct.order());
  //sym2ei.print("symmetrized two electron contribution");

  gradient.accumulate(sym2ei);
  
  // clean up
  tbint=0;
  tbints=0;
  sym2ei=0;
  
  done_gradient();
}

void
SCF::compute_hessian(const RefSymmSCMatrix& hessian)
{
}
