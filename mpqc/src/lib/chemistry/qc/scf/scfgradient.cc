
#include <math/scmat/offset.h>

#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/scf/scf.h>

//////////////////////////////////////////////////////////////////////////////

static void
nuc_repulsion(double * g, const RefMolecule& m)
{
  // handy things
  Molecule& mol = *m.pointer();

  for (int x=0; x < mol.natom(); x++) {
    double xyz[3];
    mol.nuclear_repulsion_1der(x, xyz);
    for (int x1=0, x3=x*3; x1 < 3; x1++,x3++)
      g[x3] += xyz[x1];
  }
}

static void
ob_gradient(const RefOneBodyDerivInt& derint, double * gradient,
            const RefSymmSCMatrix& density, const RefGaussianBasisSet& gbs_)
{
  GaussianBasisSet& gbs = *gbs_.pointer();
  Molecule& mol = *gbs_->molecule().pointer();
  
  RefSCMatrixSubblockIter diter =
    density->local_blocks(SCMatrixSubblockIter::Read);

  for (diter->begin(); diter->ready(); diter->next()) {
    SCMatrixBlock *dblk = diter->block();
    double *ddata;
    int istart, iend;
    int jstart, jend;
    int tri=0;

    // dblk can either be an LTriBlock, LTriSubBlock, or RectSubBlock
    // I don't think it can be a RectBlock
    if (SCMatrixLTriBlock::castdown(dblk)) {
      SCMatrixLTriBlock *ldblk = SCMatrixLTriBlock::castdown(dblk);
      istart = ldblk->start; iend=ldblk->end;
      jstart = ldblk->start; jend=ldblk->end;
      ddata = ldblk->data;
      tri=1;
    } else if (SCMatrixLTriSubBlock::castdown(dblk)) {
      SCMatrixLTriSubBlock *ldblk = SCMatrixLTriSubBlock::castdown(dblk);
      istart = ldblk->istart; iend=ldblk->iend;
      jstart = ldblk->jstart; jend=ldblk->jend;
      ddata = ldblk->data;
      tri=1;
    } else if (SCMatrixRectSubBlock::castdown(dblk)) {
      SCMatrixRectSubBlock *ldblk = SCMatrixRectSubBlock::castdown(dblk);
      istart = ldblk->istart; iend=ldblk->iend;
      jstart = ldblk->jstart; jend=ldblk->jend;
      ddata = ldblk->data;
      tri=1;
    } else {
      fprintf(stderr,"ob_gradient: can't figure out what density block is\n");
      abort();
    }
    
    int ishstart = gbs.function_to_shell(istart);
    int ishend = (iend) ? gbs.function_to_shell(iend-1) : 0;

    int jshstart = gbs.function_to_shell(jstart);
    int jshend = (jend) ? gbs.function_to_shell(jend-1) : 0;
    
    int ilen = iend-istart;
    int jlen = jend-jstart;
    
    for (int ish=ishstart; ish <= ishend; ish++) {
      GaussianShell& gsi = gbs(ish);
      
      int ist = gbs.shell_to_function(ish);
      int ien = ist + gsi.nfunction();

      for (int jsh=jshstart; jsh <= (tri ? ish : jshend-1); jsh++) {
        GaussianShell& gsj = gbs(jsh);

        int jst = gbs.shell_to_function(jsh);
        int jen = jst + gsj.nfunction();

        for (int x=0; x < mol.natom(); x++) {
          derint->compute_shell(ish,jsh,x);
          const double *buf = derint->buffer();

          int index=0;
          double dx=0, dy=0, dz=0;
          for (int i=ist; i < ien; i++) {
            for (int j=jst; j < jen; j++) {
              if (i < istart || i >= iend || j < jstart || j >= jend) {
                index += 3;
              } else {
                int doff = (tri) ? ij_offset(i-istart,j-jstart) :
                                   (i-istart)*jlen + j-jstart;
                double denij = ddata[doff];
                dx += buf[index++] * denij;
                dy += buf[index++] * denij;
                dz += buf[index++] * denij;
              }
            }
          }

          if (ish != jsh) {
            dx *= 2.0;
            dy *= 2.0;
            dz *= 2.0;
          }

          gradient[x*3+0] += dx;
          gradient[x*3+1] += dy;
          gradient[x*3+2] += dz;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void
SCF::compute_gradient(const RefSCVector& gradient)
{
  int i;
  
  init_gradient();

  int n3=gradient.n();
  
  double *g = new double[n3];
  memset(g,0,sizeof(double)*gradient.n());

  // do the nuclear contribution
  nuc_repulsion(g, molecule());

  double *o = new double[n3];
  memset(o,0,sizeof(double)*gradient.n());

  // form overlap contribution
  RefSymmSCMatrix dens = lagrangian();
  RefOneBodyDerivInt derint = integral()->overlap_deriv();
  ob_gradient(derint, o, dens, basis());
  
  // other one electron contributions
  dens = gradient_density();
  derint = integral()->hcore_deriv();
  ob_gradient(derint, o, dens, basis());

  scf_grp_->sum(o, n3);
  for (i=0; i < n3; i++) g[i] += o[i];
  
  dens=0;
  derint=0;
  
  // now calculate two electron contribution
  memset(o,0,sizeof(double)*gradient.n());
  two_body_deriv(o);

  for (i=0; i < n3; i++) g[i] += o[i];
  
  gradient.assign(g);
  delete[] g;
  delete[] o;
  
  done_gradient();
}

//////////////////////////////////////////////////////////////////////////////

void
SCF::compute_hessian(const RefSymmSCMatrix& hessian)
{
}
