
#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>

#include <math/scmat/offset.h>
#include <math/scmat/blkiter.h>

#include <chemistry/qc/scf/scfden.h>

SCFDensity::SCFDensity(SCF *s, const RefSCMatrix&v, double o)
  : scf_(s), vec(v), occ(o)
{
}

SCFDensity::~SCFDensity()
{
}

int
SCFDensity::has_side_effects()
{
  return 1;
}

void
SCFDensity::set_occ(double o)
{
  occ=o;
}

void
SCFDensity::process(SCMatrixBlockIter& bi)
{
  int ir=0;

  RefSCMatrix vir=vec;
  
  if (BlockedSCMatrix::castdown(vec)) {
    ir=current_block();
    vir = BlockedSCMatrix::castdown(vec)->block(ir);
  }

  int nbasis=vir.ncol();

  if (!nbasis)
    return;
      
  double *ck = new double[nbasis];

  // loop over columns of the scf vector
  for (int k=0;  k < nbasis; k++) {
    double occk = scf_->occupation(ir,k);
    if (fabs(occk) < 1.0e-8)
      break;
    
    if (fabs(occk-occ) > 1.0e-8)
      continue;
    
    RefSCVector rck = vir.get_column(k);
    rck->convert(ck);

    for (bi.reset(); bi; bi++) {
      bi.set(bi.get() + ck[bi.i()]*ck[bi.j()]);
    }
  }

  delete[] ck;
}

//////////////////////////////////////////////////////////////////////////////

SCFEnergy::SCFEnergy()
  : eelec(0), deferred_(0)
{
}

SCFEnergy::~SCFEnergy()
{
}

int
SCFEnergy::has_collect()
{
  return 1;
}

void
SCFEnergy::defer_collect(int h)
{
  deferred_=h;
}

void
SCFEnergy::collect(const RefMessageGrp&grp)
{
  if (!deferred_)
    grp->sum(eelec);
}

double
SCFEnergy::result()
{
  return eelec;
}

void
SCFEnergy::reset()
{
  eelec=0.0;
}

void
SCFEnergy::process(SCMatrixBlockIter&i, SCMatrixBlockIter&j)
{
  for (i.reset(), j.reset(); i && j; i++, j++) {
    int ii=i.i(); int jj=j.j();
    eelec += (ii==jj) ? 0.5*j.get()*i.get() : i.get()*j.get();
  }
}

//////////////////////////////////////////////////////////////////////////////

LevelShift::LevelShift(SCF *s) :
  scf_(s)
{
  shift=0.0;
}

LevelShift::~LevelShift()
{
}

int
LevelShift::has_side_effects()
{
  return 1;
}

void
LevelShift::set_shift(double s)
{
  shift=s;
}

void
LevelShift::process(SCMatrixBlockIter& i)
{
  int ir=current_block();
  for (i.reset(); i; i++) {
    if (i.i() != i.j())
      continue;
    
    double occi = scf_->occupation(ir,i.i());
    
    if (occi==2.0)
      i.set(i.get()-shift);
    else if (occi==1.0)
      i.set(i.get()-0.5*shift);
  }
}

//////////////////////////////////////////////////////////////////////////////

char *
init_pmax(const RefGaussianBasisSet& gbs_, double *pmat_data)
{
  double l2inv = 1.0/log(2.0);
  double tol = pow(2.0,-126.0);
  
  GaussianBasisSet& gbs = *gbs_.pointer();
  
  char * pmax = new char[i_offset(gbs.nshell())];

  int ish, jsh, ij;
  for (ish=ij=0; ish < gbs.nshell(); ish++) {
    int istart = gbs.shell_to_function(ish);
    int iend = istart + gbs(ish).nfunction();
    
    for (jsh=0; jsh <= ish; jsh++,ij++) {
      int jstart = gbs.shell_to_function(jsh);
      int jend = jstart + gbs(jsh).nfunction();
      
      double maxp=0, tmp;

      for (int i=istart; i < iend; i++) {
        int ijoff = i_offset(i);
        for (int j=jstart; j < ((ish==jsh) ? i+1 : jend); j++,ijoff++)
          if ((tmp=fabs(pmat_data[ijoff])) > maxp)
            maxp=tmp;
      }

      if (maxp <= tol)
        maxp=tol;

      pmax[ij] = (signed char) (log(maxp)*l2inv);
    }
  }

  return pmax;
}
