
#ifdef __GNUC__
#pragma implementation
#endif

#include <util/misc/newstring.h>
#include <math/optimize/diis.h>
#include <math/optimize/scextrapmat.h>
#include <chemistry/qc/scf/clscf.h>
#include <chemistry/qc/scf/ossscf.h>
#include <chemistry/qc/scf/hsosscf.h>
#include <chemistry/qc/scf/mcscf.h>

///////////////////////////////////////////////////////////////////////////
// MCSCF

#define CLASSNAME MCSCF
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#define PARENTS public OneBodyWavefunction
#include <util/class/classi.h>
void *
MCSCF::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = OneBodyWavefunction::_castdown(cd);
  return do_castdowns(casts,cd);
}

static void
occ(PointBag_double *z, int &nd)
{
  int Z=0;
  for (Pix i=z->first(); i; z->next(i)) Z += (int) z->get(i);

  nd = Z/2;
  if (Z%2) {
    fprintf(stderr,"MCSCF::occ: Warning, there's a leftover electron.\n");
    fprintf(stderr,"  total nuclear charge = %d, %d closed shells\n",Z,nd);
    fprintf(stderr,"  total charge = %d\n\n",Z-2*nd);
  }

  nd--;
}

void
MCSCF::init()
{
  occ(_mol->charges(),_ndocc);
  _density_reset_freq = 10;
  _maxiter = 100;
  _eliminate = 1;
  ckptdir = new_string("./");
  fname = new_string("this_here_thing");

  aorb = _ndocc;
  borb = _ndocc+1;
  root=2;
}

MCSCF::MCSCF(StateIn& s) :
  OneBodyWavefunction(s)
{
  _extrap.restore_state(s);
  _data.restore_state(s);
  _error.restore_state(s);

  _accumdih.restore_state(s);
  _accumddh.restore_state(s);
  _accumeffh.restore_state(s);

  s.get(_ndocc);
  s.get(_density_reset_freq);
  s.get(_maxiter);
  s.get(_eliminate);

  s.getstring(ckptdir);
  s.getstring(fname);
}

MCSCF::MCSCF(const RefKeyVal& keyval) :
  OneBodyWavefunction(keyval)
{
  init();
  
  _extrap = keyval->describedclassvalue("extrap");
  if (_extrap.null()) {
    _extrap = new DIIS;
  }

  _accumdih = keyval->describedclassvalue("accumdih");
  if (_accumdih.null()) {
    _accumdih = new AccumHCore;
  }
  _accumdih->init(basis(),molecule());
  
  _accumddh = keyval->describedclassvalue("accumddh");
  if (_accumddh.null()) {
    _accumddh = new AccumNullDDH;
  }
  _accumddh->init(basis(),molecule());
  
  _accumeffh = keyval->describedclassvalue("accumeffh");
  if (_accumeffh.null()) {
    _accumeffh = new GSGeneralEffH;
  }
  
  if (keyval->exists("aorb"))
    aorb = keyval->intvalue("aorb");

  if (keyval->exists("borb"))
    borb = keyval->intvalue("borb");

  if (keyval->exists("root"))
    root = keyval->intvalue("root");

  if (keyval->exists("ci1")) {
    ci1 = keyval->doublevalue("ci1");
    ci2 = sqrt(1.0 - ci1*ci1);
    occa = 2.0*ci1*ci1;
    occb = 2.0*ci2*ci2;
  }

  if (keyval->exists("ndocc"))
    _ndocc = keyval->intvalue("ndocc");

  if (keyval->exists("density_reset_freq"))
    _density_reset_freq = keyval->intvalue("density_reset_freq");

  if (keyval->exists("maxiter"))
    _maxiter = keyval->intvalue("maxiter");

  if (keyval->exists("eliminate"))
    _maxiter = keyval->booleanvalue("eliminate");

  if (keyval->exists("ckpt_dir")) {
    delete[] ckptdir;
    ckptdir = keyval->pcharvalue("ckpt_dir");
  }

  if (keyval->exists("filename")) {
    delete[] fname;
    fname = keyval->pcharvalue("filename");
  }
}

MCSCF::MCSCF(const OneBodyWavefunction& obwfn) :
  OneBodyWavefunction(obwfn)
{
  init();
}

MCSCF::MCSCF(const MCSCF& mcscf) :
  OneBodyWavefunction(mcscf)
{
  _extrap = mcscf._extrap;
  _data = mcscf._data;
  _error = mcscf._error;
  _accumdih = mcscf._accumdih;
  _accumddh = mcscf._accumddh;
  _accumeffh = mcscf._accumeffh;
  _ndocc = mcscf._ndocc;
  _density_reset_freq = mcscf._density_reset_freq;
  _maxiter = mcscf._maxiter;
  _eliminate = mcscf._eliminate;

  ckptdir = new_string(mcscf.ckptdir);
  fname = new_string(mcscf.fname);
}

MCSCF::~MCSCF()
{
}

RefSCMatrix
MCSCF::eigenvectors()
{
  return _eigenvectors;
}

void
MCSCF::save_data_state(StateOut& s)
{
  _extrap.save_state(s);
  _data.save_state(s);
  _error.save_state(s);

  _accumdih.save_state(s);
  _accumddh.save_state(s);
  _accumeffh.save_state(s);

  s.put(_ndocc);
  s.put(_density_reset_freq);
  s.put(_maxiter);
  s.put(_eliminate);

  s.putstring(ckptdir);
  s.putstring(fname);
}

double
MCSCF::occupation(int i)
{
  if (i < _ndocc) return 2.0;
  if (i == aorb) return occa;
  if (i == borb) return occb;
  return 0.0;
}

int
MCSCF::value_implemented()
{
  return 1;
}

int
MCSCF::gradient_implemented()
{
  return 1;
}

int
MCSCF::hessian_implemented()
{
  return 0;
}

void
MCSCF::print(SCostream&o)
{
  OneBodyWavefunction::print(o);
}

void
MCSCF::compute()
{
  // hack!!!!  need a way to make sure that the basis geometry is the
  // same as that in the molecule, also need a member in diis to reset it
  _accumeffh->docc(0,_ndocc);
  _accumeffh->socc(_ndocc,_ndocc+2);
  
  _extrap=new DIIS;
  
  if (_hessian.needed())
    set_desired_gradient_accuracy(desired_hessian_accuracy()/100.0);

  if (_gradient.needed())
    set_desired_value_accuracy(desired_gradient_accuracy()/100.0);

  if (_energy.needed()) {
    if (_eigenvectors.result_noupdate().null()) {
      // start from core guess
      CLSCF hcwfn(*this);
      RefSCMatrix vec = hcwfn.eigenvectors();

      // schmidt orthogonalize the vector
      vec->schmidt_orthog(overlap().pointer(),_ndocc+2);

      _eigenvectors = vec;
    }

    if (_fockc.null())
      _fockc = _eigenvectors.result_noupdate()->rowdim()->create_symmmatrix();
    
    if (_focka.null())
      _focka = _fockc.clone();
    
    if (_fockb.null())
      _fockb = _fockc.clone();
    
    if (_fockab.null())
      _fockab = _fockc.clone();
    
    if (_ka.null())
      _ka = _fockc.clone();
    
    if (_kb.null())
      _kb = _fockc.clone();
    
    RefSCDimension actived = matrixkit()->dimension(basis()->nbasis()-_ndocc);

    if (_fock_evals.null())
      _fock_evals = _focka->dim()->create_diagmatrix();
    
    printf("\n  MCSCF::compute: energy accuracy = %g\n\n",
           _energy.desired_accuracy());

    double eelec,nucrep;
    do_vector(eelec,nucrep);
      
    // this will be done elsewhere eventually
    printf("  total scf energy = %20.15f\n",eelec+nucrep);

    set_energy(eelec+nucrep);
    _energy.set_actual_accuracy(_energy.desired_accuracy());
  }

  if (_gradient.needed()) {
    RefSCVector gradient = _moldim->create_vector();

    printf("\n  MCSCF::compute: gradient accuracy = %g\n\n",
           _gradient.desired_accuracy());

    do_gradient(gradient);
    gradient.print("cartesian gradient");
    set_gradient(gradient);

    _gradient.set_actual_accuracy(_gradient.desired_accuracy());
  }
  
  if (_hessian.needed()) {
    fprintf(stderr,"MCSCF::compute: gradient not implemented\n");
    abort();
  }
  
}


void
MCSCF::do_vector(double& eelec, double& nucrep)
{
  _gr_vector = _eigenvectors.result_noupdate();
  
  // allocate storage for the temp arrays
  RefSCMatrix nvector = _gr_vector.clone();
  
  _densc = _focka.clone();
  _densc.assign(0.0);
  
  _densa = _focka.clone();
  _densa.assign(0.0);
  
  _densb = _focka.clone();
  _densb.assign(0.0);
  
  _densab = _focka.clone();
  _densab.assign(0.0);
  
  _densab2 = _focka.clone();
  _densab2.assign(0.0);
  
  _gr_hcore = _focka->clone();

  // form Hcore
  _gr_hcore.assign(0.0);
  _accumdih->accum(_gr_hcore);

  // we need the overlap down below, so let's make it before we start the
  // two-electron junk
  RefSymmSCMatrix ovlp = overlap();
  
  // initialize some junk
  centers_t *centers = basis()->convert_to_centers_t();
  if (!centers) {
    fprintf(stderr,"hoot man!  no centers\n");
    abort();
  }

  int_normalize_centers(centers);
  int_initialize_offsets2(centers,centers,centers,centers);

  nucrep = int_nuclear_repulsion(centers,centers);
  
  int flags = INT_EREP|INT_NOSTRB|INT_NOSTR1|INT_NOSTR2;
  double *intbuf = 
    int_initialize_erep(flags,0,centers,centers,centers,centers);

  int_storage(1000000);

  int_init_bounds();

  eelec=0;

  int nbasis = basis()->nbasis();

  for (int iter=0; ; iter++) {
    // form the AO basis fock matrix
    double olde=eelec;
    form_ao_fock(centers,intbuf,eelec);

    if (fabs(olde-eelec) < 1.0e-13)
      break;

    printf("iter %5d energy = %20.15f delta = %15.10g\n",
           iter,eelec+nucrep,olde-eelec);

    RefSymmSCMatrix fa = _gr_hcore.copy();
    fa.accumulate(_fockc);
    fa.accumulate(_focka);

    RefSymmSCMatrix fb = _gr_hcore.copy();
    fb.accumulate(_fockc);
    fb.accumulate(_fockb);
    
    RefSymmSCMatrix f1 = fa.copy();
    f1.scale(ci1*ci1 + 0.5*ci2*ci2);

    RefSymmSCMatrix f2 = fb.copy();
    f2.scale(ci3*ci3 + 0.5*ci2*ci2);

    RefSymmSCMatrix feff = f1.copy();
    feff.accumulate(f2);

    RefSymmSCMatrix fooa = _kb.copy();
    fooa.scale(3.0);
    fooa.accumulate(_fockb);
    fooa.scale(-1.0);
    fooa.accumulate(_focka);
    fooa.accumulate(_ka);
    fooa.scale(0.25);
    
    RefSymmSCMatrix foob = _ka.copy();
    foob.scale(3.0);
    foob.accumulate(_focka);
    foob.scale(-1.0);
    foob.accumulate(_fockb);
    foob.accumulate(_kb);
    foob.scale(0.25);
    
    RefSymmSCMatrix mo1 = f1.clone();
    mo1.assign(0.0);
    mo1.accumulate_transform(_gr_vector.t(),f1);
    
    RefSymmSCMatrix mo2 = f2.clone();
    mo2.assign(0.0);
    mo2.accumulate_transform(_gr_vector.t(),f2);
    
    RefSymmSCMatrix moa = fooa.clone();
    moa.assign(0.0);
    moa.accumulate_transform(_gr_vector.t(),fooa);
    
    RefSymmSCMatrix mob = foob.clone();
    mob.assign(0.0);
    mob.accumulate_transform(_gr_vector.t(),foob);
    
    RefSymmSCMatrix koa = foob.clone();
    koa.assign(0.0);
    koa.accumulate_transform(_gr_vector.t(),_ka);
    koa.scale(2*ci1*ci3);
    
    RefSymmSCMatrix kob = foob.clone();
    kob.assign(0.0);
    kob.accumulate_transform(_gr_vector.t(),_kb);
    kob.scale(2*ci1*ci3);
    
    RefSymmSCMatrix mof = feff.clone();
    mof.assign(0.0);
    mof.accumulate_transform(_gr_vector.t(),feff);

    mo1.scale(2.0);
    mo2.scale(2.0);
    moa.scale(2.0*ci2*ci2);
    mob.scale(2.0*ci2*ci2);
    
    for (int j=0; j < _ndocc; j++) {
      mof.set_element(aorb,j,mo2.get_element(aorb,j)+
                             moa.get_element(aorb,j)-
                             kob.get_element(aorb,j)
                             );
      mof.set_element(borb,j,mo1.get_element(borb,j)+
                             mob.get_element(borb,j)-
                             koa.get_element(borb,j)
                             );
    }

    for (int j=borb+1; j < nbasis; j++) {
      mof.set_element(j,aorb,mo1.get_element(j,aorb)-
                             moa.get_element(j,aorb)+
                             kob.get_element(j,aorb)
                             );
      mof.set_element(j,borb,mo2.get_element(j,borb)-
                             mob.get_element(j,borb)+
                             koa.get_element(j,borb)
                             );
    }
    
    mof.diagonalize(_fock_evals,nvector);
    _gr_vector = _gr_vector*nvector;
    
    // and orthogonalize vector
    _gr_vector->schmidt_orthog(ovlp.pointer(),basis()->nbasis());
  
  }
      
  _gr_vector->print("converged vector");
  
  double s = (ci1+ci3)/sqrt((ci1-ci3)*(ci1-ci3) + 2*ci2*ci2);
  if (ci1-ci3 < 0)
    s=-s;
  
  printf("ci1 = %lf ci2 = %lf ci3 = %lf\n",ci1,ci2,ci3);

  double c1 = (s+1)/sqrt(2*(1+s*s));
  double c2 = (s-1)/sqrt(2*(1+s*s));

  double theta;
  double tan2theta = sqrt(2.0)*ci2;
  double denon = ci1-ci3;
  if (fabs(denon) < 1.0e-15) {
    theta = 0.25 * 3.1415292654;
  } else {
    double ttheta = atan(tan2theta/denon);
    theta = ttheta*0.5;
  }
  
  RefSCVector ca = _gr_vector.get_column(aorb);
  RefSCVector cb = _gr_vector.get_column(borb);
  
  for (int i=0; i < nbasis; i++) {
    double u = (cos(theta)*ca.get_element(i)+sin(theta)*cb.get_element(i));
    double v = (sin(theta)*ca.get_element(i)-cos(theta)*cb.get_element(i));
    ca.set_element(i,u);
    cb.set_element(i,v);
  }

  _gr_vector.assign_column(ca,aorb);
  _gr_vector.assign_column(cb,borb);
  _gr_vector.print("tcscf vector");
  
  
  double fooe;
  form_ao_fock(centers,intbuf,fooe);

  occa = 2.0*c1*c1;
  occb = 2.0*c2*c2;

  _focka.accumulate(_fockc);
  _focka.accumulate(_gr_hcore);
  _fockb.accumulate(_fockc);
  _fockb.accumulate(_gr_hcore);
  
  printf("s = %lf, ci1 = %lf, ci2 = %lf\n",s,c1,c2);
  printf("occa = %lf, occb = %lf\n",occa,occb);

  _eigenvectors = _gr_vector;
  
  int_done_erep();
  int_done_offsets2(centers,centers,centers,centers);

  free_centers(centers);
  free(centers);

  _densc = 0;
  _densa = 0;
  _densb = 0;
  _densab = 0;
  _densab2 = 0;
  _gr_hcore = 0;
  _gr_vector = 0;
  nvector = 0;
}
