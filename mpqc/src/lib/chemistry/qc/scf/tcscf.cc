
#ifdef __GNUC__
#pragma implementation
#endif

#include <util/misc/newstring.h>
#include <math/optimize/diis.h>
#include <math/optimize/scextrapmat.h>
#include <chemistry/qc/scf/clscf.h>
#include <chemistry/qc/scf/tcscf.h>

///////////////////////////////////////////////////////////////////////////
// TCSCF

#define CLASSNAME TCSCF
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#define PARENTS public OneBodyWavefunction
#include <util/class/classi.h>
void *
TCSCF::_castdown(const ClassDesc*cd)
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
    fprintf(stderr,"TCSCF::occ: Warning, there's a leftover electron.\n");
    fprintf(stderr,"  total nuclear charge = %d, %d closed shells\n",Z,nd);
    fprintf(stderr,"  total charge = %d\n\n",Z-2*nd);
  }

  nd--;
}

void
TCSCF::init()
{
  occ(_mol->charges(),_ndocc);
  _density_reset_freq = 10;
  _maxiter = 100;
  _eliminate = 1;
  ckptdir = new_string("./");
  fname = new_string("this_here_thing");

  ci1 = ci2 = 1.0/sqrt(2.0);
  occa = 1.0;
  occb = 1.0;
}

TCSCF::TCSCF(StateIn& s) :
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

TCSCF::TCSCF(const RefKeyVal& keyval) :
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
  
  if (keyval->exists("ci1")) {
    ci1 = keyval->doublevalue("ci1");
    ci2 = keyval->doublevalue("ci2");
    occa = 2.0*ci1*ci1;
    occb = 2.0*ci2*ci2;
  }

  if (keyval->exists("occa")) {
    occa = keyval->doublevalue("occa");
    ci1 = sqrt(occa/2.0);
    ci2 = sqrt(1.0 - ci1*ci1);
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

TCSCF::TCSCF(const OneBodyWavefunction& obwfn) :
  OneBodyWavefunction(obwfn)
{
  init();
}

TCSCF::TCSCF(const TCSCF& tcscf) :
  OneBodyWavefunction(tcscf)
{
  _extrap = tcscf._extrap;
  _data = tcscf._data;
  _error = tcscf._error;
  _accumdih = tcscf._accumdih;
  _accumddh = tcscf._accumddh;
  _accumeffh = tcscf._accumeffh;
  _ndocc = tcscf._ndocc;
  _density_reset_freq = tcscf._density_reset_freq;
  _maxiter = tcscf._maxiter;
  _eliminate = tcscf._eliminate;

  ckptdir = new_string(tcscf.ckptdir);
  fname = new_string(tcscf.fname);
}

TCSCF::~TCSCF()
{
}

RefSCMatrix
TCSCF::eigenvectors()
{
  return _eigenvectors;
}

void
TCSCF::save_data_state(StateOut& s)
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
TCSCF::occupation(int i)
{
  if (i < _ndocc) return 2.0;
  if (i < _ndocc + 1) return occa;
  if (i < _ndocc + 1) return occb;
  return 0.0;
}

int
TCSCF::value_implemented()
{
  return 1;
}

int
TCSCF::gradient_implemented()
{
  return 1;
}

int
TCSCF::hessian_implemented()
{
  return 0;
}

void
TCSCF::print(SCostream&o)
{
  OneBodyWavefunction::print(o);
}

void
TCSCF::compute()
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

      _eigenvectors = vec;
    }

    // schmidt orthogonalize the vector
    _eigenvectors.result_noupdate()->schmidt_orthog(overlap().pointer(),
                                                    _ndocc+2);

    if (_focka.null())
      _focka = _eigenvectors.result_noupdate()->rowdim()->create_symmmatrix();
    
    if (_fockb.null())
      _fockb = _focka.clone();
    
    if (_ka.null())
      _ka = _focka.clone();
    
    if (_kb.null())
      _kb = _focka.clone();
    
    if (_fock_evals.null())
      _fock_evals = _focka->dim()->create_diagmatrix();
    
    printf("\n  TCSCF::compute: energy accuracy = %g\n\n",
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

    printf("\n  TCSCF::compute: gradient accuracy = %g\n\n",
           _gradient.desired_accuracy());

    do_gradient(gradient);
    gradient.print("cartesian gradient");
    set_gradient(gradient);

    _gradient.set_actual_accuracy(_gradient.desired_accuracy());
  }
  
  if (_hessian.needed()) {
    fprintf(stderr,"TCSCF::compute: gradient not implemented\n");
    abort();
  }
  
}


void
TCSCF::do_vector(double& eelec, double& nucrep)
{
  _gr_vector = _eigenvectors.result_noupdate();
  
  // allocate storage for the temp arrays
  RefSCMatrix nvector = _gr_vector.clone();
  
  _gr_dens = _focka.clone();
  _gr_dens.assign(0.0);
  
  _gr_opa_dens = _focka.clone();
  _gr_opa_dens.assign(0.0);
  
  _gr_opb_dens = _focka.clone();
  _gr_opb_dens.assign(0.0);
  
  _gr_hcore = _gr_dens->clone();

  // form Hcore
  _gr_hcore.assign(0.0);
  _accumdih->accum(_gr_hcore);

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
  for (int iter=0; ; iter++) {
    // form the AO basis fock matrix
    double olde=eelec;
    form_ao_fock(centers,intbuf,eelec);

    if (fabs(olde-eelec) < 1.0e-9)
      break;

    printf("iter %5d energy = %20.15f delta = %15.10g\n",
           iter,eelec+nucrep,olde-eelec);

    // now extrapolate the fock matrix
    // first we form the error matrix which is the offdiagonal blocks of
    // the MO fock matrix
    
    _focka.scale(ci1*ci1);
    RefSymmSCMatrix mofocka = _focka.clone();
    mofocka.assign(0.0);
    mofocka.accumulate_transform(_gr_vector.t(),_focka);

    _fockb.scale(ci2*ci2);
    RefSymmSCMatrix mofockb = _fockb.clone();
    mofockb.assign(0.0);
    mofockb.accumulate_transform(_gr_vector.t(),_fockb);

    _ka.scale(ci1*ci2);
    RefSymmSCMatrix moka = _ka.clone();
    moka.assign(0.0);
    moka.accumulate_transform(_gr_vector.t(),_ka);

    _kb.scale(ci1*ci2);
    RefSymmSCMatrix mokb = _kb.clone();
    mokb.assign(0.0);
    mokb.accumulate_transform(_gr_vector.t(),_kb);

    RefSymmSCMatrix efff = mofocka.copy();
    efff.accumulate(mofockb);
    
    for (int i=0; i < _ndocc; i++) {
      efff.set_element(_ndocc,i,
           mofockb.get_element(_ndocc,i)-mokb.get_element(_ndocc,i));
      efff.set_element(_ndocc+1,i,
           mofocka.get_element(_ndocc+1,i)-moka.get_element(_ndocc+1,i));
    }
                       
    efff.set_element(_ndocc+1,_ndocc,
                     mofocka.get_element(_ndocc+1,_ndocc)-
                     mofockb.get_element(_ndocc+1,_ndocc)+
                     mokb.get_element(_ndocc+1,_ndocc)-
                     moka.get_element(_ndocc+1,_ndocc));
    
    for (int i=_ndocc+2; i < basis()->nbasis(); i++) {
      efff.set_element(i,_ndocc,
           mofocka.get_element(i,_ndocc)+mokb.get_element(i,_ndocc));
      efff.set_element(i,_ndocc+1,
           mofockb.get_element(i,_ndocc+1)+moka.get_element(i,_ndocc+1));
    }
    
    // diagonalize MO fock to get MO vector
    efff.diagonalize(_fock_evals,nvector);
    efff=0;

    // transform MO vector to AO basis
    _gr_vector = _gr_vector * nvector;
    
    // and orthogonalize vector
    _gr_vector->schmidt_orthog(overlap().pointer(),basis()->nbasis());
  }
      
  _gr_vector.print("vector");
  _fock_evals.print("evals");

  _eigenvectors = _gr_vector;
  
  int_done_bounds();
  int_done_erep();
  int_done_offsets2(centers,centers,centers,centers);
  int_done_storage();
  free_centers(centers);
  free(centers);

  _gr_dens = 0;
  _gr_opa_dens = 0;
  _gr_opb_dens = 0;
  _gr_hcore = 0;
  _gr_vector = 0;
  nvector = 0;
}
