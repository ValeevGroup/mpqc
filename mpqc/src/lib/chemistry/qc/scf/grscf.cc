
#ifdef __GNUC__
#pragma implementation
#endif

#include <util/misc/newstring.h>
#include <math/optimize/diis.h>
#include <math/optimize/scextrapmat.h>
#include <chemistry/qc/scf/grscf.h>

///////////////////////////////////////////////////////////////////////////
// GRSCF

#define CLASSNAME GRSCF
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#define PARENTS public OneBodyWavefunction
#include <util/class/classi.h>
void *
GRSCF::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = OneBodyWavefunction::_castdown(cd);
  return do_castdowns(casts,cd);
}

static void
occ(PointBag_double *z, int &nd, int &ns)
{
  int Z=0;
  for (Pix i=z->first(); i; z->next(i)) Z += (int) z->get(i);

  nd = Z/2;
  ns = Z%2;
}

void
GRSCF::init()
{
  occ(_mol->charges(),_ndocc,_nsocc);
  _density_reset_freq = 10;
  _maxiter = 100;
  _eliminate = 1;
  ckptdir = new_string("./");
  fname = new_string("this_here_thing");
}

GRSCF::GRSCF(StateIn& s) :
  OneBodyWavefunction(s)
{
  _extrap.restore_state(s);
  _data.restore_state(s);
  _error.restore_state(s);

  _accumdih.restore_state(s);
  _accumddh.restore_state(s);
  _accumeffh.restore_state(s);

  s.get(_ndocc);
  s.get(_nsocc);
  s.get(_density_reset_freq);
  s.get(_maxiter);
  s.get(_eliminate);

  s.getstring(ckptdir);
  s.getstring(fname);
}

GRSCF::GRSCF(const RefKeyVal& keyval) :
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
  
  if (keyval->exists("ndocc"))
    _ndocc = keyval->intvalue("ndocc");

  if (keyval->exists("nsocc"))
    _nsocc = keyval->intvalue("nsocc");

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

GRSCF::GRSCF(const OneBodyWavefunction& obwfn) :
  OneBodyWavefunction(obwfn)
{
  init();
}

GRSCF::GRSCF(const GRSCF& grscf) :
  OneBodyWavefunction(grscf)
{
  _extrap = grscf._extrap;
  _data = grscf._data;
  _error = grscf._error;
  _accumdih = grscf._accumdih;
  _accumddh = grscf._accumddh;
  _accumeffh = grscf._accumeffh;
  _ndocc = grscf._ndocc;
  _nsocc = grscf._nsocc;
  _density_reset_freq = grscf._density_reset_freq;
  _maxiter = grscf._maxiter;
  _eliminate = grscf._eliminate;

  ckptdir = new_string(grscf.ckptdir);
  fname = new_string(grscf.fname);
}

GRSCF::~GRSCF()
{
}

RefSCMatrix
GRSCF::eigenvectors()
{
  return _eigenvectors;
}

void
GRSCF::save_data_state(StateOut& s)
{
  _extrap.save_state(s);
  _data.save_state(s);
  _error.save_state(s);

  _accumdih.save_state(s);
  _accumddh.save_state(s);
  _accumeffh.save_state(s);

  s.put(_ndocc);
  s.put(_nsocc);
  s.put(_density_reset_freq);
  s.put(_maxiter);
  s.put(_eliminate);

  s.putstring(ckptdir);
  s.putstring(fname);
}

double
GRSCF::occupation(int i)
{
  if (i < _ndocc) return 2.0;
  if (i < _ndocc + _nsocc) return 1.0;
  return 0.0;
}

int
GRSCF::value_implemented()
{
  return 1;
}

int
GRSCF::gradient_implemented()
{
  return 1;
}

int
GRSCF::hessian_implemented()
{
  return 0;
}

void
GRSCF::print(SCostream&o)
{
  OneBodyWavefunction::print(o);
}

void
GRSCF::compute()
{
  // hack!!!!  need a way to make sure that the basis geometry is the
  // same as that in the molecule
  _extrap=new DIIS;
  
  for (int i=0; i < molecule()->natom(); i++) {
    basis()->r(i,0) = molecule()->atom(i)[0];
    basis()->r(i,1) = molecule()->atom(i)[1];
    basis()->r(i,2) = molecule()->atom(i)[2];
  }
    
  if (_energy.needed()) {
    if (_eigenvectors.result_noupdate().null()) {
      // start from core guess
      HCoreWfn hcwfn(*this);
      RefSCMatrix vec = hcwfn.eigenvectors();

      //vec->print("guess scf vector");
      
      _eigenvectors = vec;
    }

    // schmidt orthogonalize the vector
    _eigenvectors.result_noupdate()->schmidt_orthog(overlap().pointer(),
                                                    _ndocc+_nsocc);

    if (_fock.null()) {
      _fock = _eigenvectors.result_noupdate()->rowdim()->create_symmmatrix();
    }
    
    if (_nsocc && _op_fock.null()) {
      _op_fock = _fock.clone();
    }
    
    if (_fock_evals.null()) {
      _fock_evals = _fock->dim()->create_diagmatrix();
    }
    
    double eelec,nucrep;
    do_vector(eelec,nucrep);
      
    // this will be done elsewhere eventually
    printf("  total scf energy = %20.15f\n",eelec+nucrep);

    set_energy(eelec+nucrep);
    _energy.set_actual_accuracy(_energy.desired_accuracy());
  }

  if (_gradient.needed()) {
    RefSCVector gradient = _moldim->create_vector();
    do_gradient(gradient);
    gradient.print("cartesian gradient");
    set_gradient(gradient);
    _gradient.set_actual_accuracy(_gradient.desired_accuracy());
  }
  
  if (_hessian.needed()) {
    fprintf(stderr,"GRSCF::compute: gradient not implemented\n");
    abort();
  }
  
}


void
GRSCF::do_vector(double& eelec, double& nucrep)
{
  _gr_vector = _eigenvectors.result_noupdate();
  //_gr_vector.print("start vector");
  
  // allocate storage for the temp arrays
  _gr_nvector = _gr_vector.clone();
  
  _gr_dens = _fock.clone();
  _gr_dens.assign(0.0);
  
  _gr_dens_diff = _gr_dens->clone();
  _gr_dens_diff.assign(0.0);
  
  _gr_gmat = _gr_dens->clone();
  _gr_hcore = _gr_dens->clone();

  if (_nsocc) {
    _gr_op_dens = _gr_dens->clone();
    _gr_op_dens.assign(0.0);
    _gr_op_dens_diff = _gr_dens->clone();
    _gr_op_dens_diff.assign(0.0);
    _gr_op_gmat = _gr_dens->clone();
  }
  
  // form Hcore
  _gr_hcore.assign(0.0);
  _accumdih->accum(_gr_hcore);
  //_gr_hcore.print("hcore");

  for (int iter=0; iter < _maxiter; iter++) {
    // form the density from the current vector 
    form_density(_gr_vector,
                 _gr_dens,_gr_dens_diff,
                 _gr_op_dens,_gr_op_dens_diff);
    
    //_gr_dens.print("density matrix");
    
    int ij=0;
    double delta=0;
    for (int i=0; i < _gr_dens_diff->n(); i++)
      for (int j=0; j <= i; j++,ij++)
        delta += _gr_dens_diff.get_element(i,j)*_gr_dens_diff.get_element(i,j);
    delta = sqrt(delta/ij);

    if (delta < 1.0e-8) break;

    _gr_dens->scale(2.0);
    _gr_dens->scale_diagonal(0.5);
    
    form_ao_fock(nucrep);
    //_gr_gmat.print("g matrix");
    //_fock.print("ao fock matrix");

    _gr_dens->scale(0.5);
    _gr_dens->scale_diagonal(2.0);

    eelec = scf_energy();
    printf("iter %5d energy = %20.15f delta = %15.10g\n",
           iter,eelec+nucrep,delta);

    RefSymmSCMatrix _gr_error = _fock.clone();
    _gr_error.assign(0.0);
    _gr_error.accumulate_transform(_gr_vector.t(),_fock);

    for (int i=0; i < _gr_error->n(); i++) {
      double occi = occupation(i);
      for (int j=0; j <= i; j++) {
        double occj = occupation(j);
        if (occi == occj)
          _gr_error.set_element(i,j,0.0);
      }
    }
    
    _gr_gmat.assign(0.0);
    _gr_gmat.accumulate_transform(_gr_vector,_gr_error);
    _gr_error.assign(_gr_gmat);
    
    _data = new SymmSCMatrixSCExtrapData(_fock);
    _error = new SymmSCMatrixSCExtrapError(_gr_error);
    _extrap->extrapolate(_data,_error);
    _data=0;
    _error=0;
    //_fock.print("extrap fock");

    RefSymmSCMatrix _gr_mofock = _fock.clone();
    _gr_mofock.assign(0.0);
    _gr_mofock.accumulate_transform(_gr_vector.t(),_fock);

    //_gr_mofock.print("mo fock");
    
    _gr_mofock.diagonalize(_fock_evals,_gr_nvector);
    _gr_mofock=0;
    //_gr_nvector.print("evecs");

    _gr_vector = _gr_vector * _gr_nvector;
    
    _gr_vector->schmidt_orthog(overlap().pointer(),basis()->nbasis());
    //_gr_vector.print("new vector");
  }
      
  _eigenvectors = _gr_vector;
  
  _gr_dens = 0;
  _gr_dens_diff = 0;
  _gr_op_dens = 0;
  _gr_op_dens_diff = 0;
  _gr_gmat = 0;
  _gr_op_gmat = 0;
  _gr_hcore = 0;
  _gr_vector = 0;
  _gr_nvector = 0;
}

double
GRSCF::scf_energy()
{
  RefSymmSCMatrix t = _fock.copy();
  t.accumulate(_gr_hcore);

  double eelec=0;
  for (int i=0; i < t->n(); i++) {
    for (int j=0; j < i; j++) {
      eelec += _gr_dens.get_element(i,j)*t.get_element(i,j);
    }
    eelec += 0.5*_gr_dens.get_element(i,i)*t.get_element(i,i);
  }
  return eelec;
}
