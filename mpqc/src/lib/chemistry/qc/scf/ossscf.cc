
#ifdef __GNUC__
#pragma implementation
#endif

#include <util/misc/newstring.h>
#include <math/optimize/diis.h>
#include <math/optimize/scextrapmat.h>
#include <chemistry/qc/scf/clscf.h>
#include <chemistry/qc/scf/ossscf.h>

///////////////////////////////////////////////////////////////////////////
// OSSSCF

#define CLASSNAME OSSSCF
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#define PARENTS public OneBodyWavefunction
#include <util/class/classi.h>
void *
OSSSCF::_castdown(const ClassDesc*cd)
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
    fprintf(stderr,"OSSSCF::occ: Warning, there's a leftover electron.\n");
    fprintf(stderr,"  total nuclear charge = %d, %d closed shells\n",Z,nd);
    fprintf(stderr,"  total charge = %d\n\n",Z-2*nd);
  }

  nd--;
}

void
OSSSCF::init()
{
  occ(_mol->charges(),_ndocc);
  _density_reset_freq = 10;
  _maxiter = 100;
  _eliminate = 1;
  ckptdir = new_string("./");
  fname = new_string("this_here_thing");
}

OSSSCF::OSSSCF(StateIn& s) :
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

OSSSCF::OSSSCF(const RefKeyVal& keyval) :
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

OSSSCF::OSSSCF(const OneBodyWavefunction& obwfn) :
  OneBodyWavefunction(obwfn)
{
  init();
  _accumeffh = new GSGeneralEffH;
  _accumddh = new AccumNullDDH;
  _accumdih = new AccumHCore;
  _accumdih->init(basis(),molecule());
  _extrap = new DIIS;
}

OSSSCF::OSSSCF(const OSSSCF& ossscf) :
  OneBodyWavefunction(ossscf)
{
  _extrap = ossscf._extrap;
  _data = ossscf._data;
  _error = ossscf._error;
  _accumdih = ossscf._accumdih;
  _accumddh = ossscf._accumddh;
  _accumeffh = ossscf._accumeffh;
  _ndocc = ossscf._ndocc;
  _density_reset_freq = ossscf._density_reset_freq;
  _maxiter = ossscf._maxiter;
  _eliminate = ossscf._eliminate;

  ckptdir = new_string(ossscf.ckptdir);
  fname = new_string(ossscf.fname);
}

OSSSCF::~OSSSCF()
{
}

RefSCMatrix
OSSSCF::eigenvectors()
{
  return _eigenvectors;
}

void
OSSSCF::save_data_state(StateOut& s)
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
OSSSCF::occupation(int i)
{
  if (i < _ndocc) return 2.0;
  if (i < _ndocc + 2) return 1.0;
  return 0.0;
}

int
OSSSCF::value_implemented()
{
  return 1;
}

int
OSSSCF::gradient_implemented()
{
  return 1;
}

int
OSSSCF::hessian_implemented()
{
  return 0;
}

void
OSSSCF::print(ostream&o)
{
  OneBodyWavefunction::print(o);
}

void
OSSSCF::compute()
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
      // make sure we don't accidentally compute the gradient or hessian
      int gcomp = _gradient.compute(0);
      int hcomp = _hessian.compute(0);

      // start from core guess
      CLSCF hcwfn(*this);
      RefSCMatrix vec = hcwfn.eigenvectors();

      _eigenvectors = vec;
      _gradient.compute(gcomp);
      _hessian.compute(hcomp);
    }

    // schmidt orthogonalize the vector
    _eigenvectors.result_noupdate()->schmidt_orthog(overlap().pointer(),
                                                    _ndocc+2);

    if (_fock.null())
      _fock = _eigenvectors.result_noupdate()->rowdim()->create_symmmatrix();
    
    if (_opa_fock.null())
      _opa_fock = _fock.clone();
    
    if (_opb_fock.null())
      _opb_fock = _fock.clone();
    
    if (_fock_evals.null())
      _fock_evals = _fock->dim()->create_diagmatrix();
    
    printf("\n  OSSSCF::compute: energy accuracy = %g\n\n",
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

    printf("\n  OSSSCF::compute: gradient accuracy = %g\n\n",
           _gradient.desired_accuracy());

    do_gradient(gradient);
    gradient.print("cartesian gradient");
    set_gradient(gradient);

    _gradient.set_actual_accuracy(_gradient.desired_accuracy());
  }
  
  if (_hessian.needed()) {
    fprintf(stderr,"OSSSCF::compute: gradient not implemented\n");
    abort();
  }
  
}


void
OSSSCF::do_vector(double& eelec, double& nucrep)
{
  int i;

  _gr_vector = _eigenvectors.result_noupdate();
  
  // allocate storage for the temp arrays
  RefSCMatrix nvector = _gr_vector.clone();
  
  _gr_dens = _fock.clone();
  _gr_dens.assign(0.0);
  
  _gr_dens_diff = _gr_dens->clone();
  _gr_dens_diff.assign(0.0);
  
  _gr_gmat = _gr_dens->clone();
  _gr_gmat.assign(0.0);

  _gr_opa_dens = _fock.clone();
  _gr_opa_dens.assign(0.0);
  
  _gr_opa_dens_diff = _gr_dens->clone();
  _gr_opa_dens_diff.assign(0.0);
  
  _gr_opa_gmat = _gr_dens->clone();
  _gr_opa_gmat.assign(0.0);

  _gr_opb_dens = _fock.clone();
  _gr_opb_dens.assign(0.0);
  
  _gr_opb_dens_diff = _gr_dens->clone();
  _gr_opb_dens_diff.assign(0.0);
  
  _gr_opb_gmat = _gr_dens->clone();
  _gr_opb_gmat.assign(0.0);

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

  int_initialize_offsets2(centers,centers,centers,centers);

  nucrep = int_nuclear_repulsion(centers,centers);
  
  int flags = INT_EREP|INT_NOSTRB|INT_NOSTR1|INT_NOSTR2;
  double *intbuf = 
    int_initialize_erep(flags,0,centers,centers,centers,centers);

  int_storage(1000000);

  int_init_bounds();

  for (int iter=0; iter < _maxiter; iter++) {
    // form the density from the current vector 
    form_density(_gr_vector,
                 _gr_dens,_gr_dens_diff,
                 _gr_opa_dens,_gr_opa_dens_diff,
                 _gr_opb_dens,_gr_opb_dens_diff);

    // check convergence
    int ij=0;
    double delta=0;
    for (i=0; i < _gr_dens_diff->n(); i++)
      for (int j=0; j <= i; j++,ij++)
        delta += _gr_dens_diff.get_element(i,j)*_gr_dens_diff.get_element(i,j);
    delta = sqrt(delta/ij);

    if (delta < _energy.desired_accuracy()) break;

    // reset the density from time to time
    if (iter && !iter%10) {
      _gr_gmat.assign(0.0);
      _gr_dens_diff.assign(_gr_dens);
      _gr_opa_gmat.assign(0.0);
      _gr_opa_dens_diff.assign(_gr_opa_dens);
      _gr_opb_gmat.assign(0.0);
      _gr_opb_dens_diff.assign(_gr_opb_dens);
    }
      
    // scale the off-diagonal elements of the density matrix
    _gr_dens_diff->scale(2.0);
    _gr_dens_diff->scale_diagonal(0.5);
    _gr_opa_dens_diff->scale(2.0);
    _gr_opa_dens_diff->scale_diagonal(0.5);
    _gr_opb_dens_diff->scale(2.0);
    _gr_opb_dens_diff->scale_diagonal(0.5);
    
    // form the AO basis fock matrix
    form_ao_fock(centers,intbuf);

    // unscale the off-diagonal elements of the density matrix
    _gr_dens_diff->scale(0.5);
    _gr_dens_diff->scale_diagonal(2.0);
    _gr_opa_dens_diff->scale(0.5);
    _gr_opa_dens_diff->scale_diagonal(2.0);
    _gr_opb_dens_diff->scale(0.5);
    _gr_opb_dens_diff->scale_diagonal(2.0);

    // calculate the electronic energy
    eelec = scf_energy();
    printf("iter %5d energy = %20.15f delta = %15.10g\n",
           iter,eelec+nucrep,delta);

#if 0
    // now extrapolate the fock matrix
    // first we form the error matrix which is the offdiagonal blocks of
    // the MO fock matrix
    RefSymmSCMatrix mofock = _fock.clone();
    mofock.assign(0.0);
    mofock.accumulate_transform(_gr_vector.t(),_fock);

    RefSymmSCMatrix moofock = _op_fock.clone();
    moofock.assign(0.0);
    moofock.accumulate_transform(_gr_vector.t(),_op_fock);

    mofock.element_op(_accumeffh,moofock);

    for (i=0; i < mofock->n(); i++) {
      double occi = occupation(i);
      for (int j=0; j <= i; j++) {
        double occj = occupation(j);
        if (occi == occj)
          mofock.set_element(i,j,0.0);
      }
    }
    
    // now transform MO error to the AO basis
    RefSymmSCMatrix ao_error = mofock.clone();
    ao_error.assign(0.0);
    ao_error.accumulate_transform(_gr_vector,mofock);
    
    // and do the DIIS extrapolation
    _data = new SymmSCMatrix2SCExtrapData(_fock,_op_fock);
    _error = new SymmSCMatrixSCExtrapError(ao_error);
    _extrap->extrapolate(_data,_error);
    _data=0;
    _error=0;
    ao_error=0;

    // now transform extrapolated fock to MO basis
    mofock.assign(0.0);
    mofock.accumulate_transform(_gr_vector.t(),_fock);

    moofock.assign(0.0);
    moofock.accumulate_transform(_gr_vector.t(),_op_fock);

    mofock.element_op(_accumeffh,moofock);
    moofock=0;

    // diagonalize MO fock to get MO vector
    mofock.diagonalize(_fock_evals,nvector);
    mofock=0;
#else
    RefSymmSCMatrix mofock = _fock.clone();
    mofock.assign(0.0);
    mofock.accumulate_transform(_gr_vector.t(),_fock);

    RefSymmSCMatrix mooafock = _opa_fock.clone();
    mooafock.assign(0.0);
    mooafock.accumulate_transform(_gr_vector.t(),_opa_fock);

    RefSymmSCMatrix moobfock = _opb_fock.clone();
    moobfock.assign(0.0);
    moobfock.accumulate_transform(_gr_vector.t(),_opb_fock);

    int bvec=_ndocc+2-1;
    for (i=0; i < _ndocc; i++)
      mooafock.set_element(bvec,i,moobfock.get_element(bvec,i));
    
    for (i=_ndocc+2; i < basis()->nbasis(); i++) 
      mooafock.set_element(i,bvec,moobfock.get_element(i,bvec));

    mofock.element_op(_accumeffh,mooafock);

    mofock.diagonalize(_fock_evals,nvector);
#endif

    // transform MO vector to AO basis
    _gr_vector = _gr_vector * nvector;
    
    // and orthogonalize vector
    _gr_vector->schmidt_orthog(overlap().pointer(),basis()->nbasis());
  }
      
  _gr_vector.print("converged vector");
  _fock_evals.print("evals");
  
  _eigenvectors = _gr_vector;
  _eigenvectors.computed() = 1;
  
  int_done_bounds();
  int_done_erep();
  int_done_offsets2(centers,centers,centers,centers);
  int_done_storage();
  free_centers(centers);
  free(centers);

  _gr_dens = 0;
  _gr_dens_diff = 0;
  _gr_gmat = 0;
  _gr_opa_dens = 0;
  _gr_opa_dens_diff = 0;
  _gr_opa_gmat = 0;
  _gr_opb_dens = 0;
  _gr_opb_dens_diff = 0;
  _gr_opb_gmat = 0;
  _gr_hcore = 0;
  _gr_vector = 0;
  nvector = 0;
}

double
OSSSCF::scf_energy()
{
  RefSymmSCMatrix t = _fock.copy();
  t.accumulate(_gr_hcore);

  double eelec=0;
  for (int i=0; i < t->n(); i++) {
    for (int j=0; j < i; j++) {
      eelec += _gr_dens.get_element(i,j)*t.get_element(i,j)
               - _gr_opa_dens.get_element(i,j)*_gr_opa_gmat.get_element(i,j)
               - _gr_opb_dens.get_element(i,j)*_gr_opb_gmat.get_element(i,j);
    }
    eelec += 0.5*(_gr_dens.get_element(i,i)*t.get_element(i,i)
               - _gr_opa_dens.get_element(i,i)*_gr_opa_gmat.get_element(i,i)
               - _gr_opb_dens.get_element(i,i)*_gr_opb_gmat.get_element(i,i)
      );
  }
  return eelec;
}
