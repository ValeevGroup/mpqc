
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
    } else {
      // we must already have an old vector (and sab I hope)
      _ca.scale(sqrt(2.0*(1.0+sab)));
      _cb.scale(sqrt(2.0*(1.0-sab)));
      for (int i=0; i < _ca.n(); i++) {
        double a = _ca.get_element(i)+_cb.get_element(i);
        double b = _ca.get_element(i)-_cb.get_element(i);
        _ca.set_element(i,0.5*a);
        _cb.set_element(i,0.5*b);
      }
      _gr_vector = _eigenvectors.result_noupdate();
      _gr_vector.assign_column(_ca,aorb);
      _gr_vector.assign_column(_cb,borb);
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

    if (_fock_evalsc.null())
      _fock_evalsc = _focka->dim()->create_diagmatrix();
    
    if (_fock_evalsa.null())
      _fock_evalsa = actived->create_diagmatrix();
    
    if (_fock_evalsb.null())
      _fock_evalsb = actived->create_diagmatrix();
    
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
  RefSCMatrix nvectorc = _gr_vector.clone();
  RefSCMatrix nvectora =
    _fock_evalsa->dim()->create_matrix(_fock_evalsa->dim().pointer());
  RefSCMatrix nvectorb = nvectora.clone();
  
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

  _ca = _gr_vector.get_column(aorb);
  _cb = _gr_vector.get_column(borb);
  
  int nbasis = basis()->nbasis();

  for (int iter=0; ; iter++) {
    // form the AO basis fock matrix
    double olde=eelec;
    form_ao_fock(centers,intbuf,eelec);

    if (fabs(olde-eelec) < 1.0e-13)
      break;

    printf("iter %5d energy = %20.15f delta = %15.10g\n",
           iter,eelec+nucrep,olde-eelec);

    RefSymmSCMatrix sfc = _gr_hcore.copy();
    sfc.accumulate(_fockc);
    
    double alpha = 1.0/(1.0+sab*sab);
    RefSymmSCMatrix sfo = _fockab.copy();
    sfo.scale(sab);
    sfo.accumulate(_focka);
    sfo.accumulate(_fockb);
    sfo.scale(alpha*0.5);
    sfo.accumulate(sfc);

    RefSCMatrix densb = ovlp.dim()->create_matrix(ovlp.dim());
    RefSCMatrix densa = densb.clone();
    for (int i=0; i < nbasis; i++) {
      for (int j=0; j < i; j++) {
        densb.set_element(i,j,_densb.get_element(i,j));
        densb.set_element(j,i,_densb.get_element(i,j));
        densa.set_element(i,j,_densa.get_element(i,j));
        densa.set_element(j,i,_densa.get_element(i,j));
      }
      densb.set_element(i,i,_densb.get_element(i,i));
      densa.set_element(i,i,_densa.get_element(i,i));
    }
    densb->scale(0.5);
    densa->scale(0.5);
    
    RefSCMatrix fas = sfc * densa * ovlp;
    fas.accumulate(fas.t());

    RefSCMatrix fbs = sfc * densb * ovlp;
    fbs.accumulate(fbs.t());

    RefSymmSCMatrix jaka = _ka.copy();
    jaka.scale(3.0);
    jaka.accumulate(_focka);
    jaka.scale(0.5);

    RefSymmSCMatrix jbkb = _kb.copy();
    jbkb.scale(3.0);
    jbkb.accumulate(_fockb);
    jbkb.scale(0.5);

    RefSCMatrix sas = ovlp * densa * ovlp;
    sas.scale(-alpha*eop);
        
    RefSCMatrix sbs = ovlp * densb * ovlp;
    sbs.scale(-alpha*eop);
        
    RefSymmSCMatrix sfq = sfc.copy();
    RefSymmSCMatrix sfr = sfc.copy();
    sfq.accumulate(jbkb);
    sfr.accumulate(jaka);
    for (int i=0; i < nbasis; i++) {
      for (int j=0; j <= i; j++) {
        sfq.accumulate_element(i,j,sbs.get_element(i,j)+fbs.get_element(i,j));
        sfr.accumulate_element(i,j,sas.get_element(i,j)+fas.get_element(i,j));
      }
    }

    fas = sfq * densa * ovlp;
    fas.accumulate(fas.t());

    fbs = sfr * densb * ovlp;
    fbs.accumulate(fbs.t());

    fas.accumulate(fbs);
    fas.scale(-alpha/2.0);
    for (int i=0; i < nbasis; i++) {
      for (int j=0; j <= i; j++) {
        sfo.accumulate_element(i,j,fas.get_element(i,j));
      }
    }

    RefSymmSCMatrix g0 = sfo.clone();
    g0.assign(0.0);
    g0.accumulate_transform(_gr_vector.t(),sfo);

    g0.diagonalize(_fock_evalsc,nvectorc);
    RefSCMatrix vc = _gr_vector*nvectorc;

    // grab the active and virtual orbitals from vc
    RefSCMatrix fooc =
      _gr_vector->rowdim()->create_matrix(nvectora->coldim().pointer());
    for (int i=0; i < fooc->nrow(); i++)
      for (int j=0; j < fooc->ncol(); j++)
        if (_ndocc)
          fooc.set_element(i,j,vc.get_element(i,j+_ndocc));
        else
          fooc.set_element(i,j,_gr_vector.get_element(i,j));
          
    
    RefSymmSCMatrix ga = _fock_evalsa->dim()->create_symmmatrix();
    ga.assign(0.0);
    ga.accumulate_transform(fooc.t(),sfq);

    ga.diagonalize(_fock_evalsa,nvectora);
    RefSCMatrix va = fooc*nvectora;
    
    _ca = va.get_column(0);
    
    RefSymmSCMatrix gb = ga.clone();
    gb.assign(0.0);
    gb.accumulate_transform(va.t(),sfr);

    gb.diagonalize(_fock_evalsb,nvectorb);
    RefSCMatrix vb = va*nvectorb;
    
    _cb = vb.get_column(1);

    sab = 0;
    for (int i=0; i < nbasis; i++) {
      for (int j=0; j < nbasis; j++) {
        sab += _ca.get_element(i)*_cb.get_element(j)*ovlp.get_element(i,j);
      }
    }

    if (sab<0) {
      sab = -sab;
      _cb.scale(-1.0);
      vb.assign_column(_cb,1);
    }

    _gr_vector.assign(vc);

    for (int i=0; i < nbasis; i++)
      for (int j=0; j < vb->ncol(); j++)
        _gr_vector.set_element(i,j+_ndocc,vb.get_element(i,j));
                               
    // and orthogonalize vector
    _gr_vector->schmidt_orthog(ovlp.pointer(),basis()->nbasis());
  
  }
      
  sab=0;
  for (int i=0; i < nbasis; i++) {
    for (int j=0; j < nbasis; j++) {
      sab += _ca.get_element(i)*_cb.get_element(j)*ovlp.get_element(i,j);
    }
  }

  if (sab<0.0) {
    sab = -sab;
    _cb.scale(-1.0);
  }
  
  double ud = 1.0/sqrt(2.0*(1.0+sab));
  double vd = 1.0/sqrt(2.0*(1.0-sab));
  
  for (int i=0; i < nbasis; i++) {
    double u = ud*(_ca.get_element(i)+_cb.get_element(i));
    double v = vd*(_ca.get_element(i)-_cb.get_element(i));
    _ca.set_element(i,u);
    _cb.set_element(i,v);
  }
      
  double fooe;
  form_ao_fock(centers,intbuf,fooe);

  ci1 = (1.0+sab)/sqrt(2.0*(1.0+sab*sab));
  ci2 = -(1.0-sab)/sqrt(2.0*(1.0+sab*sab));

  occa = 2.0*ci1*ci1;
  occb = 2.0*ci2*ci2;

  _focka.accumulate(_fockc);
  _focka.accumulate(_gr_hcore);
  _fockb.accumulate(_fockc);
  _fockb.accumulate(_gr_hcore);
  
  printf("sab = %lf, ci1 = %lf, ci2 = %lf\n",sab,ci1,ci2);
  printf("occa = %lf, occb = %lf\n",occa,occb);
  
  _gr_vector.assign_column(_ca,aorb);
  _gr_vector.assign_column(_cb,borb);
  
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
  nvectorc = 0;
  nvectora = 0;
  nvectorb = 0;
}
