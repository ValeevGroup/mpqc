
#ifdef __GNUC__
#pragma implementation
#endif

#include <util/misc/newstring.h>
#include <math/optimize/diis.h>
#include <chemistry/qc/integral/integralv2.h>
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
  _accumdih = new AccumHCore;
  
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
  
  _accumddh = keyval->describedclassvalue("accumddh");
  if (_accumddh.null()) {
    _accumddh = new GSGeneralEffH;
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
  return 0;
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
  if (_energy.needed()) {
    if (_eigenvectors.result_noupdate().null()) {
      // start from core guess
      HCoreWfn hcwfn(*this);
      RefSCMatrix vec = hcwfn.eigenvectors();

      vec->schmidt_orthog(overlap().pointer(),_ndocc+_nocc);
      
      _eigenvectors = vec;

      set_energy(0.0);
      _energy.set_actual_accuracy(_energy.desired_accuracy());
    }
  }

  if (_gradient.needed()) {
    fprintf(stderr,"GRSCF::compute: gradient not implemented\n");
    abort();
  }
  
  if (_hessian.needed()) {
    fprintf(stderr,"GRSCF::compute: gradient not implemented\n");
    abort();
  }
  
}

