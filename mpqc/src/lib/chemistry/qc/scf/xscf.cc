//
// xscf.cc --- implementation of the excited-state open shell singlet SCF class
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifdef __GNUC__
#pragma implementation
#endif

#include <util/misc/newstring.h>
#include <math/optimize/diis.h>
#include <math/optimize/scextrapmat.h>
#include <chemistry/qc/scf/clscf.h>
#include <chemistry/qc/scf/ossscf.h>
#include <chemistry/qc/scf/hsosscf.h>
#include <chemistry/qc/scf/xscf.h>

///////////////////////////////////////////////////////////////////////////
// XSCF

#define CLASSNAME XSCF
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#define PARENTS public OneBodyWavefunction
#include <util/class/classi.h>
void *
XSCF::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = OneBodyWavefunction::_castdown(cd);
  return do_castdowns(casts,cd);
}

static void
occ(PointBag_double *z, int &nd)
{
  int Z=molecule()->nuclear_charge();

  nd = Z/2;
  if (Z%2) {
    fprintf(stderr,"XSCF::occ: Warning, there's a leftover electron.\n");
    fprintf(stderr,"  total nuclear charge = %d, %d closed shells\n",Z,nd);
    fprintf(stderr,"  total charge = %d\n\n",Z-2*nd);
  }

  nd--;
}

void
XSCF::init()
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

XSCF::XSCF(StateIn& s) :
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

XSCF::XSCF(const RefKeyVal& keyval) :
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

XSCF::XSCF(const OneBodyWavefunction& obwfn) :
  OneBodyWavefunction(obwfn)
{
  init();
}

XSCF::XSCF(const XSCF& xscf) :
  OneBodyWavefunction(xscf)
{
  _extrap = xscf._extrap;
  _data = xscf._data;
  _error = xscf._error;
  _accumdih = xscf._accumdih;
  _accumddh = xscf._accumddh;
  _accumeffh = xscf._accumeffh;
  _ndocc = xscf._ndocc;
  _density_reset_freq = xscf._density_reset_freq;
  _maxiter = xscf._maxiter;
  _eliminate = xscf._eliminate;

  ckptdir = new_string(xscf.ckptdir);
  fname = new_string(xscf.fname);
}

XSCF::~XSCF()
{
}

RefSCMatrix
XSCF::eigenvectors()
{
  return _eigenvectors;
}

void
XSCF::save_data_state(StateOut& s)
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
XSCF::occupation(int i)
{
  if (i < _ndocc) return 2.0;
  if (i == aorb) return occa;
  if (i == borb) return occb;
  return 0.0;
}

int
XSCF::value_implemented()
{
  return 1;
}

int
XSCF::gradient_implemented()
{
  return 1;
}

void
XSCF::print(ostream&o)
{
  OneBodyWavefunction::print(o);
}

void
XSCF::compute()
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

      // schmidt orthogonalize the vector
      vec->schmidt_orthog(overlap().pointer(),_ndocc+2);

      _eigenvectors = vec;
      _gradient.compute(gcomp);
      _hessian.compute(hcomp);
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
    
    printf("\n  XSCF::compute: energy accuracy = %g\n\n",
           _energy.desired_accuracy());

    double eelec,nucrep;
    do_vector(eelec,nucrep);
    double eother = 0.0;
    if (accumddh_.nonnull()) eother = accumddh_->e();
      
    // this will be done elsewhere eventually
    printf("  total scf energy = %15.10f\n",eelec+eother+nucrep);

    set_energy(eelec+eother+nucrep);
    _energy.set_actual_accuracy(_energy.desired_accuracy());
  }

  if (_gradient.needed()) {
    RefSCVector gradient = _moldim->create_vector();

    printf("\n  XSCF::compute: gradient accuracy = %g\n\n",
           _gradient.desired_accuracy());

    do_gradient(gradient);
    gradient.print("cartesian gradient");
    set_gradient(gradient);

    _gradient.set_actual_accuracy(_gradient.desired_accuracy());
  }
  
  if (_hessian.needed()) {
    fprintf(stderr,"XSCF::compute: gradient not implemented\n");
    abort();
  }
  
}


void
XSCF::do_vector(double& eelec, double& nucrep)
{
  int i;

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

  int_initialize_offsets2(centers,centers,centers,centers);

  nucrep = int_nuclear_repulsion(centers,centers);
  
  int flags = INT_EREP|INT_NOSTRB|INT_NOSTR1|INT_NOSTR2;
  double *intbuf = 
    int_initialize_erep(flags,0,centers,centers,centers,centers);

#ifdef SGI
    int_storage(12500000);
#else
    int_storage(1000000);
#endif

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

    double eother = 0.0;
    if (accumddh_.nonnull()) eother = accumddh_->e();

    printf("iter %5d energy = %15.10f delta = %15.10g\n",
           iter,eelec+eother+nucrep,olde-eelec-eother);

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
    for (i=0; i < nbasis; i++) {
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
    for (i=0; i < nbasis; i++) {
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
    for (i=0; i < nbasis; i++) {
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
    for (i=0; i < fooc->nrow(); i++)
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
    for (i=0; i < nbasis; i++) {
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

    for (i=0; i < nbasis; i++)
      for (int j=0; j < vb->ncol(); j++)
        _gr_vector.set_element(i,j+_ndocc,vb.get_element(i,j));
                               
    // and orthogonalize vector
    _gr_vector->schmidt_orthog(ovlp.pointer(),basis()->nbasis());
  
  }
      
  sab=0;
  for (i=0; i < nbasis; i++) {
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
  
  for (i=0; i < nbasis; i++) {
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
  
  _gr_vector.print("scfx vector");

  _gr_vector.assign_column(_ca,aorb);
  _gr_vector.assign_column(_cb,borb);
  _gr_vector->schmidt_orthog(ovlp.pointer(),basis()->nbasis());
  
  _gr_vector.print("ortho vector");
  _fock_evalsc.print("evalsc");
  _fock_evalsa.print("evalsa");
  _fock_evalsb.print("evalsb");
  
  _eigenvectors = _gr_vector;
  _eigenvectors.computed() = 1;
  
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

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
