
#ifdef __GNUC__
#pragma implementation
#endif

#include <util/misc/newstring.h>
#include <math/optimize/diis.h>
#include <math/optimize/scextrapmat.h>
#include <chemistry/qc/integral/integralv2.h>
#include <chemistry/qc/scf/grscf.h>

#include <chemistry/qc/intv2/int_libv2.h>

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

      // schmidt orthogonalize the guess vector
      vec->schmidt_orthog(overlap().pointer(),_ndocc+_nsocc);
      //vec->print("guess scf vector");
      
      _eigenvectors = vec;
    }

    double eelec,nucrep;
    do_vector(eelec,nucrep);
      
    // this will be done elsewhere eventually
    set_energy(eelec+nucrep);
    _energy.set_actual_accuracy(_energy.desired_accuracy());
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


void
GRSCF::do_vector(double& eelec, double& nucrep)
{
  _gr_vector = _eigenvectors.result_noupdate();
  
  // allocate storage for the temp arrays
  _gr_nvector = _gr_vector.clone();
  
  _gr_evals = _gr_vector->rowdim()->create_diagmatrix();
  _gr_dens = _gr_vector->rowdim()->create_symmmatrix();
  _gr_dens.assign(0.0);
  
  _gr_dens_diff = _gr_dens->clone();
  _gr_dens_diff.assign(0.0);
  
  _gr_gmat = _gr_dens->clone();
  _gr_hcore = _gr_dens->clone();
  _gr_fock = _gr_dens->clone();

  if (_nsocc) {
    _gr_op_dens = _gr_dens->clone();
    _gr_op_dens.assign(0.0);
    _gr_op_dens_diff = _gr_dens->clone();
    _gr_op_dens_diff.assign(0.0);
    _gr_op_gmat = _gr_dens->clone();
    _gr_op_fock = _gr_dens->clone();
  }
  
  // form Hcore
  _gr_hcore.assign(0.0);
  _accumdih->accum(_gr_hcore);

  for (int iter=0; iter < _maxiter; iter++) {
    // form the density from the current vector 
    form_density(_gr_vector,
                 _gr_dens,_gr_dens_diff,
                 _gr_op_dens,_gr_op_dens_diff);
    
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

    _gr_dens->scale(0.5);
    _gr_dens->scale_diagonal(2.0);

    eelec = scf_energy();
    printf("iter %5d energy = %20.15f delta = %15.10g\n",
           iter,eelec+nucrep,delta);

    RefSymmSCMatrix _gr_error = _gr_fock.clone();
    _gr_error.assign(0.0);
    _gr_error.accumulate_transform(_gr_vector.t(),_gr_fock);

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
    
    _data = new SymmSCMatrixSCExtrapData(_gr_fock);
    _error = new SymmSCMatrixSCExtrapError(_gr_error);
    _extrap->extrapolate(_data,_error);
    _data=0;
    _error=0;

    RefSymmSCMatrix _gr_mofock = _gr_fock.clone();
    _gr_mofock.assign(0.0);
    _gr_mofock.accumulate_transform(_gr_vector.t(),_gr_fock);
    
    _gr_mofock.diagonalize(_gr_evals,_gr_nvector);
    _gr_mofock=0;

    _gr_vector = _gr_vector * _gr_nvector;
    
    _gr_vector->schmidt_orthog(overlap().pointer(),basis()->nbasis());
  }
      
  _gr_dens = 0;
  _gr_dens_diff = 0;
  _gr_op_dens = 0;
  _gr_op_dens_diff = 0;
  _gr_gmat = 0;
  _gr_op_gmat = 0;
  _gr_hcore = 0;
  _gr_fock = 0;
  _gr_op_fock = 0;
  _gr_evals = 0;
  _gr_vector = 0;
  _gr_nvector = 0;
}

double
GRSCF::scf_energy()
{
  RefSymmSCMatrix t = _gr_fock.copy();
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

#define ioff(i) (((i)*((i)+1))>>1)
#define IOFF(a,b) (((a)>(b))?(ioff(a)+(b)):(ioff(b)+(a)))

void
GRSCF::form_ao_fock(double&nucrep)
{
  _gr_gmat.assign(0.0);

  centers_t *centers = basis()->convert_to_centers_t(molecule());
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

  char *shnfunc = new char[centers->nshell];
  for (int i=0; i < centers->nshell; i++)
    shnfunc[i] = INT_SH_NFUNC((centers),i);
  
  double tnint=0;
  
  for (int i=0; i < centers->nshell; i++) {
    for (int j=0; j <= i; j++) {
      for (int k=0; k <= i; k++) {
        for (int l=0; l <= ((k==i)?j:k); l++) {
          int s1=i, s2=j, s3=k, s4=l;

          int_erep(INT_EREP|INT_NOBCHK|INT_NOPERM,&s1,&s2,&s3,&s4);

          int n1 = shnfunc[s1];
          int n2 = shnfunc[s2];
          int n3 = shnfunc[s3];
          int n4 = shnfunc[s4];

          int e12 = (s2==s1);
          int e13e24 = (s3==s1) && (s4==s2);
          int e34 = (s4 == s3);
          int e_any = (e12||e13e24||e34);
          
          int index=0;
          
          if (e_any) {
            for (int bf1=0; bf1<=INT_MAX1(n1) ; bf1++) {
              int i1 = centers->func_num[s1] + bf1;

              for (int bf2=0; bf2<=INT_MAX2(e12,bf1,n2) ; bf2++) {
                int j1 = centers->func_num[s2] + bf2;
                int ij1=ioff(i1)+j1;

                for (int bf3=0; bf3<=INT_MAX3(e13e24,bf1,n3) ; bf3++) {
                  int k1 = centers->func_num[s3] + bf3;

                  for (int bf4=0;
                       bf4<=INT_MAX4(e13e24,e34,bf1,bf2,bf3,n4);bf4++){
                    if (INT_NONZERO(intbuf[index])) {
                      int l1 = centers->func_num[s4] + bf4;

                      int ii,jj,kk,ll;
                      if (ij1 >= ioff(k1)+l1) {
                        ii = i1; jj = j1; kk = k1; ll = l1;
                      } else {
                        ii = k1; jj = l1; kk = i1; ll = j1;
                      }

                      double pki_int = intbuf[index];
                      double pkval;

                      int lij,lkl;
                      
                      if (jj == kk) {
                        /*
                         * if i=j=k or j=k=l, then this integral contributes
                         * to J, K1, and K2 of G(ij), so
                         * pkval = (ijkl) - 0.25 * ((ikjl)-(ilkj))
                         *       = 0.5 * (ijkl)
                         */
                        if (ii == jj || kk == ll) {
                          lij=ioff(ii)+jj;
                          lkl=ioff(kk)+ll;
                          pkval = (lij==lkl) ? 0.25*pki_int: 0.5*pki_int;
                          _gr_gmat.accumulate_element(ii,jj,pkval*
                                                  _gr_dens.get_element(kk,ll));
                          _gr_gmat.accumulate_element(kk,ll,pkval*
                                                  _gr_dens.get_element(ii,jj));
                        } else {
                          /*
                           * if j=k, then this integral contributes
                           * to J and K1 of G(ij)
                           *
                           * pkval = (ijkl) - 0.25 * (ikjl)
                           *       = 0.75 * (ijkl)
                           */
                          lij=ioff(ii)+jj;
                          lkl=ioff(kk)+ll;
                          pkval = (lij==lkl) ? 0.375*pki_int: 0.75*pki_int;

                          _gr_gmat.accumulate_element(ii,jj,pkval*
                                                  _gr_dens.get_element(kk,ll));
                          _gr_gmat.accumulate_element(kk,ll,pkval*
                                                  _gr_dens.get_element(ii,jj));

                          /*
                           * this integral also contributes to K1 and K2 of
                           * G(il)
                           *
                           * pkval = -0.25 * ((ijkl)+(ikjl))
                           *       = -0.5 * (ijkl)
                           */
                          lij=ioff(ii)+ll;
                          lkl=IOFF(kk,jj);
                          pkval = (lij==lkl)? 0.25*pki_int: 0.5*pki_int;

                          _gr_gmat.accumulate_element(ii,ll,-pkval*
                                                  _gr_dens.get_element(kk,jj));
                          _gr_gmat.accumulate_element(kk,jj,-pkval*
                                                  _gr_dens.get_element(ii,ll));
                        }
                      } else if (ii == kk || jj == ll) {
                        /*
                         * if i=k or j=l, then this integral contributes
                         * to J and K2 of G(ij)
                         *
                         * pkval = (ijkl) - 0.25 * (ilkj)
                         *       = 0.75 * (ijkl)
                         */
                        lij=ioff(ii)+jj;
                        lkl=ioff(kk)+ll;
                        pkval = (lij==lkl) ? 0.375*pki_int: 0.75*pki_int;
                        _gr_gmat.accumulate_element(ii,jj,pkval*
                                                  _gr_dens.get_element(kk,ll));
                        _gr_gmat.accumulate_element(kk,ll,pkval*
                                                  _gr_dens.get_element(ii,jj));

                        /*
                         * this integral also contributes to K1 and K2 of
                         * G(ik)
                         *
                         * pkval = -0.25 * ((ijkl)+(ilkj))
                         *       = -0.5 * (ijkl)
                         */
                        lij=ioff(ii)+kk;
                        lkl=IOFF(jj,ll);
                        pkval = (lij==lkl) ? 0.25*pki_int : 0.5*pki_int;

                        _gr_gmat.accumulate_element(ii,kk,-pkval*
                                                  _gr_dens.get_element(jj,ll));
                        _gr_gmat.accumulate_element(jj,ll,-pkval*
                                                  _gr_dens.get_element(ii,kk));
                      } else {
                        /*
                         * This integral contributes to J of G(ij)
                         *
                         * pkval = (ijkl)
                         */
                        lij=ioff(ii)+jj;
                        lkl=ioff(kk)+ll;
                        pkval = (lij==lkl)? 0.5*pki_int : pki_int;

                        _gr_gmat.accumulate_element(ii,jj,pkval*
                                                  _gr_dens.get_element(kk,ll));
                        _gr_gmat.accumulate_element(kk,ll,pkval*
                                                  _gr_dens.get_element(ii,jj));

                        /*
                         * and to K1 of G(ik)
                         *
                         * pkval = -0.25 * (ijkl)
                         */
                        lij=ioff(ii)+kk;
                        lkl=IOFF(jj,ll);
                        pkval = (lij==lkl) ? 0.125*pki_int : 0.25*pki_int;

                        _gr_gmat.accumulate_element(ii,kk,-pkval*
                                                  _gr_dens.get_element(jj,ll));
                        _gr_gmat.accumulate_element(jj,ll,-pkval*
                                                  _gr_dens.get_element(ii,kk));

                        if ((ii != jj) && (kk != ll)) {
                          /*
                           * if i!=j and k!=l, then this integral also
                           * contributes to K2 of G(il)
                           *
                           * pkval = -0.25 * (ijkl)
                           *
                           * note: if we get here, then ik can't equal jl,
                           * so pkval wasn't multiplied by 0.5 above.
                           */
                          lij=ioff(ii)+ll;
                          lkl=IOFF(kk,jj);

                          _gr_gmat.accumulate_element(ii,ll,-pkval*
                                                  _gr_dens.get_element(kk,jj));
                          _gr_gmat.accumulate_element(kk,jj,-pkval*
                                                  _gr_dens.get_element(ii,ll));
                        }
                      }
                    }
                    index++;
                  }
                }
              }
            }
          } else {
            for (int i1=centers->func_num[s1], bf1=0; bf1<n1 ; bf1++, i1++) {
              for (int j1=centers->func_num[s2], bf2=0; bf2<n2 ; bf2++, j1++) {
                int ij1=ioff(i1)+j1;

                for (int k1=centers->func_num[s3],bf3=0; bf3<n3; bf3++,k1++) {
                  for (int l1=centers->func_num[s4],bf4=0;bf4<n4;bf4++,l1++) {
                    if (INT_NONZERO(intbuf[index])) {

                      int ii,jj,kk,ll;
                      if (ij1 >= ioff(k1)+l1) {
                        ii = i1; jj = j1; kk = k1; ll = l1;
                      } else {
                        ii = k1; jj = l1; kk = i1; ll = j1;
                      }

                      double pki_int = intbuf[index];
                      double pkval;
                      int lij,lkl;

                      if (jj == kk) {
                        /*
                         * if j=k, then this integral contributes
                         * to J and K1 of G(ij)
                         *
                         * pkval = (ijkl) - 0.25 * (ikjl)
                         *       = 0.75 * (ijkl)
                         */
                        lij=ioff(ii)+jj;
                        lkl=ioff(kk)+ll;
                        pkval = 0.75*pki_int;

                        _gr_gmat.accumulate_element(ii,jj,pkval*
                                                  _gr_dens.get_element(kk,ll));
                        _gr_gmat.accumulate_element(kk,ll,pkval*
                                                  _gr_dens.get_element(ii,jj));

                        /*
                         * this integral also contributes to K1 and K2 of
                         * G(il)
                         *
                         * pkval = -0.25 * ((ijkl)+(ikjl))
                         *       = -0.5 * (ijkl)
                         */
                        lij=ioff(ii)+ll;
                        lkl=IOFF(kk,jj);
                        pkval *= 0.666666666666666;
                        _gr_gmat.accumulate_element(ii,ll,-pkval*
                                                  _gr_dens.get_element(kk,jj));
                        _gr_gmat.accumulate_element(kk,jj,-pkval*
                                                  _gr_dens.get_element(ii,ll));
                      } else if (ii == kk || jj == ll) {
                        /*
                         * if i=k or j=l, then this integral contributes
                         * to J and K2 of G(ij)
                         *
                         * pkval = (ijkl) - 0.25 * (ilkj)
                         *       = 0.75 * (ijkl)
                         */
                        lij=ioff(ii)+jj;
                        lkl=ioff(kk)+ll;
                        pkval = 0.75*pki_int;

                        _gr_gmat.accumulate_element(ii,jj,pkval*
                                                  _gr_dens.get_element(kk,ll));
                        _gr_gmat.accumulate_element(kk,ll,pkval*
                                                  _gr_dens.get_element(ii,jj));

                        /*
                         * this integral also contributes to K1 and K2 of
                         * G(ik)
                         *
                         * pkval = -0.25 * ((ijkl)+(ilkj))
                         *       = -0.5 * (ijkl)
                         */
                        lij=ioff(ii)+kk;
                        lkl=IOFF(jj,ll);
                        pkval *= 0.666666666666666;

                        _gr_gmat.accumulate_element(ii,kk,-pkval*
                                                  _gr_dens.get_element(jj,ll));
                        _gr_gmat.accumulate_element(jj,ll,-pkval*
                                                  _gr_dens.get_element(ii,kk));
                      } else {
                        /*
                         * This integral contributes to J of G(ij)
                         *
                         * pkval = (ijkl)
                         */
                        lij=ioff(ii)+jj;
                        lkl=ioff(kk)+ll;
                        pkval = pki_int;

                        _gr_gmat.accumulate_element(ii,jj,pkval*
                                                  _gr_dens.get_element(kk,ll));
                        _gr_gmat.accumulate_element(kk,ll,pkval*
                                                  _gr_dens.get_element(ii,jj));

                        /*
                         * and to K1 of G(ik)
                         *
                         * pkval = -0.25 * (ijkl)
                         */
                        lij=ioff(ii)+kk;
                        lkl=IOFF(jj,ll);
                        pkval *= 0.25;

                        _gr_gmat.accumulate_element(ii,kk,-pkval*
                                                  _gr_dens.get_element(jj,ll));
                        _gr_gmat.accumulate_element(jj,ll,-pkval*
                                                  _gr_dens.get_element(ii,kk));

                        /*
                         * and to K2 of G(il)
                         *
                         * pkval = -0.25 * (ijkl)
                         */
                        lij=ioff(ii)+ll;
                        lkl=IOFF(kk,jj);

                        _gr_gmat.accumulate_element(ii,ll,-pkval*
                                                  _gr_dens.get_element(kk,jj));
                        _gr_gmat.accumulate_element(kk,jj,-pkval*
                                                  _gr_dens.get_element(ii,ll));
                      }
                    }
                    index++;
                  }
                }
              }
            }
          }
          tnint += (double) (n1*n2*n3*n4);
        }
      }
    }
  }
  
  //_gr_gmat.print("g matrix");
  
  int_done_erep();
  int_done_offsets2(centers,centers,centers,centers);
  int_done_storage();
  free_centers(centers);

  _gr_fock.assign(_gr_gmat);
  _gr_fock.accumulate(_gr_hcore);

  //_gr_fock.print("Fock matrix");
  
}
