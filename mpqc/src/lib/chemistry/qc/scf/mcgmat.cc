
#include <math/array/math_lib.h>
#include <math/scmat/local.h>
#include <chemistry/qc/integral/integralv2.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/scf/mcscf.h>

#define ioff(i) (((i)*((i)+1))>>1)
#define IOFF(a,b) (((a)>(b))?(ioff(a)+(b)):(ioff(b)+(a)))

static void
dens2(const RefSCMatrix& vec,
      const RefSymmSCMatrix& densc,
      const RefSymmSCMatrix& densa,
      const RefSymmSCMatrix& densb,
      const RefSymmSCMatrix& densab1,
      const RefSymmSCMatrix& densab2,
      const RefSCVector& _ca,
      const RefSCVector& _cb,
      int nbasis, int ndocc, int aorb, int borb)
{
  // find out what type of matrices we're dealing with
  if (LocalSCMatrix::castdown(vec.pointer())) {
    LocalSCMatrix *lvec = LocalSCMatrix::require_castdown(
      vec.pointer(), "MCSCF::form_density");
    LocalSymmSCMatrix *ldensc = LocalSymmSCMatrix::require_castdown(
      densc.pointer(), "MCSCF::form_density");
    LocalSymmSCMatrix *ldensa = LocalSymmSCMatrix::require_castdown(
      densa.pointer(), "MCSCF::form_density");
    LocalSymmSCMatrix *ldensb = LocalSymmSCMatrix::require_castdown(
      densb.pointer(), "MCSCF::form_density");
    LocalSymmSCMatrix *ldensab1 = LocalSymmSCMatrix::require_castdown(
      densab1.pointer(), "MCSCF::form_density");
    LocalSymmSCMatrix *ldensab2 = LocalSymmSCMatrix::require_castdown(
      densab2.pointer(), "MCSCF::form_density");

    for (int i=0; i < nbasis; i++) {
      for (int j=0; j <= i; j++) {
        double pt=0;
        for (int k=0; k < ndocc; k++)
          pt += lvec->get_element(i,k)*lvec->get_element(j,k);

#if 1
        double ptoa=_ca->get_element(i)*_ca->get_element(j);
        double ptob=_cb->get_element(i)*_cb->get_element(j);
#else
        double ptoa=lvec->get_element(i,aorb)*lvec->get_element(j,aorb);
        double ptob=lvec->get_element(i,borb)*lvec->get_element(j,borb);
#endif
        double ptoab=_ca->get_element(i)*_cb->get_element(j)
                    +_ca->get_element(j)*_cb->get_element(i);

        ldensc->set_element(i,j,pt);
        ldensa->set_element(i,j,ptoa);
        ldensb->set_element(i,j,ptob);
        ldensab1->set_element(i,j,ptoa+ptob);
        ldensab2->set_element(i,j,ptoab);
      }
    }
    ldensc->scale(2.0);
    ldensa->scale(2.0);
    ldensb->scale(2.0);
    ldensab1->scale(2.0);
    ldensab2->scale(2.0);
  }
}

void
MCSCF::form_ao_fock(centers_t *centers, double *intbuf, double& eelec)
{
  int inttol = int_bound_log(_energy.desired_accuracy()/100.0);

  char *shnfunc = new char[centers->nshell];
  for (int i=0; i < centers->nshell; i++)
    shnfunc[i] = INT_SH_NFUNC((centers),i);

  dens2(_gr_vector,_densc,_densa,_densb,_densab2,_densab,_ca,_cb,
        basis()->nbasis(),_ndocc,aorb,borb);
  
  _densc->scale(2.0);
  _densc->scale_diagonal(0.5);
  _densa->scale(2.0);
  _densa->scale_diagonal(0.5);
  _densb->scale(2.0);
  _densb->scale_diagonal(0.5);
  _densab->scale(2.0);
  _densab->scale_diagonal(0.5);
  
  _fockc.assign(0.0);
  _focka.assign(0.0);
  _fockb.assign(0.0);
  _fockab.assign(0.0);
  _ka.assign(0.0);
  _kb.assign(0.0);
  
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

                          _fockc.accumulate_element(ii,jj,pkval*
                                             _densc.get_element(kk,ll));
                          _fockc.accumulate_element(kk,ll,pkval*
                                             _densc.get_element(ii,jj));

                          _focka.accumulate_element(ii,jj,pkval*
                                             _densa.get_element(kk,ll));
                          _focka.accumulate_element(kk,ll,pkval*
                                             _densa.get_element(ii,jj));

                          _fockb.accumulate_element(ii,jj,pkval*
                                             _densb.get_element(kk,ll));
                          _fockb.accumulate_element(kk,ll,pkval*
                                             _densb.get_element(ii,jj));

                          _fockab.accumulate_element(ii,jj,pkval*
                                             _densab.get_element(kk,ll));
                          _fockab.accumulate_element(kk,ll,pkval*
                                             _densab.get_element(ii,jj));

                          _ka.accumulate_element(ii,jj,pkval*
                                         _densa.get_element(kk,ll));
                          _ka.accumulate_element(kk,ll,pkval*
                                         _densa.get_element(ii,jj));

                          _kb.accumulate_element(ii,jj,pkval*
                                         _densb.get_element(kk,ll));
                          _kb.accumulate_element(kk,ll,pkval*
                                         _densb.get_element(ii,jj));

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

                          _fockc.accumulate_element(ii,jj,pkval*
                                             _densc.get_element(kk,ll));
                          _fockc.accumulate_element(kk,ll,pkval*
                                             _densc.get_element(ii,jj));

                          _focka.accumulate_element(ii,jj,pkval*
                                             _densa.get_element(kk,ll));
                          _focka.accumulate_element(kk,ll,pkval*
                                             _densa.get_element(ii,jj));

                          _fockb.accumulate_element(ii,jj,pkval*
                                             _densb.get_element(kk,ll));
                          _fockb.accumulate_element(kk,ll,pkval*
                                             _densb.get_element(ii,jj));

                          _fockab.accumulate_element(ii,jj,pkval*
                                             _densab.get_element(kk,ll));
                          _fockab.accumulate_element(kk,ll,pkval*
                                             _densab.get_element(ii,jj));

                          pkval *= 0.33333333333333333;

                          _ka.accumulate_element(ii,jj,pkval*
                                         _densa.get_element(kk,ll));
                          _ka.accumulate_element(kk,ll,pkval*
                                         _densa.get_element(ii,jj));

                          _kb.accumulate_element(ii,jj,pkval*
                                         _densb.get_element(kk,ll));
                          _kb.accumulate_element(kk,ll,pkval*
                                         _densb.get_element(ii,jj));

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

                          _fockc.accumulate_element(ii,ll,-pkval*
                                             _densc.get_element(kk,jj));
                          _fockc.accumulate_element(kk,jj,-pkval*
                                             _densc.get_element(ii,ll));

                          _focka.accumulate_element(ii,ll,-pkval*
                                             _densa.get_element(kk,jj));
                          _focka.accumulate_element(kk,jj,-pkval*
                                             _densa.get_element(ii,ll));

                          _fockb.accumulate_element(ii,ll,-pkval*
                                             _densb.get_element(kk,jj));
                          _fockb.accumulate_element(kk,jj,-pkval*
                                             _densb.get_element(ii,ll));

                          _fockab.accumulate_element(ii,ll,-pkval*
                                             _densab.get_element(kk,jj));
                          _fockab.accumulate_element(kk,jj,-pkval*
                                             _densab.get_element(ii,ll));

                          _ka.accumulate_element(ii,ll,pkval*
                                         _densa.get_element(kk,jj));
                          _ka.accumulate_element(kk,jj,pkval*
                                         _densa.get_element(ii,ll));

                          _kb.accumulate_element(ii,ll,pkval*
                                         _densb.get_element(kk,jj));
                          _kb.accumulate_element(kk,jj,pkval*
                                         _densb.get_element(ii,ll));

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
                        _fockc.accumulate_element(ii,jj,pkval*
                                             _densc.get_element(kk,ll));
                        _fockc.accumulate_element(kk,ll,pkval*
                                             _densc.get_element(ii,jj));

                        _focka.accumulate_element(ii,jj,pkval*
                                             _densa.get_element(kk,ll));
                        _focka.accumulate_element(kk,ll,pkval*
                                             _densa.get_element(ii,jj));

                        _fockb.accumulate_element(ii,jj,pkval*
                                             _densb.get_element(kk,ll));
                        _fockb.accumulate_element(kk,ll,pkval*
                                             _densb.get_element(ii,jj));

                        _fockab.accumulate_element(ii,jj,pkval*
                                             _densab.get_element(kk,ll));
                        _fockab.accumulate_element(kk,ll,pkval*
                                             _densab.get_element(ii,jj));

                        pkval *= 0.3333333333333333;
                        _ka.accumulate_element(ii,jj,pkval*
                                         _densa.get_element(kk,ll));
                        _ka.accumulate_element(kk,ll,pkval*
                                         _densa.get_element(ii,jj));

                        _kb.accumulate_element(ii,jj,pkval*
                                         _densb.get_element(kk,ll));
                        _kb.accumulate_element(kk,ll,pkval*
                                         _densb.get_element(ii,jj));

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

                        _fockc.accumulate_element(ii,kk,-pkval*
                                             _densc.get_element(jj,ll));
                        _fockc.accumulate_element(jj,ll,-pkval*
                                             _densc.get_element(ii,kk));

                        _focka.accumulate_element(ii,kk,-pkval*
                                             _densa.get_element(jj,ll));
                        _focka.accumulate_element(jj,ll,-pkval*
                                             _densa.get_element(ii,kk));

                        _fockb.accumulate_element(ii,kk,-pkval*
                                             _densb.get_element(jj,ll));
                        _fockb.accumulate_element(jj,ll,-pkval*
                                             _densb.get_element(ii,kk));

                        _fockab.accumulate_element(ii,kk,-pkval*
                                             _densab.get_element(jj,ll));
                        _fockab.accumulate_element(jj,ll,-pkval*
                                             _densab.get_element(ii,kk));

                        _ka.accumulate_element(ii,kk,pkval*
                                         _densa.get_element(jj,ll));
                        _ka.accumulate_element(jj,ll,pkval*
                                         _densa.get_element(ii,kk));

                        _kb.accumulate_element(ii,kk,pkval*
                                         _densb.get_element(jj,ll));
                        _kb.accumulate_element(jj,ll,pkval*
                                         _densb.get_element(ii,kk));

                      } else {
                        /*
                         * This integral contributes to J of G(ij)
                         *
                         * pkval = (ijkl)
                         */
                        lij=ioff(ii)+jj;
                        lkl=ioff(kk)+ll;
                        pkval = (lij==lkl)? 0.5*pki_int : pki_int;

                        _fockc.accumulate_element(ii,jj,pkval*
                                             _densc.get_element(kk,ll));
                        _fockc.accumulate_element(kk,ll,pkval*
                                             _densc.get_element(ii,jj));

                        _focka.accumulate_element(ii,jj,pkval*
                                             _densa.get_element(kk,ll));
                        _focka.accumulate_element(kk,ll,pkval*
                                             _densa.get_element(ii,jj));

                        _fockb.accumulate_element(ii,jj,pkval*
                                             _densb.get_element(kk,ll));
                        _fockb.accumulate_element(kk,ll,pkval*
                                             _densb.get_element(ii,jj));

                        _fockab.accumulate_element(ii,jj,pkval*
                                             _densab.get_element(kk,ll));
                        _fockab.accumulate_element(kk,ll,pkval*
                                             _densab.get_element(ii,jj));

                        /*
                         * and to K1 of G(ik)
                         *
                         * pkval = -0.25 * (ijkl)
                         */
                        lij=ioff(ii)+kk;
                        lkl=IOFF(jj,ll);
                        pkval = (lij==lkl) ? 0.125*pki_int : 0.25*pki_int;

                        _fockc.accumulate_element(ii,kk,-pkval*
                                             _densc.get_element(jj,ll));
                        _fockc.accumulate_element(jj,ll,-pkval*
                                             _densc.get_element(ii,kk));

                        _focka.accumulate_element(ii,kk,-pkval*
                                             _densa.get_element(jj,ll));
                        _focka.accumulate_element(jj,ll,-pkval*
                                             _densa.get_element(ii,kk));

                        _fockb.accumulate_element(ii,kk,-pkval*
                                             _densb.get_element(jj,ll));
                        _fockb.accumulate_element(jj,ll,-pkval*
                                             _densb.get_element(ii,kk));

                        _fockab.accumulate_element(ii,kk,-pkval*
                                             _densab.get_element(jj,ll));
                        _fockab.accumulate_element(jj,ll,-pkval*
                                             _densab.get_element(ii,kk));

                        _ka.accumulate_element(ii,kk,pkval*
                                         _densa.get_element(jj,ll));
                        _ka.accumulate_element(jj,ll,pkval*
                                         _densa.get_element(ii,kk));

                        _kb.accumulate_element(ii,kk,pkval*
                                         _densb.get_element(jj,ll));
                        _kb.accumulate_element(jj,ll,pkval*
                                         _densb.get_element(ii,kk));

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

                          _fockc.accumulate_element(ii,ll,-pkval*
                                             _densc.get_element(kk,jj));
                          _fockc.accumulate_element(kk,jj,-pkval*
                                             _densc.get_element(ii,ll));

                          _focka.accumulate_element(ii,ll,-pkval*
                                             _densa.get_element(kk,jj));
                          _focka.accumulate_element(kk,jj,-pkval*
                                             _densa.get_element(ii,ll));

                          _fockb.accumulate_element(ii,ll,-pkval*
                                             _densb.get_element(kk,jj));
                          _fockb.accumulate_element(kk,jj,-pkval*
                                             _densb.get_element(ii,ll));

                          _fockab.accumulate_element(ii,ll,-pkval*
                                             _densab.get_element(kk,jj));
                          _fockab.accumulate_element(kk,jj,-pkval*
                                             _densab.get_element(ii,ll));

                          _ka.accumulate_element(ii,ll,pkval*
                                         _densa.get_element(kk,jj));
                          _ka.accumulate_element(kk,jj,pkval*
                                         _densa.get_element(ii,ll));

                          _kb.accumulate_element(ii,ll,pkval*
                                         _densb.get_element(kk,jj));
                          _kb.accumulate_element(kk,jj,pkval*
                                         _densb.get_element(ii,ll));

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

                        _fockc.accumulate_element(ii,jj,pkval*
                                             _densc.get_element(kk,ll));
                        _fockc.accumulate_element(kk,ll,pkval*
                                             _densc.get_element(ii,jj));

                        _focka.accumulate_element(ii,jj,pkval*
                                             _densa.get_element(kk,ll));
                        _focka.accumulate_element(kk,ll,pkval*
                                             _densa.get_element(ii,jj));

                        _fockb.accumulate_element(ii,jj,pkval*
                                             _densb.get_element(kk,ll));
                        _fockb.accumulate_element(kk,ll,pkval*
                                             _densb.get_element(ii,jj));

                        _fockab.accumulate_element(ii,jj,pkval*
                                             _densab.get_element(kk,ll));
                        _fockab.accumulate_element(kk,ll,pkval*
                                             _densab.get_element(ii,jj));

                        pkval *= 0.3333333333333333;
                        _ka.accumulate_element(ii,jj,pkval*
                                         _densa.get_element(kk,ll));
                        _ka.accumulate_element(kk,ll,pkval*
                                         _densa.get_element(ii,jj));

                        _kb.accumulate_element(ii,jj,pkval*
                                         _densb.get_element(kk,ll));
                        _kb.accumulate_element(kk,ll,pkval*
                                         _densb.get_element(ii,jj));

                        /*
                         * this integral also contributes to K1 and K2 of
                         * G(il)
                         *
                         * pkval = -0.25 * ((ijkl)+(ikjl))
                         *       = -0.5 * (ijkl)
                         */
                        lij=ioff(ii)+ll;
                        lkl=IOFF(kk,jj);
                        pkval = 0.5 * pki_int;

                        _fockc.accumulate_element(ii,ll,-pkval*
                                             _densc.get_element(kk,jj));
                        _fockc.accumulate_element(kk,jj,-pkval*
                                             _densc.get_element(ii,ll));

                        _focka.accumulate_element(ii,ll,-pkval*
                                             _densa.get_element(kk,jj));
                        _focka.accumulate_element(kk,jj,-pkval*
                                             _densa.get_element(ii,ll));

                        _fockb.accumulate_element(ii,ll,-pkval*
                                             _densb.get_element(kk,jj));
                        _fockb.accumulate_element(kk,jj,-pkval*
                                             _densb.get_element(ii,ll));

                        _fockab.accumulate_element(ii,ll,-pkval*
                                             _densab.get_element(kk,jj));
                        _fockab.accumulate_element(kk,jj,-pkval*
                                             _densab.get_element(ii,ll));

                        _ka.accumulate_element(ii,ll,pkval*
                                         _densa.get_element(kk,jj));
                        _ka.accumulate_element(kk,jj,pkval*
                                         _densa.get_element(ii,ll));

                        _kb.accumulate_element(ii,ll,pkval*
                                         _densb.get_element(kk,jj));
                        _kb.accumulate_element(kk,jj,pkval*
                                         _densb.get_element(ii,ll));

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

                        _fockc.accumulate_element(ii,jj,pkval*
                                             _densc.get_element(kk,ll));
                        _fockc.accumulate_element(kk,ll,pkval*
                                             _densc.get_element(ii,jj));

                        _focka.accumulate_element(ii,jj,pkval*
                                             _densa.get_element(kk,ll));
                        _focka.accumulate_element(kk,ll,pkval*
                                             _densa.get_element(ii,jj));

                        _fockb.accumulate_element(ii,jj,pkval*
                                             _densb.get_element(kk,ll));
                        _fockb.accumulate_element(kk,ll,pkval*
                                             _densb.get_element(ii,jj));

                        _fockab.accumulate_element(ii,jj,pkval*
                                             _densab.get_element(kk,ll));
                        _fockab.accumulate_element(kk,ll,pkval*
                                             _densab.get_element(ii,jj));

                        pkval *= 0.3333333333333333;
                        _ka.accumulate_element(ii,jj,pkval*
                                         _densa.get_element(kk,ll));
                        _ka.accumulate_element(kk,ll,pkval*
                                         _densa.get_element(ii,jj));

                        _kb.accumulate_element(ii,jj,pkval*
                                         _densb.get_element(kk,ll));
                        _kb.accumulate_element(kk,ll,pkval*
                                         _densb.get_element(ii,jj));

                        /*
                         * this integral also contributes to K1 and K2 of
                         * G(ik)
                         *
                         * pkval = -0.25 * ((ijkl)+(ilkj))
                         *       = -0.5 * (ijkl)
                         */
                        lij=ioff(ii)+kk;
                        lkl=IOFF(jj,ll);
                        pkval = 0.5 * pki_int;

                        _fockc.accumulate_element(ii,kk,-pkval*
                                             _densc.get_element(jj,ll));
                        _fockc.accumulate_element(jj,ll,-pkval*
                                             _densc.get_element(ii,kk));

                        _focka.accumulate_element(ii,kk,-pkval*
                                             _densa.get_element(jj,ll));
                        _focka.accumulate_element(jj,ll,-pkval*
                                             _densa.get_element(ii,kk));

                        _fockb.accumulate_element(ii,kk,-pkval*
                                             _densb.get_element(jj,ll));
                        _fockb.accumulate_element(jj,ll,-pkval*
                                             _densb.get_element(ii,kk));

                        _fockab.accumulate_element(ii,kk,-pkval*
                                             _densab.get_element(jj,ll));
                        _fockab.accumulate_element(jj,ll,-pkval*
                                             _densab.get_element(ii,kk));

                        _ka.accumulate_element(ii,kk,pkval*
                                         _densa.get_element(jj,ll));
                        _ka.accumulate_element(jj,ll,pkval*
                                         _densa.get_element(ii,kk));

                        _kb.accumulate_element(ii,kk,pkval*
                                         _densb.get_element(jj,ll));
                        _kb.accumulate_element(jj,ll,pkval*
                                         _densb.get_element(ii,kk));

                      } else {
                        /*
                         * This integral contributes to J of G(ij)
                         *
                         * pkval = (ijkl)
                         */
                        lij=ioff(ii)+jj;
                        lkl=ioff(kk)+ll;
                        pkval = pki_int;

                        _fockc.accumulate_element(ii,jj,pkval*
                                             _densc.get_element(kk,ll));
                        _fockc.accumulate_element(kk,ll,pkval*
                                             _densc.get_element(ii,jj));

                        _focka.accumulate_element(ii,jj,pkval*
                                             _densa.get_element(kk,ll));
                        _focka.accumulate_element(kk,ll,pkval*
                                             _densa.get_element(ii,jj));

                        _fockb.accumulate_element(ii,jj,pkval*
                                             _densb.get_element(kk,ll));
                        _fockb.accumulate_element(kk,ll,pkval*
                                             _densb.get_element(ii,jj));

                        _fockab.accumulate_element(ii,jj,pkval*
                                             _densab.get_element(kk,ll));
                        _fockab.accumulate_element(kk,ll,pkval*
                                             _densab.get_element(ii,jj));

                        /*
                         * and to K1 of G(ik)
                         *
                         * pkval = -0.25 * (ijkl)
                         */
                        lij=ioff(ii)+kk;
                        lkl=IOFF(jj,ll);
                        pkval *= 0.25;

                        _fockc.accumulate_element(ii,kk,-pkval*
                                             _densc.get_element(jj,ll));
                        _fockc.accumulate_element(jj,ll,-pkval*
                                             _densc.get_element(ii,kk));

                        _focka.accumulate_element(ii,kk,-pkval*
                                             _densa.get_element(jj,ll));
                        _focka.accumulate_element(jj,ll,-pkval*
                                             _densa.get_element(ii,kk));

                        _fockb.accumulate_element(ii,kk,-pkval*
                                             _densb.get_element(jj,ll));
                        _fockb.accumulate_element(jj,ll,-pkval*
                                             _densb.get_element(ii,kk));

                        _fockab.accumulate_element(ii,kk,-pkval*
                                             _densab.get_element(jj,ll));
                        _fockab.accumulate_element(jj,ll,-pkval*
                                             _densab.get_element(ii,kk));

                        _ka.accumulate_element(ii,kk,pkval*
                                         _densa.get_element(jj,ll));
                        _ka.accumulate_element(jj,ll,pkval*
                                         _densa.get_element(ii,kk));

                        _kb.accumulate_element(ii,kk,pkval*
                                         _densb.get_element(jj,ll));
                        _kb.accumulate_element(jj,ll,pkval*
                                         _densb.get_element(ii,kk));

                        /*
                         * and to K2 of G(il)
                         *
                         * pkval = -0.25 * (ijkl)
                         */
                        lij=ioff(ii)+ll;
                        lkl=IOFF(kk,jj);

                        _fockc.accumulate_element(ii,ll,-pkval*
                                             _densc.get_element(kk,jj));
                        _fockc.accumulate_element(kk,jj,-pkval*
                                             _densc.get_element(ii,ll));

                        _focka.accumulate_element(ii,ll,-pkval*
                                             _densa.get_element(kk,jj));
                        _focka.accumulate_element(kk,jj,-pkval*
                                             _densa.get_element(ii,ll));

                        _fockb.accumulate_element(ii,ll,-pkval*
                                             _densb.get_element(kk,jj));
                        _fockb.accumulate_element(kk,jj,-pkval*
                                             _densb.get_element(ii,ll));

                        _fockab.accumulate_element(ii,ll,-pkval*
                                             _densab.get_element(kk,jj));
                        _fockab.accumulate_element(kk,jj,-pkval*
                                             _densab.get_element(ii,ll));

                        _ka.accumulate_element(ii,ll,pkval*
                                         _densa.get_element(kk,jj));
                        _ka.accumulate_element(kk,jj,pkval*
                                         _densa.get_element(ii,ll));

                        _kb.accumulate_element(ii,ll,pkval*
                                         _densb.get_element(kk,jj));
                        _kb.accumulate_element(kk,jj,pkval*
                                         _densb.get_element(ii,ll));
                      }
                    }
                    index++;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  int_reduce_storage_threshold();

  delete[] shnfunc;

  _densc->scale(0.5);
  _densc->scale_diagonal(2.0);
  
  _densa->scale(0.5);
  _densa->scale_diagonal(2.0);
  
  _densb->scale(0.5);
  _densb->scale_diagonal(2.0);
  
  _densab->scale(0.5);
  _densab->scale_diagonal(2.0);
  
  double h11=0;
  double h21=0;
  double h22=0;
  double h31=0;
  double h32=0;
  double h33=0;
  
  RefSymmSCMatrix fa = _fockc.copy();
  fa.accumulate(_gr_hcore);
  fa.accumulate(_focka);
  
  RefSymmSCMatrix fb = _fockc.copy();
  fb.accumulate(_gr_hcore);
  fb.accumulate(_fockb);
  
  for (int i=0; i < basis()->nbasis(); i++) {
    for (int j=0; j < i; j++) {
      h11 += (_densc.get_element(i,j)+_densa.get_element(i,j))*
             (2*_gr_hcore.get_element(i,j)+_fockc.get_element(i,j)+
              _focka.get_element(i,j));

      h33 += (_densc.get_element(i,j)+_densb.get_element(i,j))*
             (2*_gr_hcore.get_element(i,j)+_fockc.get_element(i,j)+
              _fockb.get_element(i,j));

      h22 += (2*_densc.get_element(i,j)+_densab2.get_element(i,j))*
                _gr_hcore.get_element(i,j) +
              _densc.get_element(i,j)*(_fockc.get_element(i,j)+
                _focka.get_element(i,j)+_fockb.get_element(i,j)) +
             0.5*_densa.get_element(i,j)*(
                _fockb.get_element(i,j)+3*_kb.get_element(i,j));

      h21 += _densab.get_element(i,j)*fa.get_element(i,j)*2/sqrt(2.0);
      h32 += _densab.get_element(i,j)*fb.get_element(i,j)*2/sqrt(2.0);

      h31 += _densa.get_element(i,j)*_kb.get_element(i,j);
    }
    h11 += 0.5*(_densc.get_element(i,i)+_densa.get_element(i,i))*
           (2*_gr_hcore.get_element(i,i)+_fockc.get_element(i,i)+
            _focka.get_element(i,i));
    h33 += 0.5*(_densc.get_element(i,i)+_densb.get_element(i,i))*
           (2*_gr_hcore.get_element(i,i)+_fockc.get_element(i,i)+
            _fockb.get_element(i,i));

    int j=i;
    h22 += 0.5*((2*_densc.get_element(i,j)+_densab2.get_element(i,j))*
                _gr_hcore.get_element(i,j) +
              _densc.get_element(i,j)*(_fockc.get_element(i,j)+
                _focka.get_element(i,j)+_fockb.get_element(i,j)) +
             0.5*_densa.get_element(i,j)*(
                _fockb.get_element(i,j)+3*_kb.get_element(i,j)));

    h21 += _densab.get_element(i,j)*fa.get_element(i,j)/sqrt(2.0);
    h32 += _densab.get_element(i,j)*fb.get_element(i,j)/sqrt(2.0);

    h31 += 0.5*_densa.get_element(i,i)*_kb.get_element(i,i);
  }

  RefSCDimension l3 = new LocalSCDimension(3);
  RefSymmSCMatrix h = l3->create_symmmatrix();
  RefSCMatrix hv = l3->create_matrix(l3);
  RefDiagSCMatrix hl = l3->create_diagmatrix();

  h.set_element(0,0,h11);
  h.set_element(1,0,h21);
  h.set_element(1,1,h22);
  h.set_element(2,0,h31);
  h.set_element(2,1,h32);
  h.set_element(2,2,h33);

  h.print("ci matrix");
  h.diagonalize(hl,hv);
  hv.print("ci coeffs");
  hl.print("ci evals");
  
  double mci1=hv.get_element(0,0);
  double mci2=hv.get_element(1,0);
  double mci3=hv.get_element(2,0);
  
  double e = mci1*mci1*h11 + mci2*mci2*h22 + mci3*mci3*h33 +
             2*(mci1*mci2*h21 + mci1*mci3*h31 + mci2*mci3*h32);
  printf("e1 = %lf\n",e);
  
  mci1=hv.get_element(0,1);
  mci2=hv.get_element(1,1);
  mci3=hv.get_element(2,1);
  
  e = mci1*mci1*h11 + mci2*mci2*h22 + mci3*mci3*h33 +
             2*(mci1*mci2*h21 + mci1*mci3*h31 + mci2*mci3*h32);
  printf("e2 = %lf\n",e);
  
  double hab = 0;
  double jkab = 0;
  for (int i=0; i < basis()->nbasis(); i++) {
    for (int j=0; j < i; j++) {
      hab += _ca.get_element(i)*_cb.get_element(j)*
             _gr_hcore.get_element(i,j);
      hab += _ca.get_element(j)*_cb.get_element(i)*
             _gr_hcore.get_element(i,j);
      jkab += 0.5*_densa.get_element(i,j)*(_fockb.get_element(i,j)+
                                           3*_kb.get_element(i,j));
    }
    hab += _ca.get_element(i)*_cb.get_element(i)*_gr_hcore.get_element(i,i);
    jkab += _densa.get_element(i,i)*0.25*(_fockb.get_element(i,i)+
                                         3*_kb.get_element(i,i));
  }

  double alpha = 1.0/(1.0+sab*sab);
  
  double e1 = 2.0*alpha*sab*hab;
  double e2 = alpha*jkab;
  eop = 2.0*sab*hab+jkab;
  
  for (int i=0; i < basis()->nbasis(); i++) {
    for (int j=0; j < i; j++) {
      e1 += (2.0*_densc.get_element(i,j) + alpha*_densab2.get_element(i,j))*
            _gr_hcore.get_element(i,j);
      e2 += _densc.get_element(i,j)*_fockc.get_element(i,j) +
            alpha*(_densab2.get_element(i,j)*_fockc.get_element(i,j)+
                   2.0*sab*_densc.get_element(i,j)*_fockab.get_element(i,j));
      eop += _densab2.get_element(i,j)*_gr_hcore.get_element(i,j) +
             _densab2.get_element(i,j)*_fockc.get_element(i,j) +
             2.0*sab*_densc.get_element(i,j)*_fockab.get_element(i,j);
    }
    e1 += (2.0*_densc.get_element(i,i) + alpha*_densab2.get_element(i,i))*
          0.5*_gr_hcore.get_element(i,i);
    e2 += 0.5*_densc.get_element(i,i)*_fockc.get_element(i,i) +
          0.5*alpha*(_densab2.get_element(i,i)*_fockc.get_element(i,i)+
                 2.0*sab*_densc.get_element(i,i)*_fockab.get_element(i,i));
    eop += 0.5*_densab2.get_element(i,i)*_gr_hcore.get_element(i,i) +
           0.5*_densab2.get_element(i,i)*_fockc.get_element(i,i) +
           sab*_densc.get_element(i,i)*_fockab.get_element(i,i);
  }
  
  printf("sab = %lf\n",sab);
  eelec = e1 + e2;
  
}
