
#include <math/array/math_lib.h>
#include <math/scmat/local.h>
#include <chemistry/qc/integral/integralv2.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/scf/tcscf.h>

#define ioff(i) (((i)*((i)+1))>>1)
#define IOFF(a,b) (((a)>(b))?(ioff(a)+(b)):(ioff(b)+(a)))

static void
dens2(const RefSCMatrix& vec,
      const RefSymmSCMatrix& dens,
      const RefSymmSCMatrix& dens1,
      const RefSymmSCMatrix& dens2,
      int nbasis, int ndocc)
{
  int k;
  
  // find out what type of matrices we're dealing with
  if (LocalSCMatrix::castdown(vec.pointer())) {
    LocalSCMatrix *lvec = LocalSCMatrix::require_castdown(
      vec.pointer(), "TCSCF::form_density");
    LocalSymmSCMatrix *ldens = LocalSymmSCMatrix::require_castdown(
      dens.pointer(), "TCSCF::form_density");
    LocalSymmSCMatrix *ldens1 = LocalSymmSCMatrix::require_castdown(
      dens1.pointer(), "TCSCF::form_density");
    LocalSymmSCMatrix *ldens2 = LocalSymmSCMatrix::require_castdown(
      dens2.pointer(), "TCSCF::form_density");

    for (int i=0; i < nbasis; i++) {
      for (int j=0; j <= i; j++) {
        double pt=0;
        for (k=0; k < ndocc; k++)
          pt += 2.0*lvec->get_element(i,k)*lvec->get_element(j,k);

        double ptoa=0;
        for (k=ndocc; k < ndocc+1; k++)
          ptoa += 2.0*lvec->get_element(i,k)*lvec->get_element(j,k);

        double ptob=0;
        for (k=ndocc+1; k < ndocc+2; k++)
          ptob += 2.0*lvec->get_element(i,k)*lvec->get_element(j,k);

        ldens->set_element(i,j,pt);
        ldens1->set_element(i,j,ptoa);
        ldens2->set_element(i,j,ptob);
      }
    }
  }
}

void
TCSCF::form_ao_fock(centers_t *centers, double *intbuf, double& eelec)
{
  int i;
  
  int inttol = int_bound_log(_energy.desired_accuracy()/100.0);

  char *shnfunc = new char[centers->nshell];
  for (i=0; i < centers->nshell; i++)
    shnfunc[i] = INT_SH_NFUNC((centers),i);

  dens2(_gr_vector,_gr_dens,_gr_opa_dens,_gr_opb_dens,
        basis()->nbasis(),_ndocc);
  
  _gr_dens->scale(2.0);
  _gr_dens->scale_diagonal(0.5);
  _gr_opa_dens->scale(2.0);
  _gr_opa_dens->scale_diagonal(0.5);
  _gr_opb_dens->scale(2.0);
  _gr_opb_dens->scale_diagonal(0.5);
  
  RefSymmSCMatrix dca = _gr_dens.copy();
  dca.accumulate(_gr_opa_dens);
  
  RefSymmSCMatrix dcb = _gr_dens.copy();
  dcb.accumulate(_gr_opb_dens);
  
  RefSymmSCMatrix gca = _focka;
  gca.assign(0.0);
  RefSymmSCMatrix gcb = _fockb;
  gcb.assign(0.0);
  RefSymmSCMatrix ka = _ka;
  ka.assign(0.0);
  RefSymmSCMatrix kb = _kb;
  kb.assign(0.0);
  
  for (i=0; i < centers->nshell; i++) {
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

                          gca.accumulate_element(ii,jj,pkval*
                                             dca.get_element(kk,ll));
                          gca.accumulate_element(kk,ll,pkval*
                                             dca.get_element(ii,jj));

                          gcb.accumulate_element(ii,jj,pkval*
                                             dcb.get_element(kk,ll));
                          gcb.accumulate_element(kk,ll,pkval*
                                             dcb.get_element(ii,jj));

                          ka.accumulate_element(ii,jj,pkval*
                                         _gr_opa_dens.get_element(kk,ll));
                          ka.accumulate_element(kk,ll,pkval*
                                         _gr_opa_dens.get_element(ii,jj));

                          kb.accumulate_element(ii,jj,pkval*
                                         _gr_opb_dens.get_element(kk,ll));
                          kb.accumulate_element(kk,ll,pkval*
                                         _gr_opb_dens.get_element(ii,jj));

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

                          gca.accumulate_element(ii,jj,pkval*
                                             dca.get_element(kk,ll));
                          gca.accumulate_element(kk,ll,pkval*
                                             dca.get_element(ii,jj));

                          gcb.accumulate_element(ii,jj,pkval*
                                             dcb.get_element(kk,ll));
                          gcb.accumulate_element(kk,ll,pkval*
                                             dcb.get_element(ii,jj));

                          pkval *= 0.33333333333333333;

                          ka.accumulate_element(ii,jj,pkval*
                                         _gr_opa_dens.get_element(kk,ll));
                          ka.accumulate_element(kk,ll,pkval*
                                         _gr_opa_dens.get_element(ii,jj));

                          kb.accumulate_element(ii,jj,pkval*
                                         _gr_opb_dens.get_element(kk,ll));
                          kb.accumulate_element(kk,ll,pkval*
                                         _gr_opb_dens.get_element(ii,jj));

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

                          gca.accumulate_element(ii,ll,-pkval*
                                             dca.get_element(kk,jj));
                          gca.accumulate_element(kk,jj,-pkval*
                                             dca.get_element(ii,ll));

                          gcb.accumulate_element(ii,ll,-pkval*
                                             dcb.get_element(kk,jj));
                          gcb.accumulate_element(kk,jj,-pkval*
                                             dcb.get_element(ii,ll));

                          ka.accumulate_element(ii,ll,pkval*
                                         _gr_opa_dens.get_element(kk,jj));
                          ka.accumulate_element(kk,jj,pkval*
                                         _gr_opa_dens.get_element(ii,ll));

                          kb.accumulate_element(ii,ll,pkval*
                                         _gr_opb_dens.get_element(kk,jj));
                          kb.accumulate_element(kk,jj,pkval*
                                         _gr_opb_dens.get_element(ii,ll));

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
                        gca.accumulate_element(ii,jj,pkval*
                                             dca.get_element(kk,ll));
                        gca.accumulate_element(kk,ll,pkval*
                                             dca.get_element(ii,jj));

                        gcb.accumulate_element(ii,jj,pkval*
                                             dcb.get_element(kk,ll));
                        gcb.accumulate_element(kk,ll,pkval*
                                             dcb.get_element(ii,jj));

                        pkval *= 0.3333333333333333;
                        ka.accumulate_element(ii,jj,pkval*
                                         _gr_opa_dens.get_element(kk,ll));
                        ka.accumulate_element(kk,ll,pkval*
                                         _gr_opa_dens.get_element(ii,jj));

                        kb.accumulate_element(ii,jj,pkval*
                                         _gr_opb_dens.get_element(kk,ll));
                        kb.accumulate_element(kk,ll,pkval*
                                         _gr_opb_dens.get_element(ii,jj));

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

                        gca.accumulate_element(ii,kk,-pkval*
                                             dca.get_element(jj,ll));
                        gca.accumulate_element(jj,ll,-pkval*
                                             dca.get_element(ii,kk));

                        gcb.accumulate_element(ii,kk,-pkval*
                                             dcb.get_element(jj,ll));
                        gcb.accumulate_element(jj,ll,-pkval*
                                             dcb.get_element(ii,kk));

                        ka.accumulate_element(ii,kk,pkval*
                                         _gr_opa_dens.get_element(jj,ll));
                        ka.accumulate_element(jj,ll,pkval*
                                         _gr_opa_dens.get_element(ii,kk));

                        kb.accumulate_element(ii,kk,pkval*
                                         _gr_opb_dens.get_element(jj,ll));
                        kb.accumulate_element(jj,ll,pkval*
                                         _gr_opb_dens.get_element(ii,kk));

                      } else {
                        /*
                         * This integral contributes to J of G(ij)
                         *
                         * pkval = (ijkl)
                         */
                        lij=ioff(ii)+jj;
                        lkl=ioff(kk)+ll;
                        pkval = (lij==lkl)? 0.5*pki_int : pki_int;

                        gca.accumulate_element(ii,jj,pkval*
                                             dca.get_element(kk,ll));
                        gca.accumulate_element(kk,ll,pkval*
                                             dca.get_element(ii,jj));

                        gcb.accumulate_element(ii,jj,pkval*
                                             dcb.get_element(kk,ll));
                        gcb.accumulate_element(kk,ll,pkval*
                                             dcb.get_element(ii,jj));

                        /*
                         * and to K1 of G(ik)
                         *
                         * pkval = -0.25 * (ijkl)
                         */
                        lij=ioff(ii)+kk;
                        lkl=IOFF(jj,ll);
                        pkval = (lij==lkl) ? 0.125*pki_int : 0.25*pki_int;

                        gca.accumulate_element(ii,kk,-pkval*
                                             dca.get_element(jj,ll));
                        gca.accumulate_element(jj,ll,-pkval*
                                             dca.get_element(ii,kk));

                        gcb.accumulate_element(ii,kk,-pkval*
                                             dcb.get_element(jj,ll));
                        gcb.accumulate_element(jj,ll,-pkval*
                                             dcb.get_element(ii,kk));

                        ka.accumulate_element(ii,kk,pkval*
                                         _gr_opa_dens.get_element(jj,ll));
                        ka.accumulate_element(jj,ll,pkval*
                                         _gr_opa_dens.get_element(ii,kk));

                        kb.accumulate_element(ii,kk,pkval*
                                         _gr_opb_dens.get_element(jj,ll));
                        kb.accumulate_element(jj,ll,pkval*
                                         _gr_opb_dens.get_element(ii,kk));

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

                          gca.accumulate_element(ii,ll,-pkval*
                                             dca.get_element(kk,jj));
                          gca.accumulate_element(kk,jj,-pkval*
                                             dca.get_element(ii,ll));

                          gcb.accumulate_element(ii,ll,-pkval*
                                             dcb.get_element(kk,jj));
                          gcb.accumulate_element(kk,jj,-pkval*
                                             dcb.get_element(ii,ll));

                          ka.accumulate_element(ii,ll,pkval*
                                         _gr_opa_dens.get_element(kk,jj));
                          ka.accumulate_element(kk,jj,pkval*
                                         _gr_opa_dens.get_element(ii,ll));

                          kb.accumulate_element(ii,ll,pkval*
                                         _gr_opb_dens.get_element(kk,jj));
                          kb.accumulate_element(kk,jj,pkval*
                                         _gr_opb_dens.get_element(ii,ll));

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

                        gca.accumulate_element(ii,jj,pkval*
                                             dca.get_element(kk,ll));
                        gca.accumulate_element(kk,ll,pkval*
                                             dca.get_element(ii,jj));

                        gcb.accumulate_element(ii,jj,pkval*
                                             dcb.get_element(kk,ll));
                        gcb.accumulate_element(kk,ll,pkval*
                                             dcb.get_element(ii,jj));

                        pkval *= 0.3333333333333333;
                        ka.accumulate_element(ii,jj,pkval*
                                         _gr_opa_dens.get_element(kk,ll));
                        ka.accumulate_element(kk,ll,pkval*
                                         _gr_opa_dens.get_element(ii,jj));

                        kb.accumulate_element(ii,jj,pkval*
                                         _gr_opb_dens.get_element(kk,ll));
                        kb.accumulate_element(kk,ll,pkval*
                                         _gr_opb_dens.get_element(ii,jj));

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

                        gca.accumulate_element(ii,ll,-pkval*
                                             dca.get_element(kk,jj));
                        gca.accumulate_element(kk,jj,-pkval*
                                             dca.get_element(ii,ll));

                        gcb.accumulate_element(ii,ll,-pkval*
                                             dcb.get_element(kk,jj));
                        gcb.accumulate_element(kk,jj,-pkval*
                                             dcb.get_element(ii,ll));

                        ka.accumulate_element(ii,ll,pkval*
                                         _gr_opa_dens.get_element(kk,jj));
                        ka.accumulate_element(kk,jj,pkval*
                                         _gr_opa_dens.get_element(ii,ll));

                        kb.accumulate_element(ii,ll,pkval*
                                         _gr_opb_dens.get_element(kk,jj));
                        kb.accumulate_element(kk,jj,pkval*
                                         _gr_opb_dens.get_element(ii,ll));

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

                        gca.accumulate_element(ii,jj,pkval*
                                             dca.get_element(kk,ll));
                        gca.accumulate_element(kk,ll,pkval*
                                             dca.get_element(ii,jj));

                        gcb.accumulate_element(ii,jj,pkval*
                                             dcb.get_element(kk,ll));
                        gcb.accumulate_element(kk,ll,pkval*
                                             dcb.get_element(ii,jj));

                        pkval *= 0.3333333333333333;
                        ka.accumulate_element(ii,jj,pkval*
                                         _gr_opa_dens.get_element(kk,ll));
                        ka.accumulate_element(kk,ll,pkval*
                                         _gr_opa_dens.get_element(ii,jj));

                        kb.accumulate_element(ii,jj,pkval*
                                         _gr_opb_dens.get_element(kk,ll));
                        kb.accumulate_element(kk,ll,pkval*
                                         _gr_opb_dens.get_element(ii,jj));

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

                        gca.accumulate_element(ii,kk,-pkval*
                                             dca.get_element(jj,ll));
                        gca.accumulate_element(jj,ll,-pkval*
                                             dca.get_element(ii,kk));

                        gcb.accumulate_element(ii,kk,-pkval*
                                             dcb.get_element(jj,ll));
                        gcb.accumulate_element(jj,ll,-pkval*
                                             dcb.get_element(ii,kk));

                        ka.accumulate_element(ii,kk,pkval*
                                         _gr_opa_dens.get_element(jj,ll));
                        ka.accumulate_element(jj,ll,pkval*
                                         _gr_opa_dens.get_element(ii,kk));

                        kb.accumulate_element(ii,kk,pkval*
                                         _gr_opb_dens.get_element(jj,ll));
                        kb.accumulate_element(jj,ll,pkval*
                                         _gr_opb_dens.get_element(ii,kk));

                      } else {
                        /*
                         * This integral contributes to J of G(ij)
                         *
                         * pkval = (ijkl)
                         */
                        lij=ioff(ii)+jj;
                        lkl=ioff(kk)+ll;
                        pkval = pki_int;

                        gca.accumulate_element(ii,jj,pkval*
                                             dca.get_element(kk,ll));
                        gca.accumulate_element(kk,ll,pkval*
                                             dca.get_element(ii,jj));

                        gcb.accumulate_element(ii,jj,pkval*
                                             dcb.get_element(kk,ll));
                        gcb.accumulate_element(kk,ll,pkval*
                                             dcb.get_element(ii,jj));

                        /*
                         * and to K1 of G(ik)
                         *
                         * pkval = -0.25 * (ijkl)
                         */
                        lij=ioff(ii)+kk;
                        lkl=IOFF(jj,ll);
                        pkval *= 0.25;

                        gca.accumulate_element(ii,kk,-pkval*
                                             dca.get_element(jj,ll));
                        gca.accumulate_element(jj,ll,-pkval*
                                             dca.get_element(ii,kk));

                        gcb.accumulate_element(ii,kk,-pkval*
                                             dcb.get_element(jj,ll));
                        gcb.accumulate_element(jj,ll,-pkval*
                                             dcb.get_element(ii,kk));

                        ka.accumulate_element(ii,kk,pkval*
                                         _gr_opa_dens.get_element(jj,ll));
                        ka.accumulate_element(jj,ll,pkval*
                                         _gr_opa_dens.get_element(ii,kk));

                        kb.accumulate_element(ii,kk,pkval*
                                         _gr_opb_dens.get_element(jj,ll));
                        kb.accumulate_element(jj,ll,pkval*
                                         _gr_opb_dens.get_element(ii,kk));

                        /*
                         * and to K2 of G(il)
                         *
                         * pkval = -0.25 * (ijkl)
                         */
                        lij=ioff(ii)+ll;
                        lkl=IOFF(kk,jj);

                        gca.accumulate_element(ii,ll,-pkval*
                                             dca.get_element(kk,jj));
                        gca.accumulate_element(kk,jj,-pkval*
                                             dca.get_element(ii,ll));

                        gcb.accumulate_element(ii,ll,-pkval*
                                             dcb.get_element(kk,jj));
                        gcb.accumulate_element(kk,jj,-pkval*
                                             dcb.get_element(ii,ll));

                        ka.accumulate_element(ii,ll,pkval*
                                         _gr_opa_dens.get_element(kk,jj));
                        ka.accumulate_element(kk,jj,pkval*
                                         _gr_opa_dens.get_element(ii,ll));

                        kb.accumulate_element(ii,ll,pkval*
                                         _gr_opb_dens.get_element(kk,jj));
                        kb.accumulate_element(kk,jj,pkval*
                                         _gr_opb_dens.get_element(ii,ll));
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

  delete[] shnfunc;

  dca->scale(0.5);
  dca->scale_diagonal(2.0);
  
  dcb->scale(0.5);
  dcb->scale_diagonal(2.0);
  
  _gr_opa_dens->scale(0.5);
  _gr_opa_dens->scale_diagonal(2.0);
  
  _gr_opb_dens->scale(0.5);
  _gr_opb_dens->scale_diagonal(2.0);
  
  gca.accumulate(_gr_hcore);
  gcb.accumulate(_gr_hcore);
  
  double h11 = 0;
  double h22 = 0;
  double h12 = 0;
  double h21 = 0;
  for (i=0; i < basis()->nbasis(); i++) {
    for (int j=0; j < i; j++) {
      h11 += dca.get_element(i,j)*
             (_gr_hcore.get_element(i,j)+gca.get_element(i,j));
      h22 += dcb.get_element(i,j)*
             (_gr_hcore.get_element(i,j)+gcb.get_element(i,j));
      h12 += _gr_opa_dens.get_element(i,j)*kb.get_element(i,j);
      h21 += _gr_opb_dens.get_element(i,j)*ka.get_element(i,j);
    }
    h11 += 0.5*dca.get_element(i,i)*
             (_gr_hcore.get_element(i,i)+gca.get_element(i,i));
    h22 += 0.5*dcb.get_element(i,i)*
             (_gr_hcore.get_element(i,i)+gcb.get_element(i,i));
    h12 += 0.5*_gr_opa_dens.get_element(i,i)*kb.get_element(i,i);
    h21 += 0.5*_gr_opb_dens.get_element(i,i)*ka.get_element(i,i);
  }
  
  //printf("h11 = %lf  h22 = %lf  h12 = %lf h21 = %lf\n",h11,h22,h12,h21);
  
  RefSCDimension l2 = new LocalSCDimension(2);
  RefSymmSCMatrix h = l2->create_symmmatrix();
  RefSCMatrix hv = l2->create_matrix(l2);
  RefDiagSCMatrix hl = l2->create_diagmatrix();

  h.set_element(0,0,h11);
  h.set_element(1,1,h22);
  h.set_element(1,0,h12);
  //h.print("ci matrix");
  h.diagonalize(hl,hv);
  //hv.print("ci coeffs");
  //hl.print("ci evals");
  
  ci1 = hv.get_element(0,0);
  ci2 = hv.get_element(1,0);
  double c1c2 = ci1*ci2;
  
  occa = 2*ci1*ci1;
  occb = 2*ci2*ci2;
  
  printf("ci1 = %lf ci2 = %lf\n",ci1,ci2);
  eelec = ci1*ci1*h11 + ci2*ci2*h22 + 2*c1c2*h12;
  
  dca=0;
  dcb=0;
  gca=0;
  gcb=0;

}
