
#include <math/array/math_lib.h>
#include <chemistry/qc/integral/integralv2.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/scf/grscf.h>

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
  
  int_done_erep();
  int_done_offsets2(centers,centers,centers,centers);
  int_done_storage();
  free_centers(centers);
  free(centers);
  delete[] shnfunc;

  _fock.assign(_gr_gmat);
  _fock.accumulate(_gr_hcore);

}
