
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <tmpl.h>
#include <math/scmat/matrix.h>
#include <math/scmat/local.h>
#include <util/group/picl.h>
#include <math/dmt/libdmt.h>
#include <util/misc/libmisc.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>

#define ioff(i) (((i)*((i)+1))>>1)
#define IOFF(a,b) (((a)>(b))?(ioff(a)+(b)):(ioff(b)+(a)))

#include <chemistry/qc/dmtscf/scf_dmt.h>

#include <chemistry/qc/mbpt/gmat.h>

static void
form_max_dens(double *dpmat, signed char *maxp,
              centers_t *centers, int nbasis);

static int
make_g_d_nor(scf_struct_t *scf_info, centers_t *centers, RefSymmSCMatrix& Gmat,
         double *DPmat, double *mgdbuff, FILE *outfile);

/* this is a pointer to the buffer which holds the integrals */
static double *mbpt_gmat_intbuf;

/**************************************************************************
 *
 * calculate the closed shell G matrix
 * assume all matrices are held locally -- IMBN
 *
 * input:
 *   centers  = pointer to initialized centers struct
 *   Gmat     = matrix containing old G matrix
 *   DPmat    = matrix containing density diff matrix
 *   outfile  = FILE pointer to output (can be null if no output desired)
 *
 * on return:
 *   Gmat contains the new G matrix
 *
 * return 0 on success and -1 on failure
 */

int
mbpt_make_gmat(scf_struct_t *scf_info, centers_t *centers,
               RefSymmSCMatrix& Gmat,
               double *DPmat, FILE *outfile)
{
  int errcod;

 /* get some timing info */
  tim_enter("gmat");

//  if (scf_info->load_bal) {
//    errcod = make_g_d_lb(scf_info, centers, Gmat, DPmat, mbpt_gmat_intbuf,
//                             outfile);
//    } 
//  else {
    errcod = make_g_d_nor(scf_info, centers, Gmat, DPmat, mbpt_gmat_intbuf,
                          outfile);
//    }

  if (errcod != 0) {
    fprintf(stderr,"mbpt_gmat: trouble forming gmat 3\n");
    return -1;
    }

  tim_exit("gmat");

  return 0;
}

/************************************************************************
 *
 * Form the vector maxp; each element of maxp is the 2-based log of the
 * largest element (absolute value) in a block of the density matrix
 * (DPmat). The density matrix is of dimension nbasis x nbasis
 *
 ************************************************************************/

static void
form_max_dens(double *DPmat, signed char *maxp,
                   centers_t *centers, int nbasis)
{

  int i, j, k, l, ij;
  int isize, jsize, ioffset, joffset;
  double linv = 1.0/log(2.0);
  double tol = pow(2.0,-126.0);
  double ftmp, tmp;
  double *dpmat_ptr;

  for (i=0; i<centers->nshell; i++) {
    isize = INT_SH_NFUNC(centers,i);
    ioffset = centers->func_num[i];
    for (j=0; j<=i; j++) {
      jsize = INT_SH_NFUNC(centers,j);
      joffset = centers->func_num[j];
      tmp = 0.0;
      for (k=0; k<isize; k++) {
        dpmat_ptr = &DPmat[nbasis*(ioffset+k) + joffset];
        for (l=0; l<jsize; l++) {
          ftmp = fabs(*dpmat_ptr++);
          if (ftmp > tmp) tmp = ftmp;
          }
        }
      tmp = (tmp > tol) ? tmp : tol;
      ij = i*(i+1)/2 +j;
      maxp[ij] = (signed char) (log(tmp)*linv);
             /* log(tmp)/linv equals the 2-based log of tmp */
      }
    }

}

/***************************************************************************
 *
 * Given a centers struct and an scf struct, initialize some stuff needed
 * to form the G matrix;
 */

int
mbpt_init_gmat(centers_t *centers, scf_struct_t *scf_info, double *intbuf)
{

//  mbpt_gmat_intbuf =
//    int_initialize_erep(flags,0,centers,centers,centers,centers);
  mbpt_gmat_intbuf = intbuf;
  if (!mbpt_gmat_intbuf) {
    fprintf(stderr,"mbpt_init_gmat:  int_initialize_erep() failed\n");
    return -1;
  }

  int_storage(scf_info->int_store);

  return 0;
}

/*************************************************************************
 *
 * frees memory allocated for gmat construction (cf. init_gmat)
 */

void
mbpt_done_gmat(centers_t *centers, scf_struct_t *scf_info)
{
  int_done_storage();
}

//extern signed char *scf_bnd_Qvec;
extern signed char *int_Qvec;

/**************************************************************************
 *
 */

static int
make_g_d_nor(scf_struct_t *scf_info, centers_t *centers, RefSymmSCMatrix& Gmat,
         double *DPmat, double *mgdbuff, FILE *outfile)
{
  int tmax,imax,cpmax,pmaxijk;
  int pmaxik,pmaxjk,pmaxij,Qvecij;
  int i,j,k,l;
  int ij,kl;
  int n1,n2,n3,n4;
  int e12,e34,e13e24,e_any;
  int bf1,bf2,bf3,bf4;
  int i1,j1,k1,l1;
  int i2,j2,k2,l2;
  int ii,jj,kk,ll;
  int ij1;
  int lij,lkl;
  int index;
  int int_index,kindex;
  int nproc=numnodes0();
  int me=mynode0();
  int s1,s2,s3,s4;
  int inttol = (int) ((double) -(scf_info->intcut)/log10(2.0));
  int nshellt;

  double tnint=0.0;
  double pki_int,value;
  double *gtmp=NULL, *ptmp=NULL;
  double *dpmat_ptr;

  char *shnfunc=NULL;
  signed char *maxp=NULL;


  // Scale DPmat; this is necessary when using the gmat formation
  // program from scf (modified slightly), since this program assumes
  // that the off-diagonal elements have been scaled by a factor of 2.0
  dpmat_ptr = DPmat;
  for (i=0; i<scf_info->nbfao; i++) {
    for (j=0; j<scf_info->nbfao; j++) {
      if (i != j) *dpmat_ptr++ *= 2.0;
      else dpmat_ptr++;
      }
    }

  // Allocate and assign maxp
  if (scf_info->eliminate) {
    nshellt = centers->nshell*(centers->nshell+1)/2;
    maxp = (signed char*) malloc(sizeof(signed char)*nshellt);
    if (!(maxp)) {
      fprintf(stderr,"mkgdlb: could not malloc maxp\n");
      return -1;
      }
    form_max_dens(DPmat, maxp, centers, scf_info->nbfao);
    }

  // Allocate and assign ptmp (contains lower triangle of DPmat
  ptmp = (double*) malloc(sizeof(double)*scf_info->nbatri);
  if (!(ptmp)) {
    fprintf(stderr,"mkgdlb: could not malloc ptmp\n");
    return -1;
    }
  for (i=0; i<scf_info->nbfao; i++) {
    dpmat_ptr = &DPmat[i*scf_info->nbfao];
    for (j=0; j<=i; j++) {
      ptmp[i*(i+1)/2 + j] = *dpmat_ptr++;
      }
    }


  // "Unscale" DPmat to get the original DPmat
  dpmat_ptr = DPmat;
  for (i=0; i<scf_info->nbfao; i++) {
    for (j=0; j<scf_info->nbfao; j++) {
      if (i != j)  *dpmat_ptr++ *= 0.50;
      else dpmat_ptr++;
      }
    }

  // Allocate and initialize gtmp
  gtmp = (double *) malloc(sizeof(double)*scf_info->nbatri);
  for (i=0; i<scf_info->nbatri; i++) gtmp[i] = 0.0;
 
  // Allocate and assign shnfunc
  shnfunc = (char *) malloc(centers->nshell);
  if (!shnfunc) {
    fprintf(stderr,"make_g_d_lb: could not malloc shnfunc\n");
    return -1;
    }
  for (i=0; i < centers->nshell; i++) shnfunc[i]=INT_SH_NFUNC((centers),i);


  /********************************************************
   * Start the actual formation of the G matrix:          *
   * Loop over all shells, calculate a bunch of integrals *
   * from each shell quartet, and stick those integrals   *
   * where they belong                                    *
   ********************************************************/


  kindex=int_index=0;
  for (i=0; i<centers->nshell; i++) {

    for (j=0; j<=i; j++) {
      ij = ioff(i)+j;
//    Qvecij=(int)scf_bnd_Qvec[ij];
      Qvecij=(int)int_Qvec[ij];
      if(scf_info->eliminate) pmaxij=maxp[ij];

      for (k=0; k<=i; k++,kindex++) {
        if(kindex%nproc!=me) {
          continue;
          }

        kl=ioff(k);
        if(scf_info->eliminate) {
          pmaxijk=pmaxij;
          if((pmaxik=maxp[(ioff(i)+k)]-2)>pmaxijk) pmaxijk=pmaxik;
          if((pmaxjk=maxp[IOFF(j,k)]-2)>pmaxijk) pmaxijk=pmaxjk;
          }

        for (l=0; l<=(k==i?j:k); l++) {

//        imax = (int) scf_bnd_Qvec[kl]+Qvecij;
          imax = (int) int_Qvec[kl]+Qvecij;

          if(scf_info->eliminate) {
            cpmax = (maxp[kl]>pmaxijk) ? maxp[kl] : pmaxijk;
            if((tmax=maxp[(ioff(i)+l)]-2)>cpmax) cpmax=tmax;
            if((tmax=maxp[IOFF(j,l)]-2)>cpmax) cpmax=tmax;

            if(cpmax+imax < inttol) {
        //    /* If we are trying to save integrals on this node, then
        //     * int_index must be incremented now. */
        //    if (scf_info->int_store) int_index++;
              kl++;
              continue;
              }
            }

            s1 = i; s2 = j; s3 = k; s4 = l;

            int_erep(INT_EREP|INT_NOBCHK|INT_NOPERM,&s1,&s2,&s3,&s4);

            n1 = shnfunc[s1];
            n2 = shnfunc[s2];
            n3 = shnfunc[s3];
            n4 = shnfunc[s4];

           // Shell equivalence information
            e12    = (s2==s1);
            e13e24 = (s3==s1) && (s4==s2);
            e34    = (s4==s3);

            index = 0;

            e_any = (e12||e13e24||e34);
            if(e_any) {
              for (bf1=0; bf1<=INT_MAX1(n1) ; bf1++) {
                i2 = centers->func_num[s1] + bf1;

                for (bf2=0; bf2<=INT_MAX2(e12,bf1,n2) ; bf2++) {
                  j2 = centers->func_num[s2] + bf2;
                  if(i2>=j2) { i1=i2; j1=j2; }
                  else { i1=j2; j1=i2; }
                  ij1=ioff(i1)+j1;

                  for (bf3=0; bf3<=INT_MAX3(e13e24,bf1,n3) ; bf3++) {
                    k2 = centers->func_num[s3] + bf3;

                    for (bf4=0;bf4<=INT_MAX4(e13e24,e34,bf1,bf2,bf3,n4);bf4++){
                      if (INT_NONZERO(mgdbuff[index])) {
                        l2 = centers->func_num[s4] + bf4;

                        if(k2>=l2) { k1=k2; l1=l2; }
                        else { k1=l2; l1=k2; }

                        if(ij1 >= ioff(k1)+l1) {
                          ii = i1; jj = j1; kk = k1; ll = l1;
                          }
                        else {
                          ii = k1; jj = l1; kk = i1; ll = j1;
                          }

                        pki_int = mgdbuff[index];

                        if (jj == kk) {
                          if (ii == jj || kk == ll) {
                            lij=ioff(ii)+jj;
                            lkl=ioff(kk)+ll;
                            value=(lij==lkl)? 0.25*pki_int: 0.5*pki_int;
                            gtmp[lij] += ptmp[lkl]*value;
                            gtmp[lkl] += ptmp[lij]*value;
                            }
                          else {
                            lij=ioff(ii)+jj;
                            lkl=ioff(kk)+ll;
                            value=(lij==lkl)? 0.375*pki_int: 0.75*pki_int;
                            gtmp[lij] += ptmp[lkl]*value;
                            gtmp[lkl] += ptmp[lij]*value;

                            lij=ioff(ii)+ll;
                            lkl=IOFF(kk,jj);
                            value=(lij==lkl)? 0.25*pki_int: 0.5*pki_int;
                            gtmp[lij] -= ptmp[lkl]*value;
                            gtmp[lkl] -= ptmp[lij]*value;
                            }
                          }
                        else if (ii == kk || jj == ll) {
                          lij=ioff(ii)+jj;
                          lkl=ioff(kk)+ll;
                          value=(lij==lkl)? 0.375*pki_int: 0.75*pki_int;
                          gtmp[lij] += ptmp[lkl]*value;
                          gtmp[lkl] += ptmp[lij]*value;

                          lij=ioff(ii)+kk;
                          lkl=IOFF(jj,ll);
                          value=(lij==lkl)? 0.25*pki_int : 0.5*pki_int;
                          gtmp[lij] -= ptmp[lkl]*value;
                          gtmp[lkl] -= ptmp[lij]*value;
                          }
                        else {
                          lij=ioff(ii)+jj;
                          lkl=ioff(kk)+ll;
                          value=(lij==lkl)? 0.5*pki_int : pki_int;
                          gtmp[lij] += ptmp[lkl]*value;
                          gtmp[lkl] += ptmp[lij]*value;

                          lij=ioff(ii)+kk;
                          lkl=IOFF(jj,ll);
                          value=(lij==lkl)? 0.125*pki_int: 0.25*pki_int;
                          gtmp[lij] -= ptmp[lkl]*value;
                          gtmp[lkl] -= ptmp[lij]*value;

                          if((ii != jj) && (kk != ll)) {
                            lij=ioff(ii)+ll;
                            lkl=IOFF(kk,jj);
                            value=(lij==lkl)? 0.125*pki_int: 0.25*pki_int;
                            gtmp[lij] -= ptmp[lkl]*value;
                            gtmp[lkl] -= ptmp[lij]*value;
                            }
                          }
                        }
                      index++;
                      }
                    }
                  }
                }
              }
            else {
              for (bf1=0; bf1<n1 ; bf1++) {
                i2 = centers->func_num[s1] + bf1;

                for (bf2=0; bf2<n2 ; bf2++) {
                  j2 = centers->func_num[s2] + bf2;
                  if(i2>=j2) { i1=i2; j1=j2; }
                  else { i1=j2; j1=i2; }
                  ij1=ioff(i1)+j1;

                  for (bf3=0; bf3<n3 ; bf3++) {
                    k2 = centers->func_num[s3] + bf3;

                    for (bf4=0; bf4<n4; bf4++) {
                      if (INT_NONZERO(mgdbuff[index])) {
                        l2 = centers->func_num[s4] + bf4;

                        if(k2>=l2) { k1=k2; l1=l2; }
                        else { k1=l2; l1=k2; }

                        if(ij1 >= ioff(k1)+l1) {
                          ii = i1; jj = j1; kk = k1; ll = l1;
                          }
                        else {
                          ii = k1; jj = l1; kk = i1; ll = j1;
                          }

                        pki_int = mgdbuff[index];

                        if (jj == kk) {
                          lij=ioff(ii)+jj;
                          lkl=ioff(kk)+ll;
                          value=0.75*pki_int;
                          gtmp[lij] += ptmp[lkl]*value;
                          gtmp[lkl] += ptmp[lij]*value;

                          lij=ioff(ii)+ll;
                          lkl=IOFF(kk,jj);
                          value=0.5*pki_int;
                          gtmp[lij] -= ptmp[lkl]*value;
                          gtmp[lkl] -= ptmp[lij]*value;
                          }
                        else if (ii == kk || jj == ll) {
                          lij=ioff(ii)+jj;
                          lkl=ioff(kk)+ll;
                          value=0.75*pki_int;
                          gtmp[lij] += ptmp[lkl]*value;
                          gtmp[lkl] += ptmp[lij]*value;

                          lij=ioff(ii)+kk;
                          lkl=IOFF(jj,ll);
                          value=0.5*pki_int;
                          gtmp[lij] -= ptmp[lkl]*value;
                          gtmp[lkl] -= ptmp[lij]*value;
                          }
                        else {
                          lij=ioff(ii)+jj;
                          lkl=ioff(kk)+ll;
                          value=pki_int;
                          gtmp[lij] += ptmp[lkl]*value;
                          gtmp[lkl] += ptmp[lij]*value;
  
                          lij=ioff(ii)+kk;
                          lkl=IOFF(jj,ll);
                          value*=0.25;
                          gtmp[lij] -= ptmp[lkl]*value;
                          gtmp[lkl] -= ptmp[lij]*value;

                          lij=ioff(ii)+ll;
                          lkl=IOFF(kk,jj);
                          gtmp[lij] -= ptmp[lkl]*value;
                          gtmp[lkl] -= ptmp[lij]*value;
                          }
                        }
                      index++;
                      }
                    }
                  }
                }
              }
            tnint += (double) (n1*n2*n3*n4);
          kl++;
          int_index++;
          } // exit l loop
        }   // exit k loop
      }     // exit j loop
    }       // exit i loop

  if (scf_info->print_flg & 4) {
    gsum0(&tnint,1,5,mtype_get(),0);
    if (mynode0()==0) {
      fprintf(outfile,"  %8.0f integrals in make_g_d\n",tnint);
      fflush(outfile);
    }
  }


  // Sum up contributions to gtmp
  gop1(gtmp,scf_info->nbatri,ptmp,'+',mtype_get());


  // Put gtmp back into Gmat
  for (i=0; i<scf_info->nbfao; i++) {
    for (j=0; j<=i; j++) {
      ij = i*(i+1)/2 + j;
      Gmat->set_element(i,j,gtmp[ij]);
//    Gmat->set_element(j,i,gtmp[ij]);  don't do this - only lower triangle
      }
    }


  // Free up memory
  if (gtmp) free(gtmp);
  if (maxp) free(maxp);
  if (ptmp) free(ptmp);
  if (shnfunc) free(shnfunc);

  return 0;
}
