
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <tmpl.h>
#include <util/group/picl.h>
#include <math/array/math_lib.h>
#include <util/misc/libmisc.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>
#include <chemistry/qc/dmtscf/mgd.h>

#include <chemistry/qc/dmtscf/scf.h>
#include <chemistry/qc/dmtscf/scf_bnd.gbl>

#include <chemistry/qc/dmtscf/scf_loopg.gbl>
#include <chemistry/qc/dmtscf/scf_loopg.lcl>

/*************************************************************************
 *
 * form the skeleton K matrix using the loop driven algorithm.  this actually
 * forms delta(K) first, and then adds delta(K) to Kmat.
 *
 * input:
 *   centers  = pointer to initialized centers struct
 *   scf_info = pointer to initialized scf struct
 *   sym_info = pointer to initialized sym struct
 *   Kmat     = scattered dmt matrix containing skeleton G from last iteration
 *   DPmat    = scattered dmt matrix containing density difference matrix
 *   SScr1    = scattered dmt scratch matrix
 *   mgdbuff  = pointer to array used by integral routines to return 2ei's
 *   outfile  = FILE pointer to output
 *
 * on return:
 *   Kmat contains the K matrix for this iteration
 *
 * return 0 on success, -1 on failure
 */

GLOBAL_FUNCTION int
scf_make_k_l(centers,scf_info,sym_info,Kmat,DPmat,SScr1,mgdbuff,outfile)
centers_t *centers;
scf_struct_t *scf_info;
sym_struct_t *sym_info;
dmt_matrix Kmat;
dmt_matrix DPmat;
dmt_matrix SScr1;
double *mgdbuff;
FILE *outfile;
{
  int i,j;
  int *ib,*jb,*kb,*lb;
  int *isz,*jsz,*ksz,*lsz;
  int nl;
  int use_symmetry=(sym_info->g > 1)?1:0;
  int nlocal=dmt_nlocal(DPmat);
  double tnint=0.0;
  struct mgd gdb;
  double maxpkl;
  double maxpijkl;

  double *maxp = (double *) malloc(sizeof(double)*nlocal);
  loop_t *loop;

  assert(dmt_distribution(Kmat) == SCATTERED);
  assert(dmt_distribution(DPmat) == SCATTERED);
  assert(dmt_distribution(SScr1) == SCATTERED);

  check_alloc(maxp,"maxp");

  tim_enter("scf_mkkl");

 /* make sure "erep" is entered on all nodes */
  tim_enter("erep");
  tim_exit("erep");

 /* fill in maxp array */

  if (scf_info->eliminate) fill_max(DPmat,maxp);

 /* set up some pointers */
  ib = &gdb.si; jb = &gdb.sj; kb = &gdb.sk; lb = &gdb.sl;
  isz = &gdb.isz; jsz = &gdb.jsz; ksz = &gdb.ksz; lsz = &gdb.lsz;

  dmt_fill(SScr1,0.0);

 /* put scr1, maxp, and DPmat on the loop */

  loop = dmt_ngl_create("%m %armr",SScr1,maxp,DPmat);

  while (dmt_ngl_next(loop)) {
    dmt_ngl_create_inner(loop,0);
    while (dmt_ngl_next_inner_m(loop,kb,ksz,lb,lsz,&gdb.glp)) {
      dmt_ngl_find_am(loop,1,*kb,*lb,&maxpkl,&gdb.plp);

      for (nl=0; nl < nlocal; nl++) {
        dmt_get_block_dsc(DPmat,nl,ib,isz,jb,jsz,&gdb.ploc);
        dmt_get_block_dsc( Kmat,nl,ib,isz,jb,jsz,&gdb.gloc);

        maxpijkl = (maxp[nl] > maxpkl) ? maxp[nl] : maxpkl;

        if ((*kb > *ib) || (*kb == *ib && *lb > *jb)) continue;

        if (use_symmetry) if (!sym_info->p1[*ib]) continue;

        mgd_int_loop(centers,scf_info,sym_info,&gdb,mgdbuff,maxpijkl,&tnint);
      }
    }
  }

 /* we are finished with loop */
  dmt_ngl_kill(loop);

 /* sum scr1 and scr2 into gmats */
  dmt_sum(SScr1,Kmat);

 /* diagonal blocks need redundant elements, fill these in now */
  for (nl=0; nl < nlocal; nl++) {
    dmt_get_block_dsc( Kmat,nl,ib,isz,jb,jsz,&gdb.gloc);

    /* test for diagonal block */
    if (*ib == *jb) {
      /* filled in ignored elements */
      for (i=0; i < *isz ; i++) {
        for (j=0; j < i ; j++) {
          gdb.gloc[j*(*isz)+i] = gdb.gloc[i*(*isz)+j];
        }
      }
    }
  }

  if (outfile && (scf_info->print_flg & 4)) {
    gsum0(&tnint,1,5,mtype_get(),0);
    if (mynode0()==0)
      fprintf(outfile,"  %8.0f integrals in scf_make_k_l\n",tnint);
  }

  tim_exit("scf_mkkl");

  free(maxp);

  int_reduce_storage_threshold();

  return 0;
}


LOCAL_FUNCTION VOID
mgd_int_loop(centers,scf_info,sym_info,gdb,mgdbuff,maxpijkl,tnint)
centers_t *centers;
scf_struct_t *scf_info;
sym_struct_t *sym_info;
struct mgd *gdb;
double *mgdbuff;
double maxpijkl;
double *tnint;
{
  int g,m,leavel,nijkl;
  int num;
  int s1,s2,s3,s4;
  int n1,n2,n3,n4;
  int e12,e34,e13e24;
  int bf1,bf2,bf3,bf4;
  int i,j,k,l;
  int i1,j1,k1,l1;
  int i2,j2,k2,l2;
  int ii,jj,kk,ll;
  int ij,kl,ijkl;
  int lij,lkl;
  int gi,gj,gk,gl,gij,gkl,gijkl;
  int index,keql;
  int imax,cpmax,scale;
  int use_symmetry=(sym_info->g >1);
  int inttol = (int) ((double) -(scf_info->intcut)/log10(2.0));

  double qijkl=1.0;
  double pki_int;
  double tol = pow(2.0,-126.0);
  double linv = 1.0/log(2.0);

  cpmax = (maxpijkl>tol) ? (int) (log(maxpijkl)*linv) : 
                           (int) (log(tol)*linv);

  i=gdb->si; j=gdb->sj; k=gdb->sk; l=gdb->sl;

  if ((j==k) && ((i==j) || (k==l)))
    num=1;
  else if ((k==l)||(j==k)||(j==l)||(i==k)||(i==j))
    num=2;
  else
    num=3;

  for (m=0; m < num ; m++) {
    switch(m) {
    case 0:
      s1=i; s2=j; s3=k; s4=l;
      scale=-2;
      break;

    case 1:
      if(j==k) {
        s1=i; s2=l; s3=k; s4=j;
      } else {
        s1=i; s2=k; s3=j; s4=l;
      }
      scale=-2;
      break;

    case 2:
      s1=i; s2=l; s3=k; s4=j;
      scale=-2;
    }

    if (scf_info->eliminate) {
      imax = scf_erep_bound(s1,s2,s3,s4);
      if (scale+imax+cpmax < inttol) continue;
    }

    if (use_symmetry) {
      ij=IOFF(s1,s2);
      if (!sym_info->lamij[ij]) continue;

      kl=IOFF(s3,s4);
      ijkl=IOFF(ij,kl);

      nijkl=leavel=0;
      for (g=0; g < sym_info->g ; g++) {
        gi = sym_info->shell_map[s1][g];
        gj = sym_info->shell_map[s2][g];
        gk = sym_info->shell_map[s3][g];
        gl = sym_info->shell_map[s4][g];
        gij = IOFF(gi,gj);
        gkl = IOFF(gk,gl);
        gijkl = IOFF(gij,gkl);
        if (gijkl > ijkl) leavel=1;
        if (gijkl == ijkl) nijkl++;
      }
      if (leavel) continue;
      qijkl = (double) sym_info->g/nijkl;
    }

    n1 = INT_SH_NFUNC((centers),s1);
    n2 = INT_SH_NFUNC((centers),s2);
    n3 = INT_SH_NFUNC((centers),s3);
    n4 = INT_SH_NFUNC((centers),s4);

   /* Shell equivalency information. */
    e12    = (s2==s1);
    e13e24 = (s3==s1) && (s4==s2);
    e34    = (s4==s3);

    tim_enter("erep");
    int_erep(INT_EREP|INT_NOBCHK|INT_NOPERM,&s1,&s2,&s3,&s4);
    tim_exit("erep");

    index = 0;
    for (bf1=0; bf1<=INT_MAX1(n1); bf1++) {
      for (bf2=0; bf2<=INT_MAX2(e12,bf1,n2); bf2++) {
        for (bf3=0; bf3<=INT_MAX3(e13e24,bf1,n3); bf3++) {
          for (bf4=0; bf4<=INT_MAX4(e13e24,e34,bf1,bf2,bf3,n4); bf4++) {
            if (INT_NONZERO(mgdbuff[index])) {
              i1 = centers->func_num[s1] + bf1;
              j1 = centers->func_num[s2] + bf2;
              k1 = centers->func_num[s3] + bf3;
              l1 = centers->func_num[s4] + bf4;

              pki_int = (double) qijkl*mgdbuff[index];

              switch(m) {
              case 0:
                i2=i1; j2=j1; k2=k1; l2=l1;
                ii=bf1; jj=bf2; kk=bf3; ll=bf4;
                break;

              case 1:
                if(j==k) {
                  i2=i1; j2=l1; k2=k1; l2=j1;
                  ii=bf1; jj=bf4; kk=bf3; ll=bf2;
                } else {
                  i2=i1; j2=k1; k2=j1; l2=l1;
                  ii=bf1; jj=bf3; kk=bf2; ll=bf4;
                }
                break;

              case 2:
                i2=i1; j2=l1; k2=k1; l2=j1;
                ii=bf1; jj=bf4; kk=bf3; ll=bf2;
              }

              lij = ii*gdb->jsz+jj;
              lkl = kk*gdb->lsz+ll;

              if (j2==k2 && (i2==j2||k2==l2)) {
                put_it(scf_info,lij,lkl,i2,j2,k2,l2,pki_int,gdb,5);
              } else switch(num) {
              case 1:
                keql=(k==l);

                if (i2==k2 || j2==l2) {
                  put_it(scf_info,lij,lkl,i2,j2,k2,l2,pki_int,gdb,3);

                  if (keql) {
                    lij=ii*gdb->ksz+kk;
                    lkl=(MAX0(jj,ll))*gdb->lsz+(MIN0(jj,ll));
                  } else {
                    lij=(MAX0(ii,kk))*gdb->ksz+(MIN0(ii,kk));
                    lkl=jj*gdb->lsz+ll;
                  }
                  put_it(scf_info,lij,lkl,i2,k2,j2,l2,pki_int,gdb,4);

                } else if (j2==k2) {
                  put_it(scf_info,lij,lkl,i2,j2,k2,l2,pki_int,gdb,3);

                  if(keql) {
                    lij=ii*gdb->lsz+ll;
                    lkl=(MAX0(jj,kk))*gdb->ksz+(MIN0(jj,kk));
                  } else {
                    lkl=ii*gdb->lsz+ll;
                    lij=(MAX0(jj,kk))*gdb->ksz+(MIN0(jj,kk));
                  }
                  put_it(scf_info,lij,lkl,i2,l2,j2,k2,pki_int,gdb,4);

                } else if (i2==j2 || k2==l2) {
                  if (keql) {
                    lij=ii*gdb->ksz+kk;
                    lkl=(MAX0(jj,ll))*gdb->lsz+(MIN0(jj,ll));
                  } else {
                    lij=(MAX0(ii,kk))*gdb->ksz+(MIN0(ii,kk));
                    lkl=jj*gdb->lsz+ll;
                  }
                  put_it(scf_info,lij,lkl,i2,k2,j2,l2,pki_int,gdb,2);

                } else {
                  if (k==l) {
                    lij=ii*gdb->jsz+kk;
                    lkl=(MAX0(jj,ll))*gdb->lsz+(MIN0(jj,ll));
                    put_it(scf_info,lij,lkl,i2,k2,j2,l2,pki_int,gdb,2);

                    lij=ii*gdb->lsz+ll;
                    lkl=(MAX0(jj,kk))*gdb->ksz+(MIN0(jj,kk));
                    put_it(scf_info,lij,lkl,i2,l2,j2,k2,pki_int,gdb,2);
                  } else {
                    lij=(MAX0(ii,kk))*gdb->ksz+(MIN0(ii,kk));
                    lkl=jj*gdb->lsz+ll;
                    put_it(scf_info,lij,lkl,i2,k2,j2,l2,pki_int,gdb,2);

                    lkl=ii*gdb->lsz+ll;
                    lij=(MAX0(jj,kk))*gdb->ksz+(MIN0(jj,kk));
                    put_it(scf_info,lij,lkl,i2,l2,j2,k2,pki_int,gdb,2);
                  }
                }
                break;

              case 2:
                if (!m) {
                  if (i2==k2 || j2==l2 || j2==k2) {
                    put_it(scf_info,lij,lkl,i2,j2,k2,l2,pki_int,gdb,3);
                  } else {
                    if (j==k) {
                      lij=ii*gdb->ksz+kk;
                      lkl=jj*gdb->lsz+ll;
                      put_it(scf_info,lij,lkl,i2,k2,j2,l2,pki_int,gdb,2);

                    } else if (j==l) {
                      lij=ii*gdb->lsz+ll;
                      lkl=kk*gdb->jsz+jj;
                      put_it(scf_info,lij,lkl,i2,l2,k2,j2,pki_int,gdb,2);

                    } else if (i==k) {
                      lij=kk*gdb->jsz+jj;
                      lkl=ii*gdb->lsz+ll;
                      put_it(scf_info,lij,lkl,k2,j2,i2,l2,pki_int,gdb,2);
                    }
                  }
                } else {
                  if (k==l) lkl = (MAX0(kk,ll))*gdb->lsz+(MIN0(kk,ll));
                  if (i==j) lij = (MAX0(ii,jj))*gdb->jsz+(MIN0(ii,jj));

                  if (i2==k2 || j2==l2 || j2==k2) {
                    put_it(scf_info,lij,lkl,i2,j2,k2,l2,pki_int,gdb,2);
                  } else if (i2==j2 || k2==l2) {
                    put_it(scf_info,lij,lkl,i2,j2,k2,l2,pki_int,gdb,4);
                  } else {
                    put_it(scf_info,lij,lkl,i2,j2,k2,l2,pki_int,gdb,2);

                    if (j==k) {
                      lij=ii*gdb->ksz+kk;
                      lkl=jj*gdb->lsz+ll;
                      put_it(scf_info,lij,lkl,i2,k2,j2,l2,pki_int,gdb,2);

                    } else if (j==l) {
                      lij=ii*gdb->lsz+ll;
                      lkl=kk*gdb->jsz+jj;
                      put_it(scf_info,lij,lkl,i2,l2,k2,j2,pki_int,gdb,2);

                    } else if (i==k) {
                      lij=kk*gdb->jsz+jj;
                      lkl=ii*gdb->lsz+ll;
                      put_it(scf_info,lij,lkl,k2,j2,i2,l2,pki_int,gdb,2);
                    }
                  }
                }
                break;

              case 3:
                if (!m) {
                  if (i2==k2 || j2==l2 || j2==k2) {
                    put_it(scf_info,lij,lkl,i2,j2,k2,l2,pki_int,gdb,3);
                  }
                } else {
                  if (i2==j2 || k2==l2) {
                    put_it(scf_info,lij,lkl,i2,j2,k2,l2,pki_int,gdb,4);
                  } else {
                    put_it(scf_info,lij,lkl,i2,j2,k2,l2,pki_int,gdb,2);
                  }
                }
                break;
              
              default:
                fprintf(stderr,"scf_mkgl: num=%d! Huh? \n",num);
              }
            }
            index++;
          }
        }
      }
    }

    (*tnint) += (double) (n1*n2*n3*n4);
  }
}

LOCAL_FUNCTION VOID
put_it(scf_info,lij,lkl,i2,j2,k2,l2,pki_int,gdb,iab)
scf_struct_t *scf_info;
int lij;
int lkl;
int i2;
int j2;
int k2;
int l2;
double pki_int;
struct mgd *gdb;
int iab;
{
  if ((IOFF(i2,j2))==(IOFF(k2,l2))) pki_int *= 0.5;

  switch(iab) {
  case 1:
    break;

  case 2:
    gdb->gloc[lij] -= 0.25*gdb->plp[lkl]*pki_int;
    gdb->glp[lkl] -= 0.25*gdb->ploc[lij]*pki_int;
    break;

  case 3:
    gdb->gloc[lij] -= 0.25*gdb->plp[lkl]*pki_int;
    gdb->glp[lkl] -= 0.25*gdb->ploc[lij]*pki_int;
    break;

  case 4:
    gdb->gloc[lij] -= 0.5*gdb->plp[lkl]*pki_int;
    gdb->glp[lkl] -= 0.5*gdb->ploc[lij]*pki_int;
    break;

  case 5:
    gdb->gloc[lij] -= 0.5*gdb->plp[lkl]*pki_int;
    gdb->glp[lkl] -= 0.5*gdb->ploc[lij]*pki_int;
    break;

  default:
    fprintf(stderr,"scf_loopk: you shouldn't be here\n");
  }
}


LOCAL_FUNCTION VOID
fill_max(DPmat,maxp)
dmt_matrix DPmat;
double *maxp;
{
  int nlocal=dmt_nlocal(DPmat);
  int nl,ii;
  int ib,jb,isz,jsz;
  double tmp;
  double *pblk;

  for (nl=0; nl < nlocal ; nl++) {
    dmt_get_block_dsc(DPmat,nl,&ib,&isz,&jb,&jsz,&pblk);

    tmp=0.0;
    for (ii=0; ii < isz*jsz ; ii++)
      tmp = (fabs(pblk[ii]) > tmp) ? fabs(pblk[ii]) : tmp;

    maxp[nl]=tmp;
  }
}
