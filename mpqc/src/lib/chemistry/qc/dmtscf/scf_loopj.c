/*
 * scf_loopj.c
 *
 * Copyright (C) 1996 Limit Point Systems, Inc.
 *
 * Author: Edward Seidl <seidl@janed.com>
 * Maintainer: LPS
 *
 * This file is part of the SC Toolkit.
 *
 * The SC Toolkit is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Library General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * The SC Toolkit is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public License
 * along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
 * the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * The U.S. Government is granted a limited license as per AL 91-7.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
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

#include <chemistry/qc/dmtscf/scf_loopj.gbl>
#include <chemistry/qc/dmtscf/scf_loopj.lcl>

/*************************************************************************
 *
 * form the skeleton J matrix using the loop driven algorithm.  this actually
 * forms delta(J) first, and then adds delta(J) to Jmat.
 *
 * input:
 *   centers  = pointer to initialized centers struct
 *   scf_info = pointer to initialized scf struct
 *   sym_info = pointer to initialized sym struct
 *   Jmat     = scattered dmt matrix containing skeleton G from last iteration
 *   DPmat    = scattered dmt matrix containing density difference matrix
 *   SScr1    = scattered dmt scratch matrix
 *   SScr2    = scattered dmt scratch matrix
 *   mgdbuff  = pointer to array used by integral routines to return 2ei's
 *   outfile  = FILE pointer to output
 *
 * on return:
 *   Jmat contains the J matrix for this iteration
 *
 * return 0 on success, -1 on failure
 */

GLOBAL_FUNCTION int
scf_make_j_l(centers,scf_info,sym_info, Jmat,DPmat,SScr1,mgdbuff,outfile)
centers_t *centers;
scf_struct_t *scf_info;
sym_struct_t *sym_info;
dmt_matrix Jmat;
dmt_matrix DPmat;
dmt_matrix SScr1;
double *mgdbuff;
FILE *outfile;
{
  int i,j;
  int *ib,*jb,*kb,*lb;
  int *isz,*jsz,*ksz,*lsz;
  int nl;
  int use_symmetry = (sym_info->g > 1) ? 1 : 0;
  int nlocal = dmt_nlocal(DPmat);

  double maxpkl, maxpijkl;
  double tnint = 0.0;
  double *maxp = (double *) malloc(sizeof(double)*nlocal);

  struct mgd gdb;
  loop_t *loop;

  assert(dmt_distribution(Jmat) == SCATTERED);
  assert(dmt_distribution(DPmat) == SCATTERED);
  assert(dmt_distribution(SScr1) == SCATTERED);

  check_alloc(maxp,"maxp");

  tim_enter("scf_mkjl");

 /* make sure "erep" is entered on all nodes */
  tim_enter("erep");
  tim_exit("erep");

 /* fill in maxp array */
  if (scf_info->eliminate) fill_max(DPmat,maxp);

 /* set up some pointers */
  ib = &gdb.si; jb = &gdb.sj; kb = &gdb.sk; lb = &gdb.sl;
  isz = &gdb.isz; jsz = &gdb.jsz; ksz = &gdb.ksz; lsz = &gdb.lsz;

  dmt_fill(SScr1,0.0);

 /* put SScr1, maxp, and DPmat on the loop.  when done SScr1 will contain
  * delta(J)
  */
  loop = dmt_ngl_create("%m %armr",SScr1,maxp,DPmat);

  while (dmt_ngl_next(loop)) {
    dmt_ngl_create_inner(loop,0);
    while (dmt_ngl_next_inner_m(loop,kb,ksz,lb,lsz,&gdb.glp)) {
      dmt_ngl_find_am(loop,1,*kb,*lb,&maxpkl,&gdb.plp);

      for (nl=0; nl < nlocal; nl++) {
        dmt_get_block_dsc(DPmat,nl,ib,isz,jb,jsz,&gdb.ploc);
        dmt_get_block_dsc( Jmat,nl,ib,isz,jb,jsz,&gdb.gloc);

        maxpijkl=(maxp[nl] > maxpkl) ? maxp[nl]:maxpkl;

        if ((*kb > *ib) || (*kb == *ib && *lb > *jb)) continue;

        if (use_symmetry) if (!sym_info->p1[*ib]) continue;

        mgd_int_loop(centers,scf_info,sym_info,&gdb,mgdbuff,maxpijkl,&tnint);
      }
    }
  }

 /* we are finished with loop */
  dmt_ngl_kill(loop);

 /* sum scr1 and scr2 into gmats */
  dmt_sum(SScr1,Jmat);

 /* diagonal blocks need redundant elements, fill these in now */
  for (nl=0; nl < nlocal; nl++) {
    dmt_get_block_dsc( Jmat,nl,ib,isz,jb,jsz,&gdb.gloc);

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
      fprintf(outfile,"  %8.0f integrals in scf_make_j_l\n",tnint);
  }

  tim_exit("scf_mkjl");

  free(maxp);

  int_reduce_storage_threshold();

  return 0;
}


LOCAL_FUNCTION void
mgd_int_loop(centers,scf_info,sym_info,gdb,mgdbuff,maxpijkl,tnint)
centers_t *centers;
scf_struct_t *scf_info;
sym_struct_t *sym_info;
struct mgd *gdb;
double *mgdbuff;
double maxpijkl;
double *tnint;
{
  int g,leavel,nijkl;
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
  int index;
  int use_symmetry=(sym_info->g >1);
  int imax,cpmax;
  int inttol = (int) ((double) -(scf_info->intcut)/log10(2.0));

  double pki_int;
  double qijkl=1.0;
  double tol = pow(2.0,-126.0);
  double linv = 1.0/log(2.0);


  cpmax = (maxpijkl>tol) ? (int) (log(maxpijkl)*linv) : 
                           (int) (log(tol)*linv);

  i=gdb->si; j=gdb->sj; k=gdb->sk; l=gdb->sl;

  s1=i; s2=j; s3=k; s4=l;

  if (scf_info->eliminate) {
    imax = scf_erep_bound(s1,s2,s3,s4);
    if (imax+cpmax < inttol) return;
  }

  if (use_symmetry) {
    ij=IOFF(s1,s2);
    if (!sym_info->lamij[ij]) return;

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
      if(gijkl > ijkl) leavel=1;
      if(gijkl == ijkl) nijkl++;
    }
    if (leavel) return;
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

            i2=i1; j2=j1; k2=k1; l2=l1;
            ii=bf1; jj=bf2; kk=bf3; ll=bf4;

            lij = ii*gdb->jsz+jj;
            lkl = kk*gdb->lsz+ll;

            if((IOFF(i2,j2))==(IOFF(k2,l2))) pki_int *= 0.5;

            gdb->gloc[lij] += gdb->plp[lkl]*pki_int;
            gdb->glp[lkl] += gdb->ploc[lij]*pki_int;
          }
          index++;
        }
      }
    }
  }
  (*tnint)+= (double) (n1*n2*n3*n4);
}

LOCAL_FUNCTION void
fill_max(DPmat,maxp)
dmt_matrix DPmat;
double *maxp;
{
  int nl,ii;
  int ib,jb,isz,jsz;
  int nlocal=dmt_nlocal(DPmat);
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

/*************************************************************************
 * Local Variables:
 * mode: c
 * eval: (c-set-style "ETS")
 * End:
 */
