/*
 * mxmv.c
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

/* $Log$
 * Revision 1.5  1996/10/25 21:00:50  etseidl
 * add copyleft notice
 *
 * Revision 1.4  1996/03/23 02:38:08  cljanss
 * Everything can now be configured with autoconf.
 *
 * Revision 1.3  1995/03/17 01:50:09  cljanss
 * Removed -I. and -I$(SRCDIR) from the default include path in
 * GlobalMakefile to avoid name conflicts with system include files.
 * Modified files under src.lib to include all files relative to src.lib.
 * Makefiles under src.bin need to add the -I. and -I$(SRCDIR) back onto
 * INCLUDE and CXXINCLUDE or make other arrangements.
 *
 * Revision 1.2  1994/08/26  17:57:47  etseidl
 * get rid of rcs id's and fix a few warnings
 *
 * Revision 1.1.1.1  1993/12/29  12:53:30  etseidl
 * SC source tree 0.1
 *
 * Revision 1.2  1992/06/17  22:10:16  jannsen
 * clean up for saber-c
 *
 * Revision 1.1.1.1  1992/03/17  16:35:34  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:35:33  seidl
 * Initial revision
 *
 * Revision 1.1  1991/12/20  16:01:42  seidl
 * Initial revision
 *
 * Revision 1.3  91/12/03  01:28:35  etseidl
 * fix comments
 * 
 * Revision 1.2  1991/12/03  00:24:50  etseidl
 * change alloc_matrix to matrixallc, free_matrix to matrixfree
 *
 * Revision 1.1  1991/12/02  23:39:44  etseidl
 * Initial revision
 * */

#define IOFF(a,b) ((a)>(b))?(a)*((a)+1)/2+(b):(b)*((b)+1)/2+(a)

#include <stdio.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/matrix.h>
#include <math/array/matrixallc.h>
#include <math/array/matrixfree.h>

#include <math/array/mxmv.gbl>
#include <math/array/mxmv.lcl>

/*
 *                                                            
 * a reasonably fast matrix multiply (at least on the DEC3100)
 * written by ETS                                            
 *                                                          
 * AF,BF,and CF are double_vectors which hold the lower triangle of
 *   square, symmetric matrices
 *
 * AF and BF are copied to temporary matrices
 *                                                        
 * ta,tb and tc indicate whether the corresponding arrays are 
 *              to be converted to their transpose           
 *                                                          
 * nr,nl,nc are the number of rows,links,and columns in the
 *          final matrices to be multiplied together          
 *          if ta=0 AF should have the dimensions nr x nl    
 *          if ta=1 AF should have the dimensions nl x nr   
 *          if tb=0 BF should have the dimensions nl x nc  
 *          if tb=1 BF should have the dimensions nc x nl
 *          if tc=0 CF should have the dimensions nr x nc    
 *          if tc=1 CF should have the dimensions nc x nr   
 *                                                         
 * add is 1 if this matrix is to be added to the one passed
 *        in as CF, 0 otherwise                           
 */

static int keep_nrv=0;
static int keep_nlv=0;
static int keep_ncv=0;
static double_matrix_t aav,bbv;


GLOBAL_FUNCTION void
math_dvxdv_dv(AF,ta,BF,tb,CF,tc,nr,nl,nc,add)
double_vector_t *AF;
int ta;
double_vector_t *BF;
int tb;
double_vector_t *CF;
int tc;
int nr;
int nc;
int nl;
int add;
{
  int i,j,ij;
  int errcod;
  double_matrix_t scr1;
  double_matrix_t scr2;

  if(nr != nc) {
    fprintf(stderr,"dvxdv_dv:\n");
    fprintf(stderr,"nr must equal nc\n");
    exit(-1);
    }

  if(nr != nl) {
    fprintf(stderr,"dvxdv_dv:\n");
    fprintf(stderr,"nr must equal nl\n");
    exit(-1);
    }

  if(nc != nl) {
    fprintf(stderr,"dvxdv_dv:\n");
    fprintf(stderr,"nc must equal nl\n");
    exit(-1);
    }

  errcod = allocbn_double_matrix(&scr1,"n1 n2",nr,nl);
  if(errcod != 0) {
    fprintf(stderr,"dvxdv_dv:\n");
    fprintf(stderr,"trouble allocating memory for scr1 matrix\n");
    exit(-1);
    }
  errcod = allocbn_double_matrix(&scr2,"n1 n2",nc,nl);
  if(errcod != 0) {
    fprintf(stderr,"dvxdv_dv:\n");
    fprintf(stderr,"trouble allocating memory for scr2 matrix\n");
    exit(-1);
    }

  for(i=ij=0; i < nr ; i++)
    for(j=0; j <= i; j++,ij++) {
      scr1.d[i][j] = scr1.d[j][i] = AF->d[ij];
      scr2.d[i][j] = scr2.d[j][i] = BF->d[ij];
      }

  math_double_mxmv(scr1.d,0,scr2.d,1,CF->d,0,nr,nl,nc,add);

  free_double_matrix(&scr2);
  free_double_matrix(&scr1);
  }

GLOBAL_FUNCTION void
math_dmxdv_dv(AF,ta,BF,tb,CF,tc,nr,nl,nc,add)
double_matrix_t *AF;
int ta;
double_vector_t *BF;
int tb;
double_vector_t *CF;
int tc;
int nr;
int nc;
int nl;
int add;
{
  int i,j,ij;
  int errcod;
  double_matrix_t scr;

  if(nr != nc) {
    fprintf(stderr,"dmxdv_dv:\n");
    fprintf(stderr,"nr must equal nc\n");
    exit(-1);
    }

  if(nc != nl) {
    fprintf(stderr,"dmxdv_dv:\n");
    fprintf(stderr,"nc must equal nl\n");
    exit(-1);
    }

  errcod = allocbn_double_matrix(&scr,"n1 n2",nc,nl);
  if(errcod != 0) {
    fprintf(stderr,"dmxdv_dv:\n");
    fprintf(stderr,"trouble allocating memory for scr matrix\n");
    exit(-1);
    }

  for(i=ij=0; i < nr ; i++)
    for(j=0; j <= i; j++,ij++)
      scr.d[i][j] = scr.d[j][i] = BF->d[ij];

  math_double_mxmv(AF->d,ta,scr.d,1,CF->d,tc,nr,nl,nc,add);

  free_double_matrix(&scr);
  }

GLOBAL_FUNCTION void
math_dvxdm_dv(AF,ta,BF,tb,CF,tc,nr,nl,nc,add)
double_vector_t *AF;
int ta;
double_matrix_t *BF;
int tb;
double_vector_t *CF;
int tc;
int nr;
int nc;
int nl;
int add;
{
  int i,j,ij;
  int errcod;
  double_matrix_t scr;

  if(nr != nc) {
    fprintf(stderr,"dvxdm_dv:\n");
    fprintf(stderr,"nr must equal nc\n");
    exit(-1);
    }

  if(nr != nl) {
    fprintf(stderr,"dvxdm_dv:\n");
    fprintf(stderr,"nr must equal nl\n");
    exit(-1);
    }

  errcod = allocbn_double_matrix(&scr,"n1 n2",nr,nl);
  if(errcod != 0) {
    fprintf(stderr,"dvxdm_dv:\n");
    fprintf(stderr,"trouble allocating memory for scr matrix\n");
    exit(-1);
    }

  for(i=ij=0; i < nr ; i++)
    for(j=0; j <= i; j++,ij++)
      scr.d[i][j] = scr.d[j][i] = AF->d[ij];

  math_double_mxmv(scr.d,0,BF->d,tb,CF->d,tc,nr,nl,nc,add);

  free_double_matrix(&scr);
  }

GLOBAL_FUNCTION void
math_dmxdm_dv(AF,ta,BF,tb,CF,tc,nr,nl,nc,add)
double_matrix_t *AF;
int ta;
double_matrix_t *BF;
int tb;
double_vector_t *CF;
int tc;
int nr;
int nc;
int nl;
int add;
{
  math_double_mxmv(AF->d,ta,BF->d,tb,CF->d,tc,nr,nl,nc,add);
  }

GLOBAL_FUNCTION void
math_double_mxmv(AF,ta,BF,tb,CF,tc,nr,nl,nc,add)
double **AF;
int ta;
double **BF;
int tb;
double *CF;
int tc;
int nr;
int nc;
int nl;
int add;
{
  int odd_nc;
  int i,j,k,ij;
  int ij1,i1j,i1j1;
  int errcod;
  double t00,t01,t10,t11;
  double **a,**b;
  double *att,*bt;
  double *at1,*bt1;

  if(!aav.d) {
    errcod = allocbn_double_matrix(&aav,"n1 n2",nr,nl);
    if(errcod != 0) {
      fprintf(stderr,"math_double_mxmv:\n");
      fprintf(stderr,"trouble allocating memory for temporary matrix aav\n");
      exit(-1);
      }
    errcod = allocbn_double_matrix(&bbv,"n1 n2",nc,nl);
    if(errcod != 0) {
      fprintf(stderr,"math_double_mxmv:\n");
      fprintf(stderr,"trouble allocating memory for temporary matrix bbv\n");
      exit(-1);
      }
    keep_nrv = nr;
    keep_nlv = nl;
    keep_ncv = nc;
    }

  if(nl > keep_nlv) {
    free_double_matrix(&aav);
    free_double_matrix(&bbv);
    keep_nlv = nl;
    keep_nrv = (nr > keep_nrv) ? nr : keep_nrv;
    keep_ncv = (nc > keep_ncv) ? nc : keep_ncv;
    errcod = allocbn_double_matrix(&aav,"n1 n2",keep_nrv,keep_nlv);
    if(errcod != 0) {
      fprintf(stderr,"math_double_mxmv:\n");
      fprintf(stderr,"trouble allocating memory for temporary matrix aav\n");
      exit(-1);
      }
    errcod = allocbn_double_matrix(&bbv,"n1 n2",keep_ncv,keep_nlv);
    if(errcod != 0) {
      fprintf(stderr,"math_double_mxmv:\n");
      fprintf(stderr,"trouble allocating memory for temporary matrix bbv\n");
      exit(-1);
      }
    }
  if(nr > keep_nrv) {
    free_double_matrix(&aav);
    keep_nrv = nr;
    errcod = allocbn_double_matrix(&aav,"n1 n2",keep_nrv,keep_nlv);
    if(errcod != 0) {
      fprintf(stderr,"math_double_mxmv:\n");
      fprintf(stderr,"trouble allocating memory for temporary matrix aav\n");
      exit(-1);
      }
    }
  if(nc > keep_ncv) {
    free_double_matrix(&bbv);
    keep_ncv = nc;
    errcod = allocbn_double_matrix(&bbv,"n1 n2",keep_ncv,keep_nlv);
    if(errcod != 0) {
      fprintf(stderr,"math_double_mxmv:\n");
      fprintf(stderr,"trouble allocating memory for temporary matrix bbv\n");
      exit(-1);
      }
    }

  odd_nc = (nc)%2;

  a=aav.d;
  if(ta)
    for(i=0; i < nr ; i++)
      for(j=0; j < nl ; j++)
        a[i][j] = AF[j][i];
  else
    a=AF;

  b=bbv.d;
  if(tb)
    b=BF;
  else
    for(i=0; i < nc ; i++)
      for(j=0; j < nl ; j++)
        b[i][j] = BF[j][i];
      
  for(j=0; j < nc-1 ; j+=2) {
    for(i=0; i <= j ; i+=2) {
      att=a[i]; bt=b[j];
      at1=a[i+1]; bt1=b[j+1];
      ij = IOFF(i,j);
      ij1 = IOFF(i,j+1);
      i1j = IOFF(i+1,j);
      i1j1 = IOFF(i+1,j+1);
      if(add) {
        t00 = CF[ij];
        if(i!=j) t01 = CF[ij1];
        t10 = CF[i1j];
        t11 = CF[i1j1];
        }
      else t00=t01=t10=t11=0.0;
      for(k=nl; k ; k--,att++,bt++,at1++,bt1++) {
        t00 += *att * *bt;
        if(i!=j) t01 += *att * *bt1;
        t10 += *at1 * *bt;
        t11 += *at1 * *bt1;
        }
      CF[ij]=t00;
      if(i!=j) CF[ij1]=t01;
      CF[i1j]=t10;
      CF[i1j1]=t11;
      }
    }
  if(odd_nc) {
    for(i=0; i < nr-1 ; i+=2) {
      att=a[i]; bt=b[j];
      at1=a[i+1];
      ij = IOFF(i,j);
      i1j = IOFF(i+1,j);
      if(add) {
        t00 = CF[ij];
        t10 = CF[i1j];
        }
      else t00=t10=0.0;
      for(k= nl; k ; k--,att++,bt++,at1++) {
        t00 += *att * *bt;
        t10 += *at1 * *bt;
        }
      CF[ij]=t00;
      CF[i1j]=t10;
      }
    att=a[i]; bt=b[j];
    ij = IOFF(i,j);
    if(add) t00 = CF[ij];
    else t00=0.0;
    for(k=nl; k ; k--,att++,bt++) t00 += *att * *bt;
    CF[ij]=t00;
    }
  }
