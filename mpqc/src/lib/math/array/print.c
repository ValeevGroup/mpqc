/*
 * print.c
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
 * Revision 1.5  1996/10/25 21:00:51  etseidl
 * add copyleft notice
 *
 * Revision 1.4  1996/03/23 02:38:09  cljanss
 * Everything can now be configured with autoconf.
 *
 * Revision 1.3  1995/03/17 01:50:10  cljanss
 * Removed -I. and -I$(SRCDIR) from the default include path in
 * GlobalMakefile to avoid name conflicts with system include files.
 * Modified files under src.lib to include all files relative to src.lib.
 * Makefiles under src.bin need to add the -I. and -I$(SRCDIR) back onto
 * INCLUDE and CXXINCLUDE or make other arrangements.
 *
 * Revision 1.2  1994/08/26  17:57:48  etseidl
 * get rid of rcs id's and fix a few warnings
 *
 * Revision 1.1.1.1  1993/12/29  12:53:30  etseidl
 * SC source tree 0.1
 *
 * Revision 1.2  1992/06/17  22:10:18  jannsen
 * clean up for saber-c
 *
 * Revision 1.1.1.1  1992/03/17  16:35:37  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:35:36  seidl
 * Initial revision
 *
 * Revision 1.2  1991/12/20  16:04:51  seidl
 * fix bug in print_dv
 *
 * Revision 1.1  1991/12/20  16:04:13  seidl
 * Initial revision
 *
 * Revision 1.1  91/12/02  23:48:44  etseidl
 * Initial revision
 *  */

/* routines to pretty-print double matrices and double_vectors
 *
 * math_print_dv assumes that you are printing a the lower triangle
 * of a square matrix
 *
 * math_print_dmvv prints a double_matrix, and two double_vectors.
 * it is intended for printing out and scf vector, the eigenvalues, and
 * the occupation numbers
 *
 * math_print_dmv is intended for printing out the eigenvectors and
 * eigenvalues of a matrix
 */

#include <stdio.h>
#include <tmpl.h>
#include <math/array/matrix.h>

#include <math/array/print.gbl>
#include <math/array/print.lcl>

GLOBAL_FUNCTION void
math_print_dm(out,aa)
FILE *out;
double_matrix_t *aa;
{
  int ii,jj,kk,nn;
  int i,j;
  int m=aa->n1;
  int n=aa->n2;
  double **a = aa->d;

  a=aa->d;

  ii=0;jj=0;
  for(ii=jj=0;;) {
    ii++;
    jj++;
    kk=6*jj;
    nn=n;
    if (nn > kk) nn=kk;
    fprintf (out,"\n");
    for (i=ii; i <= nn; i++) fprintf(out,"%12d",i);
    fprintf (out,"\n");
    for (i=0; i < m; i++) {
      fprintf (out,"\n%5d",i+1);
      for (j=ii-1; j < nn; j++) {
        fprintf (out,"%12.7f",a[i][j]);
        }
      }
    fprintf (out,"\n");
    if (n <= kk) {
      fflush(out);
      return;
      }
    ii=kk;
    }
  }

GLOBAL_FUNCTION void
math_print_dv(out,aa)
FILE *out;
double_vector_t *aa;
{
  int ii,jj,kk,mm,nn;
  int i,j,i1,i2;
  double sqrt(double);
  double *a=aa->d;
  int m = (((int) sqrt(1.0+8.0*(double)aa->n))-1)/2;

  for(ii=jj=0;;) {
    ii++;
    jj++;
    kk=6*jj;
    nn = kk + kk*(kk-1)/2;
    mm=m;
    if (m > kk) mm=kk;
    fprintf (out,"\n");
    for (i=ii; i <= mm; i++) fprintf(out,"%12d",i);
    fprintf (out,"\n");
    for (i=ii; i <= m; i++) {
      i1=i*(i-1)/2+ii;
      i2=i+i*(i-1)/2;
      if (i2 > nn) i2 = i1+5;
      fprintf (out,"\n%5d",i);
      for (j=i1; j <= i2; j++) {
        fprintf (out,"%12.7f",a[j-1]);
        }
      }
    if (m <= kk) {
      fprintf(out,"\n");
      fflush(out);
      return;
      }
    ii=kk;
    }
  }

GLOBAL_FUNCTION void
math_print_dmvv(out,aa,bb,cc)
FILE *out;
double_matrix_t *aa;
double_vector_t *bb;
double_vector_t *cc;
{
  int ii,jj,kk,nn;
  int i,j;
  int m=aa->n1;
  int n=aa->n2;
  double **a = aa->d;
  double *b = bb->d;
  double *c = cc->d;

  for(ii=jj=0; ;) {
    ii++;
    jj++;
    kk=6*jj;
    nn=n;
    if (nn > kk) nn=kk;
    fprintf (out,"\n");
    for (i=ii; i <= nn; i++) fprintf(out,"%12d",i);
    fprintf (out,"\n");
    for (i=0; i < m; i++) {
      fprintf (out,"\n%5d",i+1);
      for (j=ii-1; j < nn; j++) {
        fprintf (out,"%12.7f",a[i][j]);
        }
      }
    fprintf (out,"\n");
    fprintf (out,"\n     ");
    for (j=ii-1; j < nn; j++) {
      fprintf(out,"%12.7f",b[j]);
      }
    fprintf (out,"\n");
    fprintf (out,"\n     ");
    for (j=ii-1; j < nn; j++) {
      fprintf(out,"%12.7f",c[j]);
      }
    fprintf (out,"\n");
    if (n <= kk) {
      fflush(out);
      return;
      }
    ii=kk;
    }
  }

GLOBAL_FUNCTION void
math_print_dmv(out,aa,bb)
FILE *out;
double_matrix_t *aa;
double_vector_t *bb;
{
  int ii,jj,kk,nn;
  int i,j;
  int m=aa->n1;
  int n=aa->n2;
  double **a = aa->d;
  double *b = bb->d;

  for(ii=jj=0;;) {
    ii++;
    jj++;
    kk=6*jj;
    nn=n;
    if (nn > kk) nn=kk;
    fprintf (out,"\n");
    for (i=ii; i <= nn; i++) fprintf(out,"%12d",i);
    fprintf (out,"\n");
    for (i=0; i < m; i++) {
      fprintf (out,"\n%5d",i+1);
      for (j=ii-1; j < nn; j++) {
        fprintf (out,"%12.7f",a[i][j]);
        }
      }
    fprintf (out,"\n");
    fprintf (out,"\n     ");
    for (j=ii-1; j < nn; j++) {
      fprintf(out,"%12.7f",b[j]);
      }
    fprintf (out,"\n");
    if (n <= kk) {
      fflush(out);
      return;
      }
    ii=kk;
    }
  }
