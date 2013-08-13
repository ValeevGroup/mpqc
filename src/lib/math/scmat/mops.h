//
// mops.h --- block matrix operations
//
// Copyright (C) 1997 Limit Point Systems, Inc.
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

#ifndef _math_scmat_mops_h
#define _math_scmat_mops_h

#define D1 32

// copy a chunk of rectangular matrix source into dest.  dest is D1xD1, and is
// padded with zeros

static inline void
copy_block(double **dest, double **source,
           int istart, int ni, int jstart, int nj)
{
  int ii,jj;
  
  for (ii=0; ii < ni; ii++) {
    double *di = dest[ii];
    double *si = &source[istart+ii][jstart];
    for (jj=0; jj < nj; jj++)
      di[jj] = si[jj];
    for (; jj < D1; jj++)
      di[jj] = 0;
  }

  int left=D1-ii;
  if (left)
    memset(dest[ii], 0, sizeof(double)*left*D1);
}

static inline void
copy_trans_block(double **dest, double **source,
                 int istart, int ni, int jstart, int nj)
{
  int ii,jj;
  
  memset(dest[0], 0, sizeof(double)*D1*D1);

  for (jj=0; jj < nj; jj++) {
    double *sj = &source[jstart+jj][istart];
    for (ii=0; ii < ni; ii++)
      dest[ii][jj] = sj[ii];
  }
}

// copy a chunk of symmetric matrix source into dest.  dest is D1xD1, and is
// padded with zeros
static inline void
copy_sym_block(double **dest, double **source,
               int istart, int ni, int jstart, int nj)
{
  int ii,jj;

  for (ii=0; ii < ni; ii++) {
    double *di = dest[ii];
    double *si = &source[istart+ii][jstart];
    
    if (jstart < istart)
      for (jj=0; jj < nj; jj++)
        di[jj] = si[jj];
    else if (jstart==istart)
      for (jj=0; jj <= ii; jj++)
        di[jj] = dest[jj][ii] = si[jj];
    else
      for (jj=0; jj < nj; jj++)
        di[jj] = source[jstart+jj][istart+ii];

    for (jj=nj; jj < D1; jj++)
      di[jj] = 0;
  }

  int left=D1-ii;
  if (left)
    memset(dest[ii], 0, sizeof(double)*left*D1);
}

static inline void
return_block(double **dest, double **source,
             int istart, int ni, int jstart, int nj)
{
  int ii,jj;

  for (ii=0; ii < ni; ii++)
    for (jj=0; jj < nj; jj++)
      dest[istart+ii][jstart+jj] = source[ii][jj];
}

// a, b, and c are all D1xD1 blocks
static inline void
mult_block(double **a, double **b, double **c, int ni, int nj, int nk)
{
  int ii,jj,kk;
  double t00,t10,t20,t30;
  double *a0, *a1, *a2, *a3;
  double *c0, *c1, *c2, *c3;

  for (ii=0; ii < ni; ii += 4) {
    a0=a[ii]; a1=a[ii+1]; a2=a[ii+2]; a3=a[ii+3];
    c0=c[ii]; c1=c[ii+1]; c2=c[ii+2]; c3=c[ii+3];

    for (jj=0; jj < nj; jj++) {
      double *bt = b[jj];
      t00=c0[jj]; t10=c1[jj]; t20=c2[jj]; t30=c3[jj];

      for (kk=0; kk < nk; kk += 2) {
        double b0=bt[kk], b1=bt[kk+1];
        t00 += a0[kk]*b0 + a0[kk+1]*b1;
        t10 += a1[kk]*b0 + a1[kk+1]*b1;
        t20 += a2[kk]*b0 + a2[kk+1]*b1;
        t30 += a3[kk]*b0 + a3[kk+1]*b1;
      }

      c0[jj]=t00;
      c1[jj]=t10;
      c2[jj]=t20;
      c3[jj]=t30;
    }
  }
}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
