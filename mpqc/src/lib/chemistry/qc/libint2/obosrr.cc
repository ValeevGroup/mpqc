//
// obosrr.cc
//
// Copyright (C) 2001 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
// Maintainer: EV
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

#ifdef __GNUG__
#pragma implementation
#endif

#include <chemistry/qc/libint2/int1e.h>
#include <chemistry/qc/basis/fjt.h>

using namespace std;
using namespace sc;

inline void fail()
{
  ExEnv::errn() << scprintf("failing module:\n%s",__FILE__) << endl;
  abort();
}


/* Recurrence relation are from the Obara-Saika paper - pp. 3971-3972 */

void Int1eLibint2::AI_OSrecurs_(double ***AI0, double PA[3], double PB[3],
			      double PC[3], double gamma, int iang, int jang)
{
  int a,b,m;
  int izm = 1;
  int iym = iang + 1;
  int ixm = iym * iym;
  int jzm = 1;
  int jym = jang + 1;
  int jxm = jym * jym;
  int ix,iy,iz,jx,jy,jz;
  int iind,jind;
  double pp = 1/(2*gamma);
  int mmax = iang+jang;
  double tmp = sqrt(gamma)*M_2_SQRTPI;
  double u = gamma*(PC[0]*PC[0] + PC[1]*PC[1] + PC[2]*PC[2]);
  double *F = Fm_Eval_->values(mmax,u);

	/* Computing starting integrals for recursion */

  for(m=0;m<=mmax;m++)
    AI0[0][0][m] = tmp*F[m];

	/* Upward recursion in j with i=0 */
  
  for(b=1;b<=jang;b++)
    for(jx=0;jx<=b;jx++)
    for(jy=0;jy<=b-jx;jy++) {
      jz = b-jx-jy;
      jind = jx*jxm+jy*jym+jz*jzm;
      if (jz > 0) {
        for(m=0;m<=mmax-b;m++)	/* Electrostatic potential integrals */
          AI0[0][jind][m] = PB[2]*AI0[0][jind-jzm][m] - 
                            PC[2]*AI0[0][jind-jzm][m+1];
	if (jz > 1) {
          for(m=0;m<=mmax-b;m++)
            AI0[0][jind][m] += pp*(jz-1)*(AI0[0][jind-2*jzm][m] -
                                          AI0[0][jind-2*jzm][m+1]);
        }
      }
      else 
      if (jy > 0) {
        for(m=0;m<=mmax-b;m++)
          AI0[0][jind][m] = PB[1]*AI0[0][jind-jym][m] -
                            PC[1]*AI0[0][jind-jym][m+1];
        if (jy > 1) {
          for(m=0;m<=mmax-b;m++)  
            AI0[0][jind][m] += pp*(jy-1)*(AI0[0][jind-2*jym][m] -
                                          AI0[0][jind-2*jym][m+1]);
        }
      }
      else
      if (jx > 0) {
        for(m=0;m<=mmax-b;m++)
          AI0[0][jind][m] = PB[0]*AI0[0][jind-jxm][m] -
                            PC[0]*AI0[0][jind-jxm][m+1];
        if (jx > 1) {
          for(m=0;m<=mmax-b;m++)  
            AI0[0][jind][m] += pp*(jx-1)*(AI0[0][jind-2*jxm][m] -
                                          AI0[0][jind-2*jxm][m+1]);
        }
      }
      else
	fail();
    }
 


  /* The following fragment cannot be vectorized easily, I guess :-) */
	/* Upward recursion in i with all possible j's */

  for(b=0;b<=jang;b++)
    for(jx=0;jx<=b;jx++)
    for(jy=0;jy<=b-jx;jy++) {
    jz = b-jx-jy;
    jind = jx*jxm + jy*jym + jz*jzm;
    for(a=1;a<=iang;a++)
      for(ix=0;ix<=a;ix++)
      for(iy=0;iy<=a-ix;iy++) {
        iz = a-ix-iy;
        iind = ix*ixm + iy*iym + iz*izm;
        if (iz > 0) {
          for(m=0;m<=mmax-a-b;m++)
            AI0[iind][jind][m] = PA[2]*AI0[iind-izm][jind][m] - 
                                 PC[2]*AI0[iind-izm][jind][m+1];
          if (iz > 1) {
            for(m=0;m<=mmax-a-b;m++)
              AI0[iind][jind][m] += pp*(iz-1)*
               (AI0[iind-2*izm][jind][m] - AI0[iind-2*izm][jind][m+1]);
          }
          if (jz > 0) {
            for(m=0;m<=mmax-a-b;m++)
              AI0[iind][jind][m] += pp*jz*
               (AI0[iind-izm][jind-jzm][m] - AI0[iind-izm][jind-jzm][m+1]);
          }
        }
        else
	if (iy > 0) {
          for(m=0;m<=mmax-a-b;m++)
            AI0[iind][jind][m] = PA[1]*AI0[iind-iym][jind][m] -
                                 PC[1]*AI0[iind-iym][jind][m+1];
	  if (iy > 1) {
            for(m=0;m<=mmax-a-b;m++)
              AI0[iind][jind][m] += pp*(iy-1)*
              (AI0[iind-2*iym][jind][m] - AI0[iind-2*iym][jind][m+1]);
          }
	  if (jy > 0) {
            for(m=0;m<=mmax-a-b;m++)
              AI0[iind][jind][m] += pp*jy*
               (AI0[iind-iym][jind-jym][m] - AI0[iind-iym][jind-jym][m+1]);
          }
        }
        else
	if (ix > 0) {
          for(m=0;m<=mmax-a-b;m++)
            AI0[iind][jind][m] = PA[0]*AI0[iind-ixm][jind][m] -
                                 PC[0]*AI0[iind-ixm][jind][m+1];
          if (ix > 1) {
            for(m=0;m<=mmax-a-b;m++)
              AI0[iind][jind][m] += pp*(ix-1)*
               (AI0[iind-2*ixm][jind][m] - AI0[iind-2*ixm][jind][m+1]);
          }
          if (jx > 0) {
            for(m=0;m<=mmax-a-b;m++)
              AI0[iind][jind][m] += pp*jx*
               (AI0[iind-ixm][jind-jxm][m] - AI0[iind-ixm][jind-jxm][m+1]);  
          }
        }
        else
	  fail();
      }
    }

  return;
}


void Int1eLibint2::OI_OSrecurs_(double **OIX, double **OIY, double **OIZ, double PA[3], double PB[3],
			      double gamma, int lmaxi, int lmaxj)
{
  int i,j,k;
  double pp = 1/(2*gamma);

  OIX[0][0] = OIY[0][0] = OIZ[0][0] = 1.0;

	/* Upward recursion in j for i=0 */

  OIX[0][1] = PB[0];
  OIY[0][1] = PB[1];
  OIZ[0][1] = PB[2];

  for(j=1;j<lmaxj;j++) {
    OIX[0][j+1] = PB[0]*OIX[0][j];
    OIY[0][j+1] = PB[1]*OIY[0][j];
    OIZ[0][j+1] = PB[2]*OIZ[0][j];
    OIX[0][j+1] += j*pp*OIX[0][j-1];
    OIY[0][j+1] += j*pp*OIY[0][j-1];
    OIZ[0][j+1] += j*pp*OIZ[0][j-1];
  }

	/* Upward recursion in i for all j's */

  OIX[1][0] = PA[0];
  OIY[1][0] = PA[1];
  OIZ[1][0] = PA[2];
  for(j=1;j<=lmaxj;j++) {
    OIX[1][j] = PA[0]*OIX[0][j];
    OIY[1][j] = PA[1]*OIY[0][j];
    OIZ[1][j] = PA[2]*OIZ[0][j];
    OIX[1][j] += j*pp*OIX[0][j-1];
    OIY[1][j] += j*pp*OIY[0][j-1];
    OIZ[1][j] += j*pp*OIZ[0][j-1];
  }
  for(i=1;i<lmaxi;i++) {
    OIX[i+1][0] = PA[0]*OIX[i][0];
    OIY[i+1][0] = PA[1]*OIY[i][0];
    OIZ[i+1][0] = PA[2]*OIZ[i][0];
    OIX[i+1][0] += i*pp*OIX[i-1][0];
    OIY[i+1][0] += i*pp*OIY[i-1][0];
    OIZ[i+1][0] += i*pp*OIZ[i-1][0];
    for(j=1;j<=lmaxj;j++) {
      OIX[i+1][j] = PA[0]*OIX[i][j];
      OIY[i+1][j] = PA[1]*OIY[i][j];
      OIZ[i+1][j] = PA[2]*OIZ[i][j];
      OIX[i+1][j] += i*pp*OIX[i-1][j];
      OIY[i+1][j] += i*pp*OIY[i-1][j];
      OIZ[i+1][j] += i*pp*OIZ[i-1][j];
      OIX[i+1][j] += j*pp*OIX[i][j-1];
      OIY[i+1][j] += j*pp*OIY[i][j-1];
      OIZ[i+1][j] += j*pp*OIZ[i][j-1];
    }
  }

  return;
}
