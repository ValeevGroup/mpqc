//
// obosrr.timpl.h
//
// Copyright (C) 2001 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
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

/* Recurrence relation are from the Obara-Saika paper - pp. 3971-3972 */

template <int Order>
void
Int1eLibint2::AI_OSrecurs_(double PA[3], double PB[3],
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
  int mmax = iang+jang + Order;
  double tmp = sqrt(gamma)*M_2_SQRTPI;
  double u = gamma*(PC[0]*PC[0] + PC[1]*PC[1] + PC[2]*PC[2]);
  Fm_Eval_->eval(Fm_table_,u,mmax);

	/* Computing starting integrals for recursion */

  for(m=0;m<=mmax;m++)
    AI0_[0][0][m] = tmp * Fm_table_[m];
  if (Order >= 1)
    for(m=0;m<=mmax-1;m++) {
      AIX_[0][0][m] = 2*gamma*PC[0]*AI0_[0][0][m+1];
      AIY_[0][0][m] = 2*gamma*PC[1]*AI0_[0][0][m+1];
      AIZ_[0][0][m] = 2*gamma*PC[2]*AI0_[0][0][m+1];
    }
  if (Order >= 2)
    for (m = 0; m <= mmax - 2; m++) {
      AIXX_[0][0][m] = 4 * gamma * gamma * PC[0] * PC[0] * AI0_[0][0][m + 2]
          - 2 * gamma * AI0_[0][0][m + 1];
      AIYY_[0][0][m] = 4 * gamma * gamma * PC[1] * PC[1] * AI0_[0][0][m + 2]
          - 2 * gamma * AI0_[0][0][m + 1];
      AIZZ_[0][0][m] = 4 * gamma * gamma * PC[2] * PC[2] * AI0_[0][0][m + 2]
          - 2 * gamma * AI0_[0][0][m + 1];
      AIXY_[0][0][m] = 4 * gamma * gamma * PC[0] * PC[1] * AI0_[0][0][m + 2];
      AIXZ_[0][0][m] = 4 * gamma * gamma * PC[0] * PC[2] * AI0_[0][0][m + 2];
      AIYZ_[0][0][m] = 4 * gamma * gamma * PC[1] * PC[2] * AI0_[0][0][m + 2];
    }

	/* Upward recursion in j with i=0 */
  
  for(b=1;b<=jang;b++)
    for(jx=0;jx<=b;jx++)
    for(jy=0;jy<=b-jx;jy++) {
      jz = b-jx-jy;
      jind = jx*jxm+jy*jym+jz*jzm;
      if (jz > 0) {
        for(m=0;m<=mmax-b;m++)	/* Electrostatic potential integrals */
          AI0_[0][jind][m] = PB[2]*AI0_[0][jind-jzm][m] -
                            PC[2]*AI0_[0][jind-jzm][m+1];
        if (Order >= 1)
        for(m=0;m<=mmax-b-1;m++) {  /* Electric field integrals */
          AIX_[0][jind][m] = PB[2]*AIX_[0][jind-jzm][m] -
                            PC[2]*AIX_[0][jind-jzm][m+1];
          AIY_[0][jind][m] = PB[2]*AIY_[0][jind-jzm][m] -
                            PC[2]*AIY_[0][jind-jzm][m+1];
          AIZ_[0][jind][m] = PB[2]*AIZ_[0][jind-jzm][m] -
                            PC[2]*AIZ_[0][jind-jzm][m+1] +
                                AI0_[0][jind-jzm][m+1];
        }
        if (Order >= 2)
        for(m=0;m<=mmax-b-2;m++) {  /* Gradients of the electric field */
          AIXX_[0][jind][m] = PB[2]*AIXX_[0][jind-jzm][m] -
                             PC[2]*AIXX_[0][jind-jzm][m+1];
          AIYY_[0][jind][m] = PB[2]*AIYY_[0][jind-jzm][m] -
                             PC[2]*AIYY_[0][jind-jzm][m+1];
          AIZZ_[0][jind][m] = PB[2]*AIZZ_[0][jind-jzm][m] -
                             PC[2]*AIZZ_[0][jind-jzm][m+1] +
                               2*AIZ_[0][jind-jzm][m+1];
          AIXY_[0][jind][m] = PB[2]*AIXY_[0][jind-jzm][m] -
                             PC[2]*AIXY_[0][jind-jzm][m+1];
          AIXZ_[0][jind][m] = PB[2]*AIXZ_[0][jind-jzm][m] -
                             PC[2]*AIXZ_[0][jind-jzm][m+1] +
                                 AIX_[0][jind-jzm][m+1];
          AIYZ_[0][jind][m] = PB[2]*AIYZ_[0][jind-jzm][m] -
                             PC[2]*AIYZ_[0][jind-jzm][m+1] +
                                 AIY_[0][jind-jzm][m+1];
        }
	if (jz > 1) {
          for(m=0;m<=mmax-b;m++)
            AI0_[0][jind][m] += pp*(jz-1)*(AI0_[0][jind-2*jzm][m] -
                                          AI0_[0][jind-2*jzm][m+1]);
          if (Order >= 1)
          for(m=0;m<=mmax-b-1;m++) {
            AIX_[0][jind][m] += pp*(jz-1)*(AIX_[0][jind-2*jzm][m] -
                                          AIX_[0][jind-2*jzm][m+1]);
            AIY_[0][jind][m] += pp*(jz-1)*(AIY_[0][jind-2*jzm][m] -
                                          AIY_[0][jind-2*jzm][m+1]);
            AIZ_[0][jind][m] += pp*(jz-1)*(AIZ_[0][jind-2*jzm][m] -
                                          AIZ_[0][jind-2*jzm][m+1]);
          }
          if (Order >= 2)
          for(m=0;m<=mmax-b-2;m++) {
            AIXX_[0][jind][m] += pp*(jz-1)*(AIXX_[0][jind-2*jzm][m] -
                                           AIXX_[0][jind-2*jzm][m+1]);
            AIYY_[0][jind][m] += pp*(jz-1)*(AIYY_[0][jind-2*jzm][m] -
                                           AIYY_[0][jind-2*jzm][m+1]);
            AIZZ_[0][jind][m] += pp*(jz-1)*(AIZZ_[0][jind-2*jzm][m] -
                                           AIZZ_[0][jind-2*jzm][m+1]);
            AIXY_[0][jind][m] += pp*(jz-1)*(AIXY_[0][jind-2*jzm][m] -
                                           AIXY_[0][jind-2*jzm][m+1]);
            AIXZ_[0][jind][m] += pp*(jz-1)*(AIXZ_[0][jind-2*jzm][m] -
                                           AIXZ_[0][jind-2*jzm][m+1]);
            AIYZ_[0][jind][m] += pp*(jz-1)*(AIYZ_[0][jind-2*jzm][m] -
                                           AIYZ_[0][jind-2*jzm][m+1]);
          }
        }
      }
      else 
      if (jy > 0) {
        for(m=0;m<=mmax-b;m++)
          AI0_[0][jind][m] = PB[1]*AI0_[0][jind-jym][m] -
                            PC[1]*AI0_[0][jind-jym][m+1];
        if (Order >= 1)
        for(m=0;m<=mmax-b-1;m++) {
          AIX_[0][jind][m] = PB[1]*AIX_[0][jind-jym][m] -
                            PC[1]*AIX_[0][jind-jym][m+1];
          AIY_[0][jind][m] = PB[1]*AIY_[0][jind-jym][m] -
                            PC[1]*AIY_[0][jind-jym][m+1] +
                                AI0_[0][jind-jym][m+1];
          AIZ_[0][jind][m] = PB[1]*AIZ_[0][jind-jym][m] -
                            PC[1]*AIZ_[0][jind-jym][m+1];
        }
        if (Order >= 2)
        for(m=0;m<=mmax-b-2;m++) {
          AIXX_[0][jind][m] = PB[1]*AIXX_[0][jind-jym][m] -
                             PC[1]*AIXX_[0][jind-jym][m+1];
          AIYY_[0][jind][m] = PB[1]*AIYY_[0][jind-jym][m] -
                             PC[1]*AIYY_[0][jind-jym][m+1] +
                               2*AIY_[0][jind-jym][m+1];
          AIZZ_[0][jind][m] = PB[1]*AIZZ_[0][jind-jym][m] -
                             PC[1]*AIZZ_[0][jind-jym][m+1];
          AIXY_[0][jind][m] = PB[1]*AIXY_[0][jind-jym][m] -
                             PC[1]*AIXY_[0][jind-jym][m+1] +
                                 AIX_[0][jind-jym][m+1];
          AIXZ_[0][jind][m] = PB[1]*AIXZ_[0][jind-jym][m] -
                             PC[1]*AIXZ_[0][jind-jym][m+1];
          AIYZ_[0][jind][m] = PB[1]*AIYZ_[0][jind-jym][m] -
                             PC[1]*AIYZ_[0][jind-jym][m+1] +
                                 AIZ_[0][jind-jym][m+1];
        }
        if (jy > 1) {
          for(m=0;m<=mmax-b;m++)  
            AI0_[0][jind][m] += pp*(jy-1)*(AI0_[0][jind-2*jym][m] -
                                          AI0_[0][jind-2*jym][m+1]);
          if (Order >= 1)
          for(m=0;m<=mmax-b-1;m++) {
            AIX_[0][jind][m] += pp*(jy-1)*(AIX_[0][jind-2*jym][m] -
                                          AIX_[0][jind-2*jym][m+1]);
            AIY_[0][jind][m] += pp*(jy-1)*(AIY_[0][jind-2*jym][m] -
                                          AIY_[0][jind-2*jym][m+1]);
            AIZ_[0][jind][m] += pp*(jy-1)*(AIZ_[0][jind-2*jym][m] -
                                          AIZ_[0][jind-2*jym][m+1]);
          }
          if (Order >= 2)
          for(m=0;m<=mmax-b-2;m++) {
            AIXX_[0][jind][m] += pp*(jy-1)*(AIXX_[0][jind-2*jym][m] -
                                           AIXX_[0][jind-2*jym][m+1]);
            AIYY_[0][jind][m] += pp*(jy-1)*(AIYY_[0][jind-2*jym][m] -
                                           AIYY_[0][jind-2*jym][m+1]);
            AIZZ_[0][jind][m] += pp*(jy-1)*(AIZZ_[0][jind-2*jym][m] -
                                           AIZZ_[0][jind-2*jym][m+1]);
            AIXY_[0][jind][m] += pp*(jy-1)*(AIXY_[0][jind-2*jym][m] -
                                           AIXY_[0][jind-2*jym][m+1]);
            AIXZ_[0][jind][m] += pp*(jy-1)*(AIXZ_[0][jind-2*jym][m] -
                                           AIXZ_[0][jind-2*jym][m+1]);
            AIYZ_[0][jind][m] += pp*(jy-1)*(AIYZ_[0][jind-2*jym][m] -
                                           AIYZ_[0][jind-2*jym][m+1]);
          }
        }
      }
      else
      if (jx > 0) {
        for(m=0;m<=mmax-b;m++)
          AI0_[0][jind][m] = PB[0]*AI0_[0][jind-jxm][m] -
                            PC[0]*AI0_[0][jind-jxm][m+1];
        if (Order >= 1)
        for(m=0;m<=mmax-b-1;m++) {
          AIX_[0][jind][m] = PB[0]*AIX_[0][jind-jxm][m] -
                            PC[0]*AIX_[0][jind-jxm][m+1] +
                                AI0_[0][jind-jxm][m+1];
          AIY_[0][jind][m] = PB[0]*AIY_[0][jind-jxm][m] -
                            PC[0]*AIY_[0][jind-jxm][m+1];
          AIZ_[0][jind][m] = PB[0]*AIZ_[0][jind-jxm][m] -
                            PC[0]*AIZ_[0][jind-jxm][m+1];
        }
        if (Order >= 2)
        for(m=0;m<=mmax-b-2;m++) {
          AIXX_[0][jind][m] = PB[0]*AIXX_[0][jind-jxm][m] -
                             PC[0]*AIXX_[0][jind-jxm][m+1] +
                               2*AIX_[0][jind-jxm][m+1];
          AIYY_[0][jind][m] = PB[0]*AIYY_[0][jind-jxm][m] -
                             PC[0]*AIYY_[0][jind-jxm][m+1];
          AIZZ_[0][jind][m] = PB[0]*AIZZ_[0][jind-jxm][m] -
                             PC[0]*AIZZ_[0][jind-jxm][m+1];
          AIXY_[0][jind][m] = PB[0]*AIXY_[0][jind-jxm][m] -
                             PC[0]*AIXY_[0][jind-jxm][m+1] +
                                 AIY_[0][jind-jxm][m+1];
          AIXZ_[0][jind][m] = PB[0]*AIXZ_[0][jind-jxm][m] -
                             PC[0]*AIXZ_[0][jind-jxm][m+1] +
                                 AIZ_[0][jind-jxm][m+1];
          AIYZ_[0][jind][m] = PB[0]*AIYZ_[0][jind-jxm][m] -
                             PC[0]*AIYZ_[0][jind-jxm][m+1];
        }
        if (jx > 1) {
          for(m=0;m<=mmax-b;m++)  
            AI0_[0][jind][m] += pp*(jx-1)*(AI0_[0][jind-2*jxm][m] -
                                          AI0_[0][jind-2*jxm][m+1]);
          if (Order >= 1)
          for(m=0;m<=mmax-b-1;m++) {
            AIX_[0][jind][m] += pp*(jx-1)*(AIX_[0][jind-2*jxm][m] -
                                          AIX_[0][jind-2*jxm][m+1]);
            AIY_[0][jind][m] += pp*(jx-1)*(AIY_[0][jind-2*jxm][m] -
                                          AIY_[0][jind-2*jxm][m+1]);
            AIZ_[0][jind][m] += pp*(jx-1)*(AIZ_[0][jind-2*jxm][m] -
                                          AIZ_[0][jind-2*jxm][m+1]);
          }
          if (Order >= 2)
          for(m=0;m<=mmax-b-2;m++) {
            AIXX_[0][jind][m] += pp*(jx-1)*(AIXX_[0][jind-2*jxm][m] -
                                           AIXX_[0][jind-2*jxm][m+1]);
            AIYY_[0][jind][m] += pp*(jx-1)*(AIYY_[0][jind-2*jxm][m] -
                                           AIYY_[0][jind-2*jxm][m+1]);
            AIZZ_[0][jind][m] += pp*(jx-1)*(AIZZ_[0][jind-2*jxm][m] -
                                           AIZZ_[0][jind-2*jxm][m+1]);
            AIXY_[0][jind][m] += pp*(jx-1)*(AIXY_[0][jind-2*jxm][m] -
                                           AIXY_[0][jind-2*jxm][m+1]);
            AIXZ_[0][jind][m] += pp*(jx-1)*(AIXZ_[0][jind-2*jxm][m] -
                                           AIXZ_[0][jind-2*jxm][m+1]);
            AIYZ_[0][jind][m] += pp*(jx-1)*(AIYZ_[0][jind-2*jxm][m] -
                                           AIYZ_[0][jind-2*jxm][m+1]);
          }
        }
      }
      else
        throw ProgrammingError("Logic error in Coulomb 1-e OS RR code",__FILE__,__LINE__);
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
            AI0_[iind][jind][m] = PA[2]*AI0_[iind-izm][jind][m] -
                                 PC[2]*AI0_[iind-izm][jind][m+1];
          if (Order >= 1)
          for(m=0;m<=mmax-a-b-1;m++) {  /* Electric field integrals */
            AIX_[iind][jind][m] = PA[2]*AIX_[iind-izm][jind][m] -
                                 PC[2]*AIX_[iind-izm][jind][m+1];
            AIY_[iind][jind][m] = PA[2]*AIY_[iind-izm][jind][m] -
                                 PC[2]*AIY_[iind-izm][jind][m+1];
            AIZ_[iind][jind][m] = PA[2]*AIZ_[iind-izm][jind][m] -
                                 PC[2]*AIZ_[iind-izm][jind][m+1] +
                                     AI0_[iind-izm][jind][m+1];
          }
          if (Order >= 2)
          for(m=0;m<=mmax-a-b-2;m++) {  /* Gradients of the electric field */
            AIXX_[iind][jind][m] = PA[2]*AIXX_[iind-izm][jind][m] -
                                  PC[2]*AIXX_[iind-izm][jind][m+1];
            AIYY_[iind][jind][m] = PA[2]*AIYY_[iind-izm][jind][m] -
                                  PC[2]*AIYY_[iind-izm][jind][m+1];
            AIZZ_[iind][jind][m] = PA[2]*AIZZ_[iind-izm][jind][m] -
                                  PC[2]*AIZZ_[iind-izm][jind][m+1] +
                                    2*AIZ_[iind-izm][jind][m+1];
            AIXY_[iind][jind][m] = PA[2]*AIXY_[iind-izm][jind][m] -
                                  PC[2]*AIXY_[iind-izm][jind][m+1];
            AIXZ_[iind][jind][m] = PA[2]*AIXZ_[iind-izm][jind][m] -
                                  PC[2]*AIXZ_[iind-izm][jind][m+1] +
                                      AIX_[iind-izm][jind][m+1];
            AIYZ_[iind][jind][m] = PA[2]*AIYZ_[iind-izm][jind][m] -
                                  PC[2]*AIYZ_[iind-izm][jind][m+1] +
                                      AIY_[iind-izm][jind][m+1];
          }
          if (iz > 1) {
            for(m=0;m<=mmax-a-b;m++)
              AI0_[iind][jind][m] += pp*(iz-1)*
               (AI0_[iind-2*izm][jind][m] - AI0_[iind-2*izm][jind][m+1]);
            if (Order >= 1)
            for(m=0;m<=mmax-a-b-1;m++) {
              AIX_[iind][jind][m] += pp*(iz-1)*(AIX_[iind-2*izm][jind][m] -
                                               AIX_[iind-2*izm][jind][m+1]);
              AIY_[iind][jind][m] += pp*(iz-1)*(AIY_[iind-2*izm][jind][m] -
                                               AIY_[iind-2*izm][jind][m+1]);
              AIZ_[iind][jind][m] += pp*(iz-1)*(AIZ_[iind-2*izm][jind][m] -
                                               AIZ_[iind-2*izm][jind][m+1]);
            }
            if (Order >= 2)
            for(m=0;m<=mmax-a-b-2;m++) {
              AIXX_[iind][jind][m] += pp*(iz-1)*(AIXX_[iind-2*izm][jind][m] -
                                                AIXX_[iind-2*izm][jind][m+1]);
              AIYY_[iind][jind][m] += pp*(iz-1)*(AIYY_[iind-2*izm][jind][m] -
                                                AIYY_[iind-2*izm][jind][m+1]);
              AIZZ_[iind][jind][m] += pp*(iz-1)*(AIZZ_[iind-2*izm][jind][m] -
                                                AIZZ_[iind-2*izm][jind][m+1]);
              AIXY_[iind][jind][m] += pp*(iz-1)*(AIXY_[iind-2*izm][jind][m] -
                                                AIXY_[iind-2*izm][jind][m+1]);
              AIXZ_[iind][jind][m] += pp*(iz-1)*(AIXZ_[iind-2*izm][jind][m] -
                                                AIXZ_[iind-2*izm][jind][m+1]);
              AIYZ_[iind][jind][m] += pp*(iz-1)*(AIYZ_[iind-2*izm][jind][m] -
                                                AIYZ_[iind-2*izm][jind][m+1]);
            }
          }
          if (jz > 0) {
            for(m=0;m<=mmax-a-b;m++)
              AI0_[iind][jind][m] += pp*jz*
               (AI0_[iind-izm][jind-jzm][m] - AI0_[iind-izm][jind-jzm][m+1]);
            if (Order >= 1)
            for(m=0;m<=mmax-a-b-1;m++) {
              AIX_[iind][jind][m] += pp*jz*(AIX_[iind-izm][jind-jzm][m] -
                                           AIX_[iind-izm][jind-jzm][m+1]);
              AIY_[iind][jind][m] += pp*jz*(AIY_[iind-izm][jind-jzm][m] -
                                           AIY_[iind-izm][jind-jzm][m+1]);
              AIZ_[iind][jind][m] += pp*jz*(AIZ_[iind-izm][jind-jzm][m] -
                                           AIZ_[iind-izm][jind-jzm][m+1]);
            }
            if (Order >= 2)
            for(m=0;m<=mmax-a-b-2;m++) {
              AIXX_[iind][jind][m] += pp*jz*(AIXX_[iind-izm][jind-jzm][m] -
                                            AIXX_[iind-izm][jind-jzm][m+1]);
              AIYY_[iind][jind][m] += pp*jz*(AIYY_[iind-izm][jind-jzm][m] -
                                            AIYY_[iind-izm][jind-jzm][m+1]);
              AIZZ_[iind][jind][m] += pp*jz*(AIZZ_[iind-izm][jind-jzm][m] -
                                            AIZZ_[iind-izm][jind-jzm][m+1]);
              AIXY_[iind][jind][m] += pp*jz*(AIXY_[iind-izm][jind-jzm][m] -
                                            AIXY_[iind-izm][jind-jzm][m+1]);
              AIXZ_[iind][jind][m] += pp*jz*(AIXZ_[iind-izm][jind-jzm][m] -
                                            AIXZ_[iind-izm][jind-jzm][m+1]);
              AIYZ_[iind][jind][m] += pp*jz*(AIYZ_[iind-izm][jind-jzm][m] -
                                            AIYZ_[iind-izm][jind-jzm][m+1]);
            }
          }
        }
        else
	if (iy > 0) {
          for(m=0;m<=mmax-a-b;m++)
            AI0_[iind][jind][m] = PA[1]*AI0_[iind-iym][jind][m] -
                                 PC[1]*AI0_[iind-iym][jind][m+1];
          if (Order >= 1)
          for(m=0;m<=mmax-a-b-1;m++) {
            AIX_[iind][jind][m] = PA[1]*AIX_[iind-iym][jind][m] -
                                 PC[1]*AIX_[iind-iym][jind][m+1];
            AIY_[iind][jind][m] = PA[1]*AIY_[iind-iym][jind][m] -
                                 PC[1]*AIY_[iind-iym][jind][m+1] +
                                     AI0_[iind-iym][jind][m+1];
            AIZ_[iind][jind][m] = PA[1]*AIZ_[iind-iym][jind][m] -
                                 PC[1]*AIZ_[iind-iym][jind][m+1];
          }
          if (Order >= 2)
          for(m=0;m<=mmax-a-b-2;m++) {
            AIXX_[iind][jind][m] = PA[1]*AIXX_[iind-iym][jind][m] -
                                  PC[1]*AIXX_[iind-iym][jind][m+1];
            AIYY_[iind][jind][m] = PA[1]*AIYY_[iind-iym][jind][m] -
                                  PC[1]*AIYY_[iind-iym][jind][m+1] +
                                    2*AIY_[iind-iym][jind][m+1];
            AIZZ_[iind][jind][m] = PA[1]*AIZZ_[iind-iym][jind][m] -
                                  PC[1]*AIZZ_[iind-iym][jind][m+1];
            AIXY_[iind][jind][m] = PA[1]*AIXY_[iind-iym][jind][m] -
                                  PC[1]*AIXY_[iind-iym][jind][m+1] +
                                      AIX_[iind-iym][jind][m+1];
            AIXZ_[iind][jind][m] = PA[1]*AIXZ_[iind-iym][jind][m] -
                                  PC[1]*AIXZ_[iind-iym][jind][m+1];
            AIYZ_[iind][jind][m] = PA[1]*AIYZ_[iind-iym][jind][m] -
                                  PC[1]*AIYZ_[iind-iym][jind][m+1] +
                                      AIZ_[iind-iym][jind][m+1];
          }
	  if (iy > 1) {
            for(m=0;m<=mmax-a-b;m++)
              AI0_[iind][jind][m] += pp*(iy-1)*
              (AI0_[iind-2*iym][jind][m] - AI0_[iind-2*iym][jind][m+1]);
            if (Order >= 1)
            for(m=0;m<=mmax-a-b-1;m++) {
              AIX_[iind][jind][m] += pp*(iy-1)*(AIX_[iind-2*iym][jind][m] -
                                               AIX_[iind-2*iym][jind][m+1]);
              AIY_[iind][jind][m] += pp*(iy-1)*(AIY_[iind-2*iym][jind][m] -
                                               AIY_[iind-2*iym][jind][m+1]);
              AIZ_[iind][jind][m] += pp*(iy-1)*(AIZ_[iind-2*iym][jind][m] -
                                               AIZ_[iind-2*iym][jind][m+1]);
            }
            if (Order >= 2)
            for(m=0;m<=mmax-a-b-2;m++) {
              AIXX_[iind][jind][m] += pp*(iy-1)*(AIXX_[iind-2*iym][jind][m] -
                                                AIXX_[iind-2*iym][jind][m+1]);
              AIYY_[iind][jind][m] += pp*(iy-1)*(AIYY_[iind-2*iym][jind][m] -
                                                AIYY_[iind-2*iym][jind][m+1]);
              AIZZ_[iind][jind][m] += pp*(iy-1)*(AIZZ_[iind-2*iym][jind][m] -
                                                AIZZ_[iind-2*iym][jind][m+1]);
              AIXY_[iind][jind][m] += pp*(iy-1)*(AIXY_[iind-2*iym][jind][m] -
                                                AIXY_[iind-2*iym][jind][m+1]);
              AIXZ_[iind][jind][m] += pp*(iy-1)*(AIXZ_[iind-2*iym][jind][m] -
                                                AIXZ_[iind-2*iym][jind][m+1]);
              AIYZ_[iind][jind][m] += pp*(iy-1)*(AIYZ_[iind-2*iym][jind][m] -
                                                AIYZ_[iind-2*iym][jind][m+1]);
            }
          }
	  if (jy > 0) {
            for(m=0;m<=mmax-a-b;m++)
              AI0_[iind][jind][m] += pp*jy*
               (AI0_[iind-iym][jind-jym][m] - AI0_[iind-iym][jind-jym][m+1]);
            if (Order >= 1)
            for(m=0;m<=mmax-a-b-1;m++) {
              AIX_[iind][jind][m] += pp*jy*(AIX_[iind-iym][jind-jym][m] -
                                           AIX_[iind-iym][jind-jym][m+1]);
              AIY_[iind][jind][m] += pp*jy*(AIY_[iind-iym][jind-jym][m] -
                                           AIY_[iind-iym][jind-jym][m+1]);
              AIZ_[iind][jind][m] += pp*jy*(AIZ_[iind-iym][jind-jym][m] -
                                           AIZ_[iind-iym][jind-jym][m+1]);
            }
            if (Order >= 2)
            for(m=0;m<=mmax-a-b-2;m++) {
              AIXX_[iind][jind][m] += pp*jy*(AIXX_[iind-iym][jind-jym][m] -
                                            AIXX_[iind-iym][jind-jym][m+1]);
              AIYY_[iind][jind][m] += pp*jy*(AIYY_[iind-iym][jind-jym][m] -
                                            AIYY_[iind-iym][jind-jym][m+1]);
              AIZZ_[iind][jind][m] += pp*jy*(AIZZ_[iind-iym][jind-jym][m] -
                                            AIZZ_[iind-iym][jind-jym][m+1]);
              AIXY_[iind][jind][m] += pp*jy*(AIXY_[iind-iym][jind-jym][m] -
                                            AIXY_[iind-iym][jind-jym][m+1]);
              AIXZ_[iind][jind][m] += pp*jy*(AIXZ_[iind-iym][jind-jym][m] -
                                            AIXZ_[iind-iym][jind-jym][m+1]);
              AIYZ_[iind][jind][m] += pp*jy*(AIYZ_[iind-iym][jind-jym][m] -
                                            AIYZ_[iind-iym][jind-jym][m+1]);
            }
          }
        }
        else
	if (ix > 0) {
          for(m=0;m<=mmax-a-b;m++)
            AI0_[iind][jind][m] = PA[0]*AI0_[iind-ixm][jind][m] -
                                 PC[0]*AI0_[iind-ixm][jind][m+1];
          if (Order >= 1)
          for(m=0;m<=mmax-a-b-1;m++) {  /* Electric field integrals */
            AIX_[iind][jind][m] = PA[0]*AIX_[iind-ixm][jind][m] -
                                 PC[0]*AIX_[iind-ixm][jind][m+1] +
                                     AI0_[iind-ixm][jind][m+1];
            AIY_[iind][jind][m] = PA[0]*AIY_[iind-ixm][jind][m] -
                                 PC[0]*AIY_[iind-ixm][jind][m+1];
            AIZ_[iind][jind][m] = PA[0]*AIZ_[iind-ixm][jind][m] -
                                 PC[0]*AIZ_[iind-ixm][jind][m+1];
          }
          if (Order >= 2)
          for(m=0;m<=mmax-a-b-2;m++) {  /* Gradients of the electric field */
            AIXX_[iind][jind][m] = PA[0]*AIXX_[iind-ixm][jind][m] -
                                  PC[0]*AIXX_[iind-ixm][jind][m+1] +
                                    2*AIX_[iind-ixm][jind][m+1];
            AIYY_[iind][jind][m] = PA[0]*AIYY_[iind-ixm][jind][m] -
                                  PC[0]*AIYY_[iind-ixm][jind][m+1];
            AIZZ_[iind][jind][m] = PA[0]*AIZZ_[iind-ixm][jind][m] -
                                  PC[0]*AIZZ_[iind-ixm][jind][m+1];
            AIXY_[iind][jind][m] = PA[0]*AIXY_[iind-ixm][jind][m] -
                                  PC[0]*AIXY_[iind-ixm][jind][m+1] +
                                      AIY_[iind-ixm][jind][m+1];
            AIXZ_[iind][jind][m] = PA[0]*AIXZ_[iind-ixm][jind][m] -
                                  PC[0]*AIXZ_[iind-ixm][jind][m+1] +
                                      AIZ_[iind-ixm][jind][m+1];
            AIYZ_[iind][jind][m] = PA[0]*AIYZ_[iind-ixm][jind][m] -
                                  PC[0]*AIYZ_[iind-ixm][jind][m+1];
          }
          if (ix > 1) {
            for(m=0;m<=mmax-a-b;m++)
              AI0_[iind][jind][m] += pp*(ix-1)*
               (AI0_[iind-2*ixm][jind][m] - AI0_[iind-2*ixm][jind][m+1]);
            if (Order >= 1)
            for(m=0;m<=mmax-a-b-1;m++) {
              AIX_[iind][jind][m] += pp*(ix-1)*(AIX_[iind-2*ixm][jind][m] -
                                               AIX_[iind-2*ixm][jind][m+1]);
              AIY_[iind][jind][m] += pp*(ix-1)*(AIY_[iind-2*ixm][jind][m] -
                                               AIY_[iind-2*ixm][jind][m+1]);
              AIZ_[iind][jind][m] += pp*(ix-1)*(AIZ_[iind-2*ixm][jind][m] -
                                               AIZ_[iind-2*ixm][jind][m+1]);
            }
            if (Order >= 2)
            for(m=0;m<=mmax-a-b-2;m++) {
              AIXX_[iind][jind][m] += pp*(ix-1)*(AIXX_[iind-2*ixm][jind][m] -
                                                AIXX_[iind-2*ixm][jind][m+1]);
              AIYY_[iind][jind][m] += pp*(ix-1)*(AIYY_[iind-2*ixm][jind][m] -
                                                AIYY_[iind-2*ixm][jind][m+1]);
              AIZZ_[iind][jind][m] += pp*(ix-1)*(AIZZ_[iind-2*ixm][jind][m] -
                                                AIZZ_[iind-2*ixm][jind][m+1]);
              AIXY_[iind][jind][m] += pp*(ix-1)*(AIXY_[iind-2*ixm][jind][m] -
                                                AIXY_[iind-2*ixm][jind][m+1]);
              AIXZ_[iind][jind][m] += pp*(ix-1)*(AIXZ_[iind-2*ixm][jind][m] -
                                                AIXZ_[iind-2*ixm][jind][m+1]);
              AIYZ_[iind][jind][m] += pp*(ix-1)*(AIYZ_[iind-2*ixm][jind][m] -
                                                AIYZ_[iind-2*ixm][jind][m+1]);
            }
          }
          if (jx > 0) {
            for(m=0;m<=mmax-a-b;m++)
              AI0_[iind][jind][m] += pp*jx*
               (AI0_[iind-ixm][jind-jxm][m] - AI0_[iind-ixm][jind-jxm][m+1]);
            if (Order >= 1)
            for(m=0;m<=mmax-a-b-1;m++) {
              AIX_[iind][jind][m] += pp*jx*(AIX_[iind-ixm][jind-jxm][m] -
                                           AIX_[iind-ixm][jind-jxm][m+1]);
              AIY_[iind][jind][m] += pp*jx*(AIY_[iind-ixm][jind-jxm][m] -
                                           AIY_[iind-ixm][jind-jxm][m+1]);
              AIZ_[iind][jind][m] += pp*jx*(AIZ_[iind-ixm][jind-jxm][m] -
                                           AIZ_[iind-ixm][jind-jxm][m+1]);
            }
            if (Order >= 2)
            for(m=0;m<=mmax-a-b-2;m++) {
              AIXX_[iind][jind][m] += pp*jx*(AIXX_[iind-ixm][jind-jxm][m] -
                                            AIXX_[iind-ixm][jind-jxm][m+1]);
              AIYY_[iind][jind][m] += pp*jx*(AIYY_[iind-ixm][jind-jxm][m] -
                                            AIYY_[iind-ixm][jind-jxm][m+1]);
              AIZZ_[iind][jind][m] += pp*jx*(AIZZ_[iind-ixm][jind-jxm][m] -
                                            AIZZ_[iind-ixm][jind-jxm][m+1]);
              AIXY_[iind][jind][m] += pp*jx*(AIXY_[iind-ixm][jind-jxm][m] -
                                            AIXY_[iind-ixm][jind-jxm][m+1]);
              AIXZ_[iind][jind][m] += pp*jx*(AIXZ_[iind-ixm][jind-jxm][m] -
                                            AIXZ_[iind-ixm][jind-jxm][m+1]);
              AIYZ_[iind][jind][m] += pp*jx*(AIYZ_[iind-ixm][jind-jxm][m] -
                                            AIYZ_[iind-ixm][jind-jxm][m+1]);
            }
          }
        }
        else
          throw ProgrammingError("Logic error in Coulomb 1-e OS RR code",__FILE__,__LINE__);
      }
    }

  return;
}

