//
// shift2e.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
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

#include <util/misc/formio.h>
#include <util/misc/consumableresources.h>
#include <chemistry/qc/intv3/macros.h>
#include <chemistry/qc/intv3/int2e.h>

using namespace std;
using namespace sc;

//#undef CHECK_INTEGRAL_ALGORITHM
//#define CHECK_INTEGRAL_ALGORITHM 1

static inline void
iswtch(int *i, int *j)
{
  int tmp;

  tmp = *i;
  *i = *j;
  *j = tmp;
}

/* This initializes the shift routines.  It is called by int_initialize_erep.
 * It is passed the maximum am to be found on each center.
 */
void
Int2eV3::int_init_shiftgc(int order, int am1, int am2, int am3, int am4)
{
  /* The intermediate integral arrays are allocated by the
   * build initialization routine. */

  used_storage_shift_ = 0;

  /* Convert the am1-4 to their canonical ordering. */
  if (am2>am1) {
    iswtch(&am1,&am2);
    }
  if (am4>am3) {
    iswtch(&am3,&am4);
    }
  if ((am3 > am1)||((am3 == am1)&&(am4 > am2))) {
    iswtch(&am1,&am3);
    iswtch(&am2,&am4);
    }

  /* If the center permutation 1<->3 and 2<->4 is performed, then
   * we may need the am for center 2 to be as big as for center 4. */
  if (am4 > am2) am2 = am4;

  /* If derivatives are needed am1 will need to be larger. */
  if (order==1) am1++;
  /* For derivative integral bounds am3 will need to be larger. */
  if (order==1 && int_derivative_bounds) am3++;

  // Set up the new intermediate arrays.
  int e, c, d;
  int ndata34_e = 0;
  for (e=am1; e<=am1+am2; e++) {
    int size_e = INT_NCART(e);
    ndata34_e += size_e;
    }
  int ndata34_f = 0;
  for (d=1; d<=am4; d++) {
    int size_d = INT_NCART(d);
    int size_dm1 = INT_NCART(d-1);
    int off_cp1_dm1 = INT_NCART(am3) * size_dm1;
    int off_c_d = 0;
    for (c=am3; c<=am3+am4-d; c++) {
      int size_c = INT_NCART(c);
      int size_cp1 = INT_NCART(c+1);
      off_c_d += size_c * size_d;
      off_cp1_dm1 += size_cp1 * size_dm1;
      }
    if (off_c_d > ndata34_f) ndata34_f = off_c_d;
    if (off_cp1_dm1 > ndata34_f) ndata34_f = off_cp1_dm1;
    }
  int ndata34 = ndata34_e * ndata34_f;

  int ndata12 = 0;
  int a, b;
  int size_c_d = INT_NCART(am3)*INT_NCART(am4);
  for (b=1; b<=am2; b++) {
    int size_b = INT_NCART(b);
    int size_bm1 = INT_NCART(b-1);
    int off_a_b = 0;
    int off_ap1_bm1 = INT_NCART(am1) * size_bm1 * size_c_d;
    for (a=am1; a<=am1+am2-b; a++) {
      int size_a = INT_NCART(a);
      int size_ap1 = INT_NCART(a+1);
      off_a_b += size_a * size_b * size_c_d;
      off_ap1_bm1 += size_ap1 * size_bm1 * size_c_d;
      }
    if (off_a_b > ndata12) ndata12 = off_a_b;
    if (off_ap1_bm1 > ndata12) ndata12 = off_ap1_bm1;
    }
  int ndatamax = (ndata12>ndata34?ndata12:ndata34);
  buf34 = allocate<double>(ndata34);
  buf12 = allocate<double>(ndata12);
  bufshared = allocate<double>(ndatamax);

  used_storage_shift_ += sizeof(double)*(ndata34+ndata12+ndatamax);

  used_storage_ += used_storage_shift_;
  }

void
Int2eV3::int_done_shiftgc()
{
  used_storage_ -= used_storage_shift_;
  deallocate(buf12);
  deallocate(buf34);
  deallocate(bufshared);
  }

/* This is the principle entry point for the am shifting routines.
 * tam1-4 is the target angular momentum on centers 1-4
 * sh1-4 are the shell numbers on centers 1-4
 */
double *
Int2eV3::int_shiftgcam(int gc1, int gc2, int gc3, int gc4,
                       int tam1, int tam2, int tam3, int tam4, int peAB)
{
  int am1,am2,am3,am4;

  /* Copy the gc{1,2,3,4} into g{1,2,3,4} (static globals). */
  g1 = gc1;
  g2 = gc2;
  g3 = gc3;
  g4 = gc4;

  /* Compute the angular momentum quartet. */
  am1 = tam1;
  am2 = tam2;
  am3 = tam3;
  am4 = tam4;

  // (a0|b0) does need shifting
  if (am2==0 && am4==0) {
    return e0f0_con_ints_array[g1][g2][g3][g4](am1,am3);
    }

  /* Copy the A B equivalency info into a static global variable. */
  eAB = peAB;

  /* Compute the intermediates. */
  AmB[0] =  build.int_v_r10 - build.int_v_r20;
  AmB[1] =  build.int_v_r11 - build.int_v_r21;
  AmB[2] =  build.int_v_r12 - build.int_v_r22;
  CmD[0] =  build.int_v_r30 - build.int_v_r40;
  CmD[1] =  build.int_v_r31 - build.int_v_r41;
  CmD[2] =  build.int_v_r32 - build.int_v_r42;

#if CHECK_INTEGRAL_ALGORITHM > 1
  ExEnv::outn() << "generating ("
       << am1 << "," << am2 << "," << am3 << "," << am4 << ")"
       << ":" << endl;
#endif

  // the (e0|f0) integrals have been initialized
  IntV3Arraydoublep2 &e0f0 = e0f0_con_ints_array[g1][g2][g3][g4];

  // generate (e0|cd) for each needed e
  int e, c, d;
  int off_e = 0;
  int size34 = INT_NCART(am3)*INT_NCART(am4);
  double *buf34_1 = buf34;
  double *buf34_2 = bufshared;
  for (e=am1; e<=am1+am2; e++) {
    int size_e = INT_NCART(e);
    for (d=1; d<=am4; d++) {
      int size_d = INT_NCART(d);
      int size_dm1 = INT_NCART(d-1);
      int off_c_dm1 = 0;
      int off_cp1_dm1 = size_e * INT_NCART(am3) * size_dm1;
      int off_c_d = 0;
      for (c=am3; c<=am3+am4-d; c++) {
        int size_c = INT_NCART(c);
        int size_cp1 = INT_NCART(c+1);
        double *I0001, *I0010, *I0000;
        if (d==am4) {
          I0001 = &buf12[off_e];
          }
        else I0001 = &buf34_1[off_c_d];
        if (d==1) {
          I0010 = e0f0(e,c+1);
          I0000 = e0f0(e,c);
          }
        else {
          I0010 = &buf34_2[off_cp1_dm1];
          I0000 = &buf34_2[off_c_dm1];
          }
        shiftam_34(I0001,I0010,I0000,e,0,c,d);
        off_c_d += size_e * size_c * size_d;
        off_c_dm1 = off_cp1_dm1;
        off_cp1_dm1 += size_e * size_cp1 * size_dm1;
        }
      // swap the buffers.
      double *tmp = buf34_1;
      buf34_1 = buf34_2;
      buf34_2 = tmp;
      }
    off_e += size_e * size34;
    }

  // generate (ab|cd)
  int a, b;
  int size_c_d = size34;
  double *buf12_1 = bufshared;
  double *buf12_2 = buf12;
  for (b=1; b<=am2; b++) {
    int size_b = INT_NCART(b);
    int size_bm1 = INT_NCART(b-1);
    int off_a_b = 0;
    int off_ap1_bm1 = INT_NCART(am1) * size_bm1 * size_c_d;
    int off_a_bm1 = 0;
    for (a=am1; a<=am1+am2-b; a++) {
      int size_a = INT_NCART(a);
      int size_ap1 = INT_NCART(a+1);
      double *I0100 = &buf12_1[off_a_b];
      double *I1000;
      double *I0000;
      if (b==1 && am4 == 0) {
        I1000 = e0f0(a+1,am3);
        if (eAB) I0000 = 0;
        else I0000 = e0f0(a,am3);
        }
      else {
        I1000 = &buf12_2[off_ap1_bm1];
        if (eAB) I0000 = 0;
        else I0000 = &buf12_2[off_a_bm1];
        }
      if (eAB) shiftam_12eAB(I0100,I1000,I0000,a,b,am3,am4);
      else shiftam_12(I0100,I1000,I0000,a,b,am3,am4);
      off_a_b += size_a * size_b * size_c_d;
      off_a_bm1 = off_ap1_bm1;
      off_ap1_bm1 += size_ap1 * size_bm1 * size_c_d;
      }
      // swap the buffers.
      double *tmp = buf12_1;
      buf12_1 = buf12_2;
      buf12_2 = tmp;
    }

  /* Construct the target integrals. */
  return buf12_2;
  }

/* Shift angular momentum from center 1 to center 2.
 * I0100 are the target integrals.
 * am1-4 is the angular momentum on each of the centers in the target set.
 */
void
Int2eV3::shiftam_12(double *I0100, double *I1000, double *I0000,
                    int am1, int am2, int am3, int am4)
{
  int i;
  int i1,k1;
  int size2, size2m134, size34;

#if CHECK_INTEGRAL_ALGORITHM > 1
  ExEnv::outn() << "(" << am1 << "," << am2 << "," << am3 << "," << am4 << ")"
       << " <- "
       << "(" << am1+1 << "," << am2-1 << "," << am3 << "," << am4 << ")"
       << "(" << am1 << "," << am2-1 << "," << am3 << "," << am4 << ")"
       << endl;
#endif

  size2m134 = INT_NCART(am2-1)*INT_NCART(am3)*INT_NCART(am4);
  size34    = INT_NCART(am3)*INT_NCART(am4);
  size2     = INT_NCART(am2);

  int size_zcontrib = am2*size34;
  int size_xcontrib = (size2-(am2+1))*size34;

  double AmB0 = AmB[0];
  double AmB1 = AmB[1];
  double AmB2 = AmB[2];

  /* Loop over the target integrals. */
  double *restrictxx I0100i=I0100;
  int cartindex1 = 0;
  for (i1=0; i1<=am1; i1++) {
    for (k1=0; k1<=am1-i1; k1++) {
      //int j1 = am1 - i1 - k1;
      int ci1x1 = (cartindex1 + am1 + 2) * size2m134;
      int ci1y1 = (cartindex1 + i1) * size2m134;
      int ci1z1 = (cartindex1 + i1 + 1)  * size2m134;
      //note:
      //ci1x1 = INT_CARTINDEX(am1+1,i1+1,j1) * size2m134;
      //ci1y1 = INT_CARTINDEX(am1+1,i1,j1+1) * size2m134;
      //ci1z1 = INT_CARTINDEX(am1+1,i1,j1) * size2m134;
      int ci1 = cartindex1 * size2m134;
      // i2 == 0, k2 == 0, j2 == am2 (>0)
      double *I1000i=&I1000[ci1y1];
      double *I0000i=&I0000[ci1];
      for (i=0; i<size34; i++) {
        I0100i[i] = I1000i[i] + I0000i[i] * AmB1;
        }
      I0100i=&I0100i[size34];
      // i2 == 0, k2 > 0
      I1000i=&I1000[ci1z1];
      I0000i=&I0000[ci1];
      for (i=0; i<size_zcontrib; i++) {
        I0100i[i] = I1000i[i] + I0000i[i] * AmB2;
        }
      I0100i=&I0100i[size_zcontrib];
      // i2 >= 1
      I1000i=&I1000[ci1x1];
      I0000i=&I0000[ci1];
      for (i=0; i<size_xcontrib; i++) {
        I0100i[i] = I1000i[i] + I0000i[i] * AmB0;
        }
      I0100i=&I0100i[size_xcontrib];

      cartindex1++;
      }
    }
  }


/* Shift angular momentum from center 1 to center 2 when centers
 * one and two are the same.
 * I0100 are the target integrals.
 * am1-4 is the angular momentum on each of the centers in the target set.
 */
void
Int2eV3::shiftam_12eAB(double *I0100, double *I1000, double *I0000,
                       int am1, int am2, int am3, int am4)
{
  int i;
  int i1,k1;
  int size2, size2m134, size34;

#if CHECK_INTEGRAL_ALGORITHM > 1
  ExEnv::outn() << "(" << am1 << "," << am2 << "," << am3 << "," << am4 << ")"
       << " <- "
       << "(" << am1+1 << "," << am2-1 << "," << am3 << "," << am4 << ")"
       << "(" << am1 << "," << am2-1 << "," << am3 << "," << am4 << ")"
       << endl;
#endif

  size2m134 = INT_NCART(am2-1)*INT_NCART(am3)*INT_NCART(am4);
  size34    = INT_NCART(am3)*INT_NCART(am4);
  size2     = INT_NCART(am2);

  int size_zcontrib = am2*size34;
  int size_xcontrib = (size2-(am2+1))*size34;

  /* Loop over the target integrals. */
  double *restrictxx I0100i=I0100;
  int cartindex1 = 0;
  for (i1=0; i1<=am1; i1++) {
    for (k1=0; k1<=am1-i1; k1++) {
      //int j1 = am1 - i1 - k1;
      int ci1x1 = (cartindex1 + am1 + 2) * size2m134;
      int ci1y1 = (cartindex1 + i1) * size2m134;
      int ci1z1 = (cartindex1 + i1 + 1)  * size2m134;
      //note:
      //ci1x1 = INT_CARTINDEX(am1+1,i1+1,j1) * size2m134;
      //ci1y1 = INT_CARTINDEX(am1+1,i1,j1+1) * size2m134;
      //ci1z1 = INT_CARTINDEX(am1+1,i1,j1) * size2m134;
      // i2 == 0, k2 == 0, j2 == am2 (>0)
      double *I1000i=&I1000[ci1y1];
      for (i=0; i<size34; i++) {
        I0100i[i] = I1000i[i];
        }
      I0100i=&I0100i[size34];
      // i2 == 0, k2 > 0
      I1000i=&I1000[ci1z1];
      for (i=0; i<size_zcontrib; i++) {
        I0100i[i] = I1000i[i];
        }
      I0100i=&I0100i[size_zcontrib];
      // i2 >= 1
      I1000i=&I1000[ci1x1];
      for (i=0; i<size_xcontrib; i++) {
        I0100i[i] = I1000i[i];
        }
      I0100i=&I0100i[size_xcontrib];

      cartindex1++;
      }
    }
  }

void
Int2eV3::shiftam_34(double *restrictxx I0001, double *I0010, double *I0000,
                    int am1, int am2, int am3, int am4)
{
  int i1,k1,cartindex1;
  int i2,k2,cartindex2;
  int i3,k3,cartindex3;
  int i4,k4,cartindex4;
  int cartindex1234;
  int size23p14m1,size3p14m1,size4m1,size234m1,size34m1;

#if CHECK_INTEGRAL_ALGORITHM > 1
  ExEnv::outn() << "(" << am1 << "," << am2 << "," << am3 << "," << am4 << ")"
       << " <- "
       << "(" << am1 << "," << am2 << "," << am3+1 << "," << am4-1 << ")"
       << "(" << am1 << "," << am2 << "," << am3 << "," << am4-1 << ")"
       << endl;
#endif

  size23p14m1 = INT_NCART(am2)*INT_NCART(am3+1)*INT_NCART(am4-1);
  size3p14m1 = INT_NCART(am3+1)*INT_NCART(am4-1);
  size4m1 = INT_NCART(am4-1);

  size234m1 = INT_NCART(am2)*INT_NCART(am3)*INT_NCART(am4-1);
  size34m1 = INT_NCART(am3)*INT_NCART(am4-1);

  double CmD0 = CmD[0];
  double CmD1 = CmD[1];
  double CmD2 = CmD[2];

  /* Loop over the target integrals. */
  cartindex1234 = 0;
  cartindex1 = 0;
  for (i1=0; i1<=am1; i1++) {
    for (k1=0; k1<=am1-i1; k1++) {
      //int j1 = am1 - i1 - k1;
      int ci1_I0010 = cartindex1 * size23p14m1;
      int ci1_I0000 = cartindex1 * size234m1;
      cartindex2 = 0;
      for (i2=0; i2<=am2; i2++) {
        for (k2=0; k2<=am2-i2; k2++) {
          //int j2 = am2 - i2 - k2;
          int ci2_I0010 = ci1_I0010 + cartindex2 * size3p14m1;
          int ci2_I0000 = ci1_I0000 + cartindex2 * size34m1;
          cartindex3 = 0;
          for (i3=0; i3<=am3; i3++) {
            for (k3=0; k3<=am3-i3; k3++) {
              //int j3 = am3 - i3 - k3;
              //note: cartindex3 + am3 + 2 = INT_CARTINDEX(am3+1,i3+1,j3)
              int ci3_I0010 = ci2_I0010
                              + (cartindex3 + am3 + 2)*size4m1;
              int ci3_I0000 = ci2_I0000 + cartindex3*size4m1;
              //cartindex4 = 0;
              // this routine called only when am4 > 0
              ///// CASE 1: i4 = 0 k4 = 0 j4 = am4; shift on y
              //note: j4 = am4;
              //note: cartindex4 - i4 = INT_CARTINDEX(am4-1,i4,j4-1)
              //note: cartindex3 - i3 = INT_CARTINDEX(am3+1,i3,j3+1)
              int ci3 = cartindex3 + i3;
              I0001[cartindex1234]
                =   I0010[ci2_I0010 + ci3 * size4m1]
                + I0000[ci3_I0000]
                * CmD1;
              cartindex1234++;
              //cartindex4++;
              ///// CASE 2: i4 = 0 k4 > 0; shift on z
              ci3++;
              for (int ci4=0; ci4<am4; ci4++) {
                //note: j4 = am4 - i4 - k4;
                //note: cartindex4 - i4 - 1 = INT_CARTINDEX(am4-1,i4,j4)
                //note: ci4 = cartindex4 - i4 - 1;
                //note: cartindex3 - i3 - 1 = INT_CARTINDEX(am3+1,i3,j3)
                I0001[cartindex1234]
                  =   I0010[ci2_I0010 + ci3 * size4m1 + ci4 ]
                  + I0000[ci3_I0000 + ci4 ]
                  * CmD2;
                cartindex1234++;
                //cartindex4++;
                }
              ///// CASE 3: i4 > 0; shift on x
              int ncart_remain = INT_NCART(am4) - (am4+1);
              for (int ci4=0; ci4<ncart_remain; ci4++) {
                  //note: j4 = am4 - i4 - k4;
                  //note: cartindex4 - am4 - 1 = INT_CARTINDEX(am4-1,i4-1,j4)
                  //note: ci4 = cartindex4 - am4 - 1;
                  I0001[cartindex1234]
                    =   I0010[ci3_I0010 + ci4]
                    + I0000[ci3_I0000 + ci4]
                    * CmD0;
                  cartindex1234++;
                  //cartindex4++;
                  }
              cartindex3++;
              }
            }
          cartindex2++;
          }
        }
      cartindex1++;
      }
    }
  }


/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
