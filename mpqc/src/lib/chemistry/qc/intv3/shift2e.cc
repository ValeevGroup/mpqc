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
#include <chemistry/qc/intv3/macros.h>
#include <chemistry/qc/intv3/int2e.h>

static inline void
iswtch(int *i, int *j)
{
  int tmp;

  tmp = *i;
  *i = *j;
  *j = tmp;
}


static void
fail()
{
  cerr << scprintf("failing module:\n%s",__FILE__) << endl;
  abort();
  }

ShiftIntermediates::ShiftIntermediates()
{
  maxused_ = 0;
  data_ = 0;
  ndata_ = 0;
  nused_ = 0;
  l1_ = 0;
  l2_ = 0;
  l3_ = 0;
  l4_ = 0;
}

ShiftIntermediates::~ShiftIntermediates()
{
#if CHECK_INTEGRAL_ALGORITHM
  cout << "ShiftIntermediates: DTOR: wasted "
       << ndata_ - maxused_
       << " of " << ndata_
       << endl;
#endif
  delete[] data_;
}

void
ShiftIntermediates::out_of_memory(int i,int j,int k,int l)
{
  int size = INT_NCART(i)*INT_NCART(j)*INT_NCART(k)*INT_NCART(l);
  cerr << "ShiftIntermediates: didn't preallocate enough memory" << endl;
  cerr << "  nused_ = " << nused_ << " size = " << size << endl;
  cerr << "  ndata_ = " << ndata_ << endl;
  abort();
}

void
ShiftIntermediates::clear(int am1,int am2,int am3,int am4)
{
  nused_ = 0;
  double *****data = shell_.data();
  for (int i=0; i<=am1+am2; i++) {
    double ****datai = data[i];
    for (int j=0; j<=am2; j++) {
      double ***dataij = datai[j];
      for (int k=0; k<=am3+am4; k++) {
        double **dataijk = dataij[k];
        for (int l=0; l<=am4; l++) {
          dataijk[l] = 0;
          }
        }
      }
    }
}

double *
ShiftIntermediates::allocate(int i,int j,int k,int l)
{
  int size = INT_NCART(i)*INT_NCART(j)*INT_NCART(k)*INT_NCART(l);
  //cout << "ShiftIntermediates::allocate:"
  //     << " size = " << size
  //     << " nused_ = " << nused_
  //     << " ndata_ = " << ndata_
  //     << " data_ = " << (void*)data_
  //     << endl;
  double *data = &data_[nused_];
  shell_(i,j,k,l) = data;
  if (nused_ + size > ndata_) out_of_memory(i,j,k,l);
  nused_ += size;
  if (nused_ > maxused_) maxused_ = nused_;
  return data;
}

void
ShiftIntermediates::set_l(int l1,int l2,int l3,int l4)
{
#if CHECK_INTEGRAL_ALGORITHM
  cout << "ShiftIntermediates: set_l: wasted "
       << ndata_ - maxused_
       << " of " << ndata_
       << endl;
#endif

  l1_ = l1;
  l2_ = l2;
  l3_ = l3;
  l4_ = l4;
  shell_.set_dim(l1+l2+1,l2+1,l3+l4+1,l4+1);

  delete[] data_;
  ndata_ = 0;
  nused_ = 0;
  maxused_ = 0;
  int i;
  for (i=l1; i<=l1+l2; i++) {
    int szi = INT_NCART(i);
    for (int k=l3; k<=l3+l4; k++) {
      int szk = INT_NCART(k);
      for (int l=1; l<=k && l+k<=l3+l4 && l<=l4; l++) {
        int szl = INT_NCART(l);
        ndata_ += szi*szk*szl;
        }
      }
    }
  int szkl = INT_NCART(l3)*INT_NCART(l4);
  for (i=l1; i<=l1+l2; i++) {
    int szi = INT_NCART(i);
    for (int j=1; j<=i && i+j<=l1+l2 && j<=l2; j++) {
      int szj = INT_NCART(j);
      ndata_ += szi*szj*szkl;
      }
    }
  data_ = new double[ndata_];
#if CHECK_INTEGRAL_ALGORITHM
  cout << "ShiftIntermediates: allocated: " << ndata_ << endl;
#endif
}

int
ShiftIntermediates::nbyte()
{
  return ndata_ * sizeof(double) + shell_.nbyte();
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

  // Set up the intermediate array.
  shiftinter_.set_l(am1,am2,am3,am4);
  used_storage_shift_ += shiftinter_.nbyte();

#if CHECK_INTEGRAL_ALGORITHM
  cout << "shiftinter: " << shiftinter_.nbyte() << endl;
#endif

  used_storage_ += used_storage_shift_;
  }

void
Int2eV3::int_done_shiftgc()
{
  used_storage_ -= used_storage_shift_;
  shiftinter_.set_l(0,0,0,0);
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

  /* Copy the A B equivalency info into a static global variable. */
  eAB = peAB;

  /* Compute the intermediates. */
  AmB[0] =  build.int_v_r10 - build.int_v_r20;
  AmB[1] =  build.int_v_r11 - build.int_v_r21;
  AmB[2] =  build.int_v_r12 - build.int_v_r22;
  CmD[0] =  build.int_v_r30 - build.int_v_r40;
  CmD[1] =  build.int_v_r31 - build.int_v_r41;
  CmD[2] =  build.int_v_r32 - build.int_v_r42;

  /* Mark all of the intermediates as being noncomputed. */
  shiftinter_.clear(am1,am2,am3,am4);

#if CHECK_INTEGRAL_ALGORITHM
  cout << "generating ("
       << am1 << "," << am2 << "," << am3 << "," << am4 << ")"
       << ":" << endl;
#endif

  /* Construct the target integrals. */
  return shiftint(am1,am2,am3,am4);
  }

double *
Int2eV3::shiftint(int am1, int am2, int am3, int am4)
{
  double *buffer;

#if 0
  cout << scprintf(" S[%d(%d),%d(%d),%d(%d),%d(%d)]",
         am1,g1,
         am2,g2,
         am3,g3,
         am4,g4);
#endif

  if (am2==0 && am4==0) {
    return e0f0_con_ints_array[g1][g2][g3][g4](am1,am3);
    }
  
  /* If the integral is known, then return the pointer to its buffer. */
  if (buffer = shiftinter_(am1,am2,am3,am4)) return buffer;

  /* Find the preallocated storage for the target integrals. */
  buffer = shiftinter_.allocate(am1,am2,am3,am4);

  if (!buffer) {
    cerr << "Int2eV3::shift2e: internal error: "
         << "storage not allocated for target" << endl;
    abort();
    }

  /* Should we shift to 2 or to 4? */
  if (choose_shift(am1,am2,am3,am4) == 2) {
#if CHECK_INTEGRAL_ALGORITHM
    shiftam_12(buffer,am1,am2,am3,am4);
#else
    if (eAB) shiftam_12eAB(buffer,am1,am2,am3,am4);
    else     shiftam_12(buffer,am1,am2,am3,am4);
#endif
    }
  else {
    shiftam_34(buffer,am1,am2,am3,am4);
    }

  return buffer;
  }

/* Returns 2 if we should shift am from 1 to 2 next and
 * returns 4 if we should shift am from 3 to 4 next.
 * No attempt has been made to optimize the choice.
 */
int
Int2eV3::choose_shift(int am1, int am2, int am3, int am4)
{
  if (am2 == 0) {
    if (am4 == 0) {
      cerr << scprintf("shift: build routines missed (%d,0,%d,0)\n",am1,am3);
      fail();
      }
    return 4;
    }
  if (am4 == 0) return 2;

#if 1
  // always build on one electron first since the storage allocation
  // now computes an upper bound for intermediate storage based on this
  return 2;
#else
  int nneed2 = 0;
  int nneed4 = 0;
  if (!shiftinter_(am1+1,am2-1,am3,am4)) {
    nneed2 += INT_NCART(am1+1)*INT_NCART(am2-1)*INT_NCART(am3)*INT_NCART(am4);
    }
  if (!shiftinter_(am1,am2-1,am3,am4)) {
    nneed2 += INT_NCART(am1)*INT_NCART(am2-1)*INT_NCART(am3)*INT_NCART(am4);
    }
  if (!shiftinter_(am1,am2,am3+1,am4-1)) {
    nneed4 += INT_NCART(am1)*INT_NCART(am2)*INT_NCART(am3+1)*INT_NCART(am4-1);
    }
  if (!shiftinter_(am1,am2,am3,am4-1)) {
    nneed4 += INT_NCART(am1)*INT_NCART(am2)*INT_NCART(am3)*INT_NCART(am4-1);
    }
  if (nneed2 <= nneed4) return 2;
  return 4;
#endif
  }

/* Shift angular momentum from center 1 to center 2.
 * I0100 are the target integrals.
 * am1-4 is the angular momentum on each of the centers in the target set.
 */
void
Int2eV3::shiftam_12(double *I0100, int am1, int am2, int am3, int am4)
{
  int i;
  double *I1000;
  double *I0000;
  int i1,j1,k1;
  int i2,j2,k2;
  int size2, size2m134, size34;

#if CHECK_INTEGRAL_ALGORITHM
  cout << "(" << am1 << "," << am2 << "," << am3 << "," << am4 << ")"
       << " <- "
       << "(" << am1+1 << "," << am2-1 << "," << am3 << "," << am4 << ")"
       << "(" << am1 << "," << am2-1 << "," << am3 << "," << am4 << ")"
       << endl;
#endif

  I1000 = shiftint(am1+1,am2-1,am3,am4);
  I0000 = shiftint(am1,am2-1,am3,am4);

  size2m134 = INT_NCART(am2-1)*INT_NCART(am3)*INT_NCART(am4);
  size34    = INT_NCART(am3)*INT_NCART(am4);
  size2     = INT_NCART(am2);

  int size_zcontrib = am2*size34;
  int size_xcontrib = (size2-(am2+1))*size34;

  double AmB0 = AmB[0];
  double AmB1 = AmB[1];
  double AmB2 = AmB[2];

  /* Loop over the target integrals. */
  double *I0100i=I0100;
  int cartindex1 = 0;
  for (i1=0; i1<=am1; i1++) {
    for (k1=0; k1<=am1-i1; k1++) {
      j1 = am1 - i1 - k1;
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
Int2eV3::shiftam_12eAB(double *I0100, int am1, int am2, int am3, int am4)
{
  int i;
  double *I1000;
  int i1,j1,k1;
  int i2,j2,k2;
  int size2, size2m134, size34;

#if CHECK_INTEGRAL_ALGORITHM
  cout << "(" << am1 << "," << am2 << "," << am3 << "," << am4 << ")"
       << " <- "
       << "(" << am1+1 << "," << am2-1 << "," << am3 << "," << am4 << ")"
       << "(" << am1 << "," << am2-1 << "," << am3 << "," << am4 << ")"
       << endl;
#endif

  I1000 = shiftint(am1+1,am2-1,am3,am4);

  size2m134 = INT_NCART(am2-1)*INT_NCART(am3)*INT_NCART(am4);
  size34    = INT_NCART(am3)*INT_NCART(am4);
  size2     = INT_NCART(am2);

  int size_zcontrib = am2*size34;
  int size_xcontrib = (size2-(am2+1))*size34;

  double AmB0 = AmB[0];
  double AmB1 = AmB[1];
  double AmB2 = AmB[2];

  /* Loop over the target integrals. */
  double *I0100i=I0100;
  int cartindex1 = 0;
  for (i1=0; i1<=am1; i1++) {
    for (k1=0; k1<=am1-i1; k1++) {
      j1 = am1 - i1 - k1;
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
Int2eV3::shiftam_34(double *I0001, int am1, int am2, int am3, int am4)
{
  double *I0010;
  double *I0000;
  int i1,j1,k1,cartindex1;
  int i2,j2,k2,cartindex2;
  int i3,j3,k3,cartindex3;
  int i4,j4,k4,cartindex4;
  int cartindex1234;
  int size23p14m1,size3p14m1,size4m1,size234m1,size34m1;

#if CHECK_INTEGRAL_ALGORITHM
  cout << "(" << am1 << "," << am2 << "," << am3 << "," << am4 << ")"
       << " <- "
       << "(" << am1 << "," << am2 << "," << am3+1 << "," << am4-1 << ")"
       << "(" << am1 << "," << am2 << "," << am3 << "," << am4-1 << ")"
       << endl;
#endif

  I0010 = shiftint(am1,am2,am3+1,am4-1);
  I0000 = shiftint(am1,am2,am3,am4-1);

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
      j1 = am1 - i1 - k1;
      int ci1_I0010 = cartindex1 * size23p14m1;
      int ci1_I0000 = cartindex1 * size234m1;
      cartindex2 = 0;
      for (i2=0; i2<=am2; i2++) {
        for (k2=0; k2<=am2-i2; k2++) {
          j2 = am2 - i2 - k2;
          int ci2_I0010 = ci1_I0010 + cartindex2 * size3p14m1;
          int ci2_I0000 = ci1_I0000 + cartindex2 * size34m1;
          cartindex3 = 0;
          for (i3=0; i3<=am3; i3++) {
            for (k3=0; k3<=am3-i3; k3++) {
              j3 = am3 - i3 - k3;
              //note: cartindex3 + am3 + 2 = INT_CARTINDEX(am3+1,i3+1,j3)
              int ci3_I0010 = ci2_I0010
                              + (cartindex3 + am3 + 2)*size4m1;
              int ci3_I0000 = ci2_I0000 + cartindex3*size4m1;
              cartindex4 = 0;
              for (i4=0; i4<=am4; i4++) {
                for (k4=0; k4<=am4-i4; k4++) {
                  j4 = am4 - i4 - k4;

                  if (i4) {
                    //note: cartindex4 - am4 - 1 = INT_CARTINDEX(am4-1,i4-1,j4)
                    int ci4 = cartindex4 - am4 - 1;
                    I0001[cartindex1234]
                      =   I0010[ci3_I0010 + ci4]
                      + I0000[ci3_I0000 + ci4]
                      * CmD0;
                    }
                  else if (j4) {
                    //note: cartindex4 - i4 = INT_CARTINDEX(am4-1,i4,j4-1)
                    int ci4 = cartindex4 - i4;
                    //note: cartindex3 - i3 = INT_CARTINDEX(am3+1,i3,j3+1)
                    int ci3 = cartindex3 + i3;
                    I0001[cartindex1234]
                     =   I0010[ci2_I0010
                               + ci3 * size4m1
                               + ci4 ]
                       + I0000[ci3_I0000
                               + ci4 ]
                         * CmD1;
                    }
                  else if (k4) {
                    //note: cartindex4 - i4 - 1 = INT_CARTINDEX(am4-1,i4,j4)
                    int ci4 = cartindex4 - i4 - 1;
                    //note: cartindex3 - i3 - 1 = INT_CARTINDEX(am3+1,i3,j3)
                    int ci3 = cartindex3 + i3 + 1;
                    I0001[cartindex1234]
                     =   I0010[ci2_I0010
                               + ci3 * size4m1
                               + ci4 ]
                       + I0000[ci3_I0000
                               + ci4 ]
                         * CmD2;
#if 0
   if (cartindex1234 == 4) {
     cout << scprintf(" building with % f + % f * % f ",
                         I0010[  INT_CARTINDEX(am1,i1,j1) * size23p14m1
                               + INT_CARTINDEX(am2,i2,j2) * size3p14m1
                               + INT_CARTINDEX(am3+1,i3,j3) * size4m1
                               + INT_CARTINDEX(am4-1,i4,j4) ],
                         I0000[  INT_CARTINDEX(am1,i1,j1) * size234m1
                               + INT_CARTINDEX(am2,i2,j2) * size34m1
                               + INT_CARTINDEX(am3,i3,j3) * size4m1
                               + INT_CARTINDEX(am4-1,i4,j4) ],
                         CmD[2]);
      }
#endif
                    }

#if 0
                  if ((!am1) == (!am2) == am3 == am4) {
                    cout << scprintf("assigned I0001[%d] = % f\n",
                           cartindex1234,
                           I0001[cartindex1234]);
                    }
#endif

                  cartindex1234++;
                  cartindex4++;
                  }
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
