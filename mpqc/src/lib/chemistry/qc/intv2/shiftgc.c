
/* This file contains routines that are responsible for
 * shifting angular momentum from center 1 to 2 and from
 * center 3 to 4. */

/* $Log$
 * Revision 1.3  1994/08/26 22:45:54  etseidl
 * fix a bunch of warnings, get rid of rcs id's, get rid of bread/bwrite and
 * fread/fwrite modules
 *
 * Revision 1.2  1993/12/30  13:33:09  etseidl
 * mostly rcs id stuff
 *
 * Revision 1.4  1992/06/17  22:05:20  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.3  1992/05/26  20:25:23  jannsen
 * make derivative bounds checking optional
 * add code to allow bound intermediates computable in shell blocks
 *
 * Revision 1.2  1992/05/13  18:29:44  jannsen
 * added bounds checking for derivative integrals
 *
 * Revision 1.1.1.1  1992/03/17  16:33:54  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:33:52  seidl
 * Initial revision
 *
 * Revision 1.1  1992/02/03  15:49:49  seidl
 * Initial revision
 *
 * Revision 1.2  1992/01/30  01:31:10  cljanss
 * added int_done_shiftgc to free allocated data when done
 *
 * Revision 1.1  1991/12/14  00:19:36  cljanss
 * Initial revision
 * */
/*
 * Revision 1.5  91/09/28  19:26:55  cljanss
 * Switch to new file naming scheme
 * 
 * Revision 1.4  91/09/28  18:17:20  cljanss
 * Intermediates are no longer stored, if requested with flags.
 * 
 * Revision 1.3  91/06/19  15:19:38  janssen
 * First stab at two electron derivative integrals--might work might not
 * 
 * Revision 1.2  1991/06/16  19:32:22  janssen
 * improved performance when center A == center B
 *
 * Revision 1.1  1991/06/16  16:40:07  janssen
 * Initial revision
 * */

#include <stdio.h>
#include <tmpl.h>
#include <math/array/math_lib.h>

#include "atoms.h"
#include "int_macros.h"

#include "inter.h"
#include "buildinter.h"

#include "shiftgc.gbl"
#include "shiftgc.lcl"

/* The general contraction numbers. */
static int g1,g2,g3,g4;

/* A[] - B[] */
static double AmB[3];

/* C[] - D[] */
static double CmD[3];

/* Boolean array which gives whether or not a set of integrals has been
 * computed. */
static int_array4_t shiftinthave;

static int eAB;

/* This initializes the shift routines.  It is called by int_initialize_erep.
 * It is passed the maximum am to be found on each center.
 */
GLOBAL_FUNCTION VOID
int_init_shiftgc(order,am1,am2,am3,am4)
int order;
int am1;
int am2;
int am3;
int am4;
{
  /* The intermediate integral arrays are allocated by the
   * build initialization routine. */

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

  /* Allocate the array giving what has already been computed. */
  if (allocbn_int_array4(&shiftinthave,"n1 n2 n3 n4",
                         am1+am2+1,am2+1,am3+am4+1,am4+1)) {
    fprintf(stderr,"problem allocating shiftinthave");
    fail();
    }
  }

GLOBAL_FUNCTION VOID
int_done_shiftgc()
{
  free_int_array4(&shiftinthave);
  }

LOCAL_FUNCTION VOID
init_shiftinthave(am1,am2,am3,am4)
int am1;
int am2;
int am3;
int am4;
{
  int i,j,k,l;

  for (i=0; i<=am1+am2; i++) {
    for (j=0; j<=am2; j++) {
      for (k=0; k<=am3+am4; k++) {
        for (l=0; l<=am4; l++) {
          /* The integrals for j==0 and l==0 have been precomputed
           * by the build routine. */
          if ((j==0) && (l==0)) shiftinthave.i[i][j][k][l] = 1;
          else                  shiftinthave.i[i][j][k][l] = 0;
          }
        }
      }
    }
  }

/* This is the principle entry point for the am shifting routines.
 * tam1-4 is the target angular momentum on centers 1-4
 * sh1-4 are the shell numbers on centers 1-4
 */
GLOBAL_FUNCTION VOID
int_shiftgcam(gc1,gc2,gc3,gc4, tam1,tam2,tam3,tam4, peAB)
int gc1;
int gc2;
int gc3;
int gc4;
int tam1;
int tam2;
int tam3;
int tam4;
int peAB;
{
  int am1,am2,am3,am4;
  double *buffer;

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
  AmB[0] =  int_v_r10 - int_v_r20;
  AmB[1] =  int_v_r11 - int_v_r21;
  AmB[2] =  int_v_r12 - int_v_r22;
  CmD[0] =  int_v_r30 - int_v_r40;
  CmD[1] =  int_v_r31 - int_v_r41;
  CmD[2] =  int_v_r32 - int_v_r42;

  /* Mark all of the intermediates as being noncomputed. */
  init_shiftinthave(am1,am2,am3,am4);

  /* Construct the target integrals. */
  buffer = shiftint(am1,am2,am3,am4);

  }

LOCAL_FUNCTION double *
shiftint(am1,am2,am3,am4)
int am1;
int am2;
int am3;
int am4;
{
  double *buffer;

#if 0
  printf(" S[%d(%d),%d(%d),%d(%d),%d(%d)]",
         am1,g1,
         am2,g2,
         am3,g3,
         am4,g4);
#endif

  /* If the integral is known, then return the pointer to its buffer. */
  if (shiftinthave.i[am1][am2][am3][am4])
    return int_con_ints_array[g1][g2][g3][g4].dp[am1][am2][am3][am4];

  /* Find the preallocated storage for the target integrals. */
  buffer = int_con_ints_array[g1][g2][g3][g4].dp[am1][am2][am3][am4];

  /* Should we shift to 2 or to 4? */
  if (choose_shift(am1,am2,am3,am4) == 2) {
    if (eAB) shiftam_12eAB(buffer,am1,am2,am3,am4);
    else     shiftam_12(buffer,am1,am2,am3,am4);
    }
  else {
    shiftam_34(buffer,am1,am2,am3,am4);
    }

  /* Put the integrals in the list of precomputed integrals. */
  shiftinthave.i[am1][am2][am3][am4] = 1;

  return buffer;
  }

/* Returns 2 if we should shift am from 1 to 2 next and
 * returns 4 if we should shift am from 3 to 4 next.
 * No attempt has been made to optimize the choice.
 */
LOCAL_FUNCTION int
choose_shift(am1,am2,am3,am4)
int am1;
int am2;
int am3;
int am4;
{
  int nneed2 = 0;
  int nneed4 = 0;

  if (am2 == 0) {
    if (am4 == 0) {
      fprintf(stderr,"shift: build routines missed (%d,0,%d,0)\n",am1,am3);
      fail();
      }
    return 4;
    }
  if (am4 == 0) return 2;

  if (!shiftinthave.i[am1+1][am2-1][am3][am4]) {
    nneed2 += INT_NCART(am1+1)*INT_NCART(am2-1)*INT_NCART(am3)*INT_NCART(am4);
    }
  if (!shiftinthave.i[am1][am2-1][am3][am4]) {
    nneed2 += INT_NCART(am1)*INT_NCART(am2-1)*INT_NCART(am3)*INT_NCART(am4);
    }
  if (!shiftinthave.i[am1][am2][am3+1][am4-1]) {
    nneed4 += INT_NCART(am1)*INT_NCART(am2)*INT_NCART(am3+1)*INT_NCART(am4-1);
    }
  if (!shiftinthave.i[am1][am2-1][am3][am4-1]) {
    nneed4 += INT_NCART(am1)*INT_NCART(am2-1)*INT_NCART(am3)*INT_NCART(am4-1);
    }
  if (nneed2 <= nneed4) return 2;
  return 4;
  }

/* Shift angular momentum from center 1 to center 2.
 * I0100 are the target integrals.
 * am1-4 is the angular momentum on each of the centers in the target set.
 */
LOCAL_FUNCTION VOID
shiftam_12(I0100,am1,am2,am3,am4)
double *I0100;
int am1;
int am2;
int am3;
int am4;
{
  double *I1000;
  double *I0000;
  int i1,j1,k1;
  int i2,j2,k2;
  int cartindex34;
  int cartindex1234;
  int size2m134, size34;

  I1000 = shiftint(am1+1,am2-1,am3,am4);
  I0000 = shiftint(am1,am2-1,am3,am4);

  size2m134 = INT_NCART(am2-1)*INT_NCART(am3)*INT_NCART(am4);
  size34    = INT_NCART(am3)*INT_NCART(am4);

  /* Loop over the target integrals. */
  cartindex1234 = 0;
  for (i1=0; i1<=am1; i1++) {
    for (k1=0; k1<=am1-i1; k1++) {
      j1 = am1 - i1 - k1;
      for (i2=0; i2<=am2; i2++) {
        for (k2=0; k2<=am2-i2; k2++) {
          j2 = am2 - i2 - k2;
          for (cartindex34=0; cartindex34<size34; cartindex34++) {

            if (i2) {
              I0100[cartindex1234]
               =   I1000[  INT_CARTINDEX(am1+1,i1+1,j1) * size2m134
                         + INT_CARTINDEX(am2-1,i2-1,j2) * size34
                         + cartindex34 ]
                 + I0000[  INT_CARTINDEX(am1,i1,j1) * size2m134
                         + INT_CARTINDEX(am2-1,i2-1,j2) * size34
                         + cartindex34 ]
                   * AmB[0];
              }
            else if (j2) {
              I0100[cartindex1234]
               =   I1000[  INT_CARTINDEX(am1+1,i1,j1+1) * size2m134
                         + INT_CARTINDEX(am2-1,i2,j2-1) * size34
                         + cartindex34 ]
                 + I0000[  INT_CARTINDEX(am1,i1,j1) * size2m134
                         + INT_CARTINDEX(am2-1,i2,j2-1) * size34
                         + cartindex34 ]
                   * AmB[1];
              }
            else {
              I0100[cartindex1234]
               =   I1000[  INT_CARTINDEX(am1+1,i1,j1) * size2m134
                         + INT_CARTINDEX(am2-1,i2,j2) * size34
                         + cartindex34 ]
                 + I0000[  INT_CARTINDEX(am1,i1,j1) * size2m134
                         + INT_CARTINDEX(am2-1,i2,j2) * size34
                         + cartindex34 ]
                   * AmB[2];
              }

            cartindex1234++;
            }
          }
        }
      }
    }
  }


/* Shift angular momentum from center 1 to center 2 when centers
 * one and two are the same.
 * I0100 are the target integrals.
 * am1-4 is the angular momentum on each of the centers in the target set.
 */
LOCAL_FUNCTION VOID
shiftam_12eAB(I0100,am1,am2,am3,am4)
double *I0100;
int am1;
int am2;
int am3;
int am4;
{
  double *I1000;
  int i1,j1,k1;
  int i2,j2,k2;
  int cartindex34;
  int cartindex1234;
  int size2m134, size34;

  I1000 = shiftint(am1+1,am2-1,am3,am4);

  size2m134 = INT_NCART(am2-1)*INT_NCART(am3)*INT_NCART(am4);
  size34    = INT_NCART(am3)*INT_NCART(am4);

  /* Loop over the target integrals. */
  cartindex1234 = 0;
  for (i1=0; i1<=am1; i1++) {
    for (k1=0; k1<=am1-i1; k1++) {
      j1 = am1 - i1 - k1;
      for (i2=0; i2<=am2; i2++) {
        for (k2=0; k2<=am2-i2; k2++) {
          j2 = am2 - i2 - k2;
          for (cartindex34=0; cartindex34<size34; cartindex34++) {

            if (i2) {
              I0100[cartindex1234]
               =   I1000[  INT_CARTINDEX(am1+1,i1+1,j1) * size2m134
                         + INT_CARTINDEX(am2-1,i2-1,j2) * size34
                         + cartindex34 ];
              }
            else if (j2) {
              I0100[cartindex1234]
               =   I1000[  INT_CARTINDEX(am1+1,i1,j1+1) * size2m134
                         + INT_CARTINDEX(am2-1,i2,j2-1) * size34
                         + cartindex34 ];
              }
            else {
              I0100[cartindex1234]
               =   I1000[  INT_CARTINDEX(am1+1,i1,j1) * size2m134
                         + INT_CARTINDEX(am2-1,i2,j2) * size34
                         + cartindex34 ];
              }

            cartindex1234++;
            }
          }
        }
      }
    }
  }

LOCAL_FUNCTION VOID
shiftam_34(I0001,am1,am2,am3,am4)
double *I0001;
int am1;
int am2;
int am3;
int am4;
{
  double *I0010;
  double *I0000;
  int i1,j1,k1,cartindex1;
  int i2,j2,k2,cartindex2;
  int i3,j3,k3,cartindex3;
  int i4,j4,k4,cartindex4;
  int cartindex1234;
  int size23p14m1,size3p14m1,size4m1,size234m1,size34m1;

  I0010 = shiftint(am1,am2,am3+1,am4-1);
  I0000 = shiftint(am1,am2,am3,am4-1);

  size23p14m1 = INT_NCART(am2)*INT_NCART(am3+1)*INT_NCART(am4-1);
  size3p14m1 = INT_NCART(am3+1)*INT_NCART(am4-1);
  size4m1 = INT_NCART(am4-1);

  size234m1 = INT_NCART(am2)*INT_NCART(am3)*INT_NCART(am4-1);
  size34m1 = INT_NCART(am3)*INT_NCART(am4-1);

  /* Loop over the target integrals. */
  cartindex1 = 0;
  cartindex1234 = 0;
  for (i1=0; i1<=am1; i1++) {
    for (k1=0; k1<=am1-i1; k1++) {
      j1 = am1 - i1 - k1;
      cartindex2 = 0;
      for (i2=0; i2<=am2; i2++) {
        for (k2=0; k2<=am2-i2; k2++) {
          j2 = am2 - i2 - k2;
          cartindex3 = 0;
          for (i3=0; i3<=am3; i3++) {
            for (k3=0; k3<=am3-i3; k3++) {
              j3 = am3 - i3 - k3;
              cartindex4 = 0;
              for (i4=0; i4<=am4; i4++) {
                for (k4=0; k4<=am4-i4; k4++) {
                  j4 = am4 - i4 - k4;

                  if (i4) {
                    I0001[cartindex1234]
                     =   I0010[  INT_CARTINDEX(am1,i1,j1) * size23p14m1
                               + INT_CARTINDEX(am2,i2,j2) * size3p14m1
                               + INT_CARTINDEX(am3+1,i3+1,j3) * size4m1
                               + INT_CARTINDEX(am4-1,i4-1,j4) ]
                       + I0000[  INT_CARTINDEX(am1,i1,j1) * size234m1
                               + INT_CARTINDEX(am2,i2,j2) * size34m1
                               + INT_CARTINDEX(am3,i3,j3) * size4m1
                               + INT_CARTINDEX(am4-1,i4-1,j4) ]
                         * CmD[0];
                    }
                  else if (j4) {
                    I0001[cartindex1234]
                     =   I0010[  INT_CARTINDEX(am1,i1,j1) * size23p14m1
                               + INT_CARTINDEX(am2,i2,j2) * size3p14m1
                               + INT_CARTINDEX(am3+1,i3,j3+1) * size4m1
                               + INT_CARTINDEX(am4-1,i4,j4-1) ]
                       + I0000[  INT_CARTINDEX(am1,i1,j1) * size234m1
                               + INT_CARTINDEX(am2,i2,j2) * size34m1
                               + INT_CARTINDEX(am3,i3,j3) * size4m1
                               + INT_CARTINDEX(am4-1,i4,j4-1) ]
                         * CmD[1];
                    }
                  else if (k4) {
                    I0001[cartindex1234]
                     =   I0010[  INT_CARTINDEX(am1,i1,j1) * size23p14m1
                               + INT_CARTINDEX(am2,i2,j2) * size3p14m1
                               + INT_CARTINDEX(am3+1,i3,j3) * size4m1
                               + INT_CARTINDEX(am4-1,i4,j4) ]
                       + I0000[  INT_CARTINDEX(am1,i1,j1) * size234m1
                               + INT_CARTINDEX(am2,i2,j2) * size34m1
                               + INT_CARTINDEX(am3,i3,j3) * size4m1
                               + INT_CARTINDEX(am4-1,i4,j4) ]
                         * CmD[2];
#if 0
   if (cartindex1234 == 4) {
     printf(" building with % f + % f * % f ",
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
                    printf("assigned I0001[%d] = % f\n",
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

LOCAL_FUNCTION VOID
iswtch(i,j)
int *i;
int *j;
{
  int tmp;

  tmp = *i;
  *i = *j;
  *j = tmp;
  }

LOCAL_FUNCTION VOID
fail()
{
  fprintf(stderr,"failing module:\n%s\n",__FILE__);
  exit(1);
  }

