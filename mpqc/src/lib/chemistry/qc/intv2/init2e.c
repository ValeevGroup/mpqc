
/* $Log$
 * Revision 1.4  1994/08/26 22:45:31  etseidl
 * fix a bunch of warnings, get rid of rcs id's, get rid of bread/bwrite and
 * fread/fwrite modules
 *
 * Revision 1.3  1994/05/27  23:51:21  cljanss
 * Added support for 2 and 3 center 2 electron integrals.  Added a test porgram.
 *
 * Revision 1.2  1993/12/30  13:32:50  etseidl
 * mostly rcs id stuff
 *
 * Revision 1.5  1992/06/17  22:04:43  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.4  1992/05/26  20:25:17  jannsen
 * make derivative bounds checking optional
 * add code to allow bound intermediates computable in shell blocks
 *
 * Revision 1.3  1992/05/13  18:29:38  jannsen
 * added bounds checking for derivative integrals
 *
 * Revision 1.2  1992/03/31  01:22:09  jannsen
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.4  1992/01/30  01:30:08  cljanss
 * 1. Correct freeing of data.
 * 2. added int_find_jmax_for_con
 *
 * Revision 1.3  1992/01/10  18:00:42  cljanss
 * set int_maxsize to 0 when done, so storage.c knows that an int_inititalize
 * is needed
 *
 * Revision 1.2  1992/01/10  17:57:28  cljanss
 * setup int_integral_storage and int_maxsize globals
 *
 * Revision 1.1  1991/12/14  00:19:36  cljanss
 * Initial revision
 * */
/*
 * Revision 1.7  91/11/22  17:49:22  cljanss
 * bound matrix generated is now handled by a separate flag
 * 
 * Revision 1.6  91/09/28  19:26:49  cljanss
 * Switch to new file naming scheme
 * 
 * Revision 1.5  91/09/28  18:17:16  cljanss
 * Intermediates are no longer stored, if requested with flags.
 * 
 * Revision 1.4  91/08/09  16:53:42  cljanss
 * more work on derivative integrals (int_buffer stores der ints for 3 centers)
 * 
 * Revision 1.3  1991/06/19  23:18:19  janssen
 * added computation and checking of upper bounds to shell quartet integrals
 *
 * Revision 1.2  1991/06/19  15:19:38  janssen
 * First stab at two electron derivative integrals--might work might not
 *
 * Revision 1.1  1991/06/16  16:40:07  janssen
 * Initial revision
 * */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>

#include "atoms.h"
#include "int_flags.h"
#include "int_macros.h"
#include "int_types.h"

#define ALLOC_INTERMEDIATES
#include "inter.h"

#include "init2e.gbl"
#include "init2e.lcl"

#include "buildgc.gbl"
#include "shiftgc.gbl"
#include "int_fjt.gbl"
#include "comp_erep.gbl"
#include "utils.gbl"

#include "atomsallc.h"

/* The NCUBE exp function cannot handle large negative arguments. */
#ifndef NCUBE
#define exp_cutoff exp
#else
LOCAL_FUNCTION double
exp_cutoff(exponent)
double exponent;
{
  double r;
  if (exponent < -600.0) r = 0.0;
  else r = exp(exponent);
  return r;
  }
#endif

/* Initialize the 2e integral computation routines.
 * flags = a set of or'ed flags, must be INT_EREP for now
 * order = order of derivative, must be zero or one
 * cs1 = center structure for center 1
 * cs2 = center structure for center 2
 * cs3 = center structure for center 3
 * cs4 = center structure for center 4
 * The integrals which will be computed are (cs1 cs2|cs3 cs4).
 * This function returns the pointer to the buffer where the
 * integrals are stored.
 */
GLOBAL_FUNCTION double *
int_initialize_erep(flags,order,cs1,cs2,cs3,cs4)
int flags;
int order;
centers_t *cs1;
centers_t *cs2;
centers_t *cs3;
centers_t *cs4;
{
  int nc1,nc2,nc3,nc4;
  int jmax,jmax1,jmax2,jmax3,jmax4;

  /* Reset the variables used to get two and three center integrals. */
  int_unit2 = 0;
  int_unit4 = 0;

  /* Reset the integral storage variable. */
  int_integral_storage = 0;

  /* Turn off exponent weighted contractions. */
  int_expweight1 = 0;
  int_expweight2 = 0;
  int_expweight3 = 0;
  int_expweight4 = 0;

  /* See if flags is legit. */
  if (!(flags&INT_EREP)) {
    fprintf(stderr,"int_initialize_erep cannot handle nonINT_EREP flags yet\n");
    fail();
    }

  /* See if the order of derivative needed is allowed. */
  if (order > 1) {
    fprintf(stderr,"int_initialize_erep cannot handle order>1, yet\n");
    }

  if (order > 0) {
    if (flags&INT_NODERB) {
      int_derivative_bounds = 0;
      }
    else {
      int_derivative_bounds = 1;
      }
    }

  /* A noncritical limitation for now. */
  if ((cs1 != cs2) || (cs2 != cs3) || (cs3 != cs4)) {
    fprintf(stderr,"libint: because the int_compute_erep routine\n");
    fprintf(stderr,"might permute centers around, different centers\n");
    fprintf(stderr,"cannot be given (but this can be easily fixed)\n");
    fail();
    }

  /* Put the center pointers into the global centers pointers. */
  int_cs1 = cs1;
  int_cs2 = cs2;
  int_cs3 = cs3;
  int_cs4 = cs4;

  /* See if the intermediates are to be computed and set global variables
   * accordingly. */
  int_store1 = !(flags&INT_NOSTR1);
  int_store2 = !(flags&INT_NOSTR2);
  int_store_bounds = !(flags&INT_NOSTRB);

  /* Allocate storage for the intermediates. */
  alloc_inter(cs4->prim_offset + cs4->nprim, cs4->shell_offset + cs4->nshell);

  /* Set up the one shell intermediates, block by block. */
  if (int_store1) {
    compute_shell_1(cs1);
    if (cs2 != cs1) compute_shell_1(cs2);
    if (cs3 != cs2 && cs3 != cs1) compute_shell_1(cs3);
    if (cs4 != cs3 && cs4 != cs2 && cs4 != cs1) compute_shell_1(cs4);
  
    /* Set up the one primitive intermediates, block by block. */
    compute_prim_1(cs1);
    if (cs2 != cs1) compute_prim_1(cs2);
    if (cs3 != cs2 && cs3 != cs1) compute_prim_1(cs3);
    if (cs4 != cs3 && cs4 != cs2 && cs4 != cs1) compute_prim_1(cs4);
    }

  /* Compute the two shell intermediates, block by block. */
  if (int_store2) {
    compute_shell_2(cs1,cs1);
    if (cs2 != cs1) {
      compute_shell_2(cs1,cs2);
      compute_shell_2(cs2,cs1);
      compute_shell_2(cs2,cs2);
      }
    if (cs3 != cs2 && cs3 != cs1) {
      compute_shell_2(cs1,cs3);
      compute_shell_2(cs3,cs1);
      compute_shell_2(cs2,cs3);
      compute_shell_2(cs3,cs2);
      compute_shell_2(cs3,cs3);
      }
    if (cs4 != cs3 && cs4 != cs2 && cs4 != cs1) {
      compute_shell_2(cs1,cs4);
      compute_shell_2(cs4,cs1);
      compute_shell_2(cs2,cs4);
      compute_shell_2(cs4,cs2);
      compute_shell_2(cs3,cs4);
      compute_shell_2(cs4,cs3);
      compute_shell_2(cs4,cs4);
      }
  
    /* Compute the two primitive intermediates, block by block. */
    compute_prim_2(cs1,cs1);
    if (cs2 != cs1) {
      compute_prim_2(cs1,cs2);
      compute_prim_2(cs2,cs1);
      compute_prim_2(cs2,cs2);
      }
    if (cs3 != cs2 && cs3 != cs1) {
      compute_prim_2(cs1,cs3);
      compute_prim_2(cs3,cs1);
      compute_prim_2(cs2,cs3);
      compute_prim_2(cs3,cs2);
      compute_prim_2(cs3,cs3);
      }
    if (cs4 != cs3 && cs4 != cs2 && cs4 != cs1) {
      compute_prim_2(cs1,cs4);
      compute_prim_2(cs4,cs1);
      compute_prim_2(cs2,cs4);
      compute_prim_2(cs4,cs2);
      compute_prim_2(cs3,cs4);
      compute_prim_2(cs4,cs3);
      compute_prim_2(cs4,cs4);
      }
    }

  /* Find the max angular momentum on each center. */
  jmax1 = int_find_jmax(cs1);
  jmax2 = int_find_jmax(cs2);
  jmax3 = int_find_jmax(cs3);
  jmax4 = int_find_jmax(cs4);

  /* Find the maximum number of contractions in a shell on each center. */
  nc1 = int_find_nconmax(cs1);
  nc2 = int_find_nconmax(cs2);
  nc3 = int_find_nconmax(cs3);
  nc4 = int_find_nconmax(cs4);

  /* Initialize the Fj(T) routine. */
  jmax = jmax1+jmax2+jmax3+jmax4;
  if (int_derivative_bounds) {
    int_initialize_fjt(jmax + 2 * order); /* The 2 is for bounds checking */
    }
  else {
    int_initialize_fjt(jmax + order);
    }

  /* Initialize the build and shift routines. */
  int_init_buildgc(order,jmax1,jmax2,jmax3,jmax4,nc1,nc2,nc3,nc4);
  int_init_shiftgc(order,jmax1,jmax2,jmax3,jmax4);

  /* Allocate storage for the integral buffer. */
  int_maxsize = int_find_nfuncmax(cs1)
               *int_find_nfuncmax(cs2)
               *int_find_nfuncmax(cs3)
               *int_find_nfuncmax(cs4);
  if (order==0) {
    int_buffer = (double *) malloc(sizeof(double) * int_maxsize);
    int_derint_buffer = NULL;
    }
  else if (order==1) {
    int nderint;
    nderint = int_find_nfuncmax_aminc(cs1,1)
             *int_find_nfuncmax_aminc(cs2,1)
             *int_find_nfuncmax_aminc(cs3,1)
             *int_find_nfuncmax_aminc(cs4,1);
 
    /* Allocate the integral buffers. */
    int_buffer = (double *) malloc(sizeof(double) * 9*int_maxsize);
    int_derint_buffer = (double *) malloc(sizeof(double) * nderint);
    if (!int_derint_buffer) {
      fprintf(stderr,"couldn't malloc intermed storage for derivative ints\n");
      fail();
      }
    }

  if (!int_buffer) {
    fprintf(stderr,"couldn't allocate integrals\n");
    fail();
    }

  /* Compute the Q(M,N) which are used to compute the maximum
   * value of a integrals in shell quartets. */
  if (int_store_bounds) {
    compute_Q(cs1,cs1);
    if (cs2 != cs1) {
      compute_Q(cs1,cs2);
      compute_Q(cs2,cs1);
      compute_Q(cs2,cs2);
      }
    if (cs3 != cs2 && cs3 != cs1) {
      compute_Q(cs1,cs3);
      compute_Q(cs3,cs1);
      compute_Q(cs2,cs3);
      compute_Q(cs3,cs2);
      compute_Q(cs3,cs3);
      }
    if (cs4 != cs3 && cs4 != cs2 && cs4 != cs1) {
      compute_Q(cs1,cs4);
      compute_Q(cs4,cs1);
      compute_Q(cs2,cs4);
      compute_Q(cs4,cs2);
      compute_Q(cs3,cs4);
      compute_Q(cs4,cs3);
      compute_Q(cs4,cs4);
      }
    }

  return int_buffer;
  }

/* This is called when no more 2 electron integrals are needed.
 * It will free the intermediates. */
GLOBAL_FUNCTION VOID
int_done_erep()
{
  if (int_derint_buffer) free(int_derint_buffer);
  free(int_buffer);
  if (int_store1) {
    free_doublep_vector(&int_shell_r);
    free_int_vector(&int_shell_to_prim);
    }
  if (int_store2) {
    free_double_matrix(&int_prim_zeta);
    free_double_matrix(&int_prim_oo2zeta);
    free_double_array3(&int_prim_p);
    free_double_matrix(&int_prim_k);
    }
  if (int_store_bounds) {
    free_double_matrix(&int_shell_Q);
    }
  int_maxsize = 0;
  int_done_buildgc();
  int_done_shiftgc();
  int_done_fjt();
  }

/* Allocates storage for the intermediates.  The arguments are the
 * total number of unique primitive and shells. */
LOCAL_FUNCTION VOID
alloc_inter(nprim,nshell)
int nprim;
int nshell;
{
  if (int_store1) {
    if (  allocbn_doublep_vector(&int_shell_r, "n", nshell)
        ||allocbn_int_vector(&int_shell_to_prim, "n", nshell)
        ) {
      fprintf(stderr,"problem allocating O(n) integral intermediates for");
      fprintf(stderr," %d shells and %d primitives\n",nshell,nprim);
      fail();
      }
    }
  if (int_store2) {
    if (  allocbn_double_matrix(&int_prim_zeta, "n1 n2", nprim, nprim)
        ||allocbn_double_matrix(&int_prim_oo2zeta, "n1 n2", nprim, nprim)
        ||allocbn_double_array3(&int_prim_p, "n1 n2 n3", nprim, nprim, 3)
        ||allocbn_double_matrix(&int_prim_k, "n1 n2", nprim, nprim)
        ) {
      fprintf(stderr,"problem allocating O(n^2) integral intermediates for");
      fprintf(stderr," %d shells and %d primitives\n",nshell,nprim);
      fail();
      }
    }
  if (int_store_bounds) {
    if (allocbn_double_matrix(&int_shell_Q, "n1 n2", nshell, nshell)) {
      fprintf(stderr,"problem allocating O(n^2) integral bounds array");
      fprintf(stderr," %d shells\n",nshell);
      fail();
      }
    }
  }

LOCAL_FUNCTION VOID
compute_shell_1(cs)
centers_t *cs;
{
  int i,j;
  int offset;
  int iprim;

  offset = cs->shell_offset;
  iprim = cs->prim_offset;
  for (i=0; i<cs->n; i++) {
    for (j=0; j<cs->center[i].basis.n; j++) {

      /* The offset shell geometry vectors. */
      int_shell_r.dp[offset] = cs->center[i].r;

      /* The number of the first offset primitive in a offset shell. */
      int_shell_to_prim.i[offset] = iprim;

      offset++;
      iprim += cs->center[i].basis.shell[j].nprim;
      }
    }
  }

LOCAL_FUNCTION VOID
compute_prim_1(cs1)
centers_t *cs1;
{
#if 0 /* There are no 1 primitive intermediates. */
  int offset1;
  int i1,j1,k1;
  shell_t *shell1;

  offset1 = cs1->prim_offset;
  for (i1=0; i1<cs1->n; i1++) {
    for (j1=0; j1<cs1->center[i1].basis.n; j1++) {
      shell1 = &cs1->center[i1].basis.shell[j1];
      for (k1=0; k1<cs1->center[i1].basis.shell[j1].nprim; k1++) {

        /* The primitive normalization constants. */
        /* This is the part of the normalization constant
         * which doesn't depend on the exponents of
         * x, y, and z for each basis function. */
        int_prim_norm.d[offset1] =
            0.712705470354990 /* (2/pi)^(3/4) */
          * pow(4.0 * shell1->exp[k1], 0.5 * shell1->type.am)
          * pow(shell1->exp[k1], 0.75);

        offset1++;
        }
      }
    }
#endif
  }

LOCAL_FUNCTION VOID
compute_shell_2(cs1,cs2)
centers_t *cs1;
centers_t *cs2;
{
  /* There are no 2 shell intermediates. */
  }

/* The 2 primitive intermediates. */
LOCAL_FUNCTION VOID
compute_prim_2(cs1,cs2)
centers_t *cs1;
centers_t *cs2;
{
  int offset1, offset2;
  int i1,j1,k1,i2,j2,k2;
  shell_t *shell1,*shell2;
  int i;
  /* This is 2^(1/2) * pi^(5/4) */
  CONST double sqrt2pi54 = 5.9149671727956129;
  double AmB,AmB2;

  offset1 = cs1->prim_offset;
  for (i1=0; i1<cs1->n; i1++) {
    for (j1=0; j1<cs1->center[i1].basis.n; j1++) {
      shell1 = &cs1->center[i1].basis.shell[j1];
      for (k1=0; k1<cs1->center[i1].basis.shell[j1].nprim; k1++) {
        offset2 = cs2->prim_offset;
        for (i2=0; i2<cs2->n; i2++) {
          for (j2=0; j2<cs2->center[i2].basis.n; j2++) {
            shell2 = &cs2->center[i2].basis.shell[j2];
            for (k2=0; k2<cs2->center[i2].basis.shell[j2].nprim; k2++) {

              /* The zeta = alpha + beta intermediate. */
              int_prim_zeta.d[offset1][offset2] =
                shell1->exp[k1] + shell2->exp[k2];

              /* The 1/(2 zeta) intermediate times 2.0. */
              int_prim_oo2zeta.d[offset1][offset2] =
                1.0/int_prim_zeta.d[offset1][offset2];

              /* The p = (alpha A + beta B) / zeta */
              for (i=0; i<3; i++) {
                int_prim_p.d[offset1][offset2][i] =
                  (  shell1->exp[k1] * cs1->center[i1].r[i]
                   + shell2->exp[k2] * cs2->center[i2].r[i])
                  *  int_prim_oo2zeta.d[offset1][offset2];
                }

              /* Compute AmB^2 */
              AmB2 = 0.0;
              for (i=0; i<3; i++) {
                AmB = cs2->center[i2].r[i] - cs1->center[i1].r[i];
                AmB2 += AmB*AmB;
                }

              /* Compute the K intermediate. */
              int_prim_k.d[offset1][offset2] =
                   sqrt2pi54
                 * int_prim_oo2zeta.d[offset1][offset2]
                 * exp_cutoff( -   shell1->exp[k1] * shell2->exp[k2]
                          * int_prim_oo2zeta.d[offset1][offset2]
                          * AmB2 );

              /* Finish the 1/(2 zeta) intermediate. */
              int_prim_oo2zeta.d[offset1][offset2] =
                0.5 * int_prim_oo2zeta.d[offset1][offset2];

              offset2++;
              }
            }
          }
        offset1++;
        }
      }
    }
  }

LOCAL_FUNCTION VOID
fail()
{
  fprintf(stderr,"failing module:\n%s\n",__FILE__);
  exit(1);
  }

LOCAL_FUNCTION VOID
compute_Q(cs1,cs2)
centers_t *cs1;
centers_t *cs2;
{
  int i1,j1;
  int i2,j2;
  int offset1,offset2;
  int shell1,shell2;
  int sh1,sh2,sh3,sh4;
  int nfunc1,nfunc2;
  double max;
  double integral;
  int i,j;

  offset1 = cs1->shell_offset;
  shell1 = 0;
  for (i1=0; i1<cs1->n; i1++) {
    for (j1=0; j1<cs1->center[i1].basis.n; j1++) {
      nfunc1 = cs1->center[i1].basis.shell[j1].nfunc;

      offset2 = cs2->shell_offset;
      shell2 = 0;
      for (i2=0; i2<cs2->n; i2++) {
        for (j2=0; j2<cs2->center[i2].basis.n; j2++) {
          nfunc2 = cs2->center[i2].basis.shell[j2].nfunc;

          sh1 = shell1;
          sh2 = shell2;
          sh3 = shell1;
          sh4 = shell2;

          int_erep(INT_EREP|INT_REDUND|INT_NOPERM|INT_NOBCHK,
                   &sh1,&sh2,&sh3,&sh4);

          /* Find the biggest (ij|ij) integral. */
          max = 0.0;
          for (i=0; i<nfunc1; i++) {
            for (j=0; j<nfunc2; j++) {
              int index = i * nfunc2 * nfunc1 * nfunc2
                        + j * nfunc1 * nfunc2
                        + i * nfunc2
                        + j;
              integral = int_buffer[ index ];
              if (fabs(integral) > max) max = fabs(integral);
              }
            }

          /* Compute the Q value. */
          int_shell_Q.d[offset1][offset2] = sqrt(max);

          offset2++;
          shell2++;
          }
        }

      offset1++;
      shell1++;
      }
    }
  }
