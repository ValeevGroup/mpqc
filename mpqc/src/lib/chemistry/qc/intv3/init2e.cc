
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>

#include <chemistry/qc/intv2/atoms.h>
#include <chemistry/qc/intv3/flags.h>
#include <chemistry/qc/intv3/macros.h>
#include <chemistry/qc/intv3/types.h>
#include <chemistry/qc/intv3/int2e.h>
#include <chemistry/qc/intv3/utils.h>
#include <chemistry/qc/intv2/atomsallc.h>

static void
fail()
{
  fprintf(stderr,"failing module:\n%s\n",__FILE__);
  exit(1);
}

/* Initialize the 2e integral computation routines.
 * storage = the amount of storage available in bytes
 * order = order of derivative, must be zero or one
 * cs1 = center structure for center 1
 * cs2 = center structure for center 2
 * cs3 = center structure for center 3
 * cs4 = center structure for center 4
 * The integrals which will be computed are (cs1 cs2|cs3 cs4).
 * This function returns the pointer to the buffer where the
 * integrals are stored.
 */
double *
Int2eV3::int_initialize_erep(int storage, int order,
                             centers_t *cs1, centers_t *cs2,
                             centers_t *cs3, centers_t *cs4)
{
  int nc1,nc2,nc3,nc4;
  int jmax,jmax1,jmax2,jmax3,jmax4;

  redundant_ = 1;
  permute_ = 0;

  int_unit_shell = 0;

  /* Reset the variables used to get two and three center integrals. */
  int_unit2 = 0;
  int_unit4 = 0;

  /* Reset the integral storage variables. */
  int_integral_storage = 0;
  used_storage_ = 0;

  /* Turn off exponent weighted contractions. */
  int_expweight1 = 0;
  int_expweight2 = 0;
  int_expweight3 = 0;
  int_expweight4 = 0;

  /* See if the order of derivative needed is allowed. */
  if (order > 1) {
    fprintf(stderr,"int_initialize_erep cannot handle order>1, yet\n");
    }

  if (order > 0) {
      // a little bit less storage is used if int_derivative_bounds is 0
      // but for now it'll always be 1
      int_derivative_bounds = 1;
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

  // this size estimate is only accurate if all centers are the same
  int size_inter_1 = cs1->nshell * (sizeof(double*)+sizeof(int));
  if (storage - used_storage_ >= size_inter_1) {
      int_store1 = 1;
      used_storage_ += size_inter_1;
    }
  else int_store1 = 0;

  // this size estimate is only accurate if all centers are the same
  int size_inter_2 = cs1->nprim * cs1->nprim * (7*sizeof(double));
  if (storage - used_storage_ >= size_inter_2) {
      int_store2 = 1;
      used_storage_ += size_inter_2;
    }
  else int_store2 = 0;

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
      fjt_ = new FJT(jmax + 2*order); /* The 2 is for bounds checking */
    }
  else {
      fjt_ = new FJT(jmax + order);
    }

  /* Initialize the build and shift routines. */
  int_init_buildgc(order,jmax1,jmax2,jmax3,jmax4,nc1,nc2,nc3,nc4);
  int_init_shiftgc(order,jmax1,jmax2,jmax3,jmax4);

  /* Allocate storage for the integral buffer. */
  int maxsize = int_find_ncartmax(cs1)
               *int_find_ncartmax(cs2)
               *int_find_ncartmax(cs3)
               *int_find_ncartmax(cs4);
  if (order==0) {
    int_buffer = (double *) malloc(sizeof(double) * maxsize);
    int_derint_buffer = NULL;
    }
  else if (order==1) {
    int nderint;
    nderint = int_find_ncartmax_aminc(cs1,1)
             *int_find_ncartmax_aminc(cs2,1)
             *int_find_ncartmax_aminc(cs3,1)
             *int_find_ncartmax_aminc(cs4,1);
 
    /* Allocate the integral buffers. */
    int_buffer = (double *) malloc(sizeof(double) * 9*maxsize);
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

  return int_buffer;
  }

/* This is called when no more 2 electron integrals are needed.
 * It will free the intermediates. */
void
Int2eV3::int_done_erep()
{
  if (int_unit_shell) delete_int_unit_shell();
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
  int_done_buildgc();
  int_done_shiftgc();
}

/* Allocates storage for the intermediates.  The arguments are the
 * total number of unique primitive and shells. */
void
Int2eV3::alloc_inter(int nprim,int nshell)
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
  }

void
Int2eV3::compute_shell_1(centers_t *cs)
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

void
Int2eV3::compute_prim_1(centers_t *cs1)
{
}

void
Int2eV3::compute_shell_2(centers_t *cs1,centers_t *cs2)
{
  /* There are no 2 shell intermediates. */
}

/* The 2 primitive intermediates. */
void
Int2eV3::compute_prim_2(centers_t *cs1,centers_t *cs2)
{
  int offset1, offset2;
  int i1,j1,k1,i2,j2,k2;
  shell_t *shell1,*shell2;
  int i;
  /* This is 2^(1/2) * pi^(5/4) */
  const double sqrt2pi54 = 5.9149671727956129;
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
                 * exp( -   shell1->exp[k1] * shell2->exp[k2]
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
