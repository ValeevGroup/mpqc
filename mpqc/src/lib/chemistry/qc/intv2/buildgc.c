
/* $Log$
 * Revision 1.1  1993/12/29 12:53:01  etseidl
 * Initial revision
 *
 * Revision 1.6  1992/06/17  22:04:29  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.5  1992/05/26  20:25:14  jannsen
 * make derivative bounds checking optional
 * add code to allow bound intermediates computable in shell blocks
 *
 * Revision 1.4  1992/05/19  20:53:15  seidl
 * in function build_using_gcs, take as much as possible out of inner loops.
 * small performance improvement
 *
 * Revision 1.3  1992/05/13  18:29:28  jannsen
 * added bounds checking for derivative integrals
 *
 * Revision 1.2  1992/03/31  01:21:30  jannsen
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.5  1992/03/10  22:04:12  cljanss
 * First pass at NCUBE V3.0 conversion and a bug fix
 *
 * Revision 1.4  1992/01/30  01:26:00  cljanss
 * 1. The build_using_non_gcs routine is used for shells without gc's
 * 2. Free memory when done.
 * 3. Don't call init_inthave if it isn't needed.
 *
 * Revision 1.3  1992/01/10  17:56:18  cljanss
 * converted to NCUBE
 *
 * Revision 1.2  1992/01/08  22:39:14  cljanss
 * added routines to deal with cases where no generalized cases are used, but
 * there isn't much performance enhancement
 *
 * Revision 1.1  1991/12/14  00:19:36  cljanss
 * Initial revision
 * */
/* Copied from build.c on 11-22-91. */
static char *rcsid = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>

#include "atoms.h"
#include "inter.h"
#include "fjttable.h"
#include "int_macros.h"

#define ALLOC_BUILDINTER
#include "buildinter.h"

#include "int_fjt.gbl"
#include "int_print.gbl"

/* these statics are needed by add_store and free_store */
#define STORAGE_CHUNK 4096
struct store_list {
  void* data[STORAGE_CHUNK];
  struct store_list* p;
  };
typedef struct store_list store_list_t;
static int n_store_last;
static store_list_t* store=NULL;

#include "buildgc.gbl"
#include "buildgc.lcl"

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

  /* Offset shell numbers. */
static int osh1, osh2, osh3, osh4;
  /* Offset primitive numbers. */
static int opr1, opr2, opr3, opr4;

  /* Boolean array which gives whether or not an array is computed. */
static int_array3_t inthave;

  /* Saved initialization parameters used to free data. */
static int saved_am12,saved_am34,saved_ncon;

  /* MG is the maximum angular momentum for which we will use
   * the generated build routines. */
#include "MG.h"
#define MINA(x) (((x)<MG)?(x):MG)
typedef int (*intfunc)();
static intfunc build_routine[4][4][4][4][2];
#define DECLARE_BUILD(i,j,k,l) \
  extern int int2v ## i ## j ## k ## l (),\
             int2v ## i ## j ## k ## l ## eAB ();
#if defined(NCUBE_V2)
/* Explicit declarations are required on the NCUBE since ansi cpp isn't there */
extern int int2v0100 (), int2v0100eAB (); 
extern int int2v0101 (), int2v0101eAB (); 
extern int int2v0111 (), int2v0111eAB (); 
extern int int2v0200 (), int2v0200eAB (); 
extern int int2v0201 (), int2v0201eAB (); 
extern int int2v0202 (), int2v0202eAB (); 
extern int int2v0211 (), int2v0211eAB (); 
extern int int2v0212 (), int2v0212eAB (); 
extern int int2v0222 (), int2v0222eAB (); 
extern int int2v1100 (), int2v1100eAB (); 
extern int int2v1111 (), int2v1111eAB (); 
extern int int2v1200 (), int2v1200eAB (); 
extern int int2v1201 (), int2v1201eAB (); 
extern int int2v1211 (), int2v1211eAB (); 
extern int int2v1212 (), int2v1212eAB (); 
extern int int2v1222 (), int2v1222eAB (); 
extern int int2v2200 (), int2v2200eAB (); 
extern int int2v2201 (), int2v2201eAB (); 
extern int int2v2211 (), int2v2211eAB (); 
extern int int2v2222 (), int2v2222eAB (); 
#else
#if (MG == 1) || (MG == 2) || (MG == 3) || (MG == 4)
DECLARE_BUILD(0,1,0,0)
DECLARE_BUILD(0,1,0,1)
DECLARE_BUILD(0,1,1,1)
DECLARE_BUILD(1,1,0,0)
DECLARE_BUILD(1,1,1,1)
#endif

#if (MG == 2) || (MG == 3) || (MG == 4)
DECLARE_BUILD(0,2,0,0)
DECLARE_BUILD(0,2,0,1)
DECLARE_BUILD(0,2,0,2)
DECLARE_BUILD(0,2,1,1)
DECLARE_BUILD(0,2,1,2)
DECLARE_BUILD(0,2,2,2)
DECLARE_BUILD(1,2,0,0)
DECLARE_BUILD(1,2,0,1)
DECLARE_BUILD(1,2,1,1)
DECLARE_BUILD(1,2,1,2)
DECLARE_BUILD(1,2,2,2)
DECLARE_BUILD(2,2,0,0)
DECLARE_BUILD(2,2,0,1)
DECLARE_BUILD(2,2,1,1)
DECLARE_BUILD(2,2,2,2)
#endif

#if (MG == 3) || (MG == 4)
DECLARE_BUILD(0,3,0,0)
DECLARE_BUILD(0,3,0,1)
DECLARE_BUILD(0,3,0,2)
DECLARE_BUILD(0,3,0,3)
DECLARE_BUILD(0,3,1,1)
DECLARE_BUILD(0,3,1,2)
DECLARE_BUILD(0,3,1,3)
DECLARE_BUILD(0,3,2,2)
DECLARE_BUILD(0,3,2,3)
DECLARE_BUILD(0,3,3,3)
DECLARE_BUILD(1,3,0,0)
DECLARE_BUILD(1,3,0,1)
DECLARE_BUILD(1,3,0,2)
DECLARE_BUILD(1,3,1,1)
DECLARE_BUILD(1,3,1,2)
DECLARE_BUILD(1,3,1,3)
DECLARE_BUILD(1,3,2,2)
DECLARE_BUILD(1,3,2,3)
DECLARE_BUILD(1,3,3,3)
DECLARE_BUILD(2,3,0,0)
DECLARE_BUILD(2,3,0,1)
DECLARE_BUILD(2,3,0,2)
DECLARE_BUILD(2,3,1,1)
DECLARE_BUILD(2,3,1,2)
DECLARE_BUILD(2,3,2,2)
DECLARE_BUILD(2,3,2,3)
DECLARE_BUILD(2,3,3,3)
DECLARE_BUILD(3,3,0,0)
DECLARE_BUILD(3,3,0,1)
DECLARE_BUILD(3,3,0,2)
DECLARE_BUILD(3,3,1,1)
DECLARE_BUILD(3,3,1,2)
DECLARE_BUILD(3,3,2,2)
DECLARE_BUILD(3,3,3,3)
#endif
#endif


/* This initializes the build routines.  It is called from
 * int_initialize_erep.  This allocates storage for the
 * intermediate integrals. */
GLOBAL_FUNCTION VOID
int_init_buildgc(order,am1,am2,am3,am4,nc1,nc2,nc3,nc4)
int order;
int am1;
int am2;
int am3;
int am4;
int nc1;
int nc2;
int nc3;
int nc4;
{
  int *jmax_for_con;
  int am12;
  int am34;
  int am;
  int i,j,k,l,m;
  int ci,cj,ck,cl;

  /* Convert the am1-4 to their canonical ordering. */
  if (am2>am1) {
    iswtch(&am1,&am2);
    iswtch(&nc1,&nc2);
    }
  if (am4>am3) {
    iswtch(&am3,&am4);
    iswtch(&nc3,&nc4);
    }
  if ((am3 > am1)||((am3 == am1)&&(am4 > am2))) {
    iswtch(&am1,&am3);
    iswtch(&nc1,&nc3);
    iswtch(&am2,&am4);
    iswtch(&nc2,&nc4);
    }

  /* If the center permutation 1<->3 and 2<->4 is performed, then
   * we may need the am for center 2 to be as big as for center 4. */
  if (am4 > am2) am2 = am4;

  /* As far as this routine knows the biggest nc can end up anywhere. */
  if (nc2>nc1) nc1 = nc2;
  if (nc3>nc1) nc1 = nc3;
  if (nc4>nc1) nc1 = nc4;
  nc2 = nc3 = nc4 = nc1;

  jmax_for_con = (int *) malloc(sizeof(int)*nc1);
  for (i=0; i<nc1; i++) {
    int tmp;
    jmax_for_con[i] = int_find_jmax_for_con(int_cs1,i);
    if (  (int_cs2 != int_cs1)
        &&((tmp=int_find_jmax_for_con(int_cs2,i))>jmax_for_con[i]))
      jmax_for_con[i] = tmp;
    if (  (int_cs3 != int_cs1) && (int_cs3 != int_cs2)
        &&((tmp=int_find_jmax_for_con(int_cs3,i))>jmax_for_con[i]))
      jmax_for_con[i] = tmp;
    if (  (int_cs4 != int_cs1) && (int_cs4 != int_cs2) && (int_cs4 != int_cs3)
        &&((tmp=int_find_jmax_for_con(int_cs4,i))>jmax_for_con[i]))
      jmax_for_con[i] = tmp;
    }

  /* If derivatives are needed, then am1 can be bigger. */
  if (order==1) am1++;
  /* To compute derivative integral bounds, am3 can be bigger also. */
  if (order==1 && int_derivative_bounds) am3++;

  am12 = am1 + am2;
  am34 = am3 + am4;
  am = am12 + am34;

  /* Allocate the intlist. */
  if (  allocbn_int_array3(&inthave,"n1 n2 n3",am12+1,am34+1,am+1)
      ||allocbn_doublep_array3(&int_v_list,"n1 n2 n3",am12+1,am34+1,am+1)) {
    fprintf(stderr,"problem allocating integral intermediates for");
    fprintf(stderr," am12 = %d, am34 = %d, and am = %d \n",am12,am34,am);
    fail();
    }

  /* Set all slots to NULL */
  for (i=0; i<=am12; i++) {
    for (j=0; j<=am34; j++) {
      for (k=0; k<=am12+am34; k++) {
        int_v_list.dp[i][j][k] = NULL;
        }
      }
    }

  /* Allocate storage for the needed slots. */
  for (i=0; i<=am12; i++) {
    for (j=0; j<=am34; j++) {
      for (k=0; k<=am12+am34-i-j; k++) {
        int_v_list.dp[i][j][k] =
            (double *) malloc(sizeof(double)*INT_NCART(i)*INT_NCART(j));
        if (!int_v_list.dp[i][j][k]) {
          fprintf(stderr,"couldn't allocate all integral intermediates\n");
          fail();
          }
        }
      }
    }


  /* Allocate storage for the contracted integrals (these are the output
   * of the build routines). */
  /* The ci, etc, indices refer to which of the pair of contraction
   * coefficients we are using. */
  int_con_ints_array =
    (doublep_array4_t ****)malloc(sizeof(doublep_array4_t ***)*nc1);
  for (ci=0; ci<nc1; ci++) {
    int_con_ints_array[ci] =
      (doublep_array4_t ***)malloc(sizeof(doublep_array4_t **)*nc2);
    for (cj=0; cj<nc2; cj++) {
      int_con_ints_array[ci][cj] =
        (doublep_array4_t **)malloc(sizeof(doublep_array4_t *)*nc3);
      for (ck=0; ck<nc3; ck++) {
        int_con_ints_array[ci][cj][ck] =
          (doublep_array4_t *)malloc(sizeof(doublep_array4_t )*nc4);
        for (cl=0; cl<nc4; cl++) {
  int am12_for_con;
  int am34_for_con;

  am12_for_con = jmax_for_con[ci] + jmax_for_con[cj] + order;
  if ((jmax_for_con[ck]!=am3)||(jmax_for_con[cl]!=am4)) {
    am34_for_con = jmax_for_con[ck] + jmax_for_con[cl] + order;
    }
  else {
    am34_for_con = jmax_for_con[ck] + jmax_for_con[cl];
    }

  /* Find the biggest am for this contraction. */

  if (allocbn_doublep_array4(&int_con_ints_array[ci][cj][ck][cl],"n1 n2 n3 n4",
                             am12+1,am12+1,am34+1,am34+1)) {
    fprintf(stderr,"couldn't allocate contracted integral array\n");
    fail();
    }
  /* Allocate storage for the integrals which will be used by the shift
   * routine. */
  for (i=0; i<=am12; i++) {
    for (j=0; j<=am12; j++) {
      for (k=0; k<=am34; k++) {
        for (l=0; l<=am34; l++) {
          int_con_ints_array[ci][cj][ck][cl].dp[i][j][k][l] = NULL;
          }
        }
      }
    }
  for (i=0; i<=am12_for_con; i++) {
    for (j=0; j<=am12_for_con-i; j++) {
      for (k=0; k<=am34_for_con; k++) {
        for (l=0; l<=am34_for_con-k; l++) {
/* If there are Pople style s=p shells and the shells are ordered
 * first s and then p and there are no p or d shells on the molecule,
 * then this algorithm would will allocate a little more storage
 * than needed.  General contraction should be ordered high to
 * low angular momentum for this reason. */
#undef NO_SHARED_INT_INTER
#ifndef NO_SHARED_INT_INTER
                /* Share storage for certain cases. */
          if ((j==0)&&(l==0)) {
            int_con_ints_array[ci][cj][ck][cl].dp[i][j][k][l] =
                (double *) malloc(  sizeof(double)
                                                 * INT_NCART(i)
                                                 * INT_NCART(j)
                                                 * INT_NCART(k)
                                                 * INT_NCART(l));
             add_store(int_con_ints_array[ci][cj][ck][cl].dp[i][j][k][l]);
#if 0
             printf("mallocing %d %d %d %d %d %d %d %d\n",
                    ci,cj,ck,cl,i,j,k,l);
#endif
             }
           else if ((ci==0)&&(cj==0)&&(ck==0)&&(cl==0)) {
            int_con_ints_array[ci][cj][ck][cl].dp[i][j][k][l] =
              (double *) malloc(  sizeof(double)
                                                 * INT_NCART(i)
                                                 * INT_NCART(j)
                                                 * INT_NCART(k)
                                                 * INT_NCART(l));
             add_store(int_con_ints_array[ci][cj][ck][cl].dp[i][j][k][l]);
#if 0
             printf("mallocing %d %d %d %d %d %d %d %d\n",
                    ci,cj,ck,cl,i,j,k,l);
#endif
             }
           else if (int_con_ints_array[0][0][0][0].dp[i][j][k][l]) {
             int_con_ints_array[ci][cj][ck][cl].dp[i][j][k][l] =
               int_con_ints_array[0][0][0][0].dp[i][j][k][l];
#if 0
             printf("copying %d %d %d %d %d %d %d %d\n",
                    ci,cj,ck,cl,i,j,k,l);
#endif
             }
           else {
            int_con_ints_array[ci][cj][ck][cl].dp[i][j][k][l] =
              (double *) malloc(  sizeof(double)
                                                 * INT_NCART(i)
                                                 * INT_NCART(j)
                                                 * INT_NCART(k)
                                                 * INT_NCART(l));
             add_store(int_con_ints_array[ci][cj][ck][cl].dp[i][j][k][l]);
#if 0
             printf("mallocing %d %d %d %d %d %d %d %d\n",
                    ci,cj,ck,cl,i,j,k,l);
#endif
             }
#else
          int_con_ints_array[ci][cj][ck][cl].dp[i][j][k][l] =
                (double *) malloc(  sizeof(double)
                                                 * INT_NCART(i)
                                                 * INT_NCART(j)
                                                 * INT_NCART(k)
                                                 * INT_NCART(l));
          add_store(int_con_ints_array[ci][cj][ck][cl].dp[i][j][k][l]);
#endif
          if (!int_con_ints_array[ci][cj][ck][cl].dp[i][j][k][l]) {
            fprintf(stderr,"couldn't allocate contracted integral storage\n");
            fail();
            }
          }
        }
      }
    }
          }
        }
      }
    }

  /* This pointer is needed if int_buildam (without gc) is called. */
  int_con_ints = &int_con_ints_array[0][0][0][0];

  /* Initialize the build_routine array. */
  for (i=0; i<4; i++) {
    for (j=0; j<4; j++) {
      for (k=0; k<4; k++) {
        for (l=0; l<4; l++) {
          for (m=0; m<2; m++) {
    build_routine[i][j][k][l][m] = impossible_integral;
            }
          }
        }
      }
    }

#define ASSIGN_BUILD(i,j,k,l) \
  build_routine[i][j][k][l][0]= int2v ## i ## j ## k ## l ;\
  build_routine[i][j][k][l][1]= int2v ## i ## j ## k ## l ## eAB;

#if defined(NCUBE_V2)
  build_routine[0][1][0][0][0]= int2v0100 ; build_routine[0][1][0][0][1]= int2v0100eAB; 
  build_routine[0][1][0][1][0]= int2v0101 ; build_routine[0][1][0][1][1]= int2v0101eAB; 
  build_routine[0][1][1][1][0]= int2v0111 ; build_routine[0][1][1][1][1]= int2v0111eAB; 
  build_routine[0][2][0][0][0]= int2v0200 ; build_routine[0][2][0][0][1]= int2v0200eAB; 
  build_routine[0][2][0][1][0]= int2v0201 ; build_routine[0][2][0][1][1]= int2v0201eAB; 
  build_routine[0][2][0][2][0]= int2v0202 ; build_routine[0][2][0][2][1]= int2v0202eAB; 
  build_routine[0][2][1][1][0]= int2v0211 ; build_routine[0][2][1][1][1]= int2v0211eAB; 
  build_routine[0][2][1][2][0]= int2v0212 ; build_routine[0][2][1][2][1]= int2v0212eAB; 
  build_routine[0][2][2][2][0]= int2v0222 ; build_routine[0][2][2][2][1]= int2v0222eAB; 
  build_routine[1][1][0][0][0]= int2v1100 ; build_routine[1][1][0][0][1]= int2v1100eAB; 
  build_routine[1][1][1][1][0]= int2v1111 ; build_routine[1][1][1][1][1]= int2v1111eAB; 
  build_routine[1][2][0][0][0]= int2v1200 ; build_routine[1][2][0][0][1]= int2v1200eAB; 
  build_routine[1][2][0][1][0]= int2v1201 ; build_routine[1][2][0][1][1]= int2v1201eAB; 
  build_routine[1][2][1][1][0]= int2v1211 ; build_routine[1][2][1][1][1]= int2v1211eAB; 
  build_routine[1][2][1][2][0]= int2v1212 ; build_routine[1][2][1][2][1]= int2v1212eAB; 
  build_routine[1][2][2][2][0]= int2v1222 ; build_routine[1][2][2][2][1]= int2v1222eAB; 
  build_routine[2][2][0][0][0]= int2v2200 ; build_routine[2][2][0][0][1]= int2v2200eAB; 
  build_routine[2][2][0][1][0]= int2v2201 ; build_routine[2][2][0][1][1]= int2v2201eAB; 
  build_routine[2][2][1][1][0]= int2v2211 ; build_routine[2][2][1][1][1]= int2v2211eAB; 
  build_routine[2][2][2][2][0]= int2v2222 ; build_routine[2][2][2][2][1]= int2v2222eAB; 
#else
#if (MG == 1) || (MG == 2) || (MG == 3) || (MG == 4)
  ASSIGN_BUILD(0,1,0,0)
  ASSIGN_BUILD(0,1,0,1)
  ASSIGN_BUILD(0,1,1,1)
  ASSIGN_BUILD(1,1,0,0)
  ASSIGN_BUILD(1,1,1,1)
#endif

#if (MG == 2) || (MG == 3) || (MG == 4)
  ASSIGN_BUILD(0,2,0,0)
  ASSIGN_BUILD(0,2,0,1)
  ASSIGN_BUILD(0,2,0,2)
  ASSIGN_BUILD(0,2,1,1)
  ASSIGN_BUILD(0,2,1,2)
  ASSIGN_BUILD(0,2,2,2)
  ASSIGN_BUILD(1,2,0,0)
  ASSIGN_BUILD(1,2,0,1)
  ASSIGN_BUILD(1,2,1,1)
  ASSIGN_BUILD(1,2,1,2)
  ASSIGN_BUILD(1,2,2,2)
  ASSIGN_BUILD(2,2,0,0)
  ASSIGN_BUILD(2,2,0,1)
  ASSIGN_BUILD(2,2,1,1)
  ASSIGN_BUILD(2,2,2,2)
#endif

#if (MG == 3) || (MG == 4)
  ASSIGN_BUILD(0,3,0,0)
  ASSIGN_BUILD(0,3,0,1)
  ASSIGN_BUILD(0,3,0,2)
  ASSIGN_BUILD(0,3,0,3)
  ASSIGN_BUILD(0,3,1,1)
  ASSIGN_BUILD(0,3,1,2)
  ASSIGN_BUILD(0,3,1,3)
  ASSIGN_BUILD(0,3,2,2)
  ASSIGN_BUILD(0,3,2,3)
  ASSIGN_BUILD(0,3,3,3)
  ASSIGN_BUILD(1,3,0,0)
  ASSIGN_BUILD(1,3,0,1)
  ASSIGN_BUILD(1,3,0,2)
  ASSIGN_BUILD(1,3,1,1)
  ASSIGN_BUILD(1,3,1,2)
  ASSIGN_BUILD(1,3,1,3)
  ASSIGN_BUILD(1,3,2,2)
  ASSIGN_BUILD(1,3,2,3)
  ASSIGN_BUILD(1,3,3,3)
  ASSIGN_BUILD(2,3,0,0)
  ASSIGN_BUILD(2,3,0,1)
  ASSIGN_BUILD(2,3,0,2)
  ASSIGN_BUILD(2,3,1,1)
  ASSIGN_BUILD(2,3,1,2)
  ASSIGN_BUILD(2,3,2,2)
  ASSIGN_BUILD(2,3,2,3)
  ASSIGN_BUILD(2,3,3,3)
  ASSIGN_BUILD(3,3,0,0)
  ASSIGN_BUILD(3,3,0,1)
  ASSIGN_BUILD(3,3,0,2)
  ASSIGN_BUILD(3,3,1,1)
  ASSIGN_BUILD(3,3,1,2)
  ASSIGN_BUILD(3,3,2,2)
  ASSIGN_BUILD(3,3,3,3)
#endif
#endif

  free(jmax_for_con);
  saved_am12 = am12;
  saved_am34 = am34;
  saved_ncon = nc1;
  }

GLOBAL_FUNCTION VOID
int_done_buildgc()
{
  int i,j,k;
  int ci,cj,ck,cl;

  free_int_array3(&inthave);

  for (i=0; i<=saved_am12; i++) {
    for (j=0; j<=saved_am34; j++) {
      for (k=0; k<=saved_am12+saved_am34-i-j; k++) {
        free(int_v_list.dp[i][j][k]);
        }
      }
    }

  free_doublep_array3(&int_v_list);

#if 1
  free_store();
#else
  /* Since some storage is shared, this is a bit messy. */
  for (ci=0; ci<saved_ncon; ci++) {
    for (cj=0; cj<saved_ncon; cj++) {
      for (ck=0; ck<saved_ncon; ck++) {
        for (cl=0; cl<saved_ncon; cl++) {
          for (i=0; i<=saved_am12; i++) {
            for (j=0; j<=saved_am12; j++) {
              for (k=0; k<=saved_am34; k++) {
                for (l=0; l<=saved_am34; l++) {
  if (int_con_ints_array[ci][cj][ck][cl].dp[i][j][k][l]) {
    double *tmp = int_con_ints_array[ci][cj][ck][cl].dp[i][j][k][l];
    int ci2,cj2,ck2,cl2,i2,j2,k2,l2;

    /* Find all occurences of this buffer is set them to NULL. */

    for (ci2=0; ci2<saved_ncon; ci2++) {
      for (cj2=0; cj2<saved_ncon; cj2++) {
        for (ck2=0; ck2<saved_ncon; ck2++) {
          for (cl2=0; cl2<saved_ncon; cl2++) {
            for (i2=0; i2<=saved_am12; i2++) {
              for (j2=0; j2<=saved_am12; j2++) {
                for (k2=0; k2<=saved_am34; k2++) {
                  for (l2=0; l2<=saved_am34; l2++) {
    if (int_con_ints_array[ci2][cj2][ck2][cl2].dp[i2][j2][k2][l2] == tmp) {
      int_con_ints_array[ci2][cj2][ck2][cl2].dp[i2][j2][k2][l2] = NULL;
      }

                    }
                  }
                }
              }
            }
          }
        }
      }
    /* Now we finally free the storage. */
    free(tmp);
    }
                  }
                }
              }
            }
          }
        }
      }
    }
#endif

  for (ci=0; ci<saved_ncon; ci++) {
    for (cj=0; cj<saved_ncon; cj++) {
      for (ck=0; ck<saved_ncon; ck++) {
        for (cl=0; cl<saved_ncon; cl++) {
          free_doublep_array4(&int_con_ints_array[ci][cj][ck][cl]);
          }
        free(int_con_ints_array[ci][cj][ck]);
        }
      free(int_con_ints_array[ci][cj]);
      }
    free(int_con_ints_array[ci]);
    }
  free(int_con_ints_array);

  }

/* add_store maintains a list of free storage allocated by int_init_buildgc */
LOCAL_FUNCTION VOID
add_store(p)
void *p;
{
  if (!store) {
    store = (store_list_t*) malloc(sizeof(store_list_t));
    assert(store);
    store->p = NULL;
    n_store_last = 0;
    }
  if (n_store_last >= STORAGE_CHUNK) {
    store_list_t* tmp = (store_list_t*) malloc(sizeof(store_list_t));
    assert(tmp);
    tmp->p = store;
    store = tmp;
    n_store_last = 0;
    }
  store->data[n_store_last++] = p;
  }

/* free_store frees the memory that add_store keeps track of */
LOCAL_FUNCTION VOID
free_store()
{
  _free_store(store,n_store_last);
  store = NULL;
  }

LOCAL_FUNCTION VOID
_free_store(s,n)
store_list_t* s;
int n;
{
  int i;
  if (!s) return;
  for (i=0; i<n; i++) {
    free(s->data[i]);
    }
  _free_store(s->p,STORAGE_CHUNK);
  free(s);
  }


GLOBAL_FUNCTION VOID
int_buildgcam(minam1,minam2,minam3,minam4,
              maxam1,maxam2,maxam3,maxam4,
              dam1,dam2,dam3,dam4,
              sh1,sh2,sh3,sh4, eAB)
int minam1;
int minam2;
int minam3;
int minam4;
int maxam1;
int maxam2;
int maxam3;
int maxam4;
int dam1;
int dam2;
int dam3;
int dam4;
int sh1;
int sh2;
int sh3;
int sh4;
int eAB;
{
  int k,m,n;
  int ci,cj,ck,cl;
  int maxam12,maxam34;
  int nc1,nc2,nc3,nc4;

  if (maxam1<0 || maxam2<0 || maxam3<0 || maxam4<0) return;
  if (minam1<0) minam1=0;
  if (minam2<0) minam2=0;
  if (minam3<0) minam3=0;
  if (minam4<0) minam4=0;

  maxam12 = maxam1 + maxam2;
  maxam34 = maxam3 + maxam4;

#if 0
  fprintf(stdout,"min=(%d,%d,%d,%d) max=(%d,%d,%d,%d)\n",
          minam1,
          minam2,
          minam3,
          minam4,
          maxam1,
          maxam2,
          maxam3,
          maxam4);
#endif

  /* Compute the offset shell numbers. */
  osh1 = sh1 + int_cs1->shell_offset;
  osh2 = sh2 + int_cs2->shell_offset;
  osh3 = sh3 + int_cs3->shell_offset;
  osh4 = sh4 + int_cs4->shell_offset;

  nc1 = int_cs1->center[int_cs1->center_num[sh1]]
                 .basis.shell[int_cs1->shell_num[sh1]].ncon;
  nc2 = int_cs2->center[int_cs2->center_num[sh2]]
                 .basis.shell[int_cs2->shell_num[sh2]].ncon;
  nc3 = int_cs3->center[int_cs3->center_num[sh3]]
                 .basis.shell[int_cs3->shell_num[sh3]].ncon;
  nc4 = int_cs4->center[int_cs4->center_num[sh4]]
                 .basis.shell[int_cs4->shell_num[sh4]].ncon;

  /* Zero the target contracted integrals that the build routine
   * will accumulate into. */
  for (m=minam1; m<=maxam12; m++) {
    for (n=minam3; n<=maxam34; n++) {
  for (ci=0; ci<nc1; ci++) {
    if (m < int_shell1->type[ci].am+dam1) continue;
    for (cj=0; cj<nc2; cj++) {
      if (int_shell1->type[ci].am+dam1+int_shell2->type[cj].am+dam2 < m)
        continue;
      for (ck=0; ck<nc3; ck++) {
        if (n < int_shell3->type[ck].am+dam3) continue;
        for (cl=0; cl<nc4; cl++) {
          if (int_shell3->type[ck].am+dam3 +int_shell4->type[cl].am+dam4 < n)
            continue;
      for (k=0; k<INT_NCART(m)*INT_NCART(n); k++) {
        int_con_ints_array[ci][cj][ck][cl].dp[m][0][n][0][k] = 0.0;
        }
          }
        }
      }
    }
      }
    }

  gen_shell_intermediates(sh1,sh2,sh3,sh4);

  /* If enough of the functions come from generalized contractions
   * to make it workwhile, then don't do redundant primitives
   * at the additional cost of slower normalization computations.
   */
  if (nc1 + nc2 + nc3 + nc4 > 4)
    build_using_gcs(nc1,nc2,nc3,nc4,
                    minam1,minam3,maxam12,maxam34,dam1,dam2,dam3,dam4,eAB);
  else
    build_not_using_gcs(nc1,nc2,nc3,nc4,
                    minam1,minam3,maxam12,maxam34,dam1,dam2,dam3,dam4,eAB);

  }

LOCAL_FUNCTION VOID
build_not_using_gcs(nc1,nc2,nc3,nc4,minam1,minam3,maxam12,maxam34,dam1,dam2,dam3,dam4,eAB)
int nc1;
int nc2;
int nc3;
int nc4;
int minam1;
int minam3;
int maxam12;
int maxam34;
int dam1;
int dam2;
int dam3;
int dam4;
int eAB;
{
  int have_all_ints;
  int i,j,k,l,m,n,o;
  int ci,cj,ck,cl;
  double *bufferprim;
  double *con_ints;
  double *tmpbufferprim;

#if 0
  printf("not_gcs: %d%d%d%d\n",
         int_expweight1,
         int_expweight2,
         int_expweight3,
         int_expweight4
         );
#endif

          /* Sum thru all possible contractions. */
  for (ci=0; ci<nc1; ci++) {
    for (cj=0; cj<nc2; cj++) {
      for (ck=0; ck<nc3; ck++) {
        for (cl=0; cl<nc4; cl++) {


  /* Loop over the primitives. */
  for (i=0; i<int_shell1->nprim; i++) {
    double coef0;
    coef0 = int_shell1->coef[ci][i];
    if (int_expweight1) coef0 = coef0
                                    * int_shell1->exp[i];
    /* This factor of two comes from the derivative integral formula. */
    if (int_expweight1) coef0 *= 2.0;
    if (int_expweight2) coef0 *= 2.0;
    if (int_expweight3) coef0 *= 2.0;
    if (int_expweight4) coef0 *= 2.0;
    if (int_store1) opr1 = int_shell_to_prim.i[osh1] + i;
    for (j=0; j<int_shell2->nprim; j++) {
      double coef1;
      coef1 = int_shell2->coef[cj][j];
      if (int_expweight2) coef1 *=  coef0
                                      * int_shell2->exp[j];
      else                     coef1 *= coef0;
      if (int_store1) opr2 = int_shell_to_prim.i[osh2] + j;
      for (k=0; k<int_shell3->nprim; k++) {
        double coef2;
        coef2 = int_shell3->coef[ck][k];
        if (int_expweight3) coef2 *=  coef1
                                        * int_shell3->exp[k];
        else                     coef2 *= coef1;
        if (int_store1) opr3 = int_shell_to_prim.i[osh3] + k;
        for (l=0; l<int_shell4->nprim; l++) {
          double coef3;
          coef3 = int_shell4->coef[cl][l];
          if (int_expweight4) coef3 *=  coef2
                                          * int_shell4->exp[l];
          else                     coef3 *= coef2;
          if (int_store1) opr4 = int_shell_to_prim.i[osh4] + l;

          /* Produce the remaining intermediates. */
          gen_prim_intermediates_with_norm(i,j,k,l, maxam12+maxam34,coef3);

          /* Generate the target integrals. */
          if ((maxam12 == 0) && (maxam34 == 0)) {
            /* Do nothing: gen_prim_intermediates has set everything up. */
            have_all_ints = 1;
            }
          else if ((minam1<=MG)&&(minam3<=MG)&&(maxam12<=MG)&&(maxam34<=MG)) {
            if (build_routine[minam1]
                             [maxam12]
                             [minam3]
                             [maxam34][eAB] == impossible_integral) {
              fprintf(stderr,"trying to build with int2v%d%d%d%d (exact)\n",
                      minam1,maxam12,minam3,maxam34);
              }
            if ((*build_routine[minam1]
                               [maxam12]
                               [minam3]
                               [maxam34][eAB])()) {
              /* Mark the integrals as computed. */
              have_all_ints = 1;
              }
            else {
              have_all_ints = 0;
              }
            }
          else if (MG > 0) {
            int backminam1  = MINA(minam1);
            int backmaxam12 = MINA(maxam12);
            int backminam3  = MINA(minam3);
            int backmaxam34 = MINA(maxam34);
            if ((backmaxam12==backmaxam34)&&(backminam1>backminam3))
              backminam3 = backminam1;
            /* We cannot build everything, so build what we can. */
            if (build_routine[backminam1]
                             [backmaxam12]
                             [backminam3]
                             [backmaxam34][eAB] == impossible_integral) {
              fprintf(stderr,"trying to build with int2v%d%d%d%d\n",
                      backminam1,backmaxam12,backminam3,backmaxam34);
              }
            have_all_ints = 0;
            if ((*build_routine[backminam1]
                               [backmaxam12]
                               [backminam3]
                               [backmaxam34][eAB])()) {

              /* Mark all the intermediates as being noncomputed. */
              init_inthave(maxam12,maxam34);

              /* Mark the integrals as computed. */
              for (m=backminam1; m<=backmaxam12; m++) {
                for (n=backminam3; n<=backmaxam34; n++) {
                  inthave.i[m][n][0] = 1;
                  }
                }
              }
            else {
              fprintf(stderr,"backup build routine not available\n");
              fail();
              }
            }
          else {
            init_inthave(maxam12,maxam34);
            have_all_ints = 0;
            }

          /* Contract the primitive target integrals. */
          /* Throw out all unneeded contractions. */
          for (m=minam1; m<=maxam12; m++) {
            if (m < int_shell1->type[ci].am+dam1) continue;
            if (int_shell1->type[ci].am+dam1+int_shell2->type[cj].am+dam2 < m)
              continue;
            for (n=minam3; n<=maxam34; n++) {
              int sizemn;
              if (n < int_shell3->type[ck].am+dam3) continue;
              if (int_shell3->type[ck].am+dam3 +int_shell4->type[cl].am+dam4
                   < n)
                 continue;

              sizemn = INT_NCART(m)*INT_NCART(n);
              if (have_all_ints) bufferprim = int_v_list.dp[m][n][0];
              else bufferprim = buildprim(m, n, 0);

          con_ints = int_con_ints_array[ci][cj][ck][cl].dp[m][0][n][0];
          tmpbufferprim = bufferprim;

                /* Sum the integrals into the contracted integrals. */
                for (o=sizemn; o!=0; o--) {
                  *con_ints++ += *tmpbufferprim++;
                  }

                }
              }

          }
        }
      }
    }

          }
        }
      }
    }

#if 0
  print_int_array3(stdout,&inthave);
#endif

  }

LOCAL_FUNCTION VOID
build_using_gcs(nc1,nc2,nc3,nc4,minam1,minam3,maxam12,maxam34,dam1,dam2,dam3,dam4,eAB)
int nc1;
int nc2;
int nc3;
int nc4;
int minam1;
int minam3;
int maxam12;
int maxam34;
int dam1;
int dam2;
int dam3;
int dam4;
int eAB;
{
  int have_all_ints;
  int i,j,k,l,m,n,o;
  int ci,cj,ck,cl;
  int ist1,ist3;
  int nm3;
  int maxam1234=maxam12+maxam34;
  double coef0,coef1,coef2,coef3;
  double ishl1expi=1.0, ishl2expj=1.0, ishl3expk=1.0;
  double *bufferprim;
  double *con_ints;
  double *tmpbufferprim;
  doublep_array4_t *iciaptr;
  double c0scale;
  intfunc brptr=build_routine[minam1][maxam12][minam3][maxam34][eAB];

  /* Loop over the primitives. */
  for (i=0; i<int_shell1->nprim; i++) {
    if (int_store1) opr1 = int_shell_to_prim.i[osh1] + i;
    if (int_expweight1) ishl1expi=2.0*int_shell1->exp[i];

    for (j=0; j<int_shell2->nprim; j++) {
      if (int_store1) opr2 = int_shell_to_prim.i[osh2] + j;
      ishl2expj = (int_expweight2) ? 
                        2.0*int_shell2->exp[j]*ishl1expi : ishl1expi;

      for (k=0; k<int_shell3->nprim; k++) {
        if (int_store1) opr3 = int_shell_to_prim.i[osh3] + k;
        ishl3expk = (int_expweight3) ? 
                        2.0*int_shell3->exp[k]*ishl2expj : ishl2expj;

        for (l=0; l<int_shell4->nprim; l++) {
          if (int_store1) opr4 = int_shell_to_prim.i[osh4] + l;
          c0scale = (int_expweight4) ? 
                        2.0*int_shell4->exp[l]*ishl3expk : ishl3expk;

          /* Produce the remaining intermediates. */
          gen_prim_intermediates(i,j,k,l, maxam1234);

          /* Generate the target integrals. */
          if (!maxam1234) {
            /* Do nothing: gen_prim_intermediates has set everything up. */
            have_all_ints = 1;
            }
          else if ((minam1<=MG)&&(minam3<=MG)&&(maxam12<=MG)&&(maxam34<=MG)) {
            if (brptr == impossible_integral) {
              fprintf(stderr,"trying to build with int2v%d%d%d%d (exact)\n",
                      minam1,maxam12,minam3,maxam34);
              }
            if ((*brptr)()) {
              /* Mark the integrals as computed. */
              have_all_ints = 1;
              }
            else {
              have_all_ints = 0;
              }
            }
          else if (MG > 0) {
            int backminam1  = MINA(minam1);
            int backmaxam12 = MINA(maxam12);
            int backminam3  = MINA(minam3);
            int backmaxam34 = MINA(maxam34);
            intfunc brptr2;

            if ((backmaxam12==backmaxam34)&&(backminam1>backminam3))
              backminam3 = backminam1;
            brptr2=build_routine
                       [backminam1][backmaxam12][backminam3][backmaxam34][eAB];

            /* We cannot build everything, so build what we can. */
            if (brptr2 == impossible_integral) {
              fprintf(stderr,"trying to build with int2v%d%d%d%d\n",
                      backminam1,backmaxam12,backminam3,backmaxam34);
              }
            have_all_ints = 0;
            if ((*brptr2)()) {
              /* Mark all the intermediates as being noncomputed. */
              init_inthave(maxam12,maxam34);

              /* Mark the integrals as computed. */
              for (m=backminam1; m<=backmaxam12; m++) {
                for (n=backminam3; n<=backmaxam34; n++) {
                  inthave.i[m][n][0] = 1;
                  }
                }
              }
            else {
              fprintf(stderr,"backup build routine not available\n");
              fail();
              }
            }
          else {
            init_inthave(maxam12,maxam34);
            have_all_ints = 0;
            }

          /* Contract the primitive target integrals. */
          for (m=minam1; m<=maxam12; m++) {
            for (n=minam3; n<=maxam34; n++) {
              int sizemn = INT_NCART(m)*INT_NCART(n);
              if (have_all_ints) bufferprim = int_v_list.dp[m][n][0];
              else bufferprim = buildprim(m, n, 0);

          /* Sum thru all possible contractions.
           * Throw out all unneeded contractions. */

  for (ci=0; ci<nc1; ci++) {
    if (m < (ist1=int_shell1->type[ci].am+dam1)) continue;
    coef0 = int_shell1->coef[ci][i]*c0scale;

    for (cj=0; cj<nc2; cj++) {
      if (ist1+int_shell2->type[cj].am+dam2 < m) continue;

      coef1 = int_shell2->coef[cj][j]*coef0;

      for (ck=0; ck<nc3; ck++) {
        if (n < (ist3=int_shell3->type[ck].am+dam3)) continue;

        coef2 = int_shell3->coef[ck][k]*coef1;

        iciaptr = int_con_ints_array[ci][cj][ck];

        nm3=n-ist3-dam4;
        for (cl=0; cl<nc4; cl++,iciaptr++) {
          if (int_shell4->type[cl].am < nm3) continue;

          coef3 = int_shell4->coef[cl][l]*coef2;

          con_ints = *((*iciaptr).dp[m][0][n]);
          tmpbufferprim = bufferprim;

          /* Sum the integrals into the contracted integrals. */
          for (o=sizemn; o; o--) *con_ints++ += coef3 * *tmpbufferprim++;
          }
        }
      }
    }


              }
            }
          }
        }
      }
    }
  }

/* This routine constructs intermediates needed for each quartet of
 * primitives.  It is given the total angular momentum as the argument
 * and requires that the global primitive offsets and other global
 * constants be initialized. */
LOCAL_FUNCTION VOID
gen_prim_intermediates(pr1,pr2,pr3,pr4,am)
int pr1;
int pr2;
int pr3;
int pr4;
int am;
{
  int i;
  double T;
  double pmq,pmq2;
  double AmB,AmB2;
  /* This is 2^(1/2) * pi^(5/4) */
  CONST double sqrt2pi54 = 5.9149671727956129;
  double conv_to_s;

  if (int_store2) {
    int_v_zeta12 = int_prim_zeta.d[opr1][opr2];
    int_v_zeta34 = int_prim_zeta.d[opr3][opr4];
    int_v_oo2zeta12 = int_prim_oo2zeta.d[opr1][opr2];
    int_v_oo2zeta34 = int_prim_oo2zeta.d[opr3][opr4];
    int_v_p120 = int_prim_p.d[opr1][opr2][0];
    int_v_p121 = int_prim_p.d[opr1][opr2][1];
    int_v_p122 = int_prim_p.d[opr1][opr2][2];
    int_v_p340 = int_prim_p.d[opr3][opr4][0];
    int_v_p341 = int_prim_p.d[opr3][opr4][1];
    int_v_p342 = int_prim_p.d[opr3][opr4][2];
    int_v_k12 = int_prim_k.d[opr1][opr2];
    int_v_k34 = int_prim_k.d[opr3][opr4];
    }
  else {
    int_v_zeta12 = int_shell1->exp[pr1] + int_shell2->exp[pr2];
    int_v_zeta34 = int_shell3->exp[pr3] + int_shell4->exp[pr4];
    int_v_oo2zeta12 = 1.0/int_v_zeta12;
    int_v_oo2zeta34 = 1.0/int_v_zeta34;
    int_v_p120 = int_v_oo2zeta12 * ( int_shell1->exp[pr1] * int_v_r10 + int_shell2->exp[pr2] * int_v_r20 );
    int_v_p121 = int_v_oo2zeta12 * ( int_shell1->exp[pr1] * int_v_r11 + int_shell2->exp[pr2] * int_v_r21 );
    int_v_p122 = int_v_oo2zeta12 * ( int_shell1->exp[pr1] * int_v_r12 + int_shell2->exp[pr2] * int_v_r22 );
    int_v_p340 = int_v_oo2zeta34 * ( int_shell3->exp[pr3] * int_v_r30 + int_shell4->exp[pr4] * int_v_r40 );
    int_v_p341 = int_v_oo2zeta34 * ( int_shell3->exp[pr3] * int_v_r31 + int_shell4->exp[pr4] * int_v_r41 );
    int_v_p342 = int_v_oo2zeta34 * ( int_shell3->exp[pr3] * int_v_r32 + int_shell4->exp[pr4] * int_v_r42 );

    /* Compute AmB^2 for shell 1 and 2. */
    AmB = int_v_r20 - int_v_r10;
    AmB2 = AmB*AmB;
    AmB = int_v_r21 - int_v_r11;
    AmB2 += AmB*AmB;
    AmB = int_v_r22 - int_v_r12;
    AmB2 += AmB*AmB;

    int_v_k12 =    sqrt2pi54
                 * int_v_oo2zeta12
                 * exp_cutoff( -   int_shell1->exp[pr1] * int_shell2->exp[pr2]
                          * int_v_oo2zeta12
                          * AmB2 );

    /* Compute AmB^2 for shells 3 and 4. */
    AmB = int_v_r40 - int_v_r30;
    AmB2 = AmB*AmB;
    AmB = int_v_r41 - int_v_r31;
    AmB2 += AmB*AmB;
    AmB = int_v_r42 - int_v_r32;
    AmB2 += AmB*AmB;

    int_v_k34 =    sqrt2pi54
                 * int_v_oo2zeta34
                 * exp_cutoff( -   int_shell3->exp[pr3] * int_shell4->exp[pr4]
                          * int_v_oo2zeta34
                          * AmB2 );

    int_v_oo2zeta12 *= 0.5;
    int_v_oo2zeta34 *= 0.5;
    }

  int_v_ooze = 1.0/(int_v_zeta12 + int_v_zeta34);

  int_v_W0 = int_v_ooze*(  int_v_zeta12 * int_v_p120
                         + int_v_zeta34 * int_v_p340 );
  int_v_W1 = int_v_ooze*(  int_v_zeta12 * int_v_p121
                         + int_v_zeta34 * int_v_p341 );
  int_v_W2 = int_v_ooze*(  int_v_zeta12 * int_v_p122
                         + int_v_zeta34 * int_v_p342 );

  pmq = int_v_p120 - int_v_p340;
  pmq2 = pmq*pmq;
  pmq = int_v_p121 - int_v_p341;
  pmq2 += pmq*pmq;
  pmq = int_v_p122 - int_v_p342;
  pmq2 += pmq*pmq;

  T =   int_v_zeta12
      * int_v_zeta34
      * int_v_ooze * pmq2;

  int_fjt(am,T);

  /* Convert the fjttable produced by int_fjt into the S integrals */
  conv_to_s = sqrt(int_v_ooze) * int_v_k12 * int_v_k34;
  for (i=0; i<=am; i++) {
    int_v_list.dp[0][0][i][0] =   int_fjttable.d[i] * conv_to_s;
#if 0
    fprintf(stdout,"int_v_list.dp[0][0][%d][0] = %lf\n",
            i,int_v_list.dp[0][0][i][0]);
#endif
    }

  }

/* This is like gen_prim_intermediates, except the normalization is
 * put into the ssss integrals. */
LOCAL_FUNCTION VOID
gen_prim_intermediates_with_norm(pr1,pr2,pr3,pr4,am,norm)
int pr1;
int pr2;
int pr3;
int pr4;
int am;
double norm;
{
  int i;
  double T;
  double pmq,pmq2;
  double AmB,AmB2;
  /* This is 2^(1/2) * pi^(5/4) */
  CONST double sqrt2pi54 = 5.9149671727956129;
  double conv_to_s;

  if (int_store2) {
    int_v_zeta12 = int_prim_zeta.d[opr1][opr2];
    int_v_zeta34 = int_prim_zeta.d[opr3][opr4];
    int_v_oo2zeta12 = int_prim_oo2zeta.d[opr1][opr2];
    int_v_oo2zeta34 = int_prim_oo2zeta.d[opr3][opr4];
    int_v_p120 = int_prim_p.d[opr1][opr2][0];
    int_v_p121 = int_prim_p.d[opr1][opr2][1];
    int_v_p122 = int_prim_p.d[opr1][opr2][2];
    int_v_p340 = int_prim_p.d[opr3][opr4][0];
    int_v_p341 = int_prim_p.d[opr3][opr4][1];
    int_v_p342 = int_prim_p.d[opr3][opr4][2];
    int_v_k12 = int_prim_k.d[opr1][opr2];
    int_v_k34 = int_prim_k.d[opr3][opr4];
    }
  else {
    int_v_zeta12 = int_shell1->exp[pr1] + int_shell2->exp[pr2];
    int_v_zeta34 = int_shell3->exp[pr3] + int_shell4->exp[pr4];
    int_v_oo2zeta12 = 1.0/int_v_zeta12;
    int_v_oo2zeta34 = 1.0/int_v_zeta34;
    int_v_p120 = int_v_oo2zeta12 * ( int_shell1->exp[pr1] * int_v_r10 + int_shell2->exp[pr2] * int_v_r20 );
    int_v_p121 = int_v_oo2zeta12 * ( int_shell1->exp[pr1] * int_v_r11 + int_shell2->exp[pr2] * int_v_r21 );
    int_v_p122 = int_v_oo2zeta12 * ( int_shell1->exp[pr1] * int_v_r12 + int_shell2->exp[pr2] * int_v_r22 );
    int_v_p340 = int_v_oo2zeta34 * ( int_shell3->exp[pr3] * int_v_r30 + int_shell4->exp[pr4] * int_v_r40 );
    int_v_p341 = int_v_oo2zeta34 * ( int_shell3->exp[pr3] * int_v_r31 + int_shell4->exp[pr4] * int_v_r41 );
    int_v_p342 = int_v_oo2zeta34 * ( int_shell3->exp[pr3] * int_v_r32 + int_shell4->exp[pr4] * int_v_r42 );

    /* Compute AmB^2 for shell 1 and 2. */
    AmB = int_v_r20 - int_v_r10;
    AmB2 = AmB*AmB;
    AmB = int_v_r21 - int_v_r11;
    AmB2 += AmB*AmB;
    AmB = int_v_r22 - int_v_r12;
    AmB2 += AmB*AmB;

    int_v_k12 =    sqrt2pi54
                 * int_v_oo2zeta12
                 * exp_cutoff( -   int_shell1->exp[pr1] * int_shell2->exp[pr2]
                          * int_v_oo2zeta12
                          * AmB2 );

    /* Compute AmB^2 for shells 3 and 4. */
    AmB = int_v_r40 - int_v_r30;
    AmB2 = AmB*AmB;
    AmB = int_v_r41 - int_v_r31;
    AmB2 += AmB*AmB;
    AmB = int_v_r42 - int_v_r32;
    AmB2 += AmB*AmB;

    int_v_k34 =    sqrt2pi54
                 * int_v_oo2zeta34
                 * exp_cutoff( -   int_shell3->exp[pr3] * int_shell4->exp[pr4]
                          * int_v_oo2zeta34
                          * AmB2 );

    int_v_oo2zeta12 *= 0.5;
    int_v_oo2zeta34 *= 0.5;
    }

  int_v_ooze = 1.0/(int_v_zeta12 + int_v_zeta34);

  int_v_W0 = int_v_ooze*(  int_v_zeta12 * int_v_p120
                         + int_v_zeta34 * int_v_p340 );
  int_v_W1 = int_v_ooze*(  int_v_zeta12 * int_v_p121
                         + int_v_zeta34 * int_v_p341 );
  int_v_W2 = int_v_ooze*(  int_v_zeta12 * int_v_p122
                         + int_v_zeta34 * int_v_p342 );

  pmq = int_v_p120 - int_v_p340;
  pmq2 = pmq*pmq;
  pmq = int_v_p121 - int_v_p341;
  pmq2 += pmq*pmq;
  pmq = int_v_p122 - int_v_p342;
  pmq2 += pmq*pmq;

  T =   int_v_zeta12
      * int_v_zeta34
      * int_v_ooze * pmq2;

  int_fjt(am,T);

  /* Convert the fjttable produced by int_fjt into the S integrals */
  conv_to_s = sqrt(int_v_ooze) * int_v_k12 * int_v_k34 * norm;
  for (i=0; i<=am; i++) {
    int_v_list.dp[0][0][i][0] =   int_fjttable.d[i] * conv_to_s;
#if 0
    fprintf(stdout,"int_v_list.dp[0][0][%d][0] = %lf\n",
            i,int_v_list.dp[0][0][i][0]);
#endif
    }

  }


/* This routine computes the shell intermediates. */
LOCAL_FUNCTION VOID
gen_shell_intermediates(sh1,sh2,sh3,sh4)
int sh1;
int sh2;
int sh3;
int sh4;
{
  if (int_store1) {
    int_v_r10 = int_shell_r.dp[osh1][0];
    int_v_r11 = int_shell_r.dp[osh1][1];
    int_v_r12 = int_shell_r.dp[osh1][2];
    int_v_r20 = int_shell_r.dp[osh2][0];
    int_v_r21 = int_shell_r.dp[osh2][1];
    int_v_r22 = int_shell_r.dp[osh2][2];
    int_v_r30 = int_shell_r.dp[osh3][0];
    int_v_r31 = int_shell_r.dp[osh3][1];
    int_v_r32 = int_shell_r.dp[osh3][2];
    int_v_r40 = int_shell_r.dp[osh4][0];
    int_v_r41 = int_shell_r.dp[osh4][1];
    int_v_r42 = int_shell_r.dp[osh4][2];
    }
  else {
    int_v_r10 = int_cs1->center[int_cs1->center_num[sh1]].r[0];
    int_v_r11 = int_cs1->center[int_cs1->center_num[sh1]].r[1];
    int_v_r12 = int_cs1->center[int_cs1->center_num[sh1]].r[2];
    int_v_r20 = int_cs2->center[int_cs2->center_num[sh2]].r[0];
    int_v_r21 = int_cs2->center[int_cs2->center_num[sh2]].r[1];
    int_v_r22 = int_cs2->center[int_cs2->center_num[sh2]].r[2];
    int_v_r30 = int_cs3->center[int_cs3->center_num[sh3]].r[0];
    int_v_r31 = int_cs3->center[int_cs3->center_num[sh3]].r[1];
    int_v_r32 = int_cs3->center[int_cs3->center_num[sh3]].r[2];
    int_v_r40 = int_cs4->center[int_cs4->center_num[sh4]].r[0];
    int_v_r41 = int_cs4->center[int_cs4->center_num[sh4]].r[1];
    int_v_r42 = int_cs4->center[int_cs4->center_num[sh4]].r[2];
    }
  }

/* This builds up the primitive integrals of the type [x0|y0](m). */
LOCAL_FUNCTION double *
buildprim(am12, am34, m)
int am12;
int am34;
int m;
{
  double *buffer;

  /* Is this no integral? */
  if ((am12 < 0) || (am34 < 0)) return NULL;

  /* Is this integral on the list of computed integrals? */
  if (inthave.i[am12][am34][m]) return int_v_list.dp[am12][am34][m];

  /* Find the preallocated storage for the integrals. */
  buffer = int_v_list.dp[am12][am34][m];

  /* Should we build on center 1 or center 3. */
  if (choose_center(am12,am34,m) == 1) {
    /* Build on 1. */
    buildprim_1(buffer,am12,am34,m);
    }
  else {
    /* Build on 3. */
    buildprim_3(buffer,am12,am34,m);
    }

  /* Put the integrals in the list of computed integrals. */
  inthave.i[am12][am34][m] = 1;

#if 0
  fprintf(stdout,"buildprim: (%d 0|%d 0)(%d):\n",am12,am34,m);
  int_print_n(stdout,buffer
              ,INT_NCART(am12)
              ,1
              ,INT_NCART(am34)
              ,1
              ,0,0,0);
#endif

  return buffer;
  }

/* I00 will be made [a+1 0|b 0](m) */
LOCAL_FUNCTION VOID
buildprim_1(I00,am12,am34,m)
double *I00;
int am12;
int am34;
int m;
{
  double *I10; /* = [a0|c0](m) */
  double *I11; /* = [a0|c0](m+1) */
  double *I20; /* = [a-1 0|c0](m) */
  double *I21; /* = [a-1 0|c0](m+1) */
  double *I31; /* = [a0|c-1 0](m+1) */
  int cartindex12;
  int cartindex34;
  int cartindex1234;
  int size34,size34m1;
  int i12, j12, k12;
  int i34, j34, k34;

  /* Construct the needed intermediate integrals. */
  I10 = buildprim(am12 - 1, am34, m);
  I11 = buildprim(am12 - 1, am34, m + 1);
  I20 = buildprim(am12 - 2, am34, m);
  I21 = buildprim(am12 - 2, am34, m + 1);
  I31 = buildprim(am12 - 1, am34 - 1, m + 1);

  /* The size of the am34 group of primitives. */
  size34 = INT_NCART(am34);
  /* The size of the group of primitives with ang. mom. = am34 - 1 */
  size34m1 = INT_NCART(am34-1);

  /* Construct the new integrals. */
  cartindex12 = 0;
  cartindex1234 = 0;
  for (i12=0; i12<=am12; i12++) {
    for (k12=0; k12<=am12-i12; k12++) {
      j12 = am12 - i12 - k12;
      cartindex34 = 0;
      for (i34=0; i34<=am34; i34++) {
        for (k34=0; k34<=am34-i34; k34++) {
          j34 = am34 - i34 - k34;

          /* I10 I11 I20 I21 and I31 contrib */

          /* ------------------ Build from the x position. */
          if (i12) {
            /* I10 and I11 */
            I00[cartindex1234]
              = I10[INT_CARTINDEX(am12-1,i12-1,j12)*size34 + cartindex34]
                * ( int_v_p120 - int_v_r10)
               + I11[INT_CARTINDEX(am12-1,i12-1,j12)*size34 + cartindex34]
                 * ( int_v_W0 - int_v_p120);
            if (i12 > 1) {
              /* I20 and I21 */
              I00[cartindex1234]
               +=  (i12 - 1) * int_v_oo2zeta12
                 * (  I20[INT_CARTINDEX(am12-2,i12-2,j12)*size34 + cartindex34]
                    - I21[INT_CARTINDEX(am12-2,i12-2,j12)*size34 + cartindex34]
                      * int_v_zeta34 * int_v_ooze);
              }
            if (i34) {
              /* I31 */
              I00[cartindex1234]
               +=  i34 * 0.5 * int_v_ooze
                  * I31[  INT_CARTINDEX(am12-1,i12-1,j12)*size34m1
                        + INT_CARTINDEX(am34-1,i34-1,j34)];
              }
            }
          /* ------------------ Build from the y position. */
          else if (j12) {
            I00[cartindex1234]
              = I10[INT_CARTINDEX(am12-1,i12,j12-1)*size34 + cartindex34]
                * ( int_v_p121 - int_v_r11)
               + I11[INT_CARTINDEX(am12-1,i12,j12-1)*size34 + cartindex34]
                 * ( int_v_W1 - int_v_p121);
            if (j12 > 1) {
              I00[cartindex1234]
               +=  (j12 - 1) * int_v_oo2zeta12
                 * (  I20[INT_CARTINDEX(am12-2,i12,j12-2)*size34 + cartindex34]
                    - I21[INT_CARTINDEX(am12-2,i12,j12-2)*size34 + cartindex34]
                      * int_v_zeta34 * int_v_ooze);
              }
            if (j34) {
              I00[cartindex1234]
               +=  j34 * 0.5 * int_v_ooze
                  * I31[  INT_CARTINDEX(am12-1,i12,j12-1)*size34m1
                        + INT_CARTINDEX(am34-1,i34,j34-1)];
              }
            }
          /* ------------------ Build from the z position. */
          else if (k12) {
            I00[cartindex1234]
              = I10[INT_CARTINDEX(am12-1,i12,j12)*size34 + cartindex34]
                * ( int_v_p122 - int_v_r12)
               + I11[INT_CARTINDEX(am12-1,i12,j12)*size34 + cartindex34]
                 * ( int_v_W2 - int_v_p122);
            if (k12 > 1) {
              I00[cartindex1234]
               +=  (k12 - 1) * int_v_oo2zeta12
                 * (  I20[INT_CARTINDEX(am12-2,i12,j12)*size34 + cartindex34]
                    - I21[INT_CARTINDEX(am12-2,i12,j12)*size34 + cartindex34]
                      * int_v_zeta34 * int_v_ooze);
              }
            if (k34) {
              I00[cartindex1234]
               +=  k34 * 0.5 * int_v_ooze
                  * I31[  INT_CARTINDEX(am12-1,i12,j12)*size34m1
                        + INT_CARTINDEX(am34-1,i34,j34)];
              }
            }

          /* cartindex34 == INT_CARTINDEX(am34,i34,j34) */
          cartindex34++;
          cartindex1234++;
          }
        }
      /* cartindex12 == INT_CARTINDEX(am12,i12,j12) */
      cartindex12++;
      }
    }

  }


/* I00 will be made [a 0|b+1 0](m) */
LOCAL_FUNCTION VOID
buildprim_3(I00,am12,am34,m)
double *I00;
int am12;
int am34;
int m;
{
  double *I10; /* = [a0|c0](m) */
  double *I11; /* = [a0|c0](m+1) */
  double *I20; /* = [a0|c-1 0](m) */
  double *I21; /* = [a0|c-1 0](m+1) */
  double *I31; /* = [a-1 0|c0](m+1) */
  int cartindex12;
  int cartindex34;
  int cartindex1234;
  int size34m1,size34m2;
  int i12, j12, k12;
  int i34, j34, k34;

  /* Construct the needed intermediate integrals. */
  I10 = buildprim(am12, am34 - 1, m);
  I11 = buildprim(am12, am34 - 1, m + 1);
  I20 = buildprim(am12, am34 - 2, m);
  I21 = buildprim(am12, am34 - 2, m + 1);
  I31 = buildprim(am12 - 1, am34 - 1, m + 1);

  /* The size of the group of primitives with ang. mom. = am34 - 1 */
  size34m1 = INT_NCART(am34-1);
  size34m2 = INT_NCART(am34-2);

  /* Construct the new integrals. */
  cartindex12 = 0;
  cartindex1234 = 0;
  for (i12=0; i12<=am12; i12++) {
    for (k12=0; k12<=am12-i12; k12++) {
      j12 = am12 - i12 - k12;
      cartindex34 = 0;
      for (i34=0; i34<=am34; i34++) {
        for (k34=0; k34<=am34-i34; k34++) {
          j34 = am34 - i34 - k34;

          /* I10 I11 I20 I21 and I31 contrib */

          /* ------------------ Build from the x position. */
          if (i34) {
            /* I10 and I11 */
            I00[cartindex1234]
        /*
         *    = I10[INT_CARTINDEX(am12-1,i12-1,j12)*size34 + cartindex34]
         *      * ( int_v_p120 - int_v_r10)
         *     + I11[INT_CARTINDEX(am12-1,i12-1,j12)*size34 + cartindex34]
         *       * ( int_v_W0 - int_v_p120);
         */
              = I10[cartindex12*size34m1 + INT_CARTINDEX(am34-1,i34-1,j34)]
                * ( int_v_p340 - int_v_r30)
               + I11[cartindex12*size34m1 + INT_CARTINDEX(am34-1,i34-1,j34)]
                 * ( int_v_W0 - int_v_p340);
            if (i34 > 1) {
              /* I20 and I21 */
              I00[cartindex1234]
         /*
          *    +=  (i12 - 1) * int_v_oo2zeta12
          *      * (  I20[INT_CARTINDEX(am12-2,i12-2,j12)*size34 + cartindex34]
          *         - I21[INT_CARTINDEX(am12-2,i12-2,j12)*size34 + cartindex34]
          *           * int_v_zeta34 * int_v_ooze);
          */
               +=  (i34 - 1) * int_v_oo2zeta34
                 * (  I20[cartindex12*size34m2 +INT_CARTINDEX(am34-2,i34-2,j34)]
                    - I21[cartindex12*size34m2 +INT_CARTINDEX(am34-2,i34-2,j34)]
                      * int_v_zeta12 * int_v_ooze);
              }
            if (i12) {
              /* I31 */
              I00[cartindex1234]
        /*  
         *     +=  (i34 - 1) * 0.5 * int_v_ooze
         *        * I31[  INT_CARTINDEX(am12-1,i12-1,j12)*size34m1
         *              + INT_CARTINDEX(am34-1,i34-1,j34)];
         */
               +=  i12 * 0.5 * int_v_ooze
                  * I31[  INT_CARTINDEX(am12-1,i12-1,j12)*size34m1
                        + INT_CARTINDEX(am34-1,i34-1,j34)];
              }
            }
          /* ------------------ Build from the y position. */
          else if (j34) {
            /* I10 and I11 */
            I00[cartindex1234]
              = I10[cartindex12*size34m1 + INT_CARTINDEX(am34-1,i34,j34-1)]
                * ( int_v_p341 - int_v_r31)
               + I11[cartindex12*size34m1 + INT_CARTINDEX(am34-1,i34,j34-1)]
                 * ( int_v_W1 - int_v_p341);
            if (j34 > 1) {
              /* I20 and I21 */
              I00[cartindex1234]
               +=  (j34 - 1) * int_v_oo2zeta34
                 * (  I20[cartindex12*size34m2 +INT_CARTINDEX(am34-2,i34,j34-2)]
                    - I21[cartindex12*size34m2 +INT_CARTINDEX(am34-2,i34,j34-2)]
                      * int_v_zeta12 * int_v_ooze);
              }
            if (j12) {
              /* I31 */
              I00[cartindex1234]
               +=  j12 * 0.5 * int_v_ooze
                  * I31[  INT_CARTINDEX(am12-1,i12,j12-1)*size34m1
                        + INT_CARTINDEX(am34-1,i34,j34-1)];
              }
            }
          /* ------------------ Build from the z position. */
          else if (k34) {
            /* I10 and I11 */
            I00[cartindex1234]
              = I10[cartindex12*size34m1 + INT_CARTINDEX(am34-1,i34,j34)]
                * ( int_v_p342 - int_v_r32)
               + I11[cartindex12*size34m1 + INT_CARTINDEX(am34-1,i34,j34)]
                 * ( int_v_W2 - int_v_p342);
            if (k34 > 1) {
              /* I20 and I21 */
              I00[cartindex1234]
               +=  (k34 - 1) * int_v_oo2zeta34
                 * (  I20[cartindex12*size34m2 +INT_CARTINDEX(am34-2,i34,j34)]
                    - I21[cartindex12*size34m2 +INT_CARTINDEX(am34-2,i34,j34)]
                      * int_v_zeta12 * int_v_ooze);
              }
            if (k12) {
              /* I31 */
              I00[cartindex1234]
               +=  k12 * 0.5 * int_v_ooze
                  * I31[  INT_CARTINDEX(am12-1,i12,j12)*size34m1
                        + INT_CARTINDEX(am34-1,i34,j34)];
              }
            }

          /* cartindex34 == INT_CARTINDEX(am34,i34,j34) */
          cartindex34++;
          cartindex1234++;
          }
        }
      /* cartindex12 == INT_CARTINDEX(am12,i12,j12) */
      cartindex12++;
      }
    }

  }

/* Initialize the list of integrals which have been precomputed
 * to "not computed" (=0). */
LOCAL_FUNCTION VOID
init_inthave(am12,am34)
int am12;
int am34;
{
  int i,j,k;

#if 0 /* This cleans up inthave if the whole thing is to be printed */
  zero_int_array3(&inthave);
#endif

  for (i=0; i<=am12; i++) {
    for (j=0; j<=am34; j++) {
      for (k=0; k<=am12+am34-i-j; k++) {
        if ((i==0)&&(j==0)) inthave.i[i][j][k] = 1;
        else inthave.i[i][j][k] = 0;
        }
      }
    }
  }


/* Determines which center it is best to build upon.
 * (This can be improved, perhaps.) */
LOCAL_FUNCTION int
choose_center(am12,am34,m)
int am12;
int am34;
int m;
{
  int need1 = 0;
  int need3 = 0;

  if (am12==0) return 3;
  if (am34==0) return 1;

  if (!inthave.i[am12-1][am34][m]) {
    need1 += INT_NCART(am12-1)*INT_NCART(am34);
    }
  if (!inthave.i[am12-1][am34][m+1]) {
    need1 += INT_NCART(am12-1)*INT_NCART(am34);
    }
  if ((am12>1) && (!inthave.i[am12-2][am34][m])) {
    need1 += INT_NCART(am12-2)*INT_NCART(am34);
    }
  if ((am12>1) && (!inthave.i[am12-2][am34][m+1])) {
    need1 += INT_NCART(am12-2)*INT_NCART(am34);
    }
  if (!inthave.i[am12-1][am34-1][m+1]) {
    need1 += INT_NCART(am12-1)*INT_NCART(am34-1);
    }

  if (!inthave.i[am12][am34-1][m]) {
    need3 += INT_NCART(am12)*INT_NCART(am34-1);
    }
  if (!inthave.i[am12][am34-1][m+1]) {
    need3 += INT_NCART(am12)*INT_NCART(am34-1);
    }
  if ((am34>1) && (!inthave.i[am12][am34-2][m])) {
    need3 += INT_NCART(am12)*INT_NCART(am34-2);
    }
  if ((am34>1) && (!inthave.i[am12][am34-2][m+1])) {
    need3 += INT_NCART(am12)*INT_NCART(am34-2);
    }
  if (!inthave.i[am12-1][am34-1][m+1]) {
    need3 += INT_NCART(am12-1)*INT_NCART(am34-1);
    }

  if (need1 <= need3) return 1;
  return 3;
  }

LOCAL_FUNCTION int
impossible_integral()
{
  fprintf(stderr,"tried to build a impossible integral\n");
  fail();
  return(0);
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
  fprintf(stderr,"failing module:\n%s\n",rcsid);
  exit(1);
  }
