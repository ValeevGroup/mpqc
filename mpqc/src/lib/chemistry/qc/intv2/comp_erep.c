
/* Compute the two electron repulsion integrals.
 */

/* $Log$
 * Revision 1.2  1993/12/30 13:32:48  etseidl
 * mostly rcs id stuff
 *
 * Revision 1.5  1992/06/17  22:04:37  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.4  1992/05/13  18:29:33  jannsen
 * added bounds checking for derivative integrals
 *
 * Revision 1.3  1992/04/16  17:03:40  jannsen
 * If compiled with -DTIMING the time needed to compute different classes
 * of integrals is tim_enter'ed.
 *
 * Revision 1.2  1992/03/31  01:21:54  jannsen
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.3  1992/01/30  01:29:21  cljanss
 * Permutation information is now stored with stored integrals.
 *
 * Revision 1.2  1992/01/10  17:56:36  cljanss
 * store integrals to be reused, if requested
 *
 * Revision 1.1  1991/12/14  00:19:36  cljanss
 * Initial revision
 * */
/*
 * Revision 1.14  91/11/22  17:48:51  cljanss
 * bound matrix generated is now handled by a separate flag
 * 
 * Revision 1.13  91/11/20  20:04:42  cljanss
 * Added the int_erep_all1der_v interface.
 * 
 * Revision 1.12  91/10/31  14:45:31  cljanss
 * Turned off bounds checking and fixed a bug (omit wasn't initialized)
 * 
 * Revision 1.11  91/09/30  17:05:28  cljanss
 * fixed bug for case where integrals are not stored
 * 
 * Revision 1.10  91/09/28  19:26:48  cljanss
 * Switch to new file naming scheme
 * 
 * Revision 1.9  91/09/28  18:17:14  cljanss
 * Intermediates are no longer stored, if requested with flags.
 * 
 * Revision 1.8  91/09/26  15:53:27  cljanss
 * fixed bug in normalization of derivative integrals
 * 
 * Revision 1.7  1991/09/10  19:34:45  cljanss
 * finished putting in first derivatives
 *
 * Revision 1.6  1991/08/09  17:10:59  cljanss
 * added the int_erep_v interface to int_erep.
 *
 * Revision 1.5  1991/08/09  16:52:52  cljanss
 * more work on derivative integrals done
 *
 * Revision 1.4  1991/06/19  23:18:19  janssen
 * added computation and checking of upper bounds to shell quartet integrals
 *
 * Revision 1.3  1991/06/19  15:22:47  janssen
 * first stab at derivative integrals -- might work might not
 *
 * Revision 1.2  1991/06/16  19:32:22  janssen
 * improved performance when center A == center B
 *
 * Revision 1.1  1991/06/16  16:40:07  janssen
 * Initial revision
 * */
static char *rcsid = "$Id$";

#include <stdio.h>
#include <stdarg.h>
#include <tmpl.h>
#include <math/array/math_lib.h>

#include "atoms.h"
#include "int_macros.h"
#include "int_flags.h"
#include "int_types.h"

#include "inter.h"

#include "comp_erep.gbl"
#include "comp_erep.lcl"

#include "buildgc.gbl"
#include "shiftgc.gbl"
#include "storage.gbl"

/* This returns a quick upperbound for an integral in the given
 * shell quartet. */
GLOBAL_FUNCTION double
int_erep_bound(flags,sh1,sh2,sh3,sh4)
int flags;
int sh1;
int sh2;
int sh3;
int sh4;
{
  int osh1,osh2,osh3,osh4;

  if (!int_store_bounds) {
    fprintf(stderr,"int_erep_bound: bounds not available\n");
    fail();
    }

  /* Compute the offset shell numbers. */
  osh1 = sh1 + int_cs1->shell_offset;
  osh2 = sh2 + int_cs2->shell_offset;
  osh3 = sh3 + int_cs3->shell_offset;
  osh4 = sh4 + int_cs4->shell_offset;

  return int_shell_Q.d[osh1][osh2] * int_shell_Q.d[osh3][osh4];
  }


/* This computes the 2erep integrals for a shell quartet
 * specified by psh1, psh2, psh3, psh4.
 * The routine int_initialize_2erep must be called before
 * any integrals can be computed.
 * This routine may decide to change the shell ordering.
 * The new ordering is placed in *psh1,4 on exit.
 * for the derivatives.
 */
GLOBAL_FUNCTION VOID
int_erep(flags,psh1,psh2,psh3,psh4)
int flags;
int *psh1;
int *psh2;
int *psh3;
int *psh4;
{
  compute_erep(flags,1,psh1,psh2,psh3,psh4,0,0,0,0);
  }

/* This is an alternate interface to int_erep.  It takes
 * as arguments the flags, an integer vector of shell numbers
 * and an integer vector which will be filled in with size
 * information, if it is non-NULL. */
GLOBAL_FUNCTION VOID
int_erep_v(flags,shells,sizes)
int flags;
int *shells;
int  *sizes;
{
  int_erep(flags,&(shells[0]),&(shells[1]),&(shells[2]),&(shells[3]));
  if (sizes) {
    sizes[0] = INT_SH(int_cs1,shells[0]).nfunc;
    sizes[1] = INT_SH(int_cs2,shells[1]).nfunc;
    sizes[2] = INT_SH(int_cs3,shells[2]).nfunc;
    sizes[3] = INT_SH(int_cs4,shells[3]).nfunc;
    }
  }

/* If we need a computation with adjusted angular momentum, then
 * this lower level routine can be called instead of int_erep.
 * The dam{1,2,3,4} arguments given the amount by which the
 * angular momentum is to adjusted.  This differs from libint version
 * 1 in that it used total angular momentum here.
 */
LOCAL_FUNCTION VOID
compute_erep(flags,normalize,psh1,psh2,psh3,psh4,dam1,dam2,dam3,dam4)
int flags;
int normalize;
int *psh1;
int *psh2;
int *psh3;
int *psh4;
int dam1;
int dam2;
int dam3;
int dam4;
{
#ifdef TIMING
  char section[30];
#endif
  centers_t *pcs1=int_cs1,*pcs2=int_cs2,*pcs3=int_cs3,*pcs4=int_cs4;
  int size;
  int ii;
  int size1, size2, size3, size4;
  int tam1,tam2,tam3,tam4;
  int i,j,k,l;
  int ogc1,ogc2,ogc3,ogc4;
  int osh1,osh2,osh3,osh4;
  int sh1,sh2,sh3,sh4;
  int am1,am2,am3,am4,am12, am34;
  int minam1,minam2,minam3,minam4;
  int redundant_index;
  int e12,e13e24,e34;
  int i1,j1,k1;
  int i2,j2,k2;
  int i3,j3,k3;
  int i4,j4,k4;
  int p12,p34,p13p24;
  int eAB;

#if 0
  fprintf(stdout,"compute_erep: dam: (%d %d|%d %d), normalize: %d\n",
          dam1,dam2,dam3,dam4,normalize);
#endif

  /* Compute the offset shell numbers. */
  osh1 = *psh1 + int_cs1->shell_offset;
  osh2 = *psh2 + int_cs2->shell_offset;
  osh3 = *psh3 + int_cs3->shell_offset;
  osh4 = *psh4 + int_cs4->shell_offset;

  sh1 = *psh1;
  sh2 = *psh2;
  sh3 = *psh3;
  sh4 = *psh4;

  /* Test the arguments to make sure that they are sensible. */
  if (   sh1 < 0 || sh1 >= int_cs1->nshell
      || sh2 < 0 || sh2 >= int_cs2->nshell
      || sh3 < 0 || sh3 >= int_cs3->nshell
      || sh4 < 0 || sh4 >= int_cs4->nshell) {
    fprintf(stderr,"compute_erep has been incorrectly used\n");
    fprintf(stderr,"shells (bounds): %d (%d), %d (%d), %d (%d), %d (%d)\n",
            sh1,int_cs1->nshell-1,
            sh2,int_cs2->nshell-1,
            sh3,int_cs3->nshell-1,
            sh4,int_cs4->nshell-1);
    fail();
    }

  /* If bounds checking is not turned off, get an upper bound for
   * an integral from the shell quartet and zero out the integral
   * buffer and return if this bound is small. */
  if (!(flags|INT_NOBCHK)) {
#if defined(NCUBE_V2)
    printf("bounds checking is turned off on the NCUBE_V2\n");
    printf("because a compiler limitation won't allow a\n");
    printf("tiny bit of code to compile correctly\n");
    printf("exiting");
    exit(1);
#else
    double bound = int_erep_bound(flags,sh1,sh2,sh3,sh4);
    if (!INT_NONZERO(bound)) {
      int i,bufsize;
      bufsize = INT_SH(int_cs1,sh1).nfunc
              * INT_SH(int_cs2,sh2).nfunc
              * INT_SH(int_cs3,sh3).nfunc
              * INT_SH(int_cs4,sh4).nfunc;
      for (i=0; i<bufsize; i++) int_buffer[i] = 0.0;
      return;
      }
#endif
    }

  /* Set up pointers to the current shells. */
  int_shell1 = &int_cs1->center[int_cs1->center_num[sh1]]
                  .basis.shell[int_cs1->shell_num[sh1]];
  int_shell2 = &int_cs2->center[int_cs2->center_num[sh2]]
                  .basis.shell[int_cs2->shell_num[sh2]];
  int_shell3 = &int_cs3->center[int_cs3->center_num[sh3]]
                  .basis.shell[int_cs3->shell_num[sh3]];
  int_shell4 = &int_cs4->center[int_cs4->center_num[sh4]]
                  .basis.shell[int_cs4->shell_num[sh4]];


  /* Compute the maximum angular momentum on each centers to
   * determine the most efficient way to invoke the building and shifting
   * routines.  The minimum angular momentum will be computed at the
   * same time. */
  am1 = am2 = am3 = am4 = 0;
  minam1 = int_shell1->type[0].am;
  minam2 = int_shell2->type[0].am;
  minam3 = int_shell3->type[0].am;
  minam4 = int_shell4->type[0].am;
  for (i=0; i<int_shell1->ncon; i++) {
    if (am1<int_shell1->type[i].am) am1 = int_shell1->type[i].am;
    if (minam1>int_shell1->type[i].am) minam1 = int_shell1->type[i].am;
    }
  for (i=0; i<int_shell2->ncon; i++) {
    if (am2<int_shell2->type[i].am) am2 = int_shell2->type[i].am;
    if (minam2>int_shell2->type[i].am) minam2 = int_shell2->type[i].am;
    }
  for (i=0; i<int_shell3->ncon; i++) {
    if (am3<int_shell3->type[i].am) am3 = int_shell3->type[i].am;
    if (minam3>int_shell3->type[i].am) minam3 = int_shell3->type[i].am;
    }
  for (i=0; i<int_shell4->ncon; i++) {
    if (am4<int_shell4->type[i].am) am4 = int_shell4->type[i].am;
    if (minam4>int_shell4->type[i].am) minam4 = int_shell4->type[i].am;
    }

  am1 += dam1; minam1 += dam1;
  am2 += dam2; minam2 += dam2;
  am3 += dam3; minam3 += dam3;
  am4 += dam4; minam4 += dam4;
  am12 = am1 + am2;
  am34 = am3 + am4;

  /* if no angular momentum remains on one of the centers return */
  if (am1 < 0 || am2 < 0 || am3 < 0 || am4 < 0) {
    return;
    }

#ifdef TIMING
  sprintf(section,"erep am=%d",am12+am34);
  tim_enter(section);
#endif

#if 0
  printf("on entry: (%d,%d,%d,%d) am=(%d,%d,%d,%d) perm = %d\n",
         *psh1, *psh2, *psh3, *psh4,
         am1,am2,am3,am4,
         !(INT_NOPERM&flags));
#endif

  /* Convert the integral to the most efficient form. */
  p12 = 0;
  p34 = 0;
  p13p24 = 0;

#undef OLD_PERMUTATION_ALGORITHM
#ifdef OLD_PERMUTATION_ALGORITHM
  if (am2 > am1) {
    p12 = 1;
    iswtch(&am1,&am2);iswtch(&sh1,&sh2);iswtch(psh1,psh2);iswtch(&osh1,&osh2);
    iswtch(&dam1,&dam2);
    iswtch(&minam1,&minam2);
    pswtch((VOID_PTR *)&int_shell1,(VOID_PTR *)&int_shell2);
    pswtch((VOID_PTR *)&pcs1,(VOID_PTR *)&pcs2);
    }
  if (am4 > am3) {
    p34 = 1;
    iswtch(&am3,&am4);iswtch(&sh3,&sh4);iswtch(psh3,psh4);iswtch(&osh3,&osh4);
    iswtch(&dam3,&dam4);
    iswtch(&minam3,&minam4);
    pswtch((VOID_PTR *)&int_shell3,(VOID_PTR *)&int_shell4);
    pswtch((VOID_PTR *)&pcs3,(VOID_PTR *)&pcs4);
    }
  if ((osh1 == osh4) && (osh2 == osh3) && (osh1 != osh2)) {
    /* Don't make the permutation unless we won't override what was
     * decided above about p34. */
    if (am4 == am3) {
      p34 = 1;
      iswtch(&am3,&am4);iswtch(&sh3,&sh4);iswtch(psh3,psh4);iswtch(&osh3,&osh4);
      iswtch(&dam3,&dam4);
      iswtch(&minam3,&minam4);
      pswtch((VOID_PTR *)&int_shell3,(VOID_PTR *)&int_shell4);
      pswtch((VOID_PTR *)&pcs3,(VOID_PTR *)&pcs4);
      }
    }
  if ((am3 > am1)||((am3 == am1)&&(am4 > am2))) {
    p13p24 = 1;
    iswtch(&am1,&am3);iswtch(&sh1,&sh3);iswtch(psh1,psh3);iswtch(&osh1,&osh3);
    iswtch(&am2,&am4);iswtch(&sh2,&sh4);iswtch(psh2,psh4);iswtch(&osh2,&osh4);
    iswtch(&am12,&am34);
    iswtch(&dam1,&dam3);
    iswtch(&minam1,&minam3);
    pswtch((VOID_PTR *)&int_shell1,(VOID_PTR *)&int_shell3);
    pswtch((VOID_PTR *)&pcs1,(VOID_PTR *)&pcs3);
    iswtch(&dam2,&dam4);
    iswtch(&minam2,&minam4);
    pswtch((VOID_PTR *)&int_shell2,(VOID_PTR *)&int_shell4);
    pswtch((VOID_PTR *)&pcs2,(VOID_PTR *)&pcs4);
    }
  /* This tries to make centers A and B equivalent, if possible. */
  else if (  (am3 == am1)
           &&(am4 == am2)
           &&(!(  (int_cs1 == int_cs2)
                &&(int_cs1->center_num[sh1]==int_cs2->center_num[sh2])))
           &&(   (int_cs3 == int_cs4)
               &&(int_cs3->center_num[sh3]==int_cs4->center_num[sh4]))) {
    p13p24 = 1;
    iswtch(&am1,&am3);iswtch(&sh1,&sh3);iswtch(psh1,psh3);iswtch(&osh1,&osh3);
    iswtch(&am2,&am4);iswtch(&sh2,&sh4);iswtch(psh2,psh4);iswtch(&osh2,&osh4);
    iswtch(&am12,&am34);
    iswtch(&dam1,&dam3);
    iswtch(&minam1,&minam3);
    pswtch((VOID_PTR *)&int_shell1,(VOID_PTR *)&int_shell3);
    pswtch((VOID_PTR *)&pcs1,(VOID_PTR *)&pcs3);
    iswtch(&dam2,&dam4);
    iswtch(&minam2,&minam4);
    pswtch((VOID_PTR *)&int_shell2,(VOID_PTR *)&int_shell4);
    pswtch((VOID_PTR *)&pcs2,(VOID_PTR *)&pcs4);
    }
#else /* OLD_PERMUTATION_ALGORITHM */
  if (am2 > am1) {
    p12 = 1;
    iswtch(&am1,&am2);iswtch(&sh1,&sh2);iswtch(psh1,psh2);iswtch(&osh1,&osh2);
    iswtch(&dam1,&dam2);
    iswtch(&minam1,&minam2);
    pswtch((VOID_PTR *)&int_shell1,(VOID_PTR *)&int_shell2);
    pswtch((VOID_PTR *)&pcs1,(VOID_PTR *)&pcs2);
    }
  if (am4 > am3) {
    p34 = 1;
    iswtch(&am3,&am4);iswtch(&sh3,&sh4);iswtch(psh3,psh4);iswtch(&osh3,&osh4);
    iswtch(&dam3,&dam4);
    iswtch(&minam3,&minam4);
    pswtch((VOID_PTR *)&int_shell3,(VOID_PTR *)&int_shell4);
    pswtch((VOID_PTR *)&pcs3,(VOID_PTR *)&pcs4);
    }
  if ((osh1 == osh4) && (osh2 == osh3) && (osh1 != osh2)) {
    /* Don't make the permutation unless we won't override what was
     * decided above about p34. */
    if (am4 == am3) {
      p34 = 1;
      iswtch(&am3,&am4);iswtch(&sh3,&sh4);iswtch(psh3,psh4);iswtch(&osh3,&osh4);
      iswtch(&dam3,&dam4);
      iswtch(&minam3,&minam4);
      pswtch((VOID_PTR *)&int_shell3,(VOID_PTR *)&int_shell4);
      pswtch((VOID_PTR *)&pcs3,(VOID_PTR *)&pcs4);
      }
    }
  if ((am34 > am12)||((am34 == am12)&&(minam1 > minam3))) {
    p13p24 = 1;
    iswtch(&am1,&am3);iswtch(&sh1,&sh3);iswtch(psh1,psh3);iswtch(&osh1,&osh3);
    iswtch(&am2,&am4);iswtch(&sh2,&sh4);iswtch(psh2,psh4);iswtch(&osh2,&osh4);
    iswtch(&am12,&am34);
    iswtch(&dam1,&dam3);
    iswtch(&minam1,&minam3);
    pswtch((VOID_PTR *)&int_shell1,(VOID_PTR *)&int_shell3);
    pswtch((VOID_PTR *)&pcs1,(VOID_PTR *)&pcs3);
    iswtch(&dam2,&dam4);
    iswtch(&minam2,&minam4);
    pswtch((VOID_PTR *)&int_shell2,(VOID_PTR *)&int_shell4);
    pswtch((VOID_PTR *)&pcs2,(VOID_PTR *)&pcs4);
    }
  /* This tries to make centers A and B equivalent, if possible. */
  else if (  (am3 == am1)
           &&(am4 == am2)
           &&(minam1 == minam3)
           &&(!(  (int_cs1 == int_cs2)
                &&(int_cs1->center_num[sh1]==int_cs2->center_num[sh2])))
           &&(   (int_cs3 == int_cs4)
               &&(int_cs3->center_num[sh3]==int_cs4->center_num[sh4]))) {
    p13p24 = 1;
    iswtch(&am1,&am3);iswtch(&sh1,&sh3);iswtch(psh1,psh3);iswtch(&osh1,&osh3);
    iswtch(&am2,&am4);iswtch(&sh2,&sh4);iswtch(psh2,psh4);iswtch(&osh2,&osh4);
    iswtch(&am12,&am34);
    iswtch(&dam1,&dam3);
    iswtch(&minam1,&minam3);
    pswtch((VOID_PTR *)&int_shell1,(VOID_PTR *)&int_shell3);
    pswtch((VOID_PTR *)&pcs1,(VOID_PTR *)&pcs3);
    iswtch(&dam2,&dam4);
    iswtch(&minam2,&minam4);
    pswtch((VOID_PTR *)&int_shell2,(VOID_PTR *)&int_shell4);
    pswtch((VOID_PTR *)&pcs2,(VOID_PTR *)&pcs4);
    }
#endif /* OLD_PERMUTATION_ALGORITHM */

  if (  (pcs1 == pcs2)
      &&(pcs1->center_num[sh1]==pcs2->center_num[sh2])) {
    eAB = 1;
    }
  else {
    eAB = 0;
    }

#if 0
  printf("perm info: p12 = %d, p34 = %d, p13p24 = %d, eAB = %d\n",
         p12,p34,p13p24,eAB);
#endif

  /* If the centers were permuted, then the int_expweighted variable may
   * need to be changed. */
  if (p12) {
    iswtch(&int_expweight1,&int_expweight2);
    }
  if (p34) {
    iswtch(&int_expweight3,&int_expweight4);
    }
  if (p13p24) {
    iswtch(&int_expweight1,&int_expweight3);
    iswtch(&int_expweight2,&int_expweight4);
    }

  /* Compute the shell sizes. */
  for (ii=size1=0; ii<int_shell1->ncon; ii++) 
    size1 += INT_NCART(int_shell1->type[ii].am+dam1);
  for (ii=size2=0; ii<int_shell2->ncon; ii++) 
    size2 += INT_NCART(int_shell2->type[ii].am+dam2);
  for (ii=size3=0; ii<int_shell3->ncon; ii++) 
    size3 += INT_NCART(int_shell3->type[ii].am+dam3);
  for (ii=size4=0; ii<int_shell4->ncon; ii++) 
    size4 += INT_NCART(int_shell4->type[ii].am+dam4);
  size = size1*size2*size3*size4;

  if (int_integral_storage && (size >= int_storage_threshold)) {
    if (dam1 || dam2 || dam3 || dam4) {
      fprintf(stderr,"cannot use integral storage and dam\n");
      fail();
      }
    if (int_have_stored_integral(sh1,sh2,sh3,sh4,p12,p34,p13p24))
      goto post_computation;
    }

  /* Buildam up on center 1 and 3. */
#if 0
  printf("C:");
#endif
  int_buildgcam(minam1,minam2,minam3,minam4,
                am1,am2,am3,am4,
                dam1,dam2,dam3,dam4,
                sh1,sh2,sh3,sh4, eAB);

  /* Begin loop over generalized contractions. */
  ogc1 = 0;
  for (i=0; i<int_shell1->ncon; i++) {
    tam1 = int_shell1->type[i].am + dam1;
    if (tam1 < 0) continue;
    ogc2 = 0;
    for (j=0; j<int_shell2->ncon; j++) {
      tam2 = int_shell2->type[j].am + dam2;
      if (tam2 < 0) continue;
      ogc3 = 0;
      for (k=0; k<int_shell3->ncon; k++) {
        tam3 = int_shell3->type[k].am + dam3;
        if (tam3 < 0) continue;
        ogc4 = 0;
        for (l=0; l<int_shell4->ncon; l++) {
          tam4 = int_shell4->type[l].am + dam4;
          if (tam4 < 0) continue;

  /* Shift angular momentum from 1 to 2 and from 3 to 4. */
  int_shiftgcam(i,j,k,l,tam1,tam2,tam3,tam4, eAB);

#if 0
  fprintf(stdout,"Top source int_con_ints_array[%d][%d][%d][%d]:\n",
          i,j,k,l);
  fprintf(stdout,"am = (%d,%d,%d,%d)  n = (%d,%d,%d,%d)\n",
         tam1+tam2,
         0,
         tam3+tam4,
         0,
         INT_NCART(tam1 + tam2),
         1,
         INT_NCART(tam3 + tam4),
         1
         );
  int_print_n(stdout,int_con_ints_array[i][j][k][l].dp
               [tam1 +tam2][0][tam3 +tam4][0],
               INT_NCART(tam1+tam2),
               1,
               INT_NCART(tam3+tam4),
               1,
               0,0,0);
  fprintf(stdout,"Targets:\n");
  fprintf(stdout,"am = (%d,%d,%d,%d)  n = (%d,%d,%d,%d)\n",
         tam1,
         tam2,
         tam3,
         tam4,
         INT_NCART(tam1),
         INT_NCART(tam2),
         INT_NCART(tam3),
         INT_NCART(tam4));
  int_print_n(stdout,int_con_ints_array[i][j][k][l].dp
                     [tam1]
                     [tam2]
                     [tam3]
                     [tam4],
               INT_NCART(tam1),
               INT_NCART(tam2),
               INT_NCART(tam3),
               INT_NCART(tam4),
               0,0,0);
#endif

  /* Normalize the integrals if the normalize flag is set to true. */
  if (normalize) {
    normalize_erep_given_gc(int_con_ints_array[i][j][k][l].dp
                     [tam1]
                     [tam2]
                     [tam3]
                     [tam4],
                   i,j,k,l,
                   int_shell1,int_shell2,int_shell3,int_shell4);
    }

#if 0
  fprintf(stdout,"Normalized targets:\n");
  fprintf(stdout,"am = (%d,%d,%d,%d)\n",
         tam1,
         tam2,
         tam3,
         tam4);
  int_print_n(stdout,int_con_ints_array[i][j][k][l].dp
                     [tam1]
                     [tam2]
                     [tam3]
                     [tam4],
               INT_NCART(tam1),
               INT_NCART(tam2),
               INT_NCART(tam3),
               INT_NCART(tam4),
               0,0,0);
#endif

  /* Place the integrals in the integral buffer. */
  /* If INT_NOPERM is set, then repack the integrals while copying. */
  if ((flags & INT_NOPERM)&&(p12||p34||p13p24)) {
    int newindex,pam1,pam2,pam3,pam4;
    int pi1,pj1,pk1;
    int pi2,pj2,pk2;
    int pi3,pj3,pk3;
    int pi4,pj4,pk4;
    int psize234,psize34;
    int pogc1,pogc2,pogc3,pogc4;
    int psize1,psize2,psize3,psize4;
    pam1 = tam1;
    pam2 = tam2;
    pam3 = tam3;
    pam4 = tam4;
    pogc1 = ogc1;
    pogc2 = ogc2;
    pogc3 = ogc3;
    pogc4 = ogc4;
    psize1 = size1;
    psize2 = size2;
    psize3 = size3;
    psize4 = size4;
    if (p13p24) {
      iswtch(&pam1,&pam3);
      iswtch(&pam2,&pam4);
      iswtch(&pogc1,&pogc3);
      iswtch(&pogc2,&pogc4);
      iswtch(&psize1,&psize3);
      iswtch(&psize2,&psize4);
      }
    if (p34) {
      iswtch(&pam3,&pam4);
      iswtch(&pogc3,&pogc4);
      iswtch(&psize3,&psize4);
      }
    if (p12) {
      iswtch(&pam1,&pam2);
      iswtch(&pogc1,&pogc2);
      iswtch(&psize1,&psize2);
      }
    psize34 = psize4 * psize3;
    psize234 = psize34 * psize2;
    redundant_index = 0;
    FOR_CART(i1,j1,k1,tam1)
      FOR_CART(i2,j2,k2,tam2)
        FOR_CART(i3,j3,k3,tam3)
          FOR_CART(i4,j4,k4,tam4)
            pi1=i1; pj1=j1; pk1=k1;
            pi2=i2; pj2=j2; pk2=k2;
            pi3=i3; pj3=j3; pk3=k3;
            pi4=i4; pj4=j4; pk4=k4;
            if (p13p24) {
              iswtch(&pi1,&pi3);
              iswtch(&pj1,&pj3);
              iswtch(&pk1,&pk3);
              iswtch(&pi2,&pi4);
              iswtch(&pj2,&pj4);
              iswtch(&pk2,&pk4);
              }
            if (p34) {
              iswtch(&pi3,&pi4);
              iswtch(&pj3,&pj4);
              iswtch(&pk3,&pk4);
              }
            if (p12) {
              iswtch(&pi1,&pi2);
              iswtch(&pj1,&pj2);
              iswtch(&pk1,&pk2);
              }
            newindex =  (pogc1 + INT_CARTINDEX(pam1,pi1,pj1)) * psize234
                      + (pogc2 + INT_CARTINDEX(pam2,pi2,pj2)) * psize34
                      + (pogc3 + INT_CARTINDEX(pam3,pi3,pj3)) * psize4
                      + (pogc4 + INT_CARTINDEX(pam4,pi4,pj4));
#if 0
   printf("perm: int_buffer[%3d] += % f\n",newindex,
          int_con_ints_array[i][j][k][l]
                 .dp[tam1][tam2][tam3][tam4][redundant_index]);
#endif
            int_buffer[newindex]
              = int_con_ints_array[i][j][k][l]
                 .dp[tam1][tam2][tam3][tam4][redundant_index];
            redundant_index++;
            END_FOR_CART
          END_FOR_CART
        END_FOR_CART
      END_FOR_CART
    }
  else {
    int newindex;
    int size34 =  size3 * size4;
    int size234 = size2 * size34;
    redundant_index = 0;
#if 0
    printf("Using int_con_ints_array[%d][%d][%d][%d].dp[%d][%d][%d][%d]\n",
           i,j,k,l,tam1,tam2,tam3,tam4);
#endif
    FOR_CART(i1,j1,k1,tam1)
      FOR_CART(i2,j2,k2,tam2)
        FOR_CART(i3,j3,k3,tam3)
          FOR_CART(i4,j4,k4,tam4)
            newindex =  (ogc1 + INT_CARTINDEX(tam1,i1,j1)) * size234
                      + (ogc2 + INT_CARTINDEX(tam2,i2,j2)) * size34
                      + (ogc3 + INT_CARTINDEX(tam3,i3,j3)) * size4
                      + (ogc4 + INT_CARTINDEX(tam4,i4,j4));
#if 0
   printf("perm: int_buffer[%3d] += % f\n",newindex,
          int_con_ints_array[i][j][k][l]
                 .dp[tam1][tam2][tam3][tam4][redundant_index]);
#endif
            int_buffer[newindex]
              = int_con_ints_array[i][j][k][l]
                 .dp[tam1][tam2][tam3][tam4][redundant_index];
            redundant_index++;
            END_FOR_CART
          END_FOR_CART
        END_FOR_CART
      END_FOR_CART
    }

    /* End loop over generalized contractions. */
          ogc4 += INT_NCART(tam4);
          }
        ogc3 += INT_NCART(tam3);
        }
      ogc2 += INT_NCART(tam2);
      }
    ogc1 += INT_NCART(tam1);
    }
#if 0
  printf("\n");
#endif

  if (   int_integral_storage
      && (size >= int_storage_threshold)) {
    if (int_integral_storage>=size+int_used_integral_storage) {
      int_store_integral(sh1,sh2,sh3,sh4,p12,p34,p13p24,size);
      }
    }

  /* We branch here if an integral was precomputed and the int_buffer
   * has been already filled. */
  post_computation:

#if 0
  printf("before unpermute: am=(%d,%d,%d,%d)\n",am1,am2,am3,am4);
#endif

  /* Unpermute all of the permuted quantities. */
  if ((flags & INT_NOPERM)&&(p12||p34||p13p24)) {
    if (p13p24) {
      iswtch(&sh1,&sh3);iswtch(psh1,psh3);iswtch(&osh1,&osh3);
      iswtch(&sh2,&sh4);iswtch(psh2,psh4);iswtch(&osh2,&osh4);
      iswtch(&am1,&am3);
      iswtch(&am2,&am4);
      iswtch(&am12,&am34);
      pswtch((VOID_PTR *)&int_shell1,(VOID_PTR *)&int_shell3);
      pswtch((VOID_PTR *)&pcs1,(VOID_PTR *)&pcs3);
      pswtch((VOID_PTR *)&int_shell2,(VOID_PTR *)&int_shell4);
      pswtch((VOID_PTR *)&pcs2,(VOID_PTR *)&pcs4);
      iswtch(&int_expweight1,&int_expweight3);
      iswtch(&int_expweight2,&int_expweight4);
      }
    if (p34) {
      iswtch(&sh3,&sh4);iswtch(psh3,psh4);iswtch(&osh3,&osh4);
      iswtch(&am3,&am4);
      pswtch((VOID_PTR *)&int_shell3,(VOID_PTR *)&int_shell4);
      pswtch((VOID_PTR *)&pcs3,(VOID_PTR *)&pcs4);
      iswtch(&int_expweight3,&int_expweight4);
      }
    if (p12) {
      iswtch(&sh1,&sh2);iswtch(psh1,psh2);iswtch(&osh1,&osh2);
      iswtch(&am1,&am2);
      pswtch((VOID_PTR *)&int_shell1,(VOID_PTR *)&int_shell2);
      pswtch((VOID_PTR *)&pcs1,(VOID_PTR *)&pcs2);
      iswtch(&int_expweight1,&int_expweight2);
      }
    }

  /* Remove the redundant integrals, unless INT_REDUND is specified. */
  if (!(flags&INT_REDUND)) {
    int redundant_offset = 0;
    int nonredundant_offset = 0;
    if ((osh1 == osh4)&&(osh2 == osh3)&&(osh1 != osh2)) {
      fprintf(stderr,"nonredundant integrals cannot be generated\n");
      fail();
      }
    e12 = (osh1 == osh2);
    e13e24 = ((osh1 == osh3) && (osh2 == osh4));
    e34 = (osh3 == osh4);
    nonredundant_erep(int_buffer,e12,e34,e13e24,
                           int_shell1->nfunc,
                           int_shell2->nfunc,
                           int_shell3->nfunc,
                           int_shell4->nfunc,
                           &redundant_offset,
                           &nonredundant_offset);
    }
    
#ifdef TIMING
  tim_exit(section);
#endif

#if 0
  printf("on exit: am=(%d,%d,%d,%d)\n",am1,am2,am3,am4);
#endif
  }

/* This computes the one electron derivatives for all unique
 * centers in the passed shell quartet.  One center in
 * the set of unique centers is not included.  This can
 * be computed as minus the sum of the other derivatives.
 * The list of centers for which integrals were computed can
 * be determined from the contents of der_centers.
 * The results are put into the global integral buffer in the
 * format:
 * +------------------+
 * | dercenter1 +---+ |
 * |            | x | |
 * |            +---+ |
 * |            | y | |
 * |            +---+ |
 * |            | z | |
 * |            +---+ |
 * +------------------+
 * | dercenter2 +---+ |
 * |            | x | |
 * |            +---+ |
 * |            | y | |
 * |            +---+ |
 * |            | z | |
 * |            +---+ |
 * +------------------+
 * | dercenter3 +---+ |
 * |            | x | |
 * |            +---+ |
 * |            | y | |
 * |            +---+ |
 * |            | z | |
 * |            +---+ |
 * +------------------+
 */

GLOBAL_FUNCTION VOID
int_erep_all1der(flags,psh1,psh2,psh3,psh4,der_centers)
int flags;
int *psh1;
int *psh2;
int *psh3;
int *psh4;
der_centers_t *der_centers;
{
  double *current_buffer;
  int nints;
  double *user_int_buffer;
  int omit;
  centers_t *cs[4];
  int sh[4];
  int n_unique;
  int i,j;
  shell_t *shell1,*shell2,*shell3,*shell4;
  centers_t *ucs[4];  /* The centers struct for the unique centers. */
  int ush[4];         /* The shells for the unique centers. */
  int unum[4];        /* The number of times that this unique center occurs. */
  int uam[4];         /* The total angular momentum on each unique center. */
  int am[4];
  int osh[4];

  cs[0] = int_cs1;
  cs[1] = int_cs2;
  cs[2] = int_cs3;
  cs[3] = int_cs4;

  sh[0] = *psh1;
  sh[1] = *psh2;
  sh[2] = *psh3;
  sh[3] = *psh4;

  /* Set up pointers to the current shells. */
  shell1 = &int_cs1->center[int_cs1->center_num[*psh1]]
                  .basis.shell[int_cs1->shell_num[*psh1]];
  shell2 = &int_cs2->center[int_cs2->center_num[*psh2]]
                  .basis.shell[int_cs2->shell_num[*psh2]];
  shell3 = &int_cs3->center[int_cs3->center_num[*psh3]]
                  .basis.shell[int_cs3->shell_num[*psh3]];
  shell4 = &int_cs4->center[int_cs4->center_num[*psh4]]
                  .basis.shell[int_cs4->shell_num[*psh4]];

  am[0] = int_find_jmax_shell(shell1);
  am[1] = int_find_jmax_shell(shell2);
  am[2] = int_find_jmax_shell(shell3);
  am[3] = int_find_jmax_shell(shell4);

  /* Compute the offset shell numbers. */
  osh[0] = *psh1 + int_cs1->shell_offset;
  osh[1] = *psh2 + int_cs2->shell_offset;
  osh[2] = *psh3 + int_cs3->shell_offset;
  osh[3] = *psh4 + int_cs4->shell_offset;

  /* This macro returns true if two shell centers are the same. */
#define SC(cs1,sh1,cs2,sh2) (((cs1)==(cs2))&&((cs1)->center_num[sh1]==(cs1)->center_num[sh2]))

  /* Build the list of unique centers structures and shells. */
  n_unique = 0;
  for (i=0; i<4; i++) {
    int unique = 1;
    for (j=0; j<n_unique; j++) {
      if (SC(ucs[j],ush[j],cs[i],sh[i])) {
        unique = 0;
        uam[j] += am[i];
        unum[j]++;
        break;
        }
      }
    if (unique) {
      ucs[n_unique] = cs[i];
      ush[n_unique] = sh[i];
      uam[n_unique] = am[i];
      unum[n_unique] = 1;
      n_unique++;
      }
    }

  /* Choose which unique center is to be omitted from the calculation. */
  if (n_unique == 1) {
    omit = 0;
    }
  else if (n_unique == 2) {
    if (unum[0]==3) omit = 0;
    else if (unum[1]==3) omit = 1;
    else if (uam[1]>uam[0]) omit = 1;
    else omit = 0;
    }
  else if (n_unique == 3) {
    if (unum[0]==2) omit = 0;
    else if (unum[1]==2) omit = 1;
    else omit = 2;
    }
  else {
    int max = 0;
    omit = 0;
    for (i=0; i<4; i++) {
      if (uam[i]>max) {
        max = uam[i];
        omit = i;
        }
      }
    }

  /* Save the location of the int_buffer. */
  user_int_buffer = int_buffer;
  int_buffer = int_derint_buffer;

  /* Zero out the result integrals. */
  nints = shell1->nfunc * shell2->nfunc * shell3->nfunc * shell4->nfunc;
  for (i=0; i<3*(n_unique-1)*nints; i++) user_int_buffer[i] = 0.0;

  /* Loop thru the unique centers, computing the integrals and
   * skip the derivative on the unique center specified by omit. */
  der_centers->n = 0;
  current_buffer = user_int_buffer;
  for (i=0; i<n_unique; i++) {
    if (i==omit) continue;

    der_centers->cs[der_centers->n] = ucs[i];
    der_centers->num[der_centers->n] = ucs[i]->center_num[ush[i]];
    der_centers->n++;

    for (j=0; j<4; j++) {
      if (SC(ucs[i],ush[i],cs[j],sh[j])) {
        compute_erep_1der(flags|INT_NOPERM,current_buffer,
                          psh1,psh2,psh3,psh4,j);
        }
      }

    current_buffer = &current_buffer[3*nints];
    }

  /* Put the information about the omitted center into der_centers. */
  der_centers->ocs = ucs[omit];
  der_centers->onum = ucs[omit]->center_num[ush[omit]];

  /* Normalize the integrals. */
  current_buffer = user_int_buffer;
  for (i=0; i<3*der_centers->n; i++) {
    normalize_erep(current_buffer,shell1,shell2,shell3,shell4);
    current_buffer = &current_buffer[nints];
    }

  /* Eliminate redundant integrals, unless flags specifies otherwise. */
  current_buffer = user_int_buffer;
  if (!(flags&INT_REDUND)) {
    int redundant_offset = 0;
    int nonredundant_offset = 0;
    int e12,e13e24,e34;
    int i;

    if ((osh[0] == osh[3])&&(osh[1] == osh[2])&&(osh[0] != osh[1])) {
      fprintf(stderr,"nonredundant integrals cannot be generated (1der)\n");
      fail();
      }

    /* Shell equivalence information. */
    e12 = (osh[0] == osh[1]);
    e13e24 = ((osh[0] == osh[2]) && (osh[1] == osh[3]));
    e34 = (osh[2] == osh[3]);
    /* Repack x, y, and z integrals. */
    for (i=0; i<3*der_centers->n; i++) {
      nonredundant_erep(current_buffer,e12,e34,e13e24,
                             shell1->nfunc,
                             shell2->nfunc,
                             shell3->nfunc,
                             shell4->nfunc,
                             &redundant_offset,
                             &nonredundant_offset);
      }
    }

  /* Return the integral buffers to their normal state. */
  int_derint_buffer = int_buffer;
  int_buffer = user_int_buffer;
  }

/* This will call compute_erep
 * to compute the derivatives in terms of order == 0 integrals.
 * flags are the flags to the int_comp_erep routine
 * psh1-4 are pointers to the shell numbers
 * dercenter is 0, 1, 2, or 3 -- the center within the integral
 *           which we are taking the derivative with respect to.
 * The results are accumulated in buffer, which cannot be the same
 * as the current int_buffer.
 */
LOCAL_FUNCTION VOID
compute_erep_1der(flags,buffer,psh1,psh2,psh3,psh4,dercenter)
int flags;
double *buffer;
int *psh1;
int *psh2;
int *psh3;
int *psh4;
int dercenter;
{
  int oc1,oc2,oc3,oc4;
  int ii;
  int c1,c2,c3,c4;
  int osh[4];
  int i[4],j[4],k[4],am[4];
  int index;
  int sizem234,sizem34,sizem2,sizem3,sizem4;
  int sizep234,sizep34,sizep2,sizep3,sizep4;
  shell_t *shell1,*shell2,*shell3,*shell4;

  /* Compute the offset shell numbers. */
  osh[0] = *psh1 + int_cs1->shell_offset;
  osh[1] = *psh2 + int_cs2->shell_offset;
  osh[2] = *psh3 + int_cs3->shell_offset;
  osh[3] = *psh4 + int_cs4->shell_offset;

  /* Set up pointers to the current shells. */
  shell1 = &int_cs1->center[int_cs1->center_num[*psh1]]
                  .basis.shell[int_cs1->shell_num[*psh1]];
  shell2 = &int_cs2->center[int_cs2->center_num[*psh2]]
                  .basis.shell[int_cs2->shell_num[*psh2]];
  shell3 = &int_cs3->center[int_cs3->center_num[*psh3]]
                  .basis.shell[int_cs3->shell_num[*psh3]];
  shell4 = &int_cs4->center[int_cs4->center_num[*psh4]]
                  .basis.shell[int_cs4->shell_num[*psh4]];

  if ((dercenter<0) || (dercenter > 3)) {
    fprintf(stderr,"illegal derivative center -- must be 0, 1, 2, or 3\n");
    fail();
    }

#define DCTEST(n) ((dercenter==n)?1:0)
  /* Offsets for the intermediates with angular momentum decremented. */
  for (ii=sizem2=0; ii<shell2->ncon; ii++) 
    sizem2 += INT_NCART(shell2->type[ii].am-DCTEST(1));
  for (ii=sizem3=0; ii<shell3->ncon; ii++) 
    sizem3 += INT_NCART(shell3->type[ii].am-DCTEST(2));
  for (ii=sizem4=0; ii<shell4->ncon; ii++) 
    sizem4 += INT_NCART(shell4->type[ii].am-DCTEST(3));
  sizem34 = sizem4 * sizem3;
  sizem234 = sizem34 * sizem2;

  /* Offsets for the intermediates with angular momentum incremented. */
  for (ii=sizep2=0; ii<shell2->ncon; ii++) 
    sizep2 += INT_NCART(shell2->type[ii].am+DCTEST(1));
  for (ii=sizep3=0; ii<shell3->ncon; ii++) 
    sizep3 += INT_NCART(shell3->type[ii].am+DCTEST(2));
  for (ii=sizep4=0; ii<shell4->ncon; ii++) 
    sizep4 += INT_NCART(shell4->type[ii].am+DCTEST(3));
  sizep34 = sizep4 * sizep3;
  sizep234 = sizep34 * sizep2;

  compute_erep(flags|INT_NOPERM|INT_REDUND|INT_NOBCHK,0,psh1,psh2,psh3,psh4,
                   -DCTEST(0),
                   -DCTEST(1),
                   -DCTEST(2),
                   -DCTEST(3));

  /* Trouble if cpp is nonANSI. */
#define DERLOOP(index,indexp1) \
   oc##indexp1 = 0;\
   for ( c##indexp1 =0; c##indexp1 <shell##indexp1->ncon; c##indexp1 ++) {\
     am[index] = shell##indexp1->type[c##indexp1].am;\
     FOR_CART(i[index],j[index],k[index],am[index])

#define END_DERLOOP(index,indexp1,sign) \
       END_FOR_CART\
     oc##indexp1 += INT_NCART(am[index] sign DCTEST(index));\
     }

#define ALLDERLOOPS\
    DERLOOP(0,1)\
      DERLOOP(1,2)\
        DERLOOP(2,3)\
          DERLOOP(3,4)

#define END_ALLDERLOOPS(sign)\
            END_DERLOOP(3,4,sign)\
          END_DERLOOP(2,3,sign)\
        END_DERLOOP(1,2,sign)\
      END_DERLOOP(0,1,sign)

   /* Place the contributions into the user integral buffer. */
   index = 0;
   /* The d/dx integrals */
  ALLDERLOOPS
    if (i[dercenter]>0) {
      buffer[index] -= i[dercenter] * int_buffer[
        (oc1 + INT_CARTINDEX(am[0]-DCTEST(0),i[0]-DCTEST(0),j[0])) * sizem234
       +(oc2 + INT_CARTINDEX(am[1]-DCTEST(1),i[1]-DCTEST(1),j[1])) * sizem34
       +(oc3 + INT_CARTINDEX(am[2]-DCTEST(2),i[2]-DCTEST(2),j[2])) * sizem4
       +(oc4 + INT_CARTINDEX(am[3]-DCTEST(3),i[3]-DCTEST(3),j[3]))
       ];
      }
    index++;
    END_ALLDERLOOPS(-)

   /* The d/dy integrals */
  ALLDERLOOPS
    if (j[dercenter]>0) {
    buffer[index] -= j[dercenter] * int_buffer[
         (oc1 + INT_CARTINDEX(am[0]-DCTEST(0),i[0],j[0]-DCTEST(0))) * sizem234
        +(oc2 + INT_CARTINDEX(am[1]-DCTEST(1),i[1],j[1]-DCTEST(1))) * sizem34
        +(oc3 + INT_CARTINDEX(am[2]-DCTEST(2),i[2],j[2]-DCTEST(2))) * sizem4
        +(oc4 + INT_CARTINDEX(am[3]-DCTEST(3),i[3],j[3]-DCTEST(3)))
        ];
      }
    index++;
  END_ALLDERLOOPS(-)

   /* The d/dz integrals */
  ALLDERLOOPS
    if (k[dercenter]>0) {
    buffer[index] -= k[dercenter] * int_buffer[
         (oc1 + INT_CARTINDEX(am[0]-DCTEST(0),i[0],j[0])) * sizem234
        +(oc2 + INT_CARTINDEX(am[1]-DCTEST(1),i[1],j[1])) * sizem34
        +(oc3 + INT_CARTINDEX(am[2]-DCTEST(2),i[2],j[2])) * sizem4
        +(oc4 + INT_CARTINDEX(am[3]-DCTEST(3),i[3],j[3]))
        ];
      }
    index++;
  END_ALLDERLOOPS(-)

  /* Compute the next contribution to the integrals. */
  /* Tell the build routine that we need an exponent weighted contraction
   * with the exponents taken from the dercenter and adjust the
   * angular momentum of dercenter to the needed value. */
  if (dercenter==0) int_expweight1 = 1;
  else if (dercenter==1) int_expweight2 = 1;
  else if (dercenter==2) int_expweight3 = 1;
  else if (dercenter==3) int_expweight4 = 1;
  compute_erep(flags|INT_NOPERM|INT_REDUND|INT_NOBCHK,0,psh1,psh2,psh3,psh4,
                     DCTEST(0),
                     DCTEST(1),
                     DCTEST(2),
                     DCTEST(3));
  if (dercenter==0) int_expweight1 = 0;
  else if (dercenter==1) int_expweight2 = 0;
  else if (dercenter==2) int_expweight3 = 0;
  else if (dercenter==3) int_expweight4 = 0;

  /* Place the contributions into the user integral buffer. */
  index = 0;
  /* The d/dx integrals */
  ALLDERLOOPS
#if 0
    tmp =int_buffer[
             (oc1+INT_CARTINDEX(am[0]+DCTEST(0),i[0]+DCTEST(0),j[0]))*sizep234
            +(oc2+INT_CARTINDEX(am[1]+DCTEST(1),i[1]+DCTEST(1),j[1]))*sizep34
            +(oc3+INT_CARTINDEX(am[2]+DCTEST(2),i[2]+DCTEST(2),j[2]))*sizep4
            +(oc4+INT_CARTINDEX(am[3]+DCTEST(3),i[3]+DCTEST(3),j[3]))
            ];
    if (1 || INT_NONZERO(tmp)) {
      printf("x: ((%d%d%d)(%d%d%d)(%d%d%d)(%d%d%d)) += ((%d%d%d)(%d%d%d)(%d%d%d)(%d%d%d)) (%f) (from %5d)\n",
             i[0],j[0],k[0],
             i[1],j[1],k[1],
             i[2],j[2],k[2],
             i[3],j[3],k[3],
             i[0]+DCTEST(0),j[0],k[0],
             i[1]+DCTEST(1),j[1],k[1],
             i[2]+DCTEST(2),j[2],k[2],
             i[3]+DCTEST(3),j[3],k[3],
             tmp,
             (oc1+INT_CARTINDEX(am[0]+DCTEST(0),i[0]+DCTEST(0),j[0]))*sizep234
            +(oc2+INT_CARTINDEX(am[1]+DCTEST(1),i[1]+DCTEST(1),j[1]))*sizep34
            +(oc3+INT_CARTINDEX(am[2]+DCTEST(2),i[2]+DCTEST(2),j[2]))*sizep4
            +(oc4+INT_CARTINDEX(am[3]+DCTEST(3),i[3]+DCTEST(3),j[3]))
             );
       }
#endif
          buffer[index] += int_buffer[
             (oc1+INT_CARTINDEX(am[0]+DCTEST(0),i[0]+DCTEST(0),j[0]))*sizep234
            +(oc2+INT_CARTINDEX(am[1]+DCTEST(1),i[1]+DCTEST(1),j[1]))*sizep34
            +(oc3+INT_CARTINDEX(am[2]+DCTEST(2),i[2]+DCTEST(2),j[2]))*sizep4
            +(oc4+INT_CARTINDEX(am[3]+DCTEST(3),i[3]+DCTEST(3),j[3]))
            ];
    index++;
    END_ALLDERLOOPS(+)

  /* The d/dy integrals */
  ALLDERLOOPS
          buffer[index] += int_buffer[
             (oc1+INT_CARTINDEX(am[0]+DCTEST(0),i[0],j[0]+DCTEST(0)))*sizep234
            +(oc2+INT_CARTINDEX(am[1]+DCTEST(1),i[1],j[1]+DCTEST(1)))*sizep34
            +(oc3+INT_CARTINDEX(am[2]+DCTEST(2),i[2],j[2]+DCTEST(2)))*sizep4
            +(oc4+INT_CARTINDEX(am[3]+DCTEST(3),i[3],j[3]+DCTEST(3)))
            ];
          index++;
    END_ALLDERLOOPS(+)

  /* The d/dz integrals */
  ALLDERLOOPS
          buffer[index] += int_buffer[
               (oc1 + INT_CARTINDEX(am[0]+DCTEST(0),i[0],j[0])) * sizep234
              +(oc2 + INT_CARTINDEX(am[1]+DCTEST(1),i[1],j[1])) * sizep34
              +(oc3 + INT_CARTINDEX(am[2]+DCTEST(2),i[2],j[2])) * sizep4
              +(oc4 + INT_CARTINDEX(am[3]+DCTEST(3),i[3],j[3]))
              ];
          index++;
    END_ALLDERLOOPS(+)
  }

LOCAL_FUNCTION VOID
pswtch(i,j)
VOID_PTR *i;
VOID_PTR *j;
{
  VOID_PTR tmp;

  tmp = *i;
  *i = *j;
  *j = tmp;
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

LOCAL_FUNCTION VOID
nonredundant_erep(buffer,e12,e34,e13e24,n1,n2,n3,n4,red_off,nonred_off)
double *buffer;
int e12;
int e34;
int e13e24;
int n1;
int n2;
int n3;
int n4;
int *red_off;
int *nonred_off;
{
  int redundant_index;
  int nonredundant_index;
  int i,j,k,l;

#ifdef I860
  if(mynode0() == -1 ) { printf("junk\n");}
#endif
  redundant_index = *red_off;
  nonredundant_index = *nonred_off;
  for (i=0; i<n1; i++) {
    for (j=0; j<n2; j++) {
      for (k=0; k<n3; k++) {
        for (l=0; l<n4; l++) {
          if (  (j<=INT_MAX2(e12,i,n2))
              &&(k<=INT_MAX3(e13e24,i,n3))
              &&(l<=INT_MAX4(e13e24,e34,i,j,k,n4))) {
            buffer[nonredundant_index] = buffer[redundant_index];
            nonredundant_index++;
            }
          redundant_index++;
          }
        }
      }
    }
  *red_off = redundant_index;
  *nonred_off = nonredundant_index;
  }

/* Ints is an integral buffer to be normalized.
 * shell1-4 specifies the shell quartet.
 * This function requires the redundant integrals. */
LOCAL_FUNCTION VOID
normalize_erep(ints,shell1,shell2,shell3,shell4)
double *ints;
shell_t *shell1;
shell_t *shell2;
shell_t *shell3;
shell_t *shell4;
{
  int redundant_index;
  int gc1,gc2,gc3,gc4;
  int index1,index2,index3,index4;
  double norm1,norm2,norm3,norm4;
  int i1,j1,k1;
  int i2,j2,k2;
  int i3,j3,k3;
  int i4,j4,k4;

  redundant_index = 0;
  FOR_GCCART(gc1,index1,i1,j1,k1,shell1)
    norm1 = shell1->norm[gc1][index1];
    FOR_GCCART(gc2,index2,i2,j2,k2,shell2)
      norm2 = norm1 * shell2->norm[gc2][index2];
      FOR_GCCART(gc3,index3,i3,j3,k3,shell3)
        norm3 = norm2 * shell3->norm[gc3][index3];
        FOR_GCCART(gc4,index4,i4,j4,k4,shell4)
          norm4 = norm3 * shell4->norm[gc4][index4];
          ints[redundant_index] *= norm4;
          redundant_index++;
          END_FOR_GCCART(index4)
        END_FOR_GCCART(index3)
      END_FOR_GCCART(index2)
    END_FOR_GCCART(index1)
  }

/* Ints is an integral buffer to be normalized.
 * shell1-4 specifies the shell quartet.
 * This function requires the redundant integrals. */
LOCAL_FUNCTION VOID
normalize_erep_given_gc(ints,gc1,gc2,gc3,gc4,shell1,shell2,shell3,shell4)
double *ints;
int gc1;
int gc2;
int gc3;
int gc4;
shell_t *shell1;
shell_t *shell2;
shell_t *shell3;
shell_t *shell4;
{
  int redundant_index;
  int am1,am2,am3,am4;
  int index1,index2,index3,index4;
  double norm1,norm2,norm3,norm4;
  int i1,j1,k1;
  int i2,j2,k2;
  int i3,j3,k3;
  int i4,j4,k4;

  am1 = shell1->type[gc1].am;
  am2 = shell2->type[gc2].am;
  am3 = shell3->type[gc3].am;
  am4 = shell4->type[gc4].am;

  redundant_index = 0;
  index1 = 0;
  FOR_CART(i1,j1,k1,am1)
    norm1 = shell1->norm[gc1][index1];
    index2 = 0;
    FOR_CART(i2,j2,k2,am2)
      norm2 = norm1 * shell2->norm[gc2][index2];
      index3 = 0;
      FOR_CART(i3,j3,k3,am3)
        norm3 = norm2 * shell3->norm[gc3][index3];
        index4 = 0;
        FOR_CART(i4,j4,k4,am4)
          norm4 = norm3 * shell4->norm[gc4][index4];
          ints[redundant_index] *= norm4;
          redundant_index++;
          index4++;
          END_FOR_CART
        index3++;
        END_FOR_CART
      index2++;
      END_FOR_CART
    index1++;
    END_FOR_CART
  }

/* This is an alternate interface to int_erep_all1der.  It takes
 * as arguments the flags, an integer vector of shell numbers
 * and an integer vector which will be filled in with size
 * information, if it is non-NULL, and the dercenters pointer. */
GLOBAL_FUNCTION VOID
int_erep_all1der_v(flags,shells,sizes,dercenters)
int flags;
int *shells;
int  *sizes;
der_centers_t *dercenters;
{
  int_erep_all1der(flags,&(shells[0]),&(shells[1]),&(shells[2]),&(shells[3]),
                   dercenters);
  if (sizes) {
    sizes[0] = INT_SH_NFUNC(int_cs1,shells[0]);
    sizes[1] = INT_SH_NFUNC(int_cs2,shells[1]);
    sizes[2] = INT_SH_NFUNC(int_cs3,shells[2]);
    sizes[3] = INT_SH_NFUNC(int_cs4,shells[3]);
    }
  }


GLOBAL_FUNCTION VOID
int_erep_bound1der(flags,bsh1,bsh2,size)
int flags;
int bsh1;
int bsh2;
int *size;
{
  double *current_buffer;
  int nints;
  double *user_int_buffer;
  centers_t *cs[4];
  int sh[4];
  int i;
  shell_t *shell1,*shell2,*shell3,*shell4;
  int am[4];
  int osh[4];
  int sh1 = bsh1;
  int sh2 = bsh2;
  int sh3 = bsh1;
  int sh4 = bsh2;
  int *psh1 = &sh1;
  int *psh2 = &sh2;
  int *psh3 = &sh3;
  int *psh4 = &sh4;

  cs[0] = int_cs1;
  cs[1] = int_cs2;
  cs[2] = int_cs3;
  cs[3] = int_cs4;

  sh[0] = *psh1;
  sh[1] = *psh2;
  sh[2] = *psh3;
  sh[3] = *psh4;

  /* Set up pointers to the current shells. */
  shell1 = &int_cs1->center[int_cs1->center_num[*psh1]]
                  .basis.shell[int_cs1->shell_num[*psh1]];
  shell2 = &int_cs2->center[int_cs2->center_num[*psh2]]
                  .basis.shell[int_cs2->shell_num[*psh2]];
  shell3 = &int_cs3->center[int_cs3->center_num[*psh3]]
                  .basis.shell[int_cs3->shell_num[*psh3]];
  shell4 = &int_cs4->center[int_cs4->center_num[*psh4]]
                  .basis.shell[int_cs4->shell_num[*psh4]];

  am[0] = int_find_jmax_shell(shell1);
  am[1] = int_find_jmax_shell(shell2);
  am[2] = int_find_jmax_shell(shell3);
  am[3] = int_find_jmax_shell(shell4);

  /* Compute the offset shell numbers. */
  osh[0] = *psh1 + int_cs1->shell_offset;
  osh[1] = *psh2 + int_cs2->shell_offset;
  osh[2] = *psh3 + int_cs3->shell_offset;
  osh[3] = *psh4 + int_cs4->shell_offset;

  /* Save the location of the int_buffer. */
  user_int_buffer = int_buffer;
  int_buffer = int_derint_buffer;

  /* Zero out the result integrals. */
  nints = shell1->nfunc * shell2->nfunc * shell3->nfunc * shell4->nfunc;
  for (i=0; i<3*nints; i++) user_int_buffer[i] = 0.0;

  /* Set the size so it is available to the caller. */
  *size = nints * 3;

  current_buffer = user_int_buffer;
  compute_erep_bound1der(flags|INT_NOPERM,current_buffer,
                          psh1,psh2,psh3,psh4);

  /* Normalize the integrals. */
  current_buffer = user_int_buffer;
  for (i=0; i<3; i++) {
    normalize_erep(current_buffer,shell1,shell2,shell3,shell4);
    current_buffer = &current_buffer[nints];
    }

  /* Eliminate redundant integrals, unless flags specifies otherwise. */
  current_buffer = user_int_buffer;
  if (!(flags&INT_REDUND)) {
    int redundant_offset = 0;
    int nonredundant_offset = 0;
    int e12,e13e24,e34;
    int i;

    if ((osh[0] == osh[3])&&(osh[1] == osh[2])&&(osh[0] != osh[1])) {
      fprintf(stderr,"nonredundant integrals cannot be generated (1der)\n");
      fail();
      }

    /* Shell equivalence information. */
    e12 = (osh[0] == osh[1]);
    e13e24 = ((osh[0] == osh[2]) && (osh[1] == osh[3]));
    e34 = (osh[2] == osh[3]);
    /* Repack x, y, and z integrals. */
    for (i=0; i<3; i++) {
      nonredundant_erep(current_buffer,e12,e34,e13e24,
                             shell1->nfunc,
                             shell2->nfunc,
                             shell3->nfunc,
                             shell4->nfunc,
                             &redundant_offset,
                             &nonredundant_offset);
      }
    }

  /* Return the integral buffers to their normal state. */
  int_derint_buffer = int_buffer;
  int_buffer = user_int_buffer;
  }

/* This routine computes a quantity needed to compute the
 * derivative integral bounds.
 * It fills int_buffer with (sh1+i sh2|sh1+i sh2).
 */
LOCAL_FUNCTION VOID
compute_erep_bound1der(flags,buffer,psh1,psh2,psh3,psh4)
int flags;
double *buffer;
int *psh1;
int *psh2;
int *psh3;
int *psh4;
{
  int oc1,oc2,oc3,oc4;
  int ii;
  int c1,c2,c3,c4;
  int osh[4];
  int i[4],j[4],k[4],am[4];
  int index;
  int sizem234,sizem34,sizem2,sizem3,sizem4;
  int sizep234,sizep34,sizep2,sizep3,sizep4;
  int sizepm234,sizepm34,sizepm2,sizepm3,sizepm4;
  shell_t *shell1,*shell2,*shell3,*shell4;

  /* Compute the offset shell numbers. */
  osh[0] = *psh1 + int_cs1->shell_offset;
  osh[1] = *psh2 + int_cs2->shell_offset;
  osh[2] = *psh3 + int_cs3->shell_offset;
  osh[3] = *psh4 + int_cs4->shell_offset;

  /* Set up pointers to the current shells. */
  shell1 = &int_cs1->center[int_cs1->center_num[*psh1]]
                  .basis.shell[int_cs1->shell_num[*psh1]];
  shell2 = &int_cs2->center[int_cs2->center_num[*psh2]]
                  .basis.shell[int_cs2->shell_num[*psh2]];
  shell3 = &int_cs3->center[int_cs3->center_num[*psh3]]
                  .basis.shell[int_cs3->shell_num[*psh3]];
  shell4 = &int_cs4->center[int_cs4->center_num[*psh4]]
                  .basis.shell[int_cs4->shell_num[*psh4]];

#define DCT1(n) ((((n)==0)||((n)==2))?1:0)
#define DCT2(n) ((((n)==0)||((n)==2))?((((n)==0))?1:-1):0)
  /* Offsets for the intermediates with angular momentum decremented. */
  for (ii=sizem2=0; ii<shell2->ncon; ii++) 
    sizem2 += INT_NCART(shell2->type[ii].am-DCT1(1));
  for (ii=sizem3=0; ii<shell3->ncon; ii++) 
    sizem3 += INT_NCART(shell3->type[ii].am-DCT1(2));
  for (ii=sizem4=0; ii<shell4->ncon; ii++) 
    sizem4 += INT_NCART(shell4->type[ii].am-DCT1(3));
  sizem34 = sizem4 * sizem3;
  sizem234 = sizem34 * sizem2;

  /* Offsets for the intermediates with angular momentum incremented. */
  for (ii=sizep2=0; ii<shell2->ncon; ii++) 
    sizep2 += INT_NCART(shell2->type[ii].am+DCT1(1));
  for (ii=sizep3=0; ii<shell3->ncon; ii++) 
    sizep3 += INT_NCART(shell3->type[ii].am+DCT1(2));
  for (ii=sizep4=0; ii<shell4->ncon; ii++) 
    sizep4 += INT_NCART(shell4->type[ii].am+DCT1(3));
  sizep34 = sizep4 * sizep3;
  sizep234 = sizep34 * sizep2;

  /* Offsets for the intermediates with angular momentum incremented and
   * decremented. */
  for (ii=sizepm2=0; ii<shell2->ncon; ii++) 
    sizepm2 += INT_NCART(shell2->type[ii].am+DCT2(1));
  for (ii=sizepm3=0; ii<shell3->ncon; ii++) 
    sizepm3 += INT_NCART(shell3->type[ii].am+DCT2(2));
  for (ii=sizepm4=0; ii<shell4->ncon; ii++) 
    sizepm4 += INT_NCART(shell4->type[ii].am+DCT2(3));
  sizepm34 = sizepm4 * sizepm3;
  sizepm234 = sizepm34 * sizepm2;

  /* END_DERLOOP must be redefined here because it previously depended
   * on the DCTEST macro */
#undef END_DERLOOP
#define END_DERLOOP(index,indexp1,sign) \
       END_FOR_CART\
     oc##indexp1 += INT_NCART(am[index] sign DCT1(index));\
     }

#undef END_ALLDERLOOPS
#define END_ALLDERLOOPS(sign)\
            END_DERLOOP(3,4,sign)\
          END_DERLOOP(2,3,sign)\
        END_DERLOOP(1,2,sign)\
      END_DERLOOP(0,1,sign)

  compute_erep(flags|INT_NOPERM|INT_REDUND|INT_NOBCHK,0,psh1,psh2,psh3,psh4,
                   -DCT1(0),
                   -DCT1(1),
                   -DCT1(2),
                   -DCT1(3));

   /* Place the contributions into the user integral buffer. */
   index = 0;
   /* The d/dx integrals */
  ALLDERLOOPS
    if (i[0]>0 && i[2]>0) {
      buffer[index] += i[0] * i[2] * int_buffer[
        (oc1 + INT_CARTINDEX(am[0]-DCT1(0),i[0]-DCT1(0),j[0])) * sizem234
       +(oc2 + INT_CARTINDEX(am[1]-DCT1(1),i[1]-DCT1(1),j[1])) * sizem34
       +(oc3 + INT_CARTINDEX(am[2]-DCT1(2),i[2]-DCT1(2),j[2])) * sizem4
       +(oc4 + INT_CARTINDEX(am[3]-DCT1(3),i[3]-DCT1(3),j[3]))
       ];
      }
    index++;
    END_ALLDERLOOPS(+)

   /* The d/dy integrals */
  ALLDERLOOPS
    if (j[0]>0 && j[2]>0) {
    buffer[index] += j[0] * j[2] * int_buffer[
         (oc1 + INT_CARTINDEX(am[0]-DCT1(0),i[0],j[0]-DCT1(0))) * sizem234
        +(oc2 + INT_CARTINDEX(am[1]-DCT1(1),i[1],j[1]-DCT1(1))) * sizem34
        +(oc3 + INT_CARTINDEX(am[2]-DCT1(2),i[2],j[2]-DCT1(2))) * sizem4
        +(oc4 + INT_CARTINDEX(am[3]-DCT1(3),i[3],j[3]-DCT1(3)))
        ];
      }
    index++;
  END_ALLDERLOOPS(+)

   /* The d/dz integrals */
  ALLDERLOOPS
    if (k[0]>0 && k[2]>0) {
    buffer[index] += k[0] * k[2] * int_buffer[
         (oc1 + INT_CARTINDEX(am[0]-DCT1(0),i[0],j[0])) * sizem234
        +(oc2 + INT_CARTINDEX(am[1]-DCT1(1),i[1],j[1])) * sizem34
        +(oc3 + INT_CARTINDEX(am[2]-DCT1(2),i[2],j[2])) * sizem4
        +(oc4 + INT_CARTINDEX(am[3]-DCT1(3),i[3],j[3]))
        ];
      }
    index++;
  END_ALLDERLOOPS(+)

#if 0
  printf("after -DCT1 buffer[5] is %12.8f\n",buffer[5]);
#endif

  /* Compute the next contribution to the integrals. */
  /* Tell the build routine that we need an exponent weighted contraction
   * with the exponents taken from the dercenter and adjust the
   * angular momentum of dercenter to the needed value. */
  int_expweight1 = 1;
  int_expweight3 = 1;
  compute_erep(flags|INT_NOPERM|INT_REDUND|INT_NOBCHK,0,psh1,psh2,psh3,psh4,
                     DCT1(0),
                     DCT1(1),
                     DCT1(2),
                     DCT1(3));
  int_expweight1 = 0;
  int_expweight3 = 0;

  /* Place the contributions into the user integral buffer. */
  index = 0;
  /* The d/dx integrals */
  ALLDERLOOPS
          buffer[index] += int_buffer[
             (oc1+INT_CARTINDEX(am[0]+DCT1(0),i[0]+DCT1(0),j[0]))*sizep234
            +(oc2+INT_CARTINDEX(am[1]+DCT1(1),i[1]+DCT1(1),j[1]))*sizep34
            +(oc3+INT_CARTINDEX(am[2]+DCT1(2),i[2]+DCT1(2),j[2]))*sizep4
            +(oc4+INT_CARTINDEX(am[3]+DCT1(3),i[3]+DCT1(3),j[3]))
            ];
    index++;
    END_ALLDERLOOPS(+)

  /* The d/dy integrals */
  ALLDERLOOPS
          buffer[index] += int_buffer[
             (oc1+INT_CARTINDEX(am[0]+DCT1(0),i[0],j[0]+DCT1(0)))*sizep234
            +(oc2+INT_CARTINDEX(am[1]+DCT1(1),i[1],j[1]+DCT1(1)))*sizep34
            +(oc3+INT_CARTINDEX(am[2]+DCT1(2),i[2],j[2]+DCT1(2)))*sizep4
            +(oc4+INT_CARTINDEX(am[3]+DCT1(3),i[3],j[3]+DCT1(3)))
            ];
          index++;
    END_ALLDERLOOPS(+)

  /* The d/dz integrals */
  ALLDERLOOPS
          buffer[index] += int_buffer[
               (oc1 + INT_CARTINDEX(am[0]+DCT1(0),i[0],j[0])) * sizep234
              +(oc2 + INT_CARTINDEX(am[1]+DCT1(1),i[1],j[1])) * sizep34
              +(oc3 + INT_CARTINDEX(am[2]+DCT1(2),i[2],j[2])) * sizep4
              +(oc4 + INT_CARTINDEX(am[3]+DCT1(3),i[3],j[3]))
              ];
          index++;
    END_ALLDERLOOPS(+)

#if 0
  printf("after +DCT1 buffer[5] is %12.8f\n",buffer[5]);
#endif

  /* END_DERLOOP must be redefined here because it previously depended
   * on the DCT1 macro */
#undef END_DERLOOP
#define END_DERLOOP(index,indexp1,sign) \
       END_FOR_CART\
     oc##indexp1 += INT_NCART(am[index] sign DCT2(index));\
     }

#undef END_ALLDERLOOPS
#define END_ALLDERLOOPS(sign)\
            END_DERLOOP(3,4,sign)\
          END_DERLOOP(2,3,sign)\
        END_DERLOOP(1,2,sign)\
      END_DERLOOP(0,1,sign)

  /* Compute the next contribution to the integrals. */
  /* Tell the build routine that we need an exponent weighted contraction
   * with the exponents taken from the dercenter and adjust the
   * angular momentum of dercenter to the needed value. */
  int_expweight1 = 1;
  compute_erep(flags|INT_NOPERM|INT_REDUND|INT_NOBCHK,0,psh1,psh2,psh3,psh4,
                     DCT2(0),
                     DCT2(1),
                     DCT2(2),
                     DCT2(3));
  int_expweight1 = 0;

  /* Place the contributions into the user integral buffer. */
  index = 0;
  /* The d/dx integrals */
  ALLDERLOOPS
    if (i[2] > 0)  {
          buffer[index] -= 2.0 * int_buffer[
             (oc1+INT_CARTINDEX(am[0]+DCT2(0),i[0]+DCT2(0),j[0]))*sizepm234
            +(oc2+INT_CARTINDEX(am[1]+DCT2(1),i[1]+DCT2(1),j[1]))*sizepm34
            +(oc3+INT_CARTINDEX(am[2]+DCT2(2),i[2]+DCT2(2),j[2]))*sizepm4
            +(oc4+INT_CARTINDEX(am[3]+DCT2(3),i[3]+DCT2(3),j[3]))
            ];
       }
    index++;
    END_ALLDERLOOPS(+)

  /* The d/dy integrals */
  ALLDERLOOPS
    if (j[2] > 0)  {
          buffer[index] -= 2.0 * int_buffer[
             (oc1+INT_CARTINDEX(am[0]+DCT2(0),i[0],j[0]+DCT2(0)))*sizepm234
            +(oc2+INT_CARTINDEX(am[1]+DCT2(1),i[1],j[1]+DCT2(1)))*sizepm34
            +(oc3+INT_CARTINDEX(am[2]+DCT2(2),i[2],j[2]+DCT2(2)))*sizepm4
            +(oc4+INT_CARTINDEX(am[3]+DCT2(3),i[3],j[3]+DCT2(3)))
            ];
       }
    index++;
    END_ALLDERLOOPS(+)

  /* The d/dz integrals */
  ALLDERLOOPS
    if (k[2] > 0)  {
          buffer[index] -= 2.0 * int_buffer[
               (oc1 + INT_CARTINDEX(am[0]+DCT2(0),i[0],j[0])) * sizepm234
              +(oc2 + INT_CARTINDEX(am[1]+DCT2(1),i[1],j[1])) * sizepm34
              +(oc3 + INT_CARTINDEX(am[2]+DCT2(2),i[2],j[2])) * sizepm4
              +(oc4 + INT_CARTINDEX(am[3]+DCT2(3),i[3],j[3]))
              ];
       }
    index++;
    END_ALLDERLOOPS(+)

#if 0
  printf("after +DCT2 buffer[5] is %12.8f\n",buffer[5]);
#endif
  }
