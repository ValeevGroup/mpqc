
/* $Log$
 * Revision 1.8  1995/10/25 21:15:46  cljanss
 * Cleaned up and added int_bound_to_double.  Fixed a bug that made bounds
 * underestimated by at most a factor of 4.
 *
 * Revision 1.7  1995/09/21 18:19:20  ibniels
 * Added function int_erep_2bound
 *
 * Revision 1.6  1995/08/21 19:36:19  cljanss
 * 1) New integral storage scheme using AVL trees.
 * 2) Updated bounds routines so the SCF program could use them.
 * 3) Got inttest working again.
 *
 * Revision 1.5  1995/03/17  01:49:21  cljanss
 * Removed -I. and -I$(SRCDIR) from the default include path in
 * GlobalMakefile to avoid name conflicts with system include files.
 * Modified files under src.lib to include all files relative to src.lib.
 * Makefiles under src.bin need to add the -I. and -I$(SRCDIR) back onto
 * INCLUDE and CXXINCLUDE or make other arrangements.
 *
 * Revision 1.4  1994/08/26  22:45:11  etseidl
 * fix a bunch of warnings, get rid of rcs id's, get rid of bread/bwrite and
 * fread/fwrite modules
 *
 * Revision 1.3  1994/05/18  22:05:22  etseidl
 * fix bug zeroing out int_Qvec and int_Rvec
 *
 * Revision 1.2  1993/12/30  13:32:44  etseidl
 * mostly rcs id stuff
 *
 * Revision 1.3  1992/06/17  22:04:25  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.2  1992/05/26  20:25:11  jannsen
 * make derivative bounds checking optional
 * add code to allow bound intermediates computable in shell blocks
 *
 * Revision 1.1  1992/05/13  18:29:24  jannsen
 * added bounds checking for derivative integrals
 *
 * */

/* This is the log for the original libdmtscf/scf_bnd.c:
 *
 * Revision 1.5  1992/05/04  11:01:05  seidl
 * Qvec no longer static
 *
 * Revision 1.4  1992/04/06  12:34:36  seidl
 * sandia changes
 *
 * Revision 1.3  1992/04/01  18:58:59  seidl
 * don't take sqrt twice
 *
 * Revision 1.2  1992/03/31  22:24:02  seidl
 * store bounds information as char
 *
 * Revision 1.1.1.1  1992/03/17  16:25:51  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:25:50  seidl
 * Initial revision
 *
 * Revision 1.1  1992/03/04  15:56:56  seidl
 * Initial revision
 * */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>

#include <chemistry/qc/intv2/atoms.h>
#include <chemistry/qc/intv2/int_types.h>
#include <chemistry/qc/intv2/int_flags.h>
#include <chemistry/qc/intv2/inter.h>

#define ALLOC_BOUND_GLOBALS
#include <chemistry/qc/intv2/bounds.h>

#include <chemistry/qc/intv2/bounds.gbl>
#include <chemistry/qc/intv2/bounds.lcl>

#include <chemistry/qc/intv2/comp_erep.gbl>

#define COMPUTE_Q  1
#define COMPUTE_R 2

GLOBAL_FUNCTION VOID
int_init_bounds_nocomp()
{
  int i;
  int nshell=int_cs1->nshell;
  int nsht=nshell*(nshell+1)/2;

  if (int_Qvec) free(int_Qvec);
  
  int_Qvec = (int_bound_t *) malloc(sizeof(int_bound_t)*nsht);
  if(int_Qvec==NULL) {
    fprintf(stderr,"int_init_bounds_nocomp: cannot malloc int_Qvec: %d\n",
            nsht);
    exit(1);
    }

  int_Rvec = NULL;

  int_Q = -126;
  for (i=0; i<nsht; i++) int_Qvec[i] = 0;
}

GLOBAL_FUNCTION VOID
int_init_bounds()
{
  int_init_bounds_nocomp();
  compute_bounds(&int_Q,int_Qvec,COMPUTE_Q);
  }

GLOBAL_FUNCTION VOID
int_init_bounds_1der_nocomp()
{
  int i;
  int nshell=int_cs1->nshell;
  int nsht=nshell*(nshell+1)/2;

  if (!int_derivative_bounds) {
    printf("requested der bounds but space not allocated\n");
    exit(1);
    }

  if (int_Qvec) free(int_Qvec);
  if (int_Rvec) free(int_Rvec);
  
  int_Qvec = (int_bound_t *) malloc(sizeof(int_bound_t)*nsht);
  int_Rvec = (int_bound_t *) malloc(sizeof(int_bound_t)*nsht);
  if((int_Qvec==NULL) || (int_Rvec==NULL)) {
    fprintf(stderr,"int_init_bounds_1der_nocomp: cannot malloc int_{R,Q}vec: %d\n",nsht);
    exit(1);
    }

  int_Q = -126;
  int_R = -126;
  for (i=0; i<nsht; i++) int_Qvec[i] = int_Rvec[i] = 0;
  }

GLOBAL_FUNCTION VOID
int_bounds_comp(s1,s2)
int s1;
int s2;
{
  compute_bounds_shell(&int_Q,int_Qvec,COMPUTE_Q,s1,s2);
  }

GLOBAL_FUNCTION VOID
int_bounds_1der_comp(s1,s2)
int s1;
int s2;
{
  compute_bounds_shell(&int_Q,int_Qvec,COMPUTE_Q,s1,s2);
  compute_bounds_shell(&int_R,int_Rvec,COMPUTE_R,s1,s2);
  }

GLOBAL_FUNCTION VOID
int_init_bounds_1der()
{
  int_init_bounds_1der_nocomp();
  compute_bounds(&int_Q,int_Qvec,COMPUTE_Q);
  compute_bounds(&int_R,int_Rvec,COMPUTE_R);
  }

GLOBAL_FUNCTION VOID
int_done_bounds()
{
  if (int_Qvec) free(int_Qvec);
  int_Qvec = 0;
  }

GLOBAL_FUNCTION VOID
int_done_bounds_1der()
{
  if (int_Qvec) free(int_Qvec);
  if (int_Rvec) free(int_Rvec);
  int_Qvec = 0;
  int_Rvec = 0;
  }

GLOBAL_FUNCTION int
int_erep_4bound(s1,s2,s3,s4)
int s1;
int s2;
int s3;
int s4;
{
  int ij=(s1>s2) ? ((s1*(s1+1))>>1)+s2 : ((s2*(s2+1))>>1)+s1;
  int kl=(s3>s4) ? ((s3*(s3+1))>>1)+s4 : ((s4*(s4+1))>>1)+s3;

  return((int) int_Qvec[ij]+int_Qvec[kl]);
  }

GLOBAL_FUNCTION int
int_erep_2bound(s1,s2)
int s1;
int s2;
{
  int ij=(s1>s2) ? ((s1*(s1+1))>>1)+s2 : ((s2*(s2+1))>>1)+s1;

  return((int) int_Qvec[ij]);
  }

GLOBAL_FUNCTION int
int_erep_0bound_1der()
{
#if 0
  printf("int_erep_0bound_1der(): Q: %4d R: %4d\n", int_Q,int_R);
#endif
  return 1 + int_Q + int_R;
  }

GLOBAL_FUNCTION int
int_erep_2bound_1der(s1,s2)
int s1;
int s2;
{
  int ij=(s1>s2) ? ((s1*(s1+1))>>1)+s2 : ((s2*(s2+1))>>1)+s1;
  int b1 = int_Qvec[ij] + int_R;
  int b2 = int_Q + int_Rvec[ij];

#if 0
  printf("int_erep_2bound_1der(%d,%d): Q: %4d R: %4d\n",s1,s2,
        int_Qvec[ij],int_Rvec[ij]);
#endif

  /* The actual bound is Qij R + Q Rij
   * but since I'm using log base 2 I'll use
   * 2 * max (Qij R, Q Rij) -> 1 + max (Qij + R, Q + Rij)
   */

  return 1 + ((b1>b2)? b1 : b2);
  }

GLOBAL_FUNCTION int
int_erep_4bound_1der(s1,s2,s3,s4)
int s1;
int s2;
int s3;
int s4;
{
  int ij=(s1>s2) ? ((s1*(s1+1))>>1)+s2 : ((s2*(s2+1))>>1)+s1;
  int kl=(s3>s4) ? ((s3*(s3+1))>>1)+s4 : ((s4*(s4+1))>>1)+s3;
  int b1 = int_Qvec[ij] + int_Rvec[kl];
  int b2 = int_Qvec[kl] + int_Rvec[ij];

#if 0
  printf("int_erep_4bound_1der(%d,%d,%d,%d): Q: %4d %4d R: %4d %4d\n",
         s1,s2,s3,s4,
         int_Qvec[ij],int_Qvec[kl],int_Rvec[ij],int_Rvec[kl]);
#endif

  /* The actual bound is Qij Rkl + Qkl Rij
   * but since I'm using log base 2 I'll use
   * 2 * max (Qij Rkl, Qkl Rij) -> 1 + max (Qij + Rkl, Qkl + Rij)
   */

  return 1 + ((b1>b2)? b1 : b2 );
  }

/* ripped off from clj's libintv2 */
/* (add subsequently ripped back on from ets's libdmtscf) */

/* Compute the partial bound arrays, either Q or R can be computed
 * with appropiate choice of flag. */
LOCAL_FUNCTION VOID
compute_bounds(overall,vec,flag)
int_bound_t *overall;
int_bound_t *vec;
int flag;
{
  int sh1,sh2;
  centers_t *cs1 = int_cs1;

  if ((cs1 != int_cs2)&&(cs1 != int_cs3)&&(cs1 != int_cs4)) {
    fprintf(stderr,"bounds.compute_bounds: all centers must be the same\n");
    exit(1);
    }

  *overall = -126;
  for(sh1=0; sh1 < cs1->nshell ; sh1++) {
    for(sh2=0; sh2 <= sh1 ; sh2++) {
      compute_bounds_shell(overall,vec,flag,sh1,sh2);
      }
    }
  }

/* Compute the partial bound arrays, either Q or R can be computed
 * with appropiate choice of flag. */
LOCAL_FUNCTION VOID
compute_bounds_shell(overall,vec,flag,sh1,sh2)
int_bound_t *overall;
int_bound_t *vec;
int flag;
int sh1;
int sh2;
{
  int nint;
  int shellij;
  int shells[4],size[4];
  double max;
  double tol = pow(2.0,-126.0);
  double loginv = 1.0/log(2.0);
  centers_t *cs1;
  int erep_flags = INT_EREP|INT_REDUND|INT_NOPERM|INT_NOBCHK;
  int old_int_integral_storage = int_integral_storage;
  int_integral_storage = 0;

  cs1 = int_cs1;
  if ((cs1 != int_cs2)&&(cs1 != int_cs3)&&(cs1 != int_cs4)) {
    fprintf(stderr,"bounds.compute_bounds: all centers must be the same\n");
    exit(1);
    }

  if (sh1<sh2) {
    int tmp = sh1;
    sh1 = sh2;
    sh2 = tmp;
    }

  shellij= ((sh1*(sh1+1))>>1) + sh2;
    shells[0]=shells[2]=sh1;
      shells[1]=shells[3]=sh2;

      if (flag == COMPUTE_Q) {
        int_erep_v(erep_flags,shells,size);
        nint = size[0]*size[1]*size[0]*size[1];
        max = find_max(int_buffer,nint);
#if 0
        printf("max for %d %d (size %d) is %15.11f\n", sh1, sh2, nint, max);
#endif
        }
      else if (flag == COMPUTE_R) {
        double max1,max2;
        int_erep_bound1der(erep_flags,sh1,sh2,&nint);
        max1 = find_max(int_buffer,nint);
#if 0
        printf("bound(%d) for (%d,%d) is %12.8f int_buffer =",
               flag,sh1,sh2,max1);
        for (i=0; (i<nint)&&(i<27); i++) printf(" %12.8f",int_buffer[i]);
        if (nint > 27) printf(" ...");
        printf("\n");
#endif
        int_erep_bound1der(erep_flags,sh2,sh1,&nint);
        max2 = find_max(int_buffer,nint);
        max = (max1>max2)?max1:max2;
        }
      else {
        printf("bad bound flag\n"); exit(1);
        }

    /* Compute the partial bound value. */
      max = sqrt(max);
      if (max>tol) {
        vec[shellij] = (int_bound_t) (log(max)*loginv + 0.999999999);
        }
      else {
        vec[shellij] = (int_bound_t) -126;
        }

    /* Multiply R contributions by a factor of two to account for
     * fact that contributions from both centers must be accounted
     * for. */
      if (flag == COMPUTE_R) vec[shellij]++;
      if (vec[shellij]>*overall) *overall = vec[shellij];
#if 0
      printf("bound(%d) for (%d,%d) is %4d int_buffer =",
             flag,sh1,sh2,vec[shellij]);
      for (i=0; (i<nint)&&(i<27); i++) printf(" %12.8f",int_buffer[i]);
      if (nint > 27) printf(" ...");
      printf("\n");
#endif
  int_integral_storage = old_int_integral_storage;
  }

/* This function is used to convert a double to its log base 2 rep
 * for use in bound computations. */
GLOBAL_FUNCTION int
int_bound_log(value)
double value;
{
  double tol = pow(2.0,-126.0);
  double loginv = 1.0/log(2.0);
  int_bound_t res;

  if (value > tol) res = log(value)*loginv;
  else res = -126;
  return res;
  }

/* This function is used to convert a bound to a double. */
GLOBAL_FUNCTION double
int_bound_to_double(bound)
int bound;
{
  double res = pow(2.0,bound);
  return res;
  }

/* This function is used to convert a double to its log base 2 rep
 * for use in bound computations. */
GLOBAL_FUNCTION double
int_bound_double(value)
int value;
{
  return pow(2.0,value);
  }

/* find the biggest number in the buffer */
LOCAL_FUNCTION double
find_max(int_buffer,nint)
double *int_buffer;
int nint;
{
  int i;
  double max = 0.0;
  for (i=0; i<nint; i++) {
    if (fabs(int_buffer[i]) > max) max = fabs(int_buffer[i]);
    }
  return max;
  }
