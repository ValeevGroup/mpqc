
/* $Log$
 * Revision 1.9  1995/10/25 21:19:46  cljanss
 * Adding support for pure am.  Gradients don't yet work.
 *
 * Revision 1.8  1995/09/13 22:53:41  etseidl
 * add int_accum_shell_kinetic
 *
 * Revision 1.7  1995/03/17  01:49:26  cljanss
 * Removed -I. and -I$(SRCDIR) from the default include path in
 * GlobalMakefile to avoid name conflicts with system include files.
 * Modified files under src.lib to include all files relative to src.lib.
 * Makefiles under src.bin need to add the -I. and -I$(SRCDIR) back onto
 * INCLUDE and CXXINCLUDE or make other arrangements.
 *
 * Revision 1.6  1995/03/05  06:05:27  cljanss
 * Added efield integrals.  Changed the dipole moment integral interface.
 *
 * Revision 1.5  1995/02/15  20:33:33  cljanss
 * Modified the dipole integral routines to be a bit more efficient
 * and consistent with buffer usage for the efield routines.
 *
 * Revision 1.4  1995/01/17  18:11:53  cljanss
 * Added routine int_accum_shell_point_charge.
 *
 * Revision 1.3  1994/08/26  22:45:21  etseidl
 * fix a bunch of warnings, get rid of rcs id's, get rid of bread/bwrite and
 * fread/fwrite modules
 *
 * Revision 1.2  1993/12/30  13:32:47  etseidl
 * mostly rcs id stuff
 *
 * Revision 1.5  1993/04/28  00:31:15  jannsen
 * added point charge integrals
 *
 * Revision 1.4  1992/06/17  22:04:33  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.3  1992/05/12  10:41:31  seidl
 * add routines for dipole moment integrals
 *
 * Revision 1.2  1992/03/31  01:21:43  jannsen
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.2  1992/01/08  17:14:45  cljanss
 * fixed a bug in derivative 1e integrals
 *
 * Revision 1.1  1991/12/14  00:19:36  cljanss
 * Initial revision
 * */
/*
 * Revision 1.7  91/09/28  19:26:46  cljanss
 * Switch to new file naming scheme
 * 
 * Revision 1.6  91/09/26  15:53:05  cljanss
 * fixed bug in normalization of overlap derivatives
 * 
 * Revision 1.5  91/09/18  01:08:45  cljanss
 * corrected int_accum_shell_nuclear_1der and made the local accum routines
 * into global int_accum routines.
 * 
 * Revision 1.4  1991/09/16  22:56:56  cljanss
 * added int_shell_overlap_1der and modified int_overlap_1der to use it.
 *
 * Revision 1.3  1991/09/10  19:34:45  cljanss
 * finished putting in first derivatives
 *
 * Revision 1.2  1991/08/09  16:49:10  cljanss
 * added functions to computed 1e ints by shell blocks
 *
 * Revision 1.1  1991/06/16  16:40:07  janssen
 * Initial revision
 * */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/atoms.h>
#include <chemistry/qc/intv2/int_macros.h>
#include <chemistry/qc/intv2/fjttable.h>

#include <chemistry/qc/intv2/comp_1e.gbl>
#include <chemistry/qc/intv2/comp_1e.lcl>

/*#include "initialize.gbl"*/
#include <chemistry/qc/intv2/int_fjt.gbl>
#include <chemistry/qc/intv2/utils.gbl>
#include <chemistry/qc/intv2/transform.h>

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

static double oo2zeta_a;
static double oo2zeta_b;
static double sMus;
static double sTs;
static double xi;
static double *A;
static double *B;
static double *C;
static double ss;
static double PmA[3];
static double PmB[3];
static double PmC[3];
static double zeta;
static double oo2zeta;
static shell_t *shell1, *shell2;

static int exponent_weighted = -1;
static int scale_shell_result = 0;
static double result_scale_factor = 1.0;
static int three_center = 0;
static centers_t *third_centers;
static int third_centernum;

static init_order = -1;
static double *scratch = NULL;
static double *cartesianbuffer = NULL;

static int mu;

#define IN(i,j) ((i)==(j)?1:0)
#define SELECT(x1,x2,x3,s) (((s)==0)?x1:(((s)==1)?(x2):(x3)))

/* ------------ Initialization of 1e routines. ------------------- */
/* This routine returns a buffer large enough to hold a shell doublet
 * of integrals (if order == 0) or derivative integrals (if order == 1).
 */
GLOBAL_FUNCTION double *
int_initialize_1e(flags,order,cs1,cs2)
int flags;
int order;
centers_t *cs1;
centers_t *cs2;
{
  int jmax1,jmax2,jmax;
  int scratchsize,nshell2;

  /* The efield routines look like derivatives so bump up order if
   * it is zero to allow efield integrals to be computed.
   */
  if (order == 0) order = 1;

  jmax1 = int_find_jmax(cs1);
  jmax2 = int_find_jmax(cs2);
  jmax = jmax1 + jmax2;

  int_initialize_fjt(jmax + 2*order);

  nshell2 = int_find_ncartmax(cs1)*int_find_ncartmax(cs2);

  if (order == 0) {
    init_order = 0;
    scratchsize = nshell2;
    }
  else if (order == 1) {
    init_order = 1;
    scratchsize = nshell2*3;
    }
  else {
    fprintf(stderr,"int_initialize_1e: invalid order: %d\n",order);
    exit(1);
    }

#if 0
  printf("allocating %d doubles in init_1e\n",scratchsize);
#endif
  scratch = (double *) malloc(scratchsize*sizeof(double));
  cartesianbuffer = (double *) malloc(scratchsize*sizeof(double));

  /* Return the buffer in case the user has a use for a buffer
   * of such a size. */
  return scratch;
  }

GLOBAL_FUNCTION VOID
int_done_1e()
{
  int_done_fjt();
  init_order = -1;
  free(scratch);
  free(cartesianbuffer);
  scratch = NULL;
  cartesianbuffer = NULL;
  }


/* --------------------------------------------------------------- */
/* ------------- Routines for the overlap matrix ----------------- */
/* --------------------------------------------------------------- */

/* This computes the overlap integrals between functions in two shells.
 * The result is placed in the buffer.
 */
GLOBAL_FUNCTION VOID
int_shell_overlap(cs1,cs2,buff,ish,jsh)
centers_t *cs1;
centers_t *cs2;
double *buff;
int ish;
int jsh;
{
  int c1,s1,i1,j1,k1,c2,s2,i2,j2,k2;
  int gc1,gc2;
  int index,index1,index2;
  double shell1norm,shell2norm;

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);
  index = 0;
  FOR_GCCART(gc1,index1,i1,j1,k1,shell1)
    shell1norm = shell1->norm[gc1][index1];
    FOR_GCCART(gc2,index2,i2,j2,k2,shell2)
      shell2norm = shell2->norm[gc2][index2];
      cartesianbuffer[index] = shell1norm * shell2norm * comp_shell_overlap(gc1,i1,j1,k1,gc2,i2,j2,k2);
      index++;
      END_FOR_GCCART(index2)
    END_FOR_GCCART(index1)

  int_transform_1e(cartesianbuffer, buff, shell1, shell2);
  }

/* This computes the overlap ints between functions in two shells.
 * The result is placed in the buffer.
 */
GLOBAL_FUNCTION VOID
int_shell_overlap_1der(cs1,cs2,buff,ish,jsh,dercs,centernum)
centers_t *cs1;
centers_t *cs2;
double *buff;
int ish;
int jsh;
centers_t *dercs;
int centernum;
{
  int i;
  int c1,s1,c2,s2;
  int ni,nj;

  if (!(init_order >= 0)) {
    fprintf(stderr,"int_shell_overlap: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);

  ni = shell1->nfunc;
  nj = shell2->nfunc;

#if 0
  printf("zeroing %d*%d*3 elements of buff\n",ni,nj);
#endif
  for (i=0; i<ni*nj*3; i++) {
    buff[i] = 0.0;
    }

  int_accum_shell_overlap_1der(cs1,cs2,buff,ish,jsh,dercs,centernum);
  }

/* This computes the overlap derivative integrals between functions
 * in two shells.  The result is accumulated in the buffer which is ordered
 * atom 0 x, y, z, atom 1, ... .
 */
GLOBAL_FUNCTION VOID
int_accum_shell_overlap_1der(cs1,cs2,buff,ish,jsh,dercs,centernum)
centers_t *cs1;
centers_t *cs2;
double *buff;
int ish;
int jsh;
centers_t *dercs;
int centernum;
{
  accum_shell_1der(cs1,cs2,buff,ish,jsh,dercs,centernum,comp_shell_overlap);
  }

/* Compute the overlap for the shell.  The arguments are the
 * cartesian exponents for centers 1 and 2.  The shell1 and shell2
 * globals are used. */
LOCAL_FUNCTION double
comp_shell_overlap(gc1,i1,j1,k1,gc2,i2,j2,k2)
int gc1;
int i1;
int j1;
int k1;
int gc2;
int i2;
int j2;
int k2;
{
  double exp1,exp2;
  int i,j,xyz;
  double result;
  double Pi;
  double oozeta;
  double AmB,AmB2;
  double tmp;

  if ((i1<0)||(j1<0)||(k1<0)||(i2<0)||(j2<0)||(k2<0)) return 0.0;

  /* Loop over the primitives in the shells. */
  result = 0.0;
  for (i=0; i<shell1->nprim; i++) {
    for (j=0; j<shell2->nprim; j++) {

      /* Compute the intermediates. */
      exp1 = shell1->exp[i];
      exp2 = shell2->exp[j];
      oozeta = 1.0/(exp1 + exp2);
      oo2zeta = 0.5*oozeta;
      AmB2 = 0.0;
      for (xyz=0; xyz<3; xyz++) {
        Pi = oozeta*(exp1 * A[xyz] + exp2 * B[xyz]);
        PmA[xyz] = Pi - A[xyz];
        PmB[xyz] = Pi - B[xyz];
        AmB = A[xyz] - B[xyz];
        AmB2 += AmB*AmB;
        }
      ss =   pow(3.141592653589793/(exp1+exp2),1.5)
           * exp_cutoff(- oozeta * exp1 * exp2 * AmB2);
      tmp     =  shell1->coef[gc1][i] * shell2->coef[gc2][j]
               * comp_prim_overlap(i1,j1,k1,i2,j2,k2);
      if (exponent_weighted == 0) tmp *= exp1;
      else if (exponent_weighted == 1) tmp *= exp2;
      result += tmp;
      }
    }

  return result;
  }

/* Compute the overlap between two primitive functions. */
GLOBAL_FUNCTION double
int_prim_overlap(pshell1,pshell2,pA,pB,prim1,prim2,i1,j1,k1,i2,j2,k2)
shell_t *pshell1;
shell_t *pshell2;
double *pA;
double *pB;
int prim1;
int prim2;
int i1;
int j1;
int k1;
int i2;
int j2;
int k2;
{
  int xyz;
  double Pi;
  double oozeta;
  double AmB,AmB2;

  /* Compute the intermediates. */
  oozeta = 1.0/(shell1->exp[prim1] + shell2->exp[prim2]);
  oo2zeta = 0.5*oozeta;
  AmB2 = 0.0;
  for (xyz=0; xyz<3; xyz++) {
    Pi = oozeta*(shell1->exp[prim1] * A[xyz] + shell2->exp[prim2] * B[xyz]);
    PmA[xyz] = Pi - A[xyz];
    PmB[xyz] = Pi - B[xyz];
    AmB = A[xyz] - B[xyz];
    AmB2 += AmB*AmB;
    }
  ss =   pow(3.141592653589793/(shell1->exp[prim1]+shell2->exp[prim2]),1.5)
       * exp_cutoff(- oozeta * shell1->exp[prim1] * shell2->exp[prim2] * AmB2);
  return comp_prim_overlap(i1,j1,k1,i2,j2,k2);
  }

LOCAL_FUNCTION double
comp_prim_overlap(i1,j1,k1,i2,j2,k2)
int i1;
int j1;
int k1;
int i2;
int j2;
int k2;
{
  double result;

  if (i1) {
    result = PmA[0] * comp_prim_overlap(i1-1,j1,k1,i2,j2,k2);
    if (i1>1) result += oo2zeta*(i1-1) * comp_prim_overlap(i1-2,j1,k1,i2,j2,k2);
    if (i2>0) result += oo2zeta*i2 * comp_prim_overlap(i1-1,j1,k1,i2-1,j2,k2);
    return result;
    }
  if (j1) {
    result = PmA[1] * comp_prim_overlap(i1,j1-1,k1,i2,j2,k2);
    if (j1>1) result += oo2zeta*(j1-1) * comp_prim_overlap(i1,j1-2,k1,i2,j2,k2);
    if (j2>0) result += oo2zeta*j2 * comp_prim_overlap(i1,j1-1,k1,i2,j2-1,k2);
    return result;
    }
  if (k1) {
    result = PmA[2] * comp_prim_overlap(i1,j1,k1-1,i2,j2,k2);
    if (k1>1) result += oo2zeta*(k1-1) * comp_prim_overlap(i1,j1,k1-2,i2,j2,k2);
    if (k2>0) result += oo2zeta*k2 * comp_prim_overlap(i1,j1,k1-1,i2,j2,k2-1);
    return result;
    }

  if (i2) {
    result = PmB[0] * comp_prim_overlap(i1,j1,k1,i2-1,j2,k2);
    if (i2>1) result += oo2zeta*(i2-1) * comp_prim_overlap(i1,j1,k1,i2-2,j2,k2);
    if (i1>0) result += oo2zeta*i1 * comp_prim_overlap(i1-1,j1,k1,i2-1,j2,k2);
    return result;
    }
  if (j2) {
    result = PmB[1] * comp_prim_overlap(i1,j1,k1,i2,j2-1,k2);
    if (j2>1) result += oo2zeta*(j2-1) * comp_prim_overlap(i1,j1,k1,i2,j2-2,k2);
    if (j1>0) result += oo2zeta*j1 * comp_prim_overlap(i1,j1-1,k1,i2,j2-1,k2);
    return result;
    }
  if (k2) {
    result = PmB[2] * comp_prim_overlap(i1,j1,k1,i2,j2,k2-1);
    if (k2>1) result += oo2zeta*(k2-1) * comp_prim_overlap(i1,j1,k1,i2,j2,k2-2);
    if (k1>0) result += oo2zeta*k1 * comp_prim_overlap(i1,j1,k1-1,i2,j2,k2-1);
    return result;
    }

  return ss;
  }

/* --------------------------------------------------------------- */
/* ------------- Routines for the kinetic energy ----------------- */
/* --------------------------------------------------------------- */

/* This computes the kinetic energy integrals between functions in two shells.
 * The result is placed in the buffer.
 */
GLOBAL_FUNCTION VOID
int_shell_kinetic(cs1,cs2,buff,ish,jsh)
centers_t *cs1;
centers_t *cs2;
double *buff;
int ish;
int jsh;
{
  int c1,s1,i1,j1,k1,c2,s2,i2,j2,k2;
  int cart1,cart2;
  int index;
  double norm1,norm2;
  int gc1,gc2;

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);
  index = 0;
  FOR_GCCART(gc1,cart1,i1,j1,k1,shell1)
    norm1 = shell1->norm[gc1][cart1];
    FOR_GCCART(gc2,cart2,i2,j2,k2,shell2)
      norm2 = shell2->norm[gc2][cart2];
      cartesianbuffer[index] = norm1 * norm2 * comp_shell_kinetic(gc1,i1,j1,k1,gc2,i2,j2,k2);
      index++;
      END_FOR_GCCART(cart2)
    END_FOR_GCCART(cart1)

  int_transform_1e(cartesianbuffer, buff, shell1, shell2);
  }

GLOBAL_FUNCTION VOID
int_accum_shell_kinetic(cs1,cs2,buff,ish,jsh)
centers_t *cs1;
centers_t *cs2;
double *buff;
int ish;
int jsh;
{
  int c1,s1,i1,j1,k1,c2,s2,i2,j2,k2;
  int cart1,cart2;
  int index;
  double norm1,norm2;
  int gc1,gc2;

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);
  index = 0;

  FOR_GCCART(gc1,cart1,i1,j1,k1,shell1)
    norm1 = shell1->norm[gc1][cart1];
    FOR_GCCART(gc2,cart2,i2,j2,k2,shell2)
      norm2 = shell2->norm[gc2][cart2];
      cartesianbuffer[index] = norm1 * norm2 * comp_shell_kinetic(gc1,i1,j1,k1,gc2,i2,j2,k2);
      index++;
      END_FOR_GCCART(cart2)
    END_FOR_GCCART(cart1)
  int_accum_transform_1e(cartesianbuffer, buff, shell1, shell2);
  }

/* This computes the kinetic energy derivative integrals between functions
 * in two shells.  The result is accumulated in the buffer which is ordered
 * atom 0 x, y, z, atom 1, ... .
 */
GLOBAL_FUNCTION VOID
int_accum_shell_kinetic_1der(cs1,cs2,buff,ish,jsh,dercs,centernum)
centers_t *cs1;
centers_t *cs2;
double *buff;
int ish;
int jsh;
centers_t *dercs;
int centernum;
{
  accum_shell_1der(cs1,cs2,buff,ish,jsh,dercs,centernum,comp_shell_kinetic);
  }

/* This computes the basis function part of 
 * generic derivative integrals between functions
 * in two shells.  The result is accumulated in the buffer which is ordered
 * atom 0 x, y, z, atom 1, ... .
 * The function used to compute the nonderivative integrals is shell_function.
 */
LOCAL_FUNCTION VOID
accum_shell_1der(cs1,cs2,buff,ish,jsh,dercs,centernum,shell_function)
centers_t *cs1;
centers_t *cs2;
double *buff;
int ish;
int jsh;
centers_t *dercs;
int centernum;
double (*shell_function)();
{
  int i;
  int gc1,gc2;
  int c1,s1,i1,j1,k1,c2,s2,i2,j2,k2;
  int index1,index2;
  double tmp[3];
  double shell1norm,shell2norm;
  double *ctmp = cartesianbuffer;

  /* fprintf(stdout,"accum_shell_1der: working on ( %2d | %2d )",ish,jsh);
   * if (three_center) fprintf(stdout," third_centernum = %d\n",third_centernum);
   * else fprintf(stdout,"\n");
   */

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);
  FOR_GCCART(gc1,index1,i1,j1,k1,shell1)
    shell1norm = shell1->norm[gc1][index1];
#if 0
    printf("shell1norm = shell1->norm[%d][%d] = % f\n",gc1,index1,shell1norm);
#endif
    FOR_GCCART(gc2,index2,i2,j2,k2,shell2)
      shell2norm = shell2->norm[gc2][index2];
      if ((cs1==cs2)&&(c1==c2)) {
        if (    three_center
             && !((cs1==third_centers)&&(c1==third_centernum))
             && ((cs1==dercs)&&(c1==centernum))) {
          for (i=0; i<3; i++) {
            /* Derivative wrt first shell. */
            exponent_weighted = 0;
            tmp[i] = 2.0 *
               (*shell_function)(gc1,i1+IN(i,0),j1+IN(i,1),k1+IN(i,2),gc2,i2,j2,k2);
            exponent_weighted = -1;
            if (SELECT(i1,j1,k1,i)) {
              tmp[i] -= SELECT(i1,j1,k1,i) *
                (*shell_function)(gc1,i1-IN(i,0),j1-IN(i,1),k1-IN(i,2),gc2,i2,j2,k2);
              }
            /* Derviative wrt second shell. */
            exponent_weighted = 1;
            tmp[i] += 2.0 *
               (*shell_function)(gc1,i1,j1,k1,gc2,i2+IN(i,0),j2+IN(i,1),k2+IN(i,2));
            exponent_weighted = -1;
            if (SELECT(i2,j2,k2,i)) {
              tmp[i] -= SELECT(i2,j2,k2,i) *
                (*shell_function)(gc1,i1,j1,k1,gc2,i2-IN(i,0),j2-IN(i,1),k2-IN(i,2));
              }
            }
	  }
        else {
          /* If there are two centers and they are the same, then we
           * use translational invariance to get a net contrib of 0.0 */
          for (i=0; i<3; i++) tmp[i] = 0.0;
          }
        }
      else if ((cs1==dercs)&&(c1==centernum)) {
        for (i=0; i<3; i++) {
          exponent_weighted = 0;
          tmp[i] = 2.0 *
             (*shell_function)(gc1,i1+IN(i,0),j1+IN(i,1),k1+IN(i,2),gc2,i2,j2,k2);
          exponent_weighted = -1;
          if (SELECT(i1,j1,k1,i)) {
            tmp[i] -= SELECT(i1,j1,k1,i) *
              (*shell_function)(gc1,i1-IN(i,0),j1-IN(i,1),k1-IN(i,2),gc2,i2,j2,k2);
            }
          }
        }
      else if ((cs2==dercs)&&(c2==centernum)) {
        for (i=0; i<3; i++) {
          exponent_weighted = 1;
          tmp[i] = 2.0 *
             (*shell_function)(gc1,i1,j1,k1,gc2,i2+IN(i,0),j2+IN(i,1),k2+IN(i,2));
          exponent_weighted = -1;
          if (SELECT(i2,j2,k2,i)) {
            tmp[i] -= SELECT(i2,j2,k2,i) *
              (*shell_function)(gc1,i1,j1,k1,gc2,i2-IN(i,0),j2-IN(i,1),k2-IN(i,2));
            }
          }
        }
      else {
        for (i=0; i<3; i++) tmp[i] = 0.0;
        }
#if 0
      printf("tmp1 = % f % f % f\n",tmp[0],tmp[1],tmp[2]);
#endif

#if 0
printf("exp_w = %d, scale_sh = %d, result_scale_f = % f\n",
       exponent_weighted,scale_shell_result,result_scale_factor);
#endif
      if (scale_shell_result) {
        for (i=0; i<3; i++) tmp[i] *= result_scale_factor;
        }
#if 0
      printf("tmp1 = % f % f % f\n",tmp[0],tmp[1],tmp[2]);
#endif
#if 0
      printf("shell1norm = % f, shell2norm = % f\n",shell1norm,shell2norm);
#endif
      for (i=0; i<3; i++) ctmp[i] = shell1norm * shell2norm * tmp[i];

#if 0
      printf("ctmp = % f % f % f\n",ctmp[0],ctmp[1],ctmp[2]);
#endif

      /* Increment the pointer to the xyz for the next atom. */
      ctmp += 3;
      END_FOR_GCCART(index2)
    END_FOR_GCCART(index1)

  int_accum_transform_1e_xyz(cartesianbuffer, buff, shell1, shell2);
  }

/* Compute the kinetic energy for the shell.  The arguments are the
 * cartesian exponents for centers 1 and 2.  The shell1 and shell2
 * globals are used. */
LOCAL_FUNCTION double
comp_shell_kinetic(gc1,i1,j1,k1,gc2,i2,j2,k2)
int gc1;
int i1;
int j1;
int k1;
int gc2;
int i2;
int j2;
int k2;
{
  int i,j,xyz;
  double result;
  double Pi;
  double oozeta;
  double AmB,AmB2;
  double tmp;

  /* Loop over the primitives in the shells. */
  result = 0.0;
  for (i=0; i<shell1->nprim; i++) {
    for (j=0; j<shell2->nprim; j++) {

      /* Compute the intermediates. */
      oo2zeta_a = 0.5/shell1->exp[i];
      oo2zeta_b = 0.5/shell2->exp[j];
      oozeta = 1.0/(shell1->exp[i] + shell2->exp[j]);
      oo2zeta = 0.5*oozeta;
      xi = oozeta * shell1->exp[i] * shell2->exp[j];
      AmB2 = 0.0;
      for (xyz=0; xyz<3; xyz++) {
        Pi = oozeta*(shell1->exp[i] * A[xyz] + shell2->exp[j] * B[xyz]);
        PmA[xyz] = Pi - A[xyz];
        PmB[xyz] = Pi - B[xyz];
        AmB = A[xyz] - B[xyz];
        AmB2 += AmB*AmB;
        }
      /* The s integral kinetic energy. */
      ss =   pow(3.141592653589793/(shell1->exp[i]+shell2->exp[j]),1.5)
           * exp_cutoff(- xi * AmB2);
      sTs =  ss
           * xi
           * (3.0 - 2.0 * xi * AmB2);
      tmp     =  shell1->coef[gc1][i] * shell2->coef[gc2][j]
               * comp_prim_kinetic(i1,j1,k1,i2,j2,k2);
      if (exponent_weighted == 0) tmp *= shell1->exp[i];
      else if (exponent_weighted == 1) tmp *= shell2->exp[j];
      result += tmp;
      }
    }

  /* fprintf(stdout,"comp_shell_kinetic(%d,%d,%d,%d,%d,%d): result = % 12.8lf\n",
   *         i1,j1,k1,i2,j2,k2,result);
   */
  return result;
  }

LOCAL_FUNCTION double
comp_prim_kinetic(i1,j1,k1,i2,j2,k2)
int i1;
int j1;
int k1;
int i2;
int j2;
int k2;
{
  double tmp;
  double result;

  if (i1) {
    result = PmA[0] * comp_prim_kinetic(i1-1,j1,k1,i2,j2,k2);
    if (i1>1) result += oo2zeta*(i1-1)*comp_prim_kinetic(i1-2,j1,k1,i2,j2,k2);
    if (i2) result += oo2zeta*i2*comp_prim_kinetic(i1-1,j1,k1,i2-1,j2,k2);
    tmp = comp_prim_overlap(i1,j1,k1,i2,j2,k2);
    if (i1>1) tmp -= oo2zeta_a*(i1-1)*comp_prim_overlap(i1-2,j1,k1,i2,j2,k2);
    result += 2.0 * xi * tmp;
    return result;
    }
  if (j1) {
    result = PmA[1] * comp_prim_kinetic(i1,j1-1,k1,i2,j2,k2);
    if (j1>1) result += oo2zeta*(j1-1)*comp_prim_kinetic(i1,j1-2,k1,i2,j2,k2);
    if (j2) result += oo2zeta*j2*comp_prim_kinetic(i1,j1-1,k1,i2,j2-1,k2);
    tmp = comp_prim_overlap(i1,j1,k1,i2,j2,k2);
    if (j1>1) tmp -= oo2zeta_a*(j1-1)*comp_prim_overlap(i1,j1-2,k1,i2,j2,k2);
    result += 2.0 * xi * tmp;
    return result;
    }
  if (k1) {
    result = PmA[2] * comp_prim_kinetic(i1,j1,k1-1,i2,j2,k2);
    if (k1>1) result += oo2zeta*(k1-1)*comp_prim_kinetic(i1,j1,k1-2,i2,j2,k2);
    if (k2) result += oo2zeta*k2*comp_prim_kinetic(i1,j1,k1-1,i2,j2,k2-1);
    tmp = comp_prim_overlap(i1,j1,k1,i2,j2,k2);
    if (k1>1) tmp -= oo2zeta_a*(k1-1)*comp_prim_overlap(i1,j1,k1-2,i2,j2,k2);
    result += 2.0 * xi * tmp;
    return result;
    }
  if (i2) {
    result = PmB[0] * comp_prim_kinetic(i1,j1,k1,i2-1,j2,k2);
    if (i2>1) result += oo2zeta*(i2-1)*comp_prim_kinetic(i1,j1,k1,i2-2,j2,k2);
    if (i1) result += oo2zeta*i1*comp_prim_kinetic(i1-1,j1,k1,i2-1,j2,k2);
    tmp = comp_prim_overlap(i1,j1,k1,i2,j2,k2);
    if (i2>1) tmp -= oo2zeta_b*(i2-1)*comp_prim_overlap(i1,j1,k1,i2-2,j2,k2);
    result += 2.0 * xi * tmp;
    return result;
    }
  if (j2) {
    result = PmB[1] * comp_prim_kinetic(i1,j1,k1,i2,j2-1,k2);
    if (j2>1) result += oo2zeta*(j2-1)*comp_prim_kinetic(i1,j1,k1,i2,j2-2,k2);
    if (j1) result += oo2zeta*i1*comp_prim_kinetic(i1,j1-1,k1,i2,j2-1,k2);
    tmp = comp_prim_overlap(i1,j1,k1,i2,j2,k2);
    if (j2>1) tmp -= oo2zeta_b*(j2-1)*comp_prim_overlap(i1,j1,k1,i2,j2-2,k2);
    result += 2.0 * xi * tmp;
    return result;
    }
  if (k2) {
    result = PmB[2] * comp_prim_kinetic(i1,j1,k1,i2,j2,k2-1);
    if (k2>1) result += oo2zeta*(k2-1)*comp_prim_kinetic(i1,j1,k1,i2,j2,k2-2);
    if (k1) result += oo2zeta*i1*comp_prim_kinetic(i1,j1,k1-1,i2,j2,k2-1);
    tmp = comp_prim_overlap(i1,j1,k1,i2,j2,k2);
    if (k2>1) tmp -= oo2zeta_b*(k2-1)*comp_prim_overlap(i1,j1,k1,i2,j2,k2-2);
    result += 2.0 * xi * tmp;
    return result;
    }

  return sTs;
  }

/* --------------------------------------------------------------- */
/* ------------- Routines for the nuclear attraction ------------- */
/* --------------------------------------------------------------- */

/* This computes the nuclear attraction derivative integrals between functions
 * in two shells.  The result is accumulated in the buffer which is ordered
 * atom 0 x, y, z, atom 1, ... .
 */
GLOBAL_FUNCTION VOID
int_accum_shell_nuclear_1der(cs1,cs2,buff,ish,jsh,dercs,centernum)
centers_t *cs1;
centers_t *cs2;
double *buff;
int ish;
int jsh;
centers_t *dercs;
int centernum;
{
  int_accum_shell_nuclear_hf_1der(cs1,cs2,buff,ish,jsh,dercs,centernum);
  int_accum_shell_nuclear_nonhf_1der(cs1,cs2,buff,ish,jsh,dercs,centernum);
  }

/* A correction to the Hellman-Feynman part is computed which
 * is not included in the original HF routine.  This is only needed
 * if the real Hellman-Feynman forces are desired, because the sum
 * of the hf_1der and nonhf_1der forces are still correct.
 */
GLOBAL_FUNCTION VOID
int_accum_shell_nuclear_hfc_1der(cs1,cs2,buff,ish,jsh,dercs,centernum)
centers_t *cs1;
centers_t *cs2;
double *buff;
int ish;
int jsh;
centers_t *dercs;
int centernum;
{

  /* If both ish and jsh are not on the der center,
   * then there's no correction. */
  if (!(  (cs1==dercs)
        &&(cs2==dercs)
        &&(cs1->center_num[ish]==centernum)
        &&(cs2->center_num[jsh]==centernum))) {
    return;
    }

  /* Compute the nuc attr part of the nuclear derivative for three equal
   * centers. */
  scale_shell_result = 1;
  result_scale_factor = -cs1->center[centernum].charge;
  C = cs1->center[centernum].r;
  accum_shell_efield(cs1,cs2,buff,ish,jsh);
  scale_shell_result = 0;

  }

/* This computes the nuclear attraction derivative integrals between functions
 * in two shells.  The result is accumulated in the buffer which is ordered
 * atom 0 x, y, z, atom 1, ... .  Only the Hellman-Feynman part is computed.
 */
GLOBAL_FUNCTION VOID
int_accum_shell_nuclear_hf_1der(cs1,cs2,buff,ish,jsh,dercs,centernum)
centers_t *cs1;
centers_t *cs2;
double *buff;
int ish;
int jsh;
centers_t *dercs;
int centernum;
{

  /* If both ish and jsh are on the der center, then the contrib is zero. */
  if (  (cs1==dercs)
      &&(cs2==dercs)
      &&(cs1->center_num[ish]==centernum)
      &&(cs2->center_num[jsh]==centernum)) {
    return;
    }

  /* Compute the nuc attr part of the nuclear derivative. */
  if (cs1 == dercs) {
    scale_shell_result = 1;
    result_scale_factor = -cs1->center[centernum].charge;
    C = cs1->center[centernum].r;
    accum_shell_efield(cs1,cs2,buff,ish,jsh);
    scale_shell_result = 0;
    }
  else if (cs2 == dercs) {
    scale_shell_result = 1;
    result_scale_factor = -cs2->center[centernum].charge;
    C = cs2->center[centernum].r;
    accum_shell_efield(cs1,cs2,buff,ish,jsh);
    scale_shell_result = 0;
    }

  }

/* This computes the nuclear attraction derivative integrals between functions
 * in two shells.  The result is accumulated in the buffer which is ordered
 * atom 0 x, y, z, atom 1, ... .  Only the non Hellman-Feynman part is computed.
 */
GLOBAL_FUNCTION VOID
int_accum_shell_nuclear_nonhf_1der(cs1,cs2,buff,ish,jsh,dercs,centernum)
centers_t *cs1;
centers_t *cs2;
double *buff;
int ish;
int jsh;
centers_t *dercs;
int centernum;
{
  int i;

  /* Get the basis function part of the nuclear derivative. */
  three_center = 1;
  third_centers = cs1;
  for (i=0; i<cs1->n; i++) {
    third_centernum = i;
    C = cs1->center[i].r;
    scale_shell_result = 1;
    result_scale_factor = -cs1->center[i].charge;
    accum_shell_1der(cs1,cs2,buff,ish,jsh,dercs,centernum,comp_shell_nuclear);
    scale_shell_result = 0;
#if 0
    printf("int_accum_shell_nuclear_nonhf_1der: i = %d, buff[0] = %f\n",
           i,buff[0]);
#endif
    }
  if (cs2!=cs1) {
    third_centers = cs2;
    for (i=0; i<cs2->n; i++) {
      third_centernum = i;
      C = cs2->center[i].r;
      scale_shell_result = 1;
      result_scale_factor = -cs2->center[i].charge;
      accum_shell_1der(cs1,cs2,buff,ish,jsh,dercs,centernum,comp_shell_nuclear);
      scale_shell_result = 0;
      }
    }
  three_center = 0;

  }

/* This computes the efield integrals between functions in two shells.
 * The result is accumulated in the buffer in the form bf1 x y z, bf2
 * x y z, etc.
 */
GLOBAL_FUNCTION VOID
int_accum_shell_efield(cs1,cs2,buff,ish,jsh,position)
centers_t *cs1;
centers_t *cs2;
double *buff;
int ish;
int jsh;
double *position;
{
  scale_shell_result = 0;
  C = position;
  accum_shell_efield(cs1,cs2,buff,ish,jsh);
}

/* This computes the efield integrals between functions in two shells.
 * The result is accumulated in the buffer in the form bf1 x y z, bf2
 * x y z, etc.  The globals scale_shell_result, result_scale_factor,
 * and C must be set before this is called.
 */
LOCAL_FUNCTION VOID
accum_shell_efield(cs1,cs2,buff,ish,jsh)
centers_t *cs1;
centers_t *cs2;
double *buff;
int ish;
int jsh;
{
  int i;
  int c1,s1,i1,j1,k1,c2,s2,i2,j2,k2;
  double efield[3];
  int gc1,gc2;
  int index1,index2;
  double *tmp = cartesianbuffer;

  if (!(init_order >= 1)) {
    fprintf(stderr,"accum_shell_efield: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);

  FOR_GCCART(gc1,index1,i1,j1,k1,shell1)
    FOR_GCCART(gc2,index2,i2,j2,k2,shell2)
      comp_shell_efield(efield,gc1,i1,j1,k1,gc2,i2,j2,k2);
      if (scale_shell_result) {
        for (i=0; i<3; i++) efield[i] *= result_scale_factor;
        }
      for (i=0; i<3; i++) tmp[i] = efield[i];
      tmp += 3;
      END_FOR_GCCART(index2)
    END_FOR_GCCART(index1)

  int_accum_transform_1e_xyz(cartesianbuffer, buff, shell1, shell2);
  }

/* This computes the nuc rep energy integrals between functions in two shells.
 * The result is placed in the buffer.
 */
GLOBAL_FUNCTION VOID
int_shell_nuclear(cs1,cs2,buff,ish,jsh)
centers_t *cs1;
centers_t *cs2;
double *buff;
int ish;
int jsh;
{
  int i;
  int c1,s1,i1,j1,k1,c2,s2,i2,j2,k2;
  int index;
  int gc1,gc2;
  int cart1,cart2;
  double norm1,norm2;

  if (!(init_order >= 0)) {
    fprintf(stderr,"int_shell_nuclear: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);
  index = 0;

  FOR_GCCART(gc1,cart1,i1,j1,k1,shell1)
    norm1 = shell1->norm[gc1][cart1];
    FOR_GCCART(gc2,cart2,i2,j2,k2,shell2)
      norm2 = shell2->norm[gc2][cart2];
      cartesianbuffer[index] = 0.0;
      /* Loop thru the centers on cs1. */
      for (i=0; i<cs1->n; i++) {
        C = cs1->center[i].r;
        cartesianbuffer[index] -= comp_shell_nuclear(gc1,i1,j1,k1,gc2,i2,j2,k2)
                       * cs1->center[i].charge;
        }
      /* Loop thru the centers on cs2 if necessary. */
      if (cs2 != cs1) {
        for (i=0; i<cs2->n; i++) {
          C = cs2->center[i].r;
          cartesianbuffer[index]-=comp_shell_nuclear(gc1,i1,j1,k1,gc2,i2,j2,k2)
                         * cs2->center[i].charge;
          }
        }
      cartesianbuffer[index] *= norm1 * norm2;
      index++;
      END_FOR_GCCART(cart2)
    END_FOR_GCCART(cart1)

  int_transform_1e(cartesianbuffer, buff, shell1, shell2);
  }

/* This computes the integrals between functions in two shells for
 * a point charge interaction operator.
 * The result is placed in the buffer.
 */
GLOBAL_FUNCTION VOID
int_accum_shell_point_charge(cs1,cs2,buff,ish,jsh,ncharge,charge,position)
centers_t *cs1;
centers_t *cs2;
double *buff;
int ish;
int jsh;
int ncharge;
double* charge;
double** position;
{
  int i;
  int c1,s1,i1,j1,k1,c2,s2,i2,j2,k2;
  int index;
  int gc1,gc2;
  int cart1,cart2;
  double norm1,norm2;
  double tmp;

  if (!(init_order >= 0)) {
    fprintf(stderr,"int_shell_pointcharge: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);
  index = 0;

  FOR_GCCART(gc1,cart1,i1,j1,k1,shell1)
    norm1 = shell1->norm[gc1][cart1];
    FOR_GCCART(gc2,cart2,i2,j2,k2,shell2)
      norm2 = shell2->norm[gc2][cart2];
      /* Loop thru the point charges. */
      tmp = 0.0;
      for (i=0; i<ncharge; i++) {
        C = position[i];
        tmp -=  comp_shell_nuclear(gc1,i1,j1,k1,gc2,i2,j2,k2)
                       * charge[i];
        }
      cartesianbuffer[index] = tmp * norm1 * norm2;
      index++;
      END_FOR_GCCART(cart2)
    END_FOR_GCCART(cart1)

  int_accum_transform_1e(cartesianbuffer, buff, shell1, shell2);
  }

/* This computes the integrals between functions in two shells for
 * a point charge interaction operator.
 * The result is placed in the buffer.
 */
GLOBAL_FUNCTION VOID
int_shell_point_charge(cs1,cs2,buff,ish,jsh,ncharge,charge,position)
centers_t *cs1;
centers_t *cs2;
double *buff;
int ish;
int jsh;
int ncharge;
double* charge;
double** position;
{
  int i;
  int c1,s1,i1,j1,k1,c2,s2,i2,j2,k2;
  int index;
  int gc1,gc2;
  int cart1,cart2;
  double norm1,norm2;

  if (!(init_order >= 0)) {
    fprintf(stderr,"int_shell_pointcharge: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);
  index = 0;

  FOR_GCCART(gc1,cart1,i1,j1,k1,shell1)
    norm1 = shell1->norm[gc1][cart1];
    FOR_GCCART(gc2,cart2,i2,j2,k2,shell2)
      norm2 = shell2->norm[gc2][cart2];
      cartesianbuffer[index] = 0.0;
      /* Loop thru the point charges. */
      for (i=0; i<ncharge; i++) {
        C = position[i];
        cartesianbuffer[index] -= comp_shell_nuclear(gc1,i1,j1,k1,gc2,i2,j2,k2)
                                * charge[i];
        }
      cartesianbuffer[index] *= norm1 * norm2;
      index++;
      END_FOR_GCCART(cart2)
    END_FOR_GCCART(cart1)

  int_transform_1e(cartesianbuffer, buff, shell1, shell2);
  }


/* This computes the 1e Hamiltonian integrals between functions in two shells.
 * The result is placed in the buffer.
 */
GLOBAL_FUNCTION VOID
int_shell_hcore(cs1,cs2,buff,ish,jsh)
centers_t *cs1;
centers_t *cs2;
double *buff;
int ish;
int jsh;
{
  int i;
  int c1,s1,i1,j1,k1,c2,s2,i2,j2,k2;
  int index;
  int cart1,cart2;
  double norm1,norm2;
  int gc1,gc2;

  if (!(init_order >= 0)) {
    fprintf(stderr,"int_shell_hcore: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);

  index = 0;
  FOR_GCCART(gc1,cart1,i1,j1,k1,shell1)
    norm1 = shell1->norm[gc1][cart1];
    FOR_GCCART(gc2,cart2,i2,j2,k2,shell2)
      norm2 = shell2->norm[gc2][cart2];
      cartesianbuffer[index] = comp_shell_kinetic(gc1,i1,j1,k1,gc2,i2,j2,k2);
      /* Loop thru the centers on cs1. */
      for (i=0; i<cs1->n; i++) {
        C = cs1->center[i].r;
        cartesianbuffer[index] -= comp_shell_nuclear(gc1,i1,j1,k1,gc2,i2,j2,k2)
                                * cs1->center[i].charge;
        }
      /* Loop thru the centers on cs2 if necessary. */
      if (cs2 != cs1) {
        for (i=0; i<cs2->n; i++) {
          C = cs2->center[i].r;
          cartesianbuffer[index]-=comp_shell_nuclear(gc1,i1,j1,k1,gc2,i2,j2,k2)
                                * cs2->center[i].charge;
          }
        }
      cartesianbuffer[index] *= norm1 * norm2;
      index++;
      END_FOR_GCCART(cart2)
    END_FOR_GCCART(cart1)

  int_transform_1e(cartesianbuffer, buff, shell1, shell2);
  }

/* This computes the 1e Hamiltonian deriv ints between functions in two shells.
 * The result is placed in the buffer.
 */
GLOBAL_FUNCTION VOID
int_shell_hcore_1der(cs1,cs2,buff,ish,jsh,dercs,centernum)
centers_t *cs1;
centers_t *cs2;
double *buff;
int ish;
int jsh;
centers_t *dercs;
int centernum;
{
  int i;
  int c1,s1,c2,s2;
  int ni,nj;

  if (!(init_order >= 0)) {
    fprintf(stderr,"int_shell_hcore: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);

  ni = shell1->nfunc;
  nj = shell2->nfunc;

  for (i=0; i<ni*nj*3; i++) {
    buff[i] = 0.0;
    }

  int_accum_shell_nuclear_1der(cs1,cs2,buff,ish,jsh,dercs,centernum);
  int_accum_shell_kinetic_1der(cs1,cs2,buff,ish,jsh,dercs,centernum);
  }

/* This computes the kinetic deriv ints between functions in two shells.
 * The result is placed in the buffer.
 */
GLOBAL_FUNCTION VOID
int_shell_kinetic_1der(cs1,cs2,buff,ish,jsh,dercs,centernum)
centers_t *cs1;
centers_t *cs2;
double *buff;
int ish;
int jsh;
centers_t *dercs;
int centernum;
{
  int i;
  int c1,s1,c2,s2;
  int ni,nj;

  if (!(init_order >= 0)) {
    fprintf(stderr,"int_shell_kinetic: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);

  ni = shell1->nfunc;
  nj = shell2->nfunc;

  for (i=0; i<ni*nj*3; i++) {
    buff[i] = 0.0;
    }

  int_accum_shell_kinetic_1der(cs1,cs2,buff,ish,jsh,dercs,centernum);
  }

/* This computes the nuclear deriv ints between functions in two shells.
 * The result is placed in the buffer.
 */
GLOBAL_FUNCTION VOID
int_shell_nuclear_1der(cs1,cs2,buff,ish,jsh,dercs,centernum)
centers_t *cs1;
centers_t *cs2;
double *buff;
int ish;
int jsh;
centers_t *dercs;
int centernum;
{
  int i;
  int c1,s1,c2,s2;
  int ni,nj;

  if (!(init_order >= 0)) {
    fprintf(stderr,"int_shell_nuclear: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);

  ni = shell1->nfunc;
  nj = shell2->nfunc;

  for (i=0; i<ni*nj*3; i++) {
    buff[i] = 0.0;
    }

  int_accum_shell_nuclear_1der(cs1,cs2,buff,ish,jsh,dercs,centernum);
  }

/* This computes the nuclear deriv ints between functions in two shells.
 * Only the Hellman-Feynman part is computed.
 * The result is placed in the buffer.
 */
GLOBAL_FUNCTION VOID
int_shell_nuclear_hf_1der(cs1,cs2,buff,ish,jsh,dercs,centernum)
centers_t *cs1;
centers_t *cs2;
double *buff;
int ish;
int jsh;
centers_t *dercs;
int centernum;
{
  int i;
  int c1,s1,c2,s2;
  int ni,nj;

  if (!(init_order >= 0)) {
    fprintf(stderr,"int_shell_nuclear_hf_1der: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);

  ni = shell1->nfunc;
  nj = shell2->nfunc;

  for (i=0; i<ni*nj*3; i++) {
    buff[i] = 0.0;
    }

  int_accum_shell_nuclear_hf_1der(cs1,cs2,buff,ish,jsh,dercs,centernum);
  }

/* This computes the nuclear deriv ints between functions in two shells.
 * Only the non Hellman-Feynman part is computed.
 * The result is placed in the buffer.
 */
GLOBAL_FUNCTION VOID
int_shell_nuclear_nonhf_1der(cs1,cs2,buff,ish,jsh,dercs,centernum)
centers_t *cs1;
centers_t *cs2;
double *buff;
int ish;
int jsh;
centers_t *dercs;
int centernum;
{
  int i;
  int c1,s1,c2,s2;
  int ni,nj;

  if (!(init_order >= 0)) {
    fprintf(stderr,"int_shell_nuclear_nonhf_1der: one electron routines are not init'ed\n");
    exit(1);
    }

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);

  ni = shell1->nfunc;
  nj = shell2->nfunc;

#if 0
  printf("int_shell_nuclear_nonhf_1der: zeroing %d doubles in buff\n",ni*nj*3);
#endif
  for (i=0; i<ni*nj*3; i++) {
    buff[i] = 0.0;
    }

  int_accum_shell_nuclear_nonhf_1der(cs1,cs2,buff,ish,jsh,dercs,centernum);
  }

/* Compute the nuclear attraction for the shell.  The arguments are the
 * cartesian exponents for centers 1 and 2.  The shell1 and shell2
 * globals are used. */
LOCAL_FUNCTION double
comp_shell_nuclear(gc1,i1,j1,k1,gc2,i2,j2,k2)
int gc1;
int i1;
int j1;
int k1;
int gc2;
int i2;
int j2;
int k2;
{
  int i,j,k,xyz;
  double result;
  double Pi;
  double oozeta;
  double AmB,AmB2;
  double PmC2;
  double auxcoef;
  int am;
  double tmp;

  am = i1+j1+k1+i2+j2+k2;

  /* Loop over the primitives in the shells. */
  result = 0.0;
  for (i=0; i<shell1->nprim; i++) {
    for (j=0; j<shell2->nprim; j++) {

      /* Compute the intermediates. */
      zeta = shell1->exp[i] + shell2->exp[j];
      oozeta = 1.0/zeta;
      oo2zeta = 0.5*oozeta;
      AmB2 = 0.0;
      PmC2 = 0.0;
      for (xyz=0; xyz<3; xyz++) {
        Pi = oozeta*(shell1->exp[i] * A[xyz] + shell2->exp[j] * B[xyz]);
        PmA[xyz] = Pi - A[xyz];
        PmB[xyz] = Pi - B[xyz];
        PmC[xyz] = Pi - C[xyz];
        AmB = A[xyz] - B[xyz];
        AmB2 += AmB*AmB;
        PmC2 += PmC[xyz]*PmC[xyz];
        }

      /* The auxillary integral coeficients. */
      auxcoef =   2.0 * 3.141592653589793/(shell1->exp[i]+shell2->exp[j])
           * exp_cutoff(- oozeta * shell1->exp[i] * shell2->exp[j] * AmB2);

      /* The Fm(U) intermediates. */
      int_fjt(am,zeta*PmC2);

      /* Convert the Fm(U) intermediates into the auxillary
       * nuclear attraction integrals. */
      for (k=0; k<=am; k++) {
        int_fjttable.d[k] *= auxcoef;
        }

      /* Compute the nuclear attraction integral. */
      tmp     =  shell1->coef[gc1][i] * shell2->coef[gc2][j]
               * comp_prim_nuclear(i1,j1,k1,i2,j2,k2,0);

      if (exponent_weighted == 0) tmp *= shell1->exp[i];
      else if (exponent_weighted == 1) tmp *= shell2->exp[j];

      result += tmp;
      }
    }

  /* printf("comp_shell_nuclear(%d,%d,%d,%d,%d,%d): result = % 12.8lf\n",
   *         i1,j1,k1,i2,j2,k2,result);
   */
  return result;
  }

LOCAL_FUNCTION double
comp_prim_nuclear(i1,j1,k1,i2,j2,k2,m)
int i1;
int j1;
int k1;
int i2;
int j2;
int k2;
int m;
{
  double result;

  if (i1) {
    result  = PmA[0] * comp_prim_nuclear(i1-1,j1,k1,i2,j2,k2,m);
    result -= PmC[0] * comp_prim_nuclear(i1-1,j1,k1,i2,j2,k2,m+1);
    if (i1>1) result += oo2zeta * (i1-1)
                       * (  comp_prim_nuclear(i1-2,j1,k1,i2,j2,k2,m)
                          - comp_prim_nuclear(i1-2,j1,k1,i2,j2,k2,m+1));
    if (i2) result += oo2zeta * i2
                     * (  comp_prim_nuclear(i1-1,j1,k1,i2-1,j2,k2,m)
                        - comp_prim_nuclear(i1-1,j1,k1,i2-1,j2,k2,m+1));
    }
  else if (j1) {
    result  = PmA[1] * comp_prim_nuclear(i1,j1-1,k1,i2,j2,k2,m);
    result -= PmC[1] * comp_prim_nuclear(i1,j1-1,k1,i2,j2,k2,m+1);
    if (j1>1) result += oo2zeta * (j1-1)
                       * (  comp_prim_nuclear(i1,j1-2,k1,i2,j2,k2,m)
                          - comp_prim_nuclear(i1,j1-2,k1,i2,j2,k2,m+1));
    if (j2) result += oo2zeta * j2
                     * (  comp_prim_nuclear(i1,j1-1,k1,i2,j2-1,k2,m)
                        - comp_prim_nuclear(i1,j1-1,k1,i2,j2-1,k2,m+1));
    }
  else if (k1) {
    result  = PmA[2] * comp_prim_nuclear(i1,j1,k1-1,i2,j2,k2,m);
    result -= PmC[2] * comp_prim_nuclear(i1,j1,k1-1,i2,j2,k2,m+1);
    if (k1>1) result += oo2zeta * (k1-1)
                       * (  comp_prim_nuclear(i1,j1,k1-2,i2,j2,k2,m)
                          - comp_prim_nuclear(i1,j1,k1-2,i2,j2,k2,m+1));
    if (k2) result += oo2zeta * k2
                     * (  comp_prim_nuclear(i1,j1,k1-1,i2,j2,k2-1,m)
                        - comp_prim_nuclear(i1,j1,k1-1,i2,j2,k2-1,m+1));
    }
  else if (i2) {
    result  = PmB[0] * comp_prim_nuclear(i1,j1,k1,i2-1,j2,k2,m);
    result -= PmC[0] * comp_prim_nuclear(i1,j1,k1,i2-1,j2,k2,m+1);
    if (i2>1) result += oo2zeta * (i2-1)
                       * (  comp_prim_nuclear(i1,j1,k1,i2-2,j2,k2,m)
                          - comp_prim_nuclear(i1,j1,k1,i2-2,j2,k2,m+1));
    if (i1) result += oo2zeta * i1
                     * (  comp_prim_nuclear(i1-1,j1,k1,i2-1,j2,k2,m)
                        - comp_prim_nuclear(i1-1,j1,k1,i2-1,j2,k2,m+1));
    }
  else if (j2) {
    result  = PmB[1] * comp_prim_nuclear(i1,j1,k1,i2,j2-1,k2,m);
    result -= PmC[1] * comp_prim_nuclear(i1,j1,k1,i2,j2-1,k2,m+1);
    if (j2>1) result += oo2zeta * (j2-1)
                       * (  comp_prim_nuclear(i1,j1,k1,i2,j2-2,k2,m)
                          - comp_prim_nuclear(i1,j1,k1,i2,j2-2,k2,m+1));
    if (j1) result += oo2zeta * j1
                     * (  comp_prim_nuclear(i1,j1-1,k1,i2,j2-1,k2,m)
                        - comp_prim_nuclear(i1,j1-1,k1,i2,j2-1,k2,m+1));
    }
  else if (k2) {
    result  = PmB[2] * comp_prim_nuclear(i1,j1,k1,i2,j2,k2-1,m);
    result -= PmC[2] * comp_prim_nuclear(i1,j1,k1,i2,j2,k2-1,m+1);
    if (k2>1) result += oo2zeta * (k2-1)
                       * (  comp_prim_nuclear(i1,j1,k1,i2,j2,k2-2,m)
                          - comp_prim_nuclear(i1,j1,k1,i2,j2,k2-2,m+1));
    if (k1) result += oo2zeta * k1
                     * (  comp_prim_nuclear(i1,j1,k1-1,i2,j2,k2-1,m)
                        - comp_prim_nuclear(i1,j1,k1-1,i2,j2,k2-1,m+1));
    }
  else result = int_fjttable.d[m];

  /* fprintf(stdout,"  comp_prim_nuclear(%d,%d,%d,%d,%d,%d,%d) = % 12.8lf\n",
   *         i1,j1,k1,i2,j2,k2,m,result);
   */

  return result;
  }

/* Compute the electric field integral for the shell.  The arguments are the
 * the electric field vector, the
 * cartesian exponents for centers 1 and 2.  The shell1 and shell2
 * globals are used. */
LOCAL_FUNCTION VOID
comp_shell_efield(efield,gc1,i1,j1,k1,gc2,i2,j2,k2)
double *efield;
int gc1;
int i1;
int j1;
int k1;
int gc2;
int i2;
int j2;
int k2;
{
  int i,j,k,xyz;
  double result[3];
  double norm1,norm2;
  double Pi;
  double oozeta;
  double AmB,AmB2;
  double PmC2;
  double auxcoef;
  int am;

  am = i1+j1+k1+i2+j2+k2;

  if (shell1->norm) norm1 = shell1->norm[gc1][INT_CARTINDEX(i1+j1+k1,i1,j1)];
  else norm1 = 1.0;
  if (shell2->norm) norm2 = shell2->norm[gc2][INT_CARTINDEX(i2+j2+k2,i2,j2)];
  else norm2 = 1.0;

  /* Loop over the primitives in the shells. */
  for (xyz=0; xyz<3; xyz++) result[xyz] = 0.0;
  for (i=0; i<shell1->nprim; i++) {
    for (j=0; j<shell2->nprim; j++) {

      /* Compute the intermediates. */
      zeta = shell1->exp[i] + shell2->exp[j];
      oozeta = 1.0/zeta;
      oo2zeta = 0.5*oozeta;
      AmB2 = 0.0;
      PmC2 = 0.0;
      for (xyz=0; xyz<3; xyz++) {
        Pi = oozeta*(shell1->exp[i] * A[xyz] + shell2->exp[j] * B[xyz]);
        PmA[xyz] = Pi - A[xyz];
        PmB[xyz] = Pi - B[xyz];
        PmC[xyz] = Pi - C[xyz];
        AmB = A[xyz] - B[xyz];
        AmB2 += AmB*AmB;
        PmC2 += PmC[xyz]*PmC[xyz];
        }

      /* The auxillary integral coeficients. */
      auxcoef =   2.0 * 3.141592653589793/(shell1->exp[i]+shell2->exp[j])
           * exp_cutoff(- oozeta * shell1->exp[i] * shell2->exp[j] * AmB2);

      /* The Fm(U) intermediates. */
      int_fjt(am+1,zeta*PmC2);

      /* Convert the Fm(U) intermediates into the auxillary
       * nuclear attraction integrals. */
      for (k=0; k<=am+1; k++) {
        int_fjttable.d[k] *= auxcoef;
        }

      /* Compute the nuclear attraction integral. */
      for (xyz=0; xyz<3; xyz++) {
        result[xyz] +=  shell1->coef[gc1][i] * shell2->coef[gc2][j]
                      * comp_prim_efield(xyz,i1,j1,k1,i2,j2,k2,0);
        }
      }
    }

  for (xyz=0; xyz<3; xyz++) efield[xyz] = norm1 * norm2 * result[xyz];

  /* fprintf(stdout,"comp_shell_efield(%d,%d,%d,%d,%d,%d): % 12.8lf % 12.8lf % 12.8lf\n",
   *         i1,j1,k1,i2,j2,k2,efield[0],efield[1],efield[2]);
   */
  }

LOCAL_FUNCTION double
comp_prim_efield(xyz,i1,j1,k1,i2,j2,k2,m)
int xyz;
int i1;
int j1;
int k1;
int i2;
int j2;
int k2;
int m;
{
  double result;

  /* if ((xyz != 0) || (i1 != 1)) return 0.0; */

  if (i1) {
    result  = PmA[0] * comp_prim_efield(xyz,i1-1,j1,k1,i2,j2,k2,m);
    /* fprintf(stdout," % 12.8lf",result); */
    result -= PmC[0] * comp_prim_efield(xyz,i1-1,j1,k1,i2,j2,k2,m+1);
    /* fprintf(stdout," % 12.8lf",result); */
    if (i1>1) result += oo2zeta * (i1-1)
                       * (  comp_prim_efield(xyz,i1-2,j1,k1,i2,j2,k2,m)
                          - comp_prim_efield(xyz,i1-2,j1,k1,i2,j2,k2,m+1));
    if (i2) result += oo2zeta * i2
                     * (  comp_prim_efield(xyz,i1-1,j1,k1,i2-1,j2,k2,m)
                        - comp_prim_efield(xyz,i1-1,j1,k1,i2-1,j2,k2,m+1));
    if (xyz==0) result += comp_prim_nuclear(i1-1,j1,k1,i2,j2,k2,m+1);
    /* fprintf(stdout," % 12.8lf",result); */
    }
  else if (j1) {
    result  = PmA[1] * comp_prim_efield(xyz,i1,j1-1,k1,i2,j2,k2,m);
    result -= PmC[1] * comp_prim_efield(xyz,i1,j1-1,k1,i2,j2,k2,m+1);
    if (j1>1) result += oo2zeta * (j1-1)
                       * (  comp_prim_efield(xyz,i1,j1-2,k1,i2,j2,k2,m)
                          - comp_prim_efield(xyz,i1,j1-2,k1,i2,j2,k2,m+1));
    if (j2) result += oo2zeta * j2
                     * (  comp_prim_efield(xyz,i1,j1-1,k1,i2,j2-1,k2,m)
                        - comp_prim_efield(xyz,i1,j1-1,k1,i2,j2-1,k2,m+1));
    if (xyz==1) result += comp_prim_nuclear(i1,j1-1,k1,i2,j2,k2,m+1);
    }
  else if (k1) {
    result  = PmA[2] * comp_prim_efield(xyz,i1,j1,k1-1,i2,j2,k2,m);
    result -= PmC[2] * comp_prim_efield(xyz,i1,j1,k1-1,i2,j2,k2,m+1);
    if (k1>1) result += oo2zeta * (k1-1)
                       * (  comp_prim_efield(xyz,i1,j1,k1-2,i2,j2,k2,m)
                          - comp_prim_efield(xyz,i1,j1,k1-2,i2,j2,k2,m+1));
    if (k2) result += oo2zeta * k2
                     * (  comp_prim_efield(xyz,i1,j1,k1-1,i2,j2,k2-1,m)
                        - comp_prim_efield(xyz,i1,j1,k1-1,i2,j2,k2-1,m+1));
    if (xyz==2) result += comp_prim_nuclear(i1,j1,k1-1,i2,j2,k2,m+1);
    }
  else if (i2) {
    result  = PmB[0] * comp_prim_efield(xyz,i1,j1,k1,i2-1,j2,k2,m);
    result -= PmC[0] * comp_prim_efield(xyz,i1,j1,k1,i2-1,j2,k2,m+1);
    if (i2>1) result += oo2zeta * (i2-1)
                       * (  comp_prim_efield(xyz,i1,j1,k1,i2-2,j2,k2,m)
                          - comp_prim_efield(xyz,i1,j1,k1,i2-2,j2,k2,m+1));
    if (i1) result += oo2zeta * i1
                     * (  comp_prim_efield(xyz,i1-1,j1,k1,i2-1,j2,k2,m)
                        - comp_prim_efield(xyz,i1-1,j1,k1,i2-1,j2,k2,m+1));
    if (xyz==0) result += comp_prim_nuclear(i1,j1,k1,i2-1,j2,k2,m+1);
    }
  else if (j2) {
    result  = PmB[1] * comp_prim_efield(xyz,i1,j1,k1,i2,j2-1,k2,m);
    result -= PmC[1] * comp_prim_efield(xyz,i1,j1,k1,i2,j2-1,k2,m+1);
    if (j2>1) result += oo2zeta * (j2-1)
                       * (  comp_prim_efield(xyz,i1,j1,k1,i2,j2-2,k2,m)
                          - comp_prim_efield(xyz,i1,j1,k1,i2,j2-2,k2,m+1));
    if (j1) result += oo2zeta * j1
                     * (  comp_prim_efield(xyz,i1,j1-1,k1,i2,j2-1,k2,m)
                        - comp_prim_efield(xyz,i1,j1-1,k1,i2,j2-1,k2,m+1));
    if (xyz==1) result += comp_prim_nuclear(i1,j1,k1,i2,j2-1,k2,m+1);
    }
  else if (k2) {
    result  = PmB[2] * comp_prim_efield(xyz,i1,j1,k1,i2,j2,k2-1,m);
    result -= PmC[2] * comp_prim_efield(xyz,i1,j1,k1,i2,j2,k2-1,m+1);
    if (k2>1) result += oo2zeta * (k2-1)
                       * (  comp_prim_efield(xyz,i1,j1,k1,i2,j2,k2-2,m)
                          - comp_prim_efield(xyz,i1,j1,k1,i2,j2,k2-2,m+1));
    if (k1) result += oo2zeta * k1
                     * (  comp_prim_efield(xyz,i1,j1,k1-1,i2,j2,k2-1,m)
                        - comp_prim_efield(xyz,i1,j1,k1-1,i2,j2,k2-1,m+1));
    if (xyz==2) result += comp_prim_nuclear(i1,j1,k1,i2,j2,k2-1,m+1);
    }
  else {
    /* We arrive here if we have a (s| |s) type efield integral.
     * The fjttable contains the standard (s| |s) nuc attr integrals.
     */
    result = 2.0 * zeta * PmC[xyz] * int_fjttable.d[m+1];
    }

  /* fprintf(stdout," % 12.8lf\n",result); */

  return result;
  }


/* --------------------------------------------------------------- */
/* ------------- Routines for dipole moment integrals ------------ */
/* --------------------------------------------------------------- */

/* This computes the dipole integrals between functions in two shells.
 * The result is accumulated in the buffer in the form bf1 x y z, bf2
 * x y z, etc.  The last arg, com, is the origin of the coordinate
 * system used to compute the dipole moment.
 */
GLOBAL_FUNCTION VOID
int_accum_shell_dipole(cs1,cs2,buff,ish,jsh,com)
centers_t *cs1;
centers_t *cs2;
double *buff;
int ish;
int jsh;
double *com;
{
  int c1,s1,i1,j1,k1,c2,s2,i2,j2,k2;
  int im,jm,km;
  int gc1,gc2;
  int index,index1,index2;
  double shell1norm,shell2norm;
  double dipole[3];

  C = com;

  c1 = cs1->center_num[ish];
  c2 = cs2->center_num[jsh];
  s1 = cs1->shell_num[ish];
  s2 = cs2->shell_num[jsh];
  A = cs1->center[c1].r;
  B = cs2->center[c2].r;
  shell1 = &(cs1->center[c1].basis.shell[s1]);
  shell2 = &(cs2->center[c2].basis.shell[s2]);
  index = 0;
  FOR_GCCART(gc1,index1,i1,j1,k1,shell1)
    shell1norm = shell1->norm[gc1][index1];
    FOR_GCCART(gc2,index2,i2,j2,k2,shell2)
      shell2norm = shell2->norm[gc2][index2];
      comp_shell_dipole(dipole,gc1,i1,j1,k1,gc2,i2,j2,k2);
      for(mu=0; mu < 3; mu++) {
        cartesianbuffer[index] = shell1norm * shell2norm * dipole[mu];
        index++;
        }
      END_FOR_GCCART(index2)
    END_FOR_GCCART(index1)
  int_accum_transform_1e_xyz(cartesianbuffer, buff, shell1, shell2);
  }

LOCAL_FUNCTION void
comp_shell_dipole(dipole,gc1,i1,j1,k1,gc2,i2,j2,k2)
double* dipole;
int gc1;
int i1;
int j1;
int k1;
int gc2;
int i2;
int j2;
int k2;
{
  double exp1,exp2;
  int i,j,xyz;
  double Pi;
  double oozeta;
  double AmB,AmB2;
  double tmp;

  dipole[0] = dipole[1] = dipole[2] = 0.0;

  if ((i1<0)||(j1<0)||(k1<0)||(i2<0)||(j2<0)||(k2<0)) return;

  /* Loop over the primitives in the shells. */
  for (i=0; i<shell1->nprim; i++) {
    for (j=0; j<shell2->nprim; j++) {

      /* Compute the intermediates. */
      exp1 = shell1->exp[i];
      exp2 = shell2->exp[j];
      oozeta = 1.0/(exp1 + exp2);
      oo2zeta = 0.5*oozeta;
      AmB2 = 0.0;
      for (xyz=0; xyz<3; xyz++) {
        Pi = oozeta*(exp1 * A[xyz] + exp2 * B[xyz]);
        PmA[xyz] = Pi - A[xyz];
        PmB[xyz] = Pi - B[xyz];
        PmC[xyz] = Pi - C[xyz];
        AmB = A[xyz] - B[xyz];
        AmB2 += AmB*AmB;
        }
      ss =   pow(3.141592653589793/(exp1+exp2),1.5)
           * exp_cutoff(- oozeta * exp1 * exp2 * AmB2);
      sMus = ss * PmC[mu];
      tmp     =  shell1->coef[gc1][i] * shell2->coef[gc2][j];
      if (exponent_weighted == 0) tmp *= exp1;
      else if (exponent_weighted == 1) tmp *= exp2;
      dipole[0] += tmp * comp_prim_dipole(1,0,0,i1,j1,k1,i2,j2,k2);
      dipole[1] += tmp * comp_prim_dipole(0,1,0,i1,j1,k1,i2,j2,k2);
      dipole[2] += tmp * comp_prim_dipole(0,0,1,i1,j1,k1,i2,j2,k2);
      }
    }

  }

LOCAL_FUNCTION double
comp_prim_dipole(im,jm,km,i1,j1,k1,i2,j2,k2)
int im;
int jm;
int km;
int i1;
int j1;
int k1;
int i2;
int j2;
int k2;
{
  double result;

  if (i1) {
    result = PmA[0] * comp_prim_dipole(im,jm,km,i1-1,j1,k1,i2,j2,k2);
    if (i2) 
      result += oo2zeta*i2*comp_prim_dipole(im,jm,km,i1-1,j1,k1,i2-1,j2,k2);
    if (i1>1)
      result += oo2zeta*(i1-1)*comp_prim_dipole(im,jm,km,i1-2,j1,k1,i2,j2,k2);
    if(im) result += oo2zeta*comp_prim_overlap(i1-1,j1,k1,i2,j2,k2);
    return result;
    }
  if (j1) {
    result = PmA[1] * comp_prim_dipole(im,jm,km,i1,j1-1,k1,i2,j2,k2);
    if (j2) 
      result += oo2zeta*j2*comp_prim_dipole(im,jm,km,i1,j1-1,k1,i2,j2-1,k2);
    if (j1>1) 
      result += oo2zeta*(j1-1)*comp_prim_dipole(im,jm,km,i1,j1-2,k1,i2,j2,k2);
    if(jm) result += oo2zeta*comp_prim_overlap(i1,j1-1,k1,i2,j2,k2);
    return result;
    }
  if (k1) {
    result = PmA[2] * comp_prim_dipole(im,jm,km,i1,j1,k1-1,i2,j2,k2);
    if (k2) 
      result += oo2zeta*k2*comp_prim_dipole(im,jm,km,i1,j1,k1-1,i2,j2,k2-1);
    if (k1>1) 
      result += oo2zeta*(k1-1)*comp_prim_dipole(im,jm,km,i1,j1,k1-2,i2,j2,k2);
    if(km) result += oo2zeta*comp_prim_overlap(i1,j1,k1-1,i2,j2,k2);
    return result;
    }
  if (i2) {
    result = PmB[0] * comp_prim_dipole(im,jm,km,i1,j1,k1,i2-1,j2,k2);
    if (i1) 
      result += oo2zeta*i1*comp_prim_dipole(im,jm,km,i1-1,j1,k1,i2-1,j2,k2);
    if (i2>1) 
      result += oo2zeta*(i2-1)*comp_prim_dipole(im,jm,km,i1,j1,k1,i2-2,j2,k2);
    if(im) result += oo2zeta*comp_prim_overlap(i1,j1,k1,i2-1,j2,k2);
    return result;
    }
  if (j2) {
    result = PmB[1] * comp_prim_dipole(im,jm,km,i1,j1,k1,i2,j2-1,k2);
    if (j1) 
      result += oo2zeta*i1*comp_prim_dipole(im,jm,km,i1,j1-1,k1,i2,j2-1,k2);
    if (j2>1) 
      result += oo2zeta*(j2-1)*comp_prim_dipole(im,jm,km,i1,j1,k1,i2,j2-2,k2);
    if(jm) result += oo2zeta*comp_prim_overlap(i1,j1,k1,i2,j2-1,k2);
    return result;
    }
  if (k2) {
    result = PmB[2] * comp_prim_dipole(im,jm,km,i1,j1,k1,i2,j2,k2-1);
    if (k1) 
      result += oo2zeta*i1*comp_prim_dipole(im,jm,km,i1,j1,k1-1,i2,j2,k2-1);
    if (k2>1) 
      result += oo2zeta*(k2-1)*comp_prim_dipole(im,jm,km,i1,j1,k1,i2,j2,k2-2);
    if(km) result += oo2zeta*comp_prim_overlap(i1,j1,k1,i2,j2,k2-1);
    return result;
    }

  return sMus;
  }
