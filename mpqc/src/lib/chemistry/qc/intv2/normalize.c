
/* $Log$
 * Revision 1.6  1995/11/16 00:47:38  cljanss
 * Removed normalization for individual basis functions.
 *
 * Revision 1.5  1995/10/25 21:19:54  cljanss
 * Adding support for pure am.  Gradients don't yet work.
 *
 * Revision 1.4  1995/03/17  01:49:34  cljanss
 * Removed -I. and -I$(SRCDIR) from the default include path in
 * GlobalMakefile to avoid name conflicts with system include files.
 * Modified files under src.lib to include all files relative to src.lib.
 * Makefiles under src.bin need to add the -I. and -I$(SRCDIR) back onto
 * INCLUDE and CXXINCLUDE or make other arrangements.
 *
 * Revision 1.3  1994/08/26  22:45:39  etseidl
 * fix a bunch of warnings, get rid of rcs id's, get rid of bread/bwrite and
 * fread/fwrite modules
 *
 * Revision 1.2  1993/12/30  13:32:56  etseidl
 * mostly rcs id stuff
 *
 * Revision 1.2  1992/06/17  22:04:52  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.1.1.1  1992/03/17  16:33:16  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:33:15  seidl
 * Initial revision
 *
 * Revision 1.1  1992/02/03  15:49:49  seidl
 * Initial revision
 *
 * Revision 1.2  1992/01/08  17:15:36  cljanss
 * don't exit if shell normalization info exists
 *
 * Revision 1.1  1991/12/14  00:19:36  cljanss
 * Initial revision
 * */
/*
 * Revision 1.2  91/09/28  19:26:53  cljanss
 * Switch to new file naming scheme
 * 
 * Revision 1.1  91/06/16  16:40:07  janssen
 * Initial revision
 *  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/atoms.h>
#include <chemistry/qc/intv2/int_macros.h>

#include <chemistry/qc/intv2/atomsprnt.h>

#include <chemistry/qc/intv2/normalize.gbl>
#include <chemistry/qc/intv2/normalize.lcl>
#include <chemistry/qc/intv2/comp_1e.gbl>

/* This routine must be called before anything else, even the offset
 * initializing routines (because this produces the shell nfunc info). */
GLOBAL_FUNCTION VOID
int_initialize_centers(centers)
centers_t *centers;
{
  int i,j;

  for (i=0; i<centers->n; i++) {
    for (j=0; j<centers->center[i].basis.n; j++) {
      shell_t *shell = &centers->center[i].basis.shell[j];
      int k;
      shell->nfunc = 0;
      for (k=0; k<shell->ncon; k++)
        shell->nfunc += INT_NFUNC(shell->type[k].puream,shell->type[k].am);
      int_normalize_shell(shell);
      }
    }
  }

/* This converts the contraction coefficients
 * from contraction coefficients for normalized primitives
 * to contraction coefficients for unnormalized primitives.
 */
GLOBAL_FUNCTION VOID
int_normalize_shell(shell)
shell_t *shell;
{
  int gc;
  int i,j,k;
  int index;
  double bfnormxl;
  double normalization;
  double c,ss;

  for (gc=0; gc<shell->ncon; gc++) {
    /* Convert the contraction coefficients
     * from contraction coefficients for normalized primitives
     * to contraction coefficients for unnormalized primitives.
     */
    for (i=0; i<shell->nprim; i++) {
      c = 0.25/shell->exp[i];
      ss = pow(3.141592653589793/(shell->exp[i]+shell->exp[i]),1.5);

#if 0
      printf("shell->coef[%d][%d] = before: % f",gc,i,shell->coef[gc][i]);
#endif

      shell->coef[gc][i]
        *= 1.0/sqrt(norm(shell->type[gc].am,shell->type[gc].am,c,ss));

#if 0
      printf(", after: % f\n",shell->coef[gc][i]);
#endif
      }

    normalization = int_shell_normalization(shell,gc);

    /* Adjust the contraction coefficients so that they are normalized. */
    for (i=0; i<shell->nprim; i++) {
      shell->coef[gc][i] *= normalization;
      }
    }
  }

/* Return the basis function dependent part of the norm.
 * It depends only on the power of x, y, and z in the Cartesian
 * basis function. */
LOCAL_FUNCTION double
bfnorm(i,j,k)
int i;
int j;
int k;
{
  return 1.0/(sqrt((double)  factfact(2*i-1)
                           * factfact(2*j-1)
                           * factfact(2*k-1)));
  }

LOCAL_FUNCTION long
factfact(n)
int n;
{
  long result;
  int i;

  result = 1;
  for (i=3; i<=n; i+=2) {
    result *= i;
    }
  return result;
  }

/* Compute the normalization constant for a shell.
 * returns 1/sqrt(<(x^l 0 0|(x^l 0 0)>).
 * The formula is from Obara and Saika (for the basis functions within
 * the shell that have powers of x only (a and b refer to the power
 * of x):
 * (a||b) = 1/(4 alpha) * ( a (a-1||b) + b (a||b-1) )
 */
GLOBAL_FUNCTION double
int_shell_normalization(shell,gc)
shell_t *shell;
int gc;
{
  int i,j;
  double result,c,ss;

  result = 0.0;
  for (i=0; i<shell->nprim; i++) {
    for (j=0; j<shell->nprim; j++) {
      c = 0.50/(shell->exp[i] + shell->exp[j]);
      ss = pow(3.141592653589793/(shell->exp[i]+shell->exp[j]),1.5);
      result += shell->coef[gc][i] * shell->coef[gc][j] *
               norm(shell->type[gc].am,shell->type[gc].am,c,ss);
      }
    }

  return 1.0/sqrt(result);
  }

/* Compute the norm for ((x^x1)||(x^x2)).  This is slower than need be. */
LOCAL_FUNCTION double
norm(x1,x2,c,ss)
int x1;
int x2;
double c;
double ss;
{
  if (x1 < x2) return norm(x2,x1,c,ss);
  if (x1 == 1) {
    if (x2 == 1) return c * ss;
    else return 0.0;
    }
  if (x1 == 0) return ss;
  return c * ( (x1-1) * norm(x1-2,x2,c,ss) + (x2 * norm(x1-1,x2-1,c,ss)));
  }

