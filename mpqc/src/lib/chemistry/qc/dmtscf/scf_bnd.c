
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <tmpl.h>

#include <util/misc/libmisc.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <util/group/picl.h>

#include <chemistry/qc/dmtscf/scf_bnd.gbl>
#include <chemistry/qc/dmtscf/scf_bnd.lcl>

/* this is a globally accessible variable on purpose */
signed char *scf_bnd_Qvec;

/***********************************************************************
 * this initialized the bounds information.  given a centers struct which
 * has already been passed to int_initialize_erep() and the buffer used
 * to store the integrals, create scf_dmt_Qvec[ij] = log_2(max(ij|ij))
 * we use the log so that we can fit the bounds info into a char rather than
 * a double.
 */

GLOBAL_FUNCTION int
scf_init_bounds(cs1,intbuf)
centers_t *cs1;
double *intbuf;
{
  int nshell=cs1->nshell;
  int nsht=nshell*(nshell+1)/2;
  int sh1,sh2,shellij;
  int nproc = numnodes0();
  int me = mynode0();

  int_init_bounds_nocomp();
  scf_bnd_Qvec = int_Qvec;

  for (shellij=0,sh1=0; sh1 < nshell ; sh1++) {
    for (sh2=0; sh2 <= sh1 ; sh2++,shellij++) {
      if (shellij%nproc != me) continue;
      int_bounds_comp(sh1,sh2);
      }
    }

  gop0_sc(int_Qvec,nsht,'+',mtype_get());

  return 0;
}

/*************************************************************************
 *
 * free up the memory used to store the bounds
 *
 */

GLOBAL_FUNCTION VOID
scf_done_bounds()
{
  int_done_bounds();
}

/*************************************************************************
 * 
 * given shells i, j, k, and l, return the value max(ij|ij)*max(kl|kl)
 *
 * scf_init_bound() must have been called before using this.
 */

GLOBAL_FUNCTION int
scf_erep_bound(i,j,k,l)
int i;
int j;
int k;
int l;
{
  int ij=(i>j) ? i*(i+1)/2+j : j*(j+1)/2+i;
  int kl=(k>l) ? k*(k+1)/2+l : l*(l+1)/2+k;

  return((int) int_Qvec[ij]+int_Qvec[kl]);
}
