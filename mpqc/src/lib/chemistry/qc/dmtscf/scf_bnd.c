
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <tmpl.h>

#include <util/misc/libmisc.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <comm/picl/picl.h>
#include <comm/picl/ext/piclext.h>

#include "scf_bnd.gbl"
#include "scf_bnd.lcl"

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
  
  scf_bnd_Qvec = (signed char *) malloc(sizeof(signed char)*nsht);
  if (scf_bnd_Qvec==NULL) {
    fprintf(stderr,"scf_init_bounds: cannot malloc Qvec: %d\n",nsht);
    return -1;
  }

  memset(scf_bnd_Qvec,'\0',sizeof(signed char)*nsht);

  compute_Q(cs1,intbuf);

  gop0_sc(scf_bnd_Qvec,nsht,'+',mtype_get());

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
  free(scf_bnd_Qvec);
  scf_bnd_Qvec = 0;
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

  return((int) scf_bnd_Qvec[ij]+scf_bnd_Qvec[kl]);
}

/* ripped off from clj's libintv2 */

LOCAL_FUNCTION VOID
compute_Q(cs1,intbuf)
centers_t *cs1;
double *intbuf;
{
  int shellij;
  int sh1,sh2;
  int shells[4],size[4];
  int nfunc1,nfunc2;
  double max;
  double tol = pow(2.0,-126.0);
  double loginv = 1.0/log(2.0);
  double integral;
  int i,j;
  int me=mynode0();
  int nproc=numnodes0();

  shellij=0;
  for (sh1=0; sh1 < cs1->nshell ; sh1++) {
    shells[0]=shells[2]=sh1;
    for (sh2=0; sh2 <= sh1 ; sh2++,shellij++) {
      if (shellij%nproc != me) continue;

      shells[1]=shells[3]=sh2;

      int_erep_v(INT_EREP|INT_REDUND|INT_NOPERM|INT_NOBCHK,shells,size);

      nfunc1=size[0];
      nfunc2=size[1];

    /* Find the biggest (ij|ij) integral. */
      max = 0.0;
      for (i=0; i<nfunc1; i++) {
        for (j=0; j<nfunc2; j++) {
          int index = i * nfunc2 * nfunc1 * nfunc2
                    + j * nfunc1 * nfunc2
                    + i * nfunc2
                    + j;
          integral = intbuf[ index ];
          if (fabs(integral) > max) max = fabs(integral);
        }
      }

    /* Compute the Q value. */
      max = sqrt(max);
      max = (max>tol) ? max : tol;
      scf_bnd_Qvec[shellij] = (signed char) (log(max)*loginv);
    }
  }
}
