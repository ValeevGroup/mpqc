/* $Log$
 * Revision 1.2  1994/04/11 16:57:24  etseidl
 * parallelize the initialization of the bounds information
 *
 * Revision 1.1.1.1  1993/12/29  12:53:15  etseidl
 * SC source tree 0.1
 *
 * Revision 1.6  1992/06/17  21:54:03  jannsen
 * cleaned up for saber-c
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

static char rcsid[]="$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <comm/picl/picl.h>

#include "scf_bnd.gbl"
#include "scf_bnd.lcl"

signed char *Qvec;

GLOBAL_FUNCTION VOID
scf_init_bounds(cs1,intbuf)
centers_t *cs1;
double *intbuf;
{
  int nshell=cs1->nshell;
  int nsht=nshell*(nshell+1)/2;
  
  Qvec = (signed char *) malloc(sizeof(signed char)*nsht);
  if(Qvec==NULL) {
    fprintf(stderr,"scf_init_bounds: cannot malloc Qvec: %d\n",nsht);
    exit(1);
    }

  bzero(Qvec,sizeof(signed char)*nsht);

  compute_Q(cs1,intbuf);

  gsum0(Qvec,nsht,0,mtype_get(),0);
  }

GLOBAL_FUNCTION VOID
scf_done_bounds()
{
  free(Qvec);
  }

GLOBAL_FUNCTION int
scf_erep_bound(s1,s2,s3,s4)
int s1;
int s2;
int s3;
int s4;
{
  int ij=(s1>s2) ? s1*(s1+1)/2+s2 : s2*(s2+1)/2+s1;
  int kl=(s3>s4) ? s3*(s3+1)/2+s4 : s4*(s4+1)/2+s3;

  return((int) Qvec[ij]+Qvec[kl]);
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
  for(sh1=0; sh1 < cs1->nshell ; sh1++) {
    shells[0]=shells[2]=sh1;
    for(sh2=0; sh2 <= sh1 ; sh2++,shellij++) {
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
      Qvec[shellij] = (signed char) (log(max)*loginv);
      }
    }
  }
