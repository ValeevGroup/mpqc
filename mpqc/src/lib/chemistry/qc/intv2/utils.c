
/* $Log$
 * Revision 1.6  1996/03/23 02:37:51  cljanss
 * Everything can now be configured with autoconf.
 *
 * Revision 1.5  1995/10/25 21:20:01  cljanss
 * Adding support for pure am.  Gradients don't yet work.
 *
 * Revision 1.4  1995/03/17  01:49:44  cljanss
 * Removed -I. and -I$(SRCDIR) from the default include path in
 * GlobalMakefile to avoid name conflicts with system include files.
 * Modified files under src.lib to include all files relative to src.lib.
 * Makefiles under src.bin need to add the -I. and -I$(SRCDIR) back onto
 * INCLUDE and CXXINCLUDE or make other arrangements.
 *
 * Revision 1.3  1994/08/26  22:45:58  etseidl
 * fix a bunch of warnings, get rid of rcs id's, get rid of bread/bwrite and
 * fread/fwrite modules
 *
 * Revision 1.2  1993/12/30  13:33:11  etseidl
 * mostly rcs id stuff
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>

#include <chemistry/qc/intv2/atoms.h>
#include <chemistry/qc/intv2/int_flags.h>
#include <chemistry/qc/intv2/int_macros.h>
#include <chemistry/qc/intv2/int_types.h>

#include <chemistry/qc/intv2/utils.gbl>
#include <chemistry/qc/intv2/utils.lcl>

/* Compute the number of cartesian functions in a shell. */
GLOBAL_FUNCTION int
int_ncart(sh)
    shell_t *sh;
{
  int i;
  int ret = 0;
  for (i=0; i<sh->ncon; i++) {
      ret += INT_NCART(sh->type[i].am);
    }
  return ret;
}

GLOBAL_FUNCTION int
int_find_nfuncmax(cs)
centers_t *cs;
{
  return int_find_nfuncmax_aminc(cs,0);
  }

GLOBAL_FUNCTION int
int_find_ncartmax(cs)
centers_t *cs;
{
  return int_find_ncartmax_aminc(cs,0);
  }

GLOBAL_FUNCTION int
int_find_nconmax(cs)
centers_t *cs;
{
  int i,j,ncon;

  ncon = 0;
  for (i=0; i<cs->n; i++) {
    for (j=0; j<cs->center[i].basis.n; j++) {
      if (cs->center[i].basis.shell[j].ncon > ncon)
        ncon = cs->center[i].basis.shell[j].ncon;
      }
    }

  return ncon;
  }

/* This is given a centers structure and contraction number.
 * It returns the maximum angular momentum for a basis function
 * with these constraints. */
GLOBAL_FUNCTION int
int_find_jmax_for_con(cs,con)
centers_t *cs;
int con;
{
  int i,j;
  int jmax = 0;
  for (i=0; i<cs->n; i++) {
    for (j=0; j<cs->center[i].basis.n; j++) {
      if (cs->center[i].basis.shell[j].ncon > con) {
        if (jmax < cs->center[i].basis.shell[j].type[con].am) {
          jmax = cs->center[i].basis.shell[j].type[con].am;
          }
        }
      }
    }
  return jmax;
  }

GLOBAL_FUNCTION int
int_find_jmax(cs)
centers_t *cs;
{
  int i,j,k,jmax;

  jmax = 0;
  for (i=0; i<cs->n; i++) {
    for (j=0; j<cs->center[i].basis.n; j++) {
      for (k=0; k<cs->center[i].basis.shell[j].ncon; k++) {
        if (cs->center[i].basis.shell[j].type[k].am > jmax)
          jmax = cs->center[i].basis.shell[j].type[k].am;
        }
      }
    }

  return jmax;
  }

/* Finds the maximum number of functions in a shell given that the
 * angular momentum will be increment by aminc for all functions
 * in that shell. */
GLOBAL_FUNCTION int
int_find_nfuncmax_aminc(cs,aminc)
centers_t *cs;
int aminc;
{
  int i,j,k,nfuncmax,nfunc;

  nfuncmax = 0;
  for (i=0; i<cs->n; i++) {
    for (j=0; j<cs->center[i].basis.n; j++) {
      nfunc = 0;
      for (k=0; k<cs->center[i].basis.shell[j].ncon; k++) {
        nfunc += INT_NFUNC(cs->center[i].basis.shell[j].type[k].puream,
                           cs->center[i].basis.shell[j].type[k].am+aminc);
        }
      if ((aminc == 0) && (cs->center[i].basis.shell[j].nfunc != nfunc)) {
        fprintf(stderr,"error: nfunc not init'ed for a shell\n");
        fail();
        }
      if (nfunc > nfuncmax) nfuncmax = nfunc;
      }
    }

  return nfuncmax;
  }

/* Finds the maximum number of cartesian functions in a shell given that the
 * angular momentum will be increment by aminc for all functions
 * in that shell. */
GLOBAL_FUNCTION int
int_find_ncartmax_aminc(cs,aminc)
centers_t *cs;
int aminc;
{
  int i,j,k,ncartmax,ncart;

  ncartmax = 0;
  for (i=0; i<cs->n; i++) {
    for (j=0; j<cs->center[i].basis.n; j++) {
      ncart = 0;
      for (k=0; k<cs->center[i].basis.shell[j].ncon; k++) {
        ncart += INT_NCART(cs->center[i].basis.shell[j].type[k].am+aminc);
        }
      if (ncart > ncartmax) ncartmax = ncart;
      }
    }

  return ncartmax;
  }

GLOBAL_FUNCTION int
int_find_jmax_shell(shell)
shell_t *shell;
{
  int k,jmax;

  jmax = 0;
  for (k=0; k<shell->ncon; k++) {
    if (shell->type[k].am > jmax)
      jmax = shell->type[k].am;
    }
  return jmax;
  }

LOCAL_FUNCTION void
fail()
{
  fprintf(stderr,"failing module:\n%s\n",__FILE__);
  exit(1);
  }
