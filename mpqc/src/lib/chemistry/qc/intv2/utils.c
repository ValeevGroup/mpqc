
/* $Log$
 * Revision 1.2  1993/12/30 13:33:11  etseidl
 * mostly rcs id stuff
 *
 */
static char *rcsid = "$Id$";


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>

#include "atoms.h"
#include "int_flags.h"
#include "int_macros.h"
#include "int_types.h"

#include "utils.gbl"
#include "utils.lcl"


GLOBAL_FUNCTION int
int_find_nfuncmax(cs)
centers_t *cs;
{
  return int_find_nfuncmax_aminc(cs,0);
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
        nfunc += INT_NCART(cs->center[i].basis.shell[j].type[k].am+aminc);
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

LOCAL_FUNCTION VOID
fail()
{
  fprintf(stderr,"failing module:\n%s\n",rcsid);
  exit(1);
  }
