/* $Log$
 * Revision 1.2  1993/12/30 13:33:05  etseidl
 * mostly rcs id stuff
 *
 * Revision 1.3  1992/06/17  22:05:10  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.2  1992/03/31  01:23:56  jannsen
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.1  1992/01/17  22:12:11  cljanss
 * Initial revision
 *
 * Revision 1.1  1992/01/17  12:47:20  seidl
 * Initial revision
 * */

static char rcsid[]="$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tmpl.h>
#include <util/sgen/sgen.h>

#include "int_macros.h"

#include "atoms.h"

#include "atomsrcv0.h"

void
recv0_shell(_shell,_type,_from)
shell_t *_shell;
int _type;
int _from;
{
  int nfunc;
  typedef int boolean;
  typedef char * string;
  int * recv0_test_pointer();
  int i;

  recv0_int(&(_shell->nprim),_type,_from, sizeof(int)*1);
  recv0_int(&(_shell->ncon),_type,_from, sizeof(int)*1);
  recv0_int(&(_shell->nfunc),_type,_from, sizeof(int)*1);

  if(_shell->nprim!=0) {
    if(recv0_test_pointer(_type,_from,sizeof(double *))!=NULL) {
      _shell->exp = (double *) malloc(sizeof(double )*_shell->nprim);
      bzero(_shell->exp,sizeof(double )*_shell->nprim);
      recv0_double((_shell->exp),_type,_from,
        sizeof(double)*_shell->nprim);
      }
    else _shell->exp = NULL;
    }
  else _shell->exp = NULL; /* DT. */

  if(_shell->ncon!=0) {
    if(recv0_test_pointer(_type,_from,sizeof(shell_type_t *))!=NULL) {
      _shell->type = 
        (shell_type_t *) malloc(sizeof(shell_type_t )*_shell->ncon);
      bzero(_shell->type,sizeof(shell_type_t )*_shell->ncon);
      for (i=0; i<_shell->ncon; i++)  {
        recv0_shell_type(&(_shell->type[i]),_type,_from);
        }
      }
    else _shell->type = NULL;
    }
  else _shell->type = NULL; /* DT. */

  if(_shell->ncon!=0) {
    if(recv0_test_pointer(_type,_from,sizeof(double **))!=NULL) {
      _shell->coef = (double **) malloc(sizeof(double *)*_shell->ncon);
      bzero(_shell->coef,sizeof(double *)*_shell->ncon);
      for (i=0; i<_shell->ncon; i++)  {
        if(_shell->nprim!=0) {
          if(recv0_test_pointer(_type,_from,sizeof(double *))!=NULL) {
            _shell->coef[i] = (double *) malloc(sizeof(double )*_shell->nprim);
            bzero(_shell->coef[i],sizeof(double )*_shell->nprim);
            recv0_double((_shell->coef[i]),_type,_from,
              sizeof(double)*_shell->nprim);
            }
          }
        else _shell->coef[i] = NULL;
        }
      /* Skipped: else _shell->coef[i] = NULL;*/ /* DT. */
      }
    else _shell->coef = NULL;
    }
  else _shell->coef = NULL; /* DT. */

/* hand coded part for norm */
  if(_shell->ncon!=0) {
    if(recv0_test_pointer(_type,_from,sizeof(double **))!=NULL) {
      _shell->norm = (double **) malloc(sizeof(double *)*_shell->ncon);
      bzero(_shell->norm,sizeof(double *)*_shell->ncon);
      for (i=0; i<_shell->ncon; i++)  {
        if((nfunc=INT_NCART(_shell->type[i].am))!=0) {
          if(recv0_test_pointer(_type,_from,sizeof(double *))!=NULL) {
            _shell->norm[i] = (double *) malloc(sizeof(double )*nfunc);
            bzero(_shell->norm[i],sizeof(double )*nfunc);
            recv0_double((_shell->norm[i]),_type,_from,sizeof(double)*nfunc);
            }
          }
        else _shell->norm[i] = NULL;
        }
      /* Skipped: else _shell->norm[i] = NULL;*/ /* DT. */
      }
    else _shell->norm = NULL;
    }
  else _shell->norm = NULL; /* DT. */

  return;
  }
