/* $Log$
 * Revision 1.6  1995/10/25 21:19:58  cljanss
 * Adding support for pure am.  Gradients don't yet work.
 *
 * Revision 1.5  1995/03/17  01:49:38  cljanss
 * Removed -I. and -I$(SRCDIR) from the default include path in
 * GlobalMakefile to avoid name conflicts with system include files.
 * Modified files under src.lib to include all files relative to src.lib.
 * Makefiles under src.bin need to add the -I. and -I$(SRCDIR) back onto
 * INCLUDE and CXXINCLUDE or make other arrangements.
 *
 * Revision 1.4  1994/10/13  22:26:44  etseidl
 * replace bzero with memset
 *
 * Revision 1.3  1994/08/26  22:45:48  etseidl
 * fix a bunch of warnings, get rid of rcs id's, get rid of bread/bwrite and
 * fread/fwrite modules
 *
 * Revision 1.2  1993/12/30  13:33:05  etseidl
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <tmpl.h>
#include <util/sgen/sgen.h>

#include <chemistry/qc/intv2/int_macros.h>

#include <chemistry/qc/intv2/atoms.h>

#include <chemistry/qc/intv2/atomsrcv0.h>

void
recv0_shell(_shell,_type,_from)
shell_t *_shell;
int _type;
int _from;
{
  int ncart;
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
      memset(_shell->exp,'\0',sizeof(double )*_shell->nprim);
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
      memset(_shell->type,'\0',sizeof(shell_type_t )*_shell->ncon);
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
      memset(_shell->coef,'\0',sizeof(double *)*_shell->ncon);
      for (i=0; i<_shell->ncon; i++)  {
        if(_shell->nprim!=0) {
          if(recv0_test_pointer(_type,_from,sizeof(double *))!=NULL) {
            _shell->coef[i] = (double *) malloc(sizeof(double )*_shell->nprim);
            memset(_shell->coef[i],'\0',sizeof(double )*_shell->nprim);
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
      memset(_shell->norm,'\0',sizeof(double *)*_shell->ncon);
      for (i=0; i<_shell->ncon; i++)  {
        if((ncart=INT_NCART(_shell->type[i].am))!=0) {
          if(recv0_test_pointer(_type,_from,sizeof(double *))!=NULL) {
            _shell->norm[i] = (double *) malloc(sizeof(double )*ncart);
            memset(_shell->norm[i],'\0',sizeof(double )*ncart);
            recv0_double((_shell->norm[i]),_type,_from,sizeof(double)*ncart);
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
