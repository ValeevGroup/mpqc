/* $Log$
 * Revision 1.1  1993/12/29 12:53:02  etseidl
 * Initial revision
 *
 * Revision 1.3  1992/06/17  22:04:59  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.2  1992/03/31  01:23:04  jannsen
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

#include "atomsbrd.h"

int
bread_shell(_unit,_shell,_offset)
int _unit;
shell_t *_shell;
int *_offset;
{
  int nfunc;
  typedef int boolean;
  typedef char * string;
  int * bread_test_pointer();
  int i;
  int _offset_init= *_offset;

  bread_int(_unit,&(_shell->nprim),_offset, sizeof(int)*1);
  bread_int(_unit,&(_shell->ncon),_offset, sizeof(int)*1);
  bread_int(_unit,&(_shell->nfunc),_offset, sizeof(int)*1);

  if(_shell->nprim!=0) {
    if(bread_test_pointer(_unit,_offset,sizeof(double *))!=NULL) {
      _shell->exp = (double *) malloc(sizeof(double )*_shell->nprim);
      sgen_chkmalloc(_shell->exp);
      bzero(_shell->exp,sizeof(double )*_shell->nprim);
      bread_double(_unit,(_shell->exp),_offset,
        sizeof(double)*_shell->nprim);
      }
    else _shell->exp = NULL;
    }
  else _shell->exp = NULL; /* DT. */

  if(_shell->ncon!=0) {
    if(bread_test_pointer(_unit,_offset,sizeof(shell_type_t *))!=NULL) {
      _shell->type = (shell_type_t *) 
                      malloc(sizeof(shell_type_t )*_shell->ncon);
      sgen_chkmalloc(_shell->type);
      bzero(_shell->type,sizeof(shell_type_t )*_shell->ncon);
      for (i=0; i<_shell->ncon; i++)  {
        bread_shell_type(_unit,&(_shell->type[i]),_offset);
        }
      }
    else _shell->type = NULL;
    }
  else _shell->type = NULL; /* DT. */

  if(_shell->ncon!=0) {
    if(bread_test_pointer(_unit,_offset,sizeof(double **))!=NULL) {
      _shell->coef = (double **) malloc(sizeof(double *)*_shell->ncon);
      sgen_chkmalloc(_shell->coef);
      bzero(_shell->coef,sizeof(double *)*_shell->ncon);
      for (i=0; i<_shell->ncon; i++)  {
        if(_shell->nprim!=0) {
          if(bread_test_pointer(_unit,_offset,sizeof(double *))!=NULL) {
            _shell->coef[i] = (double *) malloc(sizeof(double )*_shell->nprim);
            sgen_chkmalloc(_shell->coef[i]);
            bzero(_shell->coef[i],sizeof(double )*_shell->nprim);
            bread_double(_unit,(_shell->coef[i]),_offset,
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
    if(bread_test_pointer(_unit,_offset,sizeof(double **))!=NULL) {
      _shell->norm = (double **) malloc(sizeof(double *)*_shell->ncon);
      sgen_chkmalloc(_shell->norm);
      bzero(_shell->norm,sizeof(double *)*_shell->ncon);
      for (i=0; i<_shell->ncon; i++)  {
        if((nfunc=INT_NCART(_shell->type[i].am))!=0) {
          if(bread_test_pointer(_unit,_offset,sizeof(double *))!=NULL) {
            _shell->norm[i] = (double *) malloc(sizeof(double )*nfunc);
            sgen_chkmalloc(_shell->norm[i]);
            bzero(_shell->norm[i],sizeof(double )*nfunc);
            bread_double(_unit,(_shell->norm[i]),_offset,sizeof(double)*nfunc);
            }
          }
        else _shell->norm[i] = NULL;
        }
      /* Skipped: else _shell->norm[i] = NULL;*/ /* DT. */
      }
    else _shell->norm = NULL;
    }
  else _shell->norm = NULL; /* DT. */

  return(*_offset-_offset_init);
  }
