/* $Log$
 * Revision 1.1  1993/12/29 12:53:02  etseidl
 * Initial revision
 *
 * Revision 1.3  1992/06/17  22:05:04  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.2  1992/03/31  01:23:30  jannsen
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
#include <math.h>
#include <tmpl.h>
#include <util/sgen/sgen.h>

#include "int_macros.h"

#include "atoms.h"

#include "atomsfwr.h"

int
fwrite_shell(_fp,_shell,_offset)
FILE *_fp;
shell_t *_shell;
int *_offset;
{
  int nfunc;
  typedef int boolean;
  typedef char * string;
  int i;
  int _offset_init= *_offset;

  fwrite_int(_fp,&(_shell->nprim),_offset, sizeof(int)*1);
  fwrite_int(_fp,&(_shell->ncon),_offset, sizeof(int)*1);
  fwrite_int(_fp,&(_shell->nfunc),_offset, sizeof(int)*1);

  if(_shell->nprim!=0) {
    fwrite_pointer(_fp,&(_shell->exp),_offset,sizeof(double *));
    if(_shell->exp!=NULL) {
      fwrite_double(_fp,(_shell->exp),_offset,
        sizeof(double)*_shell->nprim);
      }
    }

  if(_shell->ncon!=0) {
    fwrite_pointer(_fp,&(_shell->type),_offset,sizeof(shell_type_t *));
    if(_shell->type!=NULL) {
      for (i=0; i<_shell->ncon; i++)  {
        fwrite_shell_type(_fp,&(_shell->type[i]),_offset);
        }
      }
    }

  if(_shell->ncon!=0) {
    fwrite_pointer(_fp,&(_shell->coef),_offset,sizeof(double **));
    if(_shell->coef!=NULL) {
      for (i=0; i<_shell->ncon; i++)  {
        if(_shell->nprim!=0) {
          fwrite_pointer(_fp,&(_shell->coef[i]),_offset,sizeof(double *));
          if(_shell->coef[i]!=NULL) {
            fwrite_double(_fp,(_shell->coef[i]),_offset,
              sizeof(double)*_shell->nprim);
            }
          }
        }
      }
    }

/* hand coded part for norm */
  if(_shell->ncon!=0) {
    fwrite_pointer(_fp,&(_shell->norm),_offset,sizeof(double **));
    if(_shell->norm!=NULL) {
      for (i=0; i<_shell->ncon; i++)  {
        if((nfunc=INT_NCART(_shell->type[i].am))!=0) {
          fwrite_pointer(_fp,&(_shell->norm[i]),_offset,sizeof(double *));
          if(_shell->norm[i]!=NULL) {
            fwrite_double(_fp,(_shell->norm[i]),_offset,sizeof(double)*nfunc);
            }
          }
        }
      }
    }

  return(*_offset-_offset_init);
  }
