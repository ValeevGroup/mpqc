/* $Log$
 * Revision 1.1  1993/12/29 12:53:03  etseidl
 * Initial revision
 *
 * Revision 1.3  1992/06/17  22:05:15  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.2  1992/03/31  01:24:06  jannsen
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

#include "atomssbc0.h"

void
sbcast0_shell(_shell,_type,_root)
shell_t *_shell;
int _type;
int _root;
{
  int nfunc;
  typedef int boolean;
  typedef char * string;
  int i;

  sbcast0_int(&(_shell->nprim),_type,_root, sizeof(int)*1);
  sbcast0_int(&(_shell->ncon),_type,_root, sizeof(int)*1);
  sbcast0_int(&(_shell->nfunc),_type,_root, sizeof(int)*1);

  if(_shell->nprim!=0) {
    sbcast0_pointer(&(_shell->exp),_type,_root,sizeof(double *));
    if(_shell->exp!=NULL) {
      sbcast0_double((_shell->exp),_type,_root,
        sizeof(double)*_shell->nprim);
      }
    }

  if(_shell->ncon!=0) {
    sbcast0_pointer(&(_shell->type),_type,_root,sizeof(shell_type_t *));
    if(_shell->type!=NULL) {
      for (i=0; i<_shell->ncon; i++)  {
        sbcast0_shell_type(&(_shell->type[i]),_type,_root);
        }
      }
    }

  if(_shell->ncon!=0) {
    sbcast0_pointer(&(_shell->coef),_type,_root,sizeof(double **));
    if(_shell->coef!=NULL) {
      for (i=0; i<_shell->ncon; i++)  {
        if(_shell->nprim!=0) {
          sbcast0_pointer(&(_shell->coef[i]),_type,_root,sizeof(double *));
          if(_shell->coef[i]!=NULL) {
            sbcast0_double((_shell->coef[i]),_type,_root,
              sizeof(double)*_shell->nprim);
            }
          }
        }
      }
    }

/* hand coded part for norm */
  if(_shell->ncon!=0) {
    sbcast0_pointer(&(_shell->norm),_type,_root,sizeof(double **));
    if(_shell->norm!=NULL) {
      for (i=0; i<_shell->ncon; i++)  {
        if((nfunc=INT_NCART(_shell->type[i].am))!=0) {
          sbcast0_pointer(&(_shell->norm[i]),_type,_root,sizeof(double *));
          if(_shell->norm[i]!=NULL) {
            sbcast0_double((_shell->norm[i]),_type,_root,sizeof(double)*nfunc);
            }
          }
        }
      }
    }

  return;
  }

