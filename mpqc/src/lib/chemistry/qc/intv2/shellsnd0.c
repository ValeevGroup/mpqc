/* $Log$
 * Revision 1.3  1994/08/26 22:45:51  etseidl
 * fix a bunch of warnings, get rid of rcs id's, get rid of bread/bwrite and
 * fread/fwrite modules
 *
 * Revision 1.2  1993/12/30  13:33:07  etseidl
 * mostly rcs id stuff
 *
 * Revision 1.3  1992/06/17  22:05:17  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.2  1992/03/31  01:24:22  jannsen
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.1  1992/01/17  22:12:11  cljanss
 * Initial revision
 *
 * Revision 1.1  1992/01/17  12:47:20  seidl
 * Initial revision
 * */

#include <stdio.h>
#include <math.h>
#include <tmpl.h>
#include <util/sgen/sgen.h>

#include "int_macros.h"

#include "atoms.h"

#include "atomssnd0.h"

void
send0_shell(_shell,_type,_dest)
shell_t *_shell;
int _type;
int _dest;
{
  int nfunc;
  typedef int boolean;
  typedef char * string;
  int i;

  send0_int(&(_shell->nprim),_type,_dest, sizeof(int)*1);
  send0_int(&(_shell->ncon),_type,_dest, sizeof(int)*1);
  send0_int(&(_shell->nfunc),_type,_dest, sizeof(int)*1);

  if(_shell->nprim!=0) {
    send0_pointer(&(_shell->exp),_type,_dest,sizeof(double *));
    if(_shell->exp!=NULL) {
      send0_double((_shell->exp),_type,_dest,
        sizeof(double)*_shell->nprim);
      }
    }

  if(_shell->ncon!=0) {
    send0_pointer(&(_shell->type),_type,_dest,sizeof(shell_type_t *));
    if(_shell->type!=NULL) {
      for (i=0; i<_shell->ncon; i++)  {
        send0_shell_type(&(_shell->type[i]),_type,_dest);
        }
      }
    }

  if(_shell->ncon!=0) {
    send0_pointer(&(_shell->coef),_type,_dest,sizeof(double **));
    if(_shell->coef!=NULL) {
      for (i=0; i<_shell->ncon; i++)  {
        if(_shell->nprim!=0) {
          send0_pointer(&(_shell->coef[i]),_type,_dest,sizeof(double *));
          if(_shell->coef[i]!=NULL) {
            send0_double((_shell->coef[i]),_type,_dest,
              sizeof(double)*_shell->nprim);
            }
          }
        }
      }
    }

/* hand coded part for norm */
  if(_shell->ncon!=0) {
    send0_pointer(&(_shell->norm),_type,_dest,sizeof(double **));
    if(_shell->norm!=NULL) {
      for (i=0; i<_shell->ncon; i++)  {
        if((nfunc=INT_NCART(_shell->type[i].am))!=0) {
          send0_pointer(&(_shell->norm[i]),_type,_dest,sizeof(double *));
          if(_shell->norm[i]!=NULL) {
            send0_double((_shell->norm[i]),_type,_dest,sizeof(double)*nfunc);
            }
          }
        }
      }
    }

  return;
  }

