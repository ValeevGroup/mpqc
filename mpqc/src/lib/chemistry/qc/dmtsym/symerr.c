/* $Log$
 * Revision 1.1  1993/12/29 12:52:58  etseidl
 * Initial revision
 *
 * Revision 1.1.1.1  1992/03/17  16:27:42  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:27:39  seidl
 * Initial revision
 *
 * Revision 1.1  1992/01/27  12:53:04  seidl
 * Initial revision
 *
 * Revision 1.1  1991/12/20  15:45:23  seidl
 * Initial revision
 *
 * Revision 1.1  1991/12/03  16:49:43  etseidl
 * Initial revision
 *
 * Revision 1.1  1991/11/22  18:28:50  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#include <stdio.h>
#include <tmpl.h>

#include "symerr.gbl"
#include "symerr.lcl"

GLOBAL_FUNCTION VOID
serror(_outfile,_sfile,_serrmsg,_slineno)
FILE *_outfile;
char *_sfile;
char *_serrmsg;
int _slineno;
{
  fprintf(_outfile,"%s: line number %d\n",_sfile,_slineno);
  fprintf(_outfile,"%s\n",_serrmsg);
  fflush(_outfile);
  }
