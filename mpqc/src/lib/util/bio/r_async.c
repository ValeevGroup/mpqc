
/* $Log$
 * Revision 1.1  1993/12/29 12:53:41  etseidl
 * Initial revision
 *
 * Revision 1.3  1992/06/17  21:42:16  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.2  1992/03/30  22:44:30  jannsen
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.2  1991/09/28  21:05:34  cljanss
 * switch to new naming convention
 *
 * Revision 1.1  91/06/16  20:28:41  seidl
 * Initial revision
 *  */

static char *rcsid = "$Id$";

#include <stdio.h>
#include <tmpl.h>
#include "param.h"
#include "types.h"
#include "bioopen.gbl"

#include "r_async.gbl"
#include "r_async.lcl"

GLOBAL_FUNCTION r_async_t *
r_async_ioopen(param,unit,name)
char *param;
int unit;
char *name;
{
  fprintf(stderr,"no r_async io yet\n");
  bioabort();
  return NULL;
  }

GLOBAL_FUNCTION VOID
r_async_ioclos(ud,status)
r_async_t *ud;
int status;
{
  fprintf(stderr,"no r_async io yet\n");
  bioabort();
  }

GLOBAL_FUNCTION VOID
r_async_iordr(ud,buffer,first,length)
r_async_t *ud;
VOID_PTR buffer;
int first;
int length;
{
  fprintf(stderr,"no r_async io yet\n");
  bioabort();
  }

GLOBAL_FUNCTION VOID
r_async_iowrr(ud,buffer,first,length)
r_async_t *ud;
VOID_PTR buffer;
int first;
int length;
{
  fprintf(stderr,"no r_async io yet\n");
  bioabort();
  }
