
/* $Log$
 * Revision 1.3  1996/03/23 02:38:54  cljanss
 * Everything can now be configured with autoconf.
 *
 * Revision 1.2  1994/10/18 23:03:47  etseidl
 * fix many warnings, use memset rather than bzero
 *
 * Revision 1.1.1.1  1993/12/29  12:53:57  etseidl
 * SC source tree 0.1
 *
 * Revision 1.1.1.1  1992/03/17  18:10:14  seidl
 * Struct GENerator 2.0
 *
 * Revision 1.1  1992/03/17  18:10:14  seidl
 * Initial revision
 *
 * Revision 1.2  1992/01/15  17:19:20  seidl
 * apply clj's patches
 *
 * Revision 1.4  1991/12/14  00:12:31  cljanss
 * check __STDC__ to see how functions are to be defined
 *
 * Revision 1.3  1991/11/18  17:21:10  cljanss
 * The variable arguments error and warning routines do not work on the sun
 *
 * Revision 1.2  91/09/28  16:40:26  cljanss
 * new naming convention is used to keep names <= 14 characters.
 * 
 * Revision 1.1  91/06/15  21:13:57  janssen
 * Initial revision
 *  */

#include <stdio.h>
#if !defined(SUN)
#include <stdarg.h>
#endif
#include "types.h"
#include "global.h"
#include "scan.gbl"

#include "error.gbl"
#include "error.lcl"

#if !defined(SUN)
GLOBAL_VA_FUNCTION void
#ifndef __STDC__
error(msg)
char *msg;
#else
error(char *msg, ...)
#endif
{
  va_list args;
  va_start(args,msg);
  fprintf(stderr,"ERROR: ");
  vfprintf(stderr,msg,args);
  fprintf(stderr,"\n");
  va_end(args);
  showpos();
  exit(1);
  }

GLOBAL_VA_FUNCTION void
#ifndef __STDC__
warn(msg)
char *msg;
#else
warn(char *msg, ...)
#endif
{
  va_list args;
  va_start(args,msg);
  fprintf(stderr,"WARNING: ");
  vfprintf(stderr,msg,args);
  fprintf(stderr,"\n");
  va_end(args);
  }
#else /* ! SUN */
GLOBAL_VA_FUNCTION void
#ifndef __STDC__
error(msg)
char *msg;
#else
error(char *msg, ...)
#endif
{
  fprintf(stderr,"ERROR: %s\n",msg);
  showpos();
  exit(1);
  }

GLOBAL_VA_FUNCTION void
#ifndef __STDC__
warn(msg)
char *msg;
#else
warn(char *msg, ...)
#endif
{
  fprintf(stderr,"WARNING: %s\n",msg);
  }
#endif /* ! SUN */

