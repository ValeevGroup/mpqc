
/* $Log$
 * Revision 1.1  1993/12/29 12:53:40  etseidl
 * Initial revision
 *
 * Revision 1.2  1992/03/30  22:34:07  jannsen
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.4  1991/09/28  21:05:31  cljanss
 * switch to new naming convention
 *
 * Revision 1.3  91/06/17  16:05:24  seidl
 * change calls to ioabort to bioabort
 * 
 * Revision 1.2  1991/06/16  20:30:55  seidl
 * fix $Log$
 * fix Revision 1.1  1993/12/29 12:53:40  etseidl
 * fix Initial revision
 * fix
 * Revision 1.2  1992/03/30  22:34:07  jannsen
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.4  1991/09/28  21:05:31  cljanss
 * switch to new naming convention
 *
 * Revision 1.3  91/06/17  16:05:24  seidl
 * change calls to ioabort to bioabort
 * 
 * */

static char *rcsid = "$Id$";

#include <tmpl.h>
#include <stdio.h>
#include "param.h"
#include "types.h"
#include "bioopen.gbl"

#include "errors.gbl"
#include "errors.lcl"

GLOBAL_FUNCTION VOID
bio_no_path_given(name)
char *name;
{
  fprintf(stderr,"%s: no path given\n",name);
  bioabort();
  }

GLOBAL_FUNCTION VOID
bio_malloc_check(caller,data)
char *caller;
VOID_PTR data;
{
  if (!data) {
    fprintf(stderr,"%s: malloc failed\n",caller);
    perror("malloc");
    bioabort();
    }
  }

GLOBAL_FUNCTION VOID
bio_fopen_check(caller,path,data)
char *caller;
char *path;
char *data;
{
  if (!data) {
    fprintf(stderr,"%s: fopen failed for %s\n",caller,path);
    perror("fopen");
    bioabort();
    }
  }

GLOBAL_FUNCTION VOID
bio_read_error(caller)
char *caller;
{
  fprintf(stderr,"%s: read failed\n",caller);
  perror("read");
  bioabort();
  }

GLOBAL_FUNCTION VOID
bio_write_error(caller)
char *caller;
{
  fprintf(stderr,"%s: write failed\n",caller);
  perror("write");
  bioabort();
  }
