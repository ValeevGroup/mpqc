
/* $Log$
 * Revision 1.1  1993/12/29 12:53:40  etseidl
 * Initial revision
 *
 * Revision 1.2  1992/06/17  21:42:12  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.1  1992/04/06  12:26:29  seidl
 * Initial revision
 *
 * Revision 1.2  1992/03/30  22:33:50  jannsen
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.6  1991/12/02  17:18:54  cljanss
 * clean up use of VOID and VOID_PTR for old compilers
 *
 * Revision 1.5  91/11/18  17:34:09  cljanss
 * changed stderr to bio_error
 * 
 * Revision 1.4  91/09/28  21:05:30  cljanss
 * switch to new naming convention
 * 
 * Revision 1.3  91/09/28  20:49:52  cljanss
 * changed the global name ptr to bio_ptr
 * 
 * Revision 1.2  1991/06/17  16:05:00  seidl
 * first working version
 *
 * Revision 1.1  1991/06/16  20:28:41  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#include <tmpl.h>
#include <string.h>
#include <stdio.h>
#include "param.h"
#include "types.h"
#include "pointers.h"

#include "file_info.gbl"
#include "r_async.gbl"
#include "ram.gbl"
#include "s_async.gbl"
#include "sequential.gbl"

extern ioFILE_t bio_ud[MAX_UNIT];

extern unsigned int bio_current_unit;

extern void bioabort();

#include "bioinit.gbl"
#include "bioinit.lcl"

/* initializes binary file "unit" with name "name" */

GLOBAL_FUNCTION int
bio_init_file(unit,name)
int unit;
char *name;

{
  int i;

  if (bio_ud[unit].status != BIO_UNINITED) {
    fprintf(bio_error,"bio_init_file: unit %d is already in use\n");
    return -1;
    }

  if(unit==0) {
    i=0;
    do {
      i++;
      } while(i < MAX_UNIT && bio_ud[i].status!=BIO_UNINITED);
    unit=i;
    }

  if(unit >= MAX_UNIT) {
    fprintf(bio_error,"bio_init_file: invalid unit number %d\n",unit);
    return(-1);
    }

  bio_ud[unit].status = BIO_UNOPENED;

  bio_init_(&unit,name);

  return(unit);
  }

LOCAL_FUNCTION VOID
bio_init_(unit,name)
int *unit;
char *name;
{
  char param[MAX_STRING];
  char method[MAX_STRING]; 
  char unitch[MAX_STRING];

  bio_current_unit = *unit;

  sprintf(unitch,"%d",*unit);

  strcpy(param,"FILES:");
  strcat(param,unitch);
  strcat(param,":");

  if (get_file_info(name,"method","%s",method) != 0) {
    strcpy(method,"sequential");
    }

  if (!strcmp(method,"sequential")) {
    bio_ud[*unit].method = BIO_SEQUENTIAL;
    bio_ud[*unit].ptr.sequential = sequential_ioinit(param,*unit,name);
    }
  else {
    fprintf(bio_error,"bioinit_: invalid/unimplemented method (%s) for unit %d\n",method,*unit);
    bioabort();
    }
  }

GLOBAL_FUNCTION VOID
bio_done_file(unit)
int unit;
{
  if (bio_ud[unit].status != BIO_UNOPENED) {
    fprintf(bio_error,"bio_done_file: %d has not be closed\n");
    return;
    }

  if (bio_ud[unit].method == BIO_SEQUENTIAL) {
    sequential_free(bio_ud[unit].ptr.sequential);
    bio_ud[unit].ptr.sequential = NULL;
    bio_ud[unit].method = BIO_NOMETHOD;
    bio_ud[unit].status = BIO_UNINITED;
    }
  else {
    fprintf(bio_error,"bioinit_: invalid/unimplemented method (%d) for unit %d\n",bio_ud[unit].method,unit);
    bioabort();
    }
  }
