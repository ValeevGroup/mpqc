
/* $Log$
 * Revision 1.1  1993/12/29 12:53:40  etseidl
 * Initial revision
 *
 * Revision 1.4  1992/06/17  21:42:14  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.3  1992/04/06  12:27:56  seidl
 * merge in sandia changes
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tmpl.h>
#include "param.h"
#include "types.h"
#define ALLOC_GLOBALS
#include "pointers.h"
#undef ALLOC_GLOBALS

#include "file_info.gbl"
#include "r_async.gbl"
#include "ram.gbl"
#include "s_async.gbl"
#include "sequential.gbl"

ioFILE_t bio_ud[MAX_UNIT];

unsigned int bio_current_unit;

#include "bioopen.gbl"
#include "bioopen.lcl"

/* Specifies where error messages are written. */
GLOBAL_FUNCTION VOID
bio_errorFILE(fp)
FILE *fp;
{
  bio_error = fp;
  }


/* opens binary file "unit" with name "name" */

GLOBAL_FUNCTION int
bio_open_file(unit,name)
int unit;
char *name;

{

  if (bio_ptr.wptr == NULL) init_ptrs();

  if (!unit && !name) {
    fprintf(bio_error,"bio_open_file: one of unit and name must be != 0\n");
    return(-1);
    }

  if(unit==0 || bio_ud[unit].status == BIO_UNINITED) {
    unit = bio_init_file(unit,name);
    }

  bio_ptr.wptr[unit] = 0;
  bioopen_(&unit);

  return(unit);
  }

GLOBAL_FUNCTION VOID
bio_close_file(unit,status)
int unit;
int status;
{
  bioclos_(&unit,&status);
  }

GLOBAL_FUNCTION int
bio_write(unit,buffer,length)
int unit;
VOID_PTR buffer;
int length;
{
   int pos;
   int fword=bio_ptr.wptr[unit];

   biowrr_(&unit,buffer,&fword,&length);
   pos=fword+length;
   bio_ptr.wptr[unit]=pos;

   return(length);
   }

GLOBAL_FUNCTION int
bio_read(unit,buffer,length)
int unit;
VOID_PTR buffer;
int length;
{
   int pos;
   int fword=bio_ptr.wptr[unit];

   biordr_(&unit,buffer,&fword,&length);
   pos=fword+length;
   bio_ptr.wptr[unit]=pos;

   return(length);
   }

GLOBAL_FUNCTION VOID 
bio_rewind(unit)
int unit;

{
  bio_ptr.wptr[unit] = 0;
  }

GLOBAL_FUNCTION int
bio_get_pos(unit)
int unit;
{
  return(bio_ptr.wptr[unit]);
  }

GLOBAL_FUNCTION VOID
bio_set_pos(unit,pos)
int unit;
int pos;
{
  bio_ptr.wptr[unit]=pos;
  }

LOCAL_FUNCTION VOID
init_ptrs()
{
  bio_ptr.wptr = (int *) malloc(sizeof(int)*MAX_UNIT);

  if (bio_ptr.wptr == NULL) {
    fprintf(bio_error,"trouble allocating memory for pointers!\n");
    exit(1);
    }

  bioinit_();
  }

LOCAL_FUNCTION VOID
bioinit_()
{
  int i;

  for (i=0; i<MAX_UNIT; i++) {
    bio_ud[i].status = BIO_UNINITED;
    bio_ud[i].method = BIO_NOMETHOD;
    bio_ud[i].ptr.sequential = NULL;
    }
  }

LOCAL_FUNCTION VOID
bioopen_(unit)
int *unit;
{

  bio_current_unit = *unit;

  if (bio_ud[*unit].method == BIO_SEQUENTIAL) {
    sequential_ioopen(bio_ud[*unit].ptr.sequential,*unit);
    }
  else {
    fprintf(bio_error,"bioopen_: invalid/unimplemented method (%d) for unit %d\n",bio_ud[*unit].method,*unit);
    bioabort();
    }
  }

LOCAL_FUNCTION VOID
bioclos_(unit,status)
int *unit;
int *status;
{

  bio_current_unit = *unit;

  if (bio_ud[*unit].method == BIO_SEQUENTIAL) {
    sequential_ioclos(bio_ud[*unit].ptr.sequential, *status);
    }
  else {
    fprintf(bio_error,"bioclos_: invalid/unimplemented method (%d) for unit %d\n",
             bio_ud[*unit].method,*unit);
    bioabort();
    }
  bio_ud[*unit].status = BIO_UNOPENED;
  }

GLOBAL_FUNCTION VOID
biowrr_(unit,buffer,first,length)
int *unit;
VOID_PTR buffer;
int *first;
int *length;
{

  bio_current_unit = *unit;

  if (*first < 0) {
    fprintf(bio_error,"iowrr_: *unit=%d, *first=%d, *length=%d\n",
            *unit, *first, *length);
    bioabort();
    }

  if (bio_ud[*unit].method == BIO_SEQUENTIAL) {
    sequential_iowrr(bio_ud[*unit].ptr.sequential,buffer,*first,*length);
    }
  else if (bio_ud[*unit].method == BIO_R_ASYNC) {
    r_async_iowrr(bio_ud[*unit].ptr.r_async,buffer,*first,*length);
    }
  else if (bio_ud[*unit].method == BIO_S_ASYNC) {
    s_async_iowrr(bio_ud[*unit].ptr.s_async,buffer,*first,*length);
    }
  else if (bio_ud[*unit].method == BIO_RAM) {
    ram_iowrr(bio_ud[*unit].ptr.ram,buffer,*first,*length);
    }
  else {
    fprintf(bio_error,"iowrr_: invalid method (%d) for unit %d\n",
            bio_ud[*unit].method,*unit);
    bioabort();
    }
  }

GLOBAL_FUNCTION VOID
biordr_(unit,buffer,first,length)
int *unit;
VOID_PTR buffer;
int *first;
int *length;
{

  bio_current_unit = *unit;

  if (*first < 0) {
    fprintf(bio_error,"iordr_: *unit=%d, *first=%d, *length=%d\n",
            *unit, *first, *length);
    bioabort();
    }

  if (bio_ud[*unit].method == BIO_SEQUENTIAL) {
    sequential_iordr(bio_ud[*unit].ptr.sequential,buffer,*first,*length);
    }
  else if (bio_ud[*unit].method == BIO_R_ASYNC) {
    r_async_iordr(bio_ud[*unit].ptr.r_async,buffer,*first,*length);
    }
  else if (bio_ud[*unit].method == BIO_S_ASYNC) {
    s_async_iordr(bio_ud[*unit].ptr.s_async,buffer,*first,*length);
    }
  else if (bio_ud[*unit].method == BIO_RAM) {
    ram_iordr(bio_ud[*unit].ptr.ram,buffer,*first,*length);
    }
  else {
    fprintf(bio_error,"iordr_: invalid method (%d) for unit %d\n",
            bio_ud[*unit].method,*unit);
    bioabort();
    }
  }

GLOBAL_FUNCTION VOID
bioabort()
{
  fprintf(bio_error,"bioabort: current unit = %d\n",bio_current_unit);
  exit(1);
  }

