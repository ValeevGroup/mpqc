
/* $Log$
 * Revision 1.1  1993/12/29 12:53:41  etseidl
 * Initial revision
 *
 * Revision 1.4  1992/06/17  21:42:22  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.3  1992/04/06  12:28:30  seidl
 * merge in sandia changes
 *
 * Revision 1.2  1992/03/30  22:45:39  jannsen
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.9  1992/01/14  19:29:47  cljanss
 * setup includes for the NCUBE
 *
 * Revision 1.8  1991/12/04  23:30:59  cljanss
 * changed the things are included
 *
 * Revision 1.7  1991/12/02  17:19:47  cljanss
 * converted to DEC
 *
 * Revision 1.6  91/11/18  17:35:21  cljanss
 * changed stderr to bio_error
 * 
 * Revision 1.5  91/09/30  17:04:42  cljanss
 * updated for NCUBE
 * 
 * Revision 1.4  91/09/28  21:05:37  cljanss
 * switch to new naming convention
 * 
 * Revision 1.3  91/06/21  03:15:09  seidl
 * change ioabort after lseek to bioabort
 * 
 * Revision 1.2  1991/06/17  16:05:59  seidl
 * first working version
 *
 * Revision 1.1  1991/06/16  20:28:41  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#include <stdio.h>
#include <string.h>
#ifdef __STDC__
#  include <stdlib.h>
#  include <fcntl.h>
#  if defined(SUN)
#    include <sys/file.h>
#  endif
#else
#  include <math.h>
#  ifdef NCUBE
#    include <fcntl.h>
#    include <sys/types.h>
#  else
#    include <sys/file.h>
#  endif
#endif
#include <tmpl.h>
#if defined(I860)
#  include <sys/types.h>
#endif
#include <sys/stat.h>
#include "param.h"
#include "types.h"
#include "pointers.h"
#include "file_info.gbl"
#include "errors.gbl"

#if SUN == 1
#include "unistd.h"
#endif

#ifndef FILENAME_MAX
#define FILENAME_MAX 256
#endif

#include "sequential.gbl"
#include "sequential.lcl"

#define DEBUG 0

extern void bioabort();

static char open_name[] = "sequential_ioopen";

GLOBAL_FUNCTION VOID
sequential_ioopen(ud,unit)
sequential_t *ud;
int unit;
{
  int i;
  char path[MAX_STRING];
  struct stat statjunk;
  char name_format[3];
  /* This file contains a list of files to be deleted
   * when we are using batch. */
  FILE *cleanfilep;
  /* If we are running in batch then this is set to the $(MBATCH)
   * environment variable to indicate the job_id. */
  int job_id;
  char *mbatchc;
  char cleanfile[FILENAME_MAX];

  /* Look at the environment to see if this is a batch job. */
  if (!(mbatchc = getenv("MBATCH="))) {
    job_id = 0;
    }
  else {
    job_id = atoi(mbatchc);
    sprintf(cleanfile,"Batch_Clean.%05d",job_id);
    /* If the file has not yet been open, then open it and
     * tell it that we want it to delete itself. */
    if (stat(cleanfile,&statjunk)!=0) {
      cleanfilep = fopen(cleanfile,"a");
      if (!cleanfilep) {
        fprintf(bio_error,"sequential I/O: couldn't open %s\n",cleanfile);
        }
      else {
        fseek(cleanfilep,(long)0,SEEK_END);
        fprintf(cleanfilep,"%s\n",cleanfile);
        fclose(cleanfilep);
        }
      }
    }

  strcpy(name_format,"%s");

  for (i=0; i<ud->n; i++) {

  /* Any file for with keep is given is kept, otherwise, any
   * file which has a /tmp as the first four characters in the
   * path is keep. */
    if (job_id&&(!strncmp("/tmp",path,4))) {
      cleanfilep = fopen(cleanfile,"a");
      if (!cleanfilep) {
        fprintf(bio_error,"sequential I/O: couldn't open %s\n",cleanfile);
        }
      else {
        fseek(cleanfilep,(long)0,SEEK_END);
        fprintf(cleanfilep,"%s\n",path);
        fclose(cleanfilep);
        }
      }

    ud->v[i].stream = open(ud->v[i].path,O_RDWR|O_CREAT,0644);
    }
  ud->next = -1;
  ud->unit = unit;

  if (ud->verbose) {
    fprintf(bio_error,"SEQ_IO: opened unit %d {\n",unit);
    fprintf(bio_error,"  blocksize = %d\n",ud->blocksize);
    for (i=0; i<ud->n; i++) {
      fprintf(bio_error,"  v[%d].path = \"%s\"\n",i,ud->v[i].path);
      }
    fprintf(bio_error,"  }\n");
    }
  }

GLOBAL_FUNCTION sequential_t *
sequential_ioinit(baseparam,unit,filename)
char *baseparam;
int unit;
char *filename;
{
  int i;
  sequential_t *ud;
  char param[MAX_STRING];
  char name[MAX_STRING];
  char volpath[MAX_STRING];
  char path[MAX_STRING];
  char name_format[3];

  strcpy(name_format,"%s");

  /* Allocate memory for the unit descriptor. */
  ud = (sequential_t *) malloc(sizeof(sequential_t));
  bio_malloc_check(open_name,ud);

 /* Find out if we want extra information to be printed about this
  * unit.  name is used as a temporary buffer here. */

  if (get_file_info(filename,"verbose","%s",name)== -1) ud->verbose = 0;
  else {
    if (!strcmp(name,"on")) ud->verbose = 1;
    else if (!strcmp(name,"yes")) ud->verbose = 1;
    else if (!strcmp(name,"y")) ud->verbose = 1;
    else if (!strcmp(name,"true")) ud->verbose = 1;
    else if (!strcmp(name,"YES")) ud->verbose = 1;
    else if (!strcmp(name,"TRUE")) ud->verbose = 1;
    else if (!strcmp(name,"ON")) ud->verbose = 1;
    else if (!strcmp(name,"Y")) ud->verbose = 1;
    else if (!strcmp(name,"1")) ud->verbose = 1;
    else ud->verbose = 0;
    }

 /* Find out how many volumes to place the unit across. */
  if (get_file_info(filename,"nvolume","%d",&ud->n) == -1) ud->n = 1;
  if (ud->n == 0) ud->n = 1;

 /* Set up the block size for this file system. */
    strcat(param,"BLOCKSIZE");
  if (get_file_info(filename,"blocksize","%d",&ud->blocksize) == -1)
      ud->blocksize = 8192;
  if (ud->blocksize == 0) ud->blocksize = 8192;

  /* Find out how the files are to be named. */
  if (get_file_info(filename,"name",name_format,name) == -1)
      strcpy(name,"sequential");

  for (i=0; i<ud->n; i++) {
    sprintf(param,"%s%d","volume",i+1);
    if (get_file_info(filename,param,"%s",volpath) == -1) {
      if (ud->n > 1) bio_no_path_given(open_name);
      else volpath[0] = '\0';
      }

    sprintf(path,"%s%s%c%s",volpath,name,'.',filename);
    ud->v[i].path = (char *) malloc(strlen(path)+1);
    bio_malloc_check(open_name,ud->v[i].path);
    strcpy(ud->v[i].path,path);

    }

  if (ud->verbose) {
    fprintf(bio_error,"SEQ_IO: init'ed unit %d {\n",unit);
    fprintf(bio_error,"  blocksize = %d\n",ud->blocksize);
    for (i=0; i<ud->n; i++) {
      fprintf(bio_error,"  v[%d].path = \"%s\"\n",i,ud->v[i].path);
      }
    fprintf(bio_error,"  }\n");
    }
  return(ud);
  }

GLOBAL_FUNCTION VOID
sequential_ioclos(ud,status)
sequential_t *ud;
int status;
{
  int i;

  for (i=0; i<ud->n; i++) {
    close(ud->v[i].stream);
    if (status == BIO_UNLINK) unlink(ud->v[i].path);
    }
  if (ud->verbose) {
    fprintf(bio_error,"SEQ_IO: closed unit %d {\n",ud->unit);
    fprintf(bio_error,"  incount = %u\n",ud->incount);
    fprintf(bio_error,"  outcount = %u\n",ud->outcount);
    fprintf(bio_error,"  }\n");
    }
  }

GLOBAL_FUNCTION VOID
sequential_free(ud)
sequential_t *ud;
{
  int i;

  for (i=0; i<ud->n; i++) {
    free(ud->v[i].path);
    ud->v[i].path = NULL;
    }
  if (ud->verbose) {
    fprintf(bio_error,"SEQ_IO: free'ed unit %d {\n",ud->unit);
    fprintf(bio_error,"  incount = %u\n",ud->incount);
    fprintf(bio_error,"  outcount = %u\n",ud->outcount);
    fprintf(bio_error,"  }\n");
    }
  free(ud);
  }

GLOBAL_FUNCTION VOID
sequential_iordr(ud,buffer,first,length)
sequential_t *ud;
char *buffer;
int first;
int length;
{
  ud->incount += length;
  sequential_iordwrr("sequential_iordr",IOOP_READ,ud,buffer,first,length);
  }

GLOBAL_FUNCTION VOID
sequential_iowrr(ud,buffer,first,length)
sequential_t *ud;
char *buffer;
int first;
int length;
{
  ud->outcount += length;
  sequential_iordwrr("sequential_iowrr",IOOP_WRITE,ud,buffer,first,length);
  }

LOCAL_FUNCTION VOID
sequential_iordwrr(caller,ioop,ud,buffer,first,length)
char *caller;
int ioop;
sequential_t *ud;
char *buffer;
int first;
int length;
{
  int i;
  long firstbyte, lastbyte;
  long ncycles;
  long remainingbytes;
  long remainingseek;
  long fullvol;
  long offset1, offset2;
  long ibuf;
  int len;

  firstbyte = first;
  lastbyte = first + length - 1;

  ncycles = firstbyte/(ud->n*ud->blocksize);
  remainingseek = firstbyte - ncycles*ud->n*ud->blocksize;
  fullvol = remainingseek/ud->blocksize;

  offset2 = ncycles * ud->blocksize;
  offset1 = offset2 + ud->blocksize;

  /* Seek all volumes to the appropiate positions. */
  if ((ud->next != firstbyte)||(ud->last_ioop != ioop)) {
    for (i=0; i<ud->n; i++) {
      long offset;
      if (i < fullvol) offset = offset1;
      else if (i == fullvol) offset = offset2 + remainingseek%ud->blocksize;
      else offset = offset2;
#if DEBUG
      fprintf(stdout,"seeking volume %d to %ld\n",i,offset);
#endif
      if (lseek(ud->v[i].stream, offset, SEEK_SET)<0) {
        fprintf(bio_error,"%s: fseek: offset = %ld, vol = %d\n",caller,offset,i);
        perror(caller);
        bioabort();
        }
      }
    if (length > 0) ud->last_ioop = ioop;
    }
  ud->next = lastbyte + 1;

  /* Do the io. */
  i = fullvol;
  len = ud->blocksize - remainingseek%ud->blocksize;
  remainingbytes = lastbyte - firstbyte + 1;
#if DEBUG
  fprintf(stdout,"%s: len=%d,remainingbytes=%ld,firstbyte=%ld,lastbyte=%ld,i=%d\n",
          caller,len,remainingbytes,firstbyte,lastbyte,i);
#endif
  ibuf = 0;
  while (remainingbytes > 0) {
    if (len > remainingbytes) len = remainingbytes;
#if DEBUG
    fprintf(stdout,"       len=%d,remainingbytes=%ld,i=%d\n",
          len,remainingbytes,i);
#endif
    if (ioop == IOOP_READ) {
      if (read(ud->v[i].stream,&buffer[ibuf],len)<1) {
        fprintf(bio_error,"%s: len = %d, volume = %d\n",caller,len,i);
        bio_read_error(caller);
        }
      }
    else if (ioop == IOOP_WRITE) {
      if (write(ud->v[i].stream,&buffer[ibuf],len)!=len) {
        fprintf(bio_error,"%s: len = %d, volume = %d\n",caller,len,i);
        bio_write_error(caller);
        }
      }
    else {
      fprintf(bio_error,"%s: illegal ioop = %d\n",caller,ioop);
      bioabort();
      }
    i++;
    if (i == ud->n) i=0;
    remainingbytes -= len;
    ibuf += len;
    len = ud->blocksize;
    }
  }
