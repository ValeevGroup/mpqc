
/* This provides an implementation of gopf2 (which is used by piclgopf2)
 * and bcast0.  It is designed for a hypercube like connectivity, but
 * works if the number of nodes is not a power of two.  These routines
 * use low level picl routines.
 */

/* $Log$
 * Revision 1.1  1993/12/29 12:53:49  etseidl
 * Initial revision
 *
 * Revision 1.6  1993/04/27  21:11:17  jannsen
 * If COMMBUFSIZE is defined, then try to avoid sending messages bigger
 * than this.
 *
 * Revision 1.5  1992/06/17  22:33:10  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.4  1992/05/26  20:07:30  jannsen
 * global ops only use the passed message types
 *
 * Revision 1.3  1992/05/04  23:04:44  jannsen
 * Changed the names of sync and barrier to the proper picl names.
 *
 * Revision 1.2  1992/04/16  18:15:14  jannsen
 * Modified global operations so message type is not ignored.
 *
 * Revision 1.1.1.1  1992/03/31  02:15:34  jannsen
 * The CCE picl routines
 *
 * Revision 1.1  1992/03/31  02:15:33  jannsen
 * Initial revision
 *
 * Revision 1.1.1.1  1992/03/31  02:11:45  jannsen
 * The CCE picl routines
 *
 * Revision 1.1  1992/03/31  02:11:44  jannsen
 * Initial revision
 *
 * Revision 1.2  91/10/07  20:09:20  cljanss
 * A variable name was the same as a passed parameter name in gcomb0,
 * which caused problems on the NCUBE.
 * 
 * Revision 1.1  91/08/21  00:12:20  cljanss
 * Initial revision
 * 
 * Revision 1.5  1991/07/17  20:44:40  cljanss
 * Checked in Mike Colvin's latest version after I made the copy.
 *
 * Revision 1.4  1991/04/10  00:15:07  colvin
 * Added Fortran stub for gsum0.
 *
 * Revision 1.3  1991/03/12  00:10:37  colvin
 * Added a version of gcomb0.
 *
 * Revision 1.2  1991/03/11  22:37:55  colvin
 * Added distributed operations which have not yet been carefully tested
 *
 * Revision 1.1  1991/03/04  23:47:45  colvin
 * Initial revision
 * */
static char rcsid[]="$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <comm/picl/picl.h>
#include "assert.h"

#ifdef COMMBUFSIZE
#undef COMMBUFSIZE
#define COMMBUFSIZE 20
#endif

#ifdef COMMBUFSIZE
static void send0_block(buf,nbytes,type,dest,ready_type)
void*buf;
int nbytes;
int type;
int dest;
int ready_type;
{
  int k,junk;

  if (nbytes<=COMMBUFSIZE) {
      send0 (buf, nbytes, type, dest);
    }
  else {
      int nbytes_actual = COMMBUFSIZE;
      for (k=0; k<nbytes; k+=COMMBUFSIZE) {
          if (k+COMMBUFSIZE>nbytes) nbytes_actual = nbytes - k;
          recv0 (&junk,4,ready_type);
          send0 (&((char *)buf)[k], nbytes_actual, type, dest);
        }
    }
}
#endif

#ifdef COMMBUFSIZE
static void recv0_block(buf,nbytes,type,orig,ready_type)
void*buf;
int nbytes;
int type;
int orig;
int ready_type;
{
  int k,junk;

  if (nbytes<=COMMBUFSIZE) {
      recv0 (buf, nbytes, type);
    }
  else {
      int nbytes_actual = COMMBUFSIZE;
      for (k=0; k<nbytes; k+=COMMBUFSIZE) {
          if (k+COMMBUFSIZE>nbytes) nbytes_actual = nbytes - k;
          send0 (&junk,4,ready_type,orig);
          recv0 (&((char *)buf)[k], nbytes_actual, type);
        }
    }
}
#endif

void gopf2(data,ndata,tmp,type,root,comb)
     char *data;
     int ndata;
     char *tmp;
     int root,type;
     long (*comb)();
{
  int n, me, nproc, host;
  int bit,lowbits;
#ifdef COMMBUFSIZE
  int ready_type = mtype_get();
#endif

  /* first figure out who I am and how many processors there are */
  who0(&nproc,&me,&host);

  /* Now calculate the size of the binary hypercube to imbed nproc into */
  n=1;
  while(n < nproc) n = n<<1;
  
  /* Now loop through the log2(n) steps required */
  for (bit=1,lowbits=0; bit < n; bit = bit<<1,lowbits=((lowbits<<1)|1)) {
      /* if i have this bit then i'm a sender */
      if (me&bit) {
          int target = me-bit;
#ifndef COMMBUFSIZE
	  send0(data,ndata,type,target);
#else
	  send0_block(data,ndata,type,target,ready_type);
#endif
          break;
	}
      /* bit is off so i'm a receiver */
      else {
          /* only recv if there was a sender */
	  if (me+bit < nproc) {
#ifndef COMMBUFSIZE
	      recv0(tmp,ndata,type);
#else
	      recv0_block(tmp,ndata,type,me+bit,ready_type);
#endif
	      (*comb)(data,tmp);
	    }
	}
    }

  if (root != 0) {
      if (root == me) {
#ifndef COMMBUFSIZE
          recv0(data,ndata,type);
#else
          recv0_block(data,ndata,type,0,ready_type);
#endif
        }
      else if (me == 0) {
#ifndef COMMBUFSIZE
          send0(data,ndata,type,root);
#else
          send0_block(data,ndata,type,root,ready_type);
#endif
        }
    }
}

void bcast0 (data, ndata, type, root)
void *data;
int ndata, type, root;
{
  int n, me, nproc, host;
  int bit,highbits;
#ifdef COMMBUFSIZE
  int ready_type = mtype_get();
#endif

  /* first figure out who I am and how many processors there are */
  who0(&nproc,&me,&host);

  /* Now calculate the size of the binary hypercube to imbed nproc into */
  n=1;
  while(n < nproc) n = n<<1;

  /* make sure that data to be sent is on node 0 */
  if (root != 0) {
      if (me == root) {
#ifndef COMMBUFSIZE
          send0(data,ndata,type,0);
#else
          send0_block(data,ndata,type,0,ready_type);
#endif
        }
      if (me == 0) {
#ifndef COMMBUFSIZE
          recv0(data,ndata,type);
#else
          recv0_block(data,ndata,type,0,ready_type);
#endif
        }
    }
  
  /* Now loop through the log2(n) steps required */
  for (bit=1,highbits=n-2; bit < n; (bit=bit<<1),(highbits=highbits<<1)) {
      /* if i don't have this bit and none above then i'm a sender */
      if (!(me&bit) && !(me&highbits)) {
          int target = me + bit;
          /* only send data that will be received */
          if (target < nproc) {
#ifndef COMMBUFSIZE
              send0(data,ndata,type,target);
#else
              send0_block(data,ndata,type,target,ready_type);
#endif
	    }
	}
      /* bit is on so i'm a receiver if bits above bit are all zero */
      else if (!(me&highbits)) {
#ifndef COMMBUFSIZE
          recv0(data,ndata,type);
#else
          recv0_block(data,ndata,type,me-bit,ready_type);
#endif
        }
    }
}


