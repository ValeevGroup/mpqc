/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *  PICL source code                                               *
 *                                                                 *
 *  We welcome questions, comments, and bug reports, and request   *
 *  that you share any modifications with us.                      *
 *                                                                 *
 *  Patrick Worley                                                 *
 *  Oak Ridge National Laboratory                                  *
 *  worley@msr.epm.ornl.gov                                        *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "globals.h"
#include "envglobals.h"

/*
 * The following routines are designed to interface FORTRAN calls to C routines
 * found in the node portion of the portable library PICL.
 */

/*-------------------------------------------------------------------------*/
/*                Low Level Routines Section                               */
/*-------------------------------------------------------------------------*/

check0_(checking)
int *checking;
{
  check0(*checking);
}

double clock0_()
{
  double clock0();
  return(clock0());
}

clocksync0_()
{
  clocksync0();
}

close0_(release)
int *release;
{
  close0(*release);
}

int host0_()
{
  int host0();
  return(host0());
}

message0_(string, lth)
char *string;
int lth;
{
  int i;

  if (lth >= 80) lth = 79;
  for (i=0; i<lth; i++) envMBUF[i] = string[i];
  envMBUF[i] = '\0';
  message0(envMBUF);
}

open0_(maxprocs, me, host)
int *maxprocs, *me, *host ;
{
  open0(maxprocs, me, host);
}

int probe0_(type)
int *type;
{
  int probe0();
  return(probe0(*type));
}

recv0_(buf, bytes, type)
char *buf;
int *bytes, *type;
{
  recv0(buf, *bytes, *type);
}

recvbegin0_(buf, bytes, type)
char *buf;
int *bytes, *type;
{
  recvbegin0(buf, *bytes, *type);
}

recvend0_(type)
int *type;
{
  recvend0(*type);
}

recvinfo0_(bytes, type, node)
int *bytes, *type, *node;
{
  recvinfo0(bytes, type, node);
}

send0_(buf, bytes, type, node)
char *buf ;
int *bytes, *type, *node ;
{
  send0(buf, *bytes, *type, *node);
}

sendbegin0_(buf, bytes, type, node)
char *buf ;
int *bytes, *type, *node ;
{
  sendbegin0(buf, *bytes, *type, *node);
}

sendend0_(type)
int *type;
{
  sendend0(*type);
}

sync0_()
{
  sync0();
}

traceblock_(tracing)
int *tracing;
{
  traceblock(*tracing);
}

traceblockbegin_(type, location, parameter)
int *type, *location, *parameter;
{
  traceblockbegin(*type, *location, *parameter);
}

traceblockend_(type, location, parameter)
int *type, *location, *parameter;
{
  traceblockend(*type, *location, *parameter);
}

traceexit_()
{
  traceexit() ;
}

tracefiles_(tempfile, permfile, verbose, lth1, lth2)
char *tempfile, *permfile;
int lth1, lth2, *verbose;
{
  int i;
  char *tmpbuf1, *tmpbuf2;

  tmpbuf1 = malloc0(lth1+1);
  if (tmpbuf1 != NULL){
    for (i=0; i<lth1; i++) tmpbuf1[i] = tempfile[i];
    tmpbuf1[i] = '\0';
    };

  tmpbuf2 = malloc0(lth2+1);
  if (tmpbuf2 != NULL){
    for (i=0; i<lth2; i++) tmpbuf2[i] = permfile[i];
    tmpbuf2[i] = '\0';
    };

  tracefiles(tmpbuf1, tmpbuf2, *verbose);
  free(tmpbuf1);
  free(tmpbuf2);
}

traceflush_()
{
  traceflush() ;
}

traceinfo_(remaining, trace, compstats, commstats)
int *remaining, *trace, *compstats, *commstats ;
{
  traceinfo(remaining, trace, compstats, commstats);
}

tracelevel_(trace, compstats, commstats)
int *trace, *compstats, *commstats ;
{
  tracelevel(*trace, *compstats, *commstats);
}

tracemark_(type)
int *type ;
{
  tracemark(*type);
}

tracemarks_(markarray, size)
int *markarray, *size ;
{
  tracemarks(markarray, *size);
}

tracemsg_(string, lth)
char *string ;
int lth;
{
  int i;

  if (lth >= 80) lth = 79;
  for (i=0; i<lth; i++) envMBUF[i] = string[i];
  envMBUF[i] = '\0';
  tracemsg(envMBUF);
}

tracenode_(tracesize, flush, sync)
int *tracesize, *flush, *sync;
{
  tracenode(*tracesize, *flush, *sync);
}

who0_(numproc, me, host)
int *numproc, *me, *host;
{
  who0(numproc, me, host);
}

/*-------------------------------------------------------------------------*/
/*                     Portable Routines Section                           */
/*-------------------------------------------------------------------------*/

barrier0_()
{
  barrier0();
}

bcast0_(buf, bytes, type, root)
char *buf;
int *bytes, *type, *root;
{
  bcast0(buf, *bytes, *type, *root);
}

bcast1_(buf, bytes, type, root)
char *buf ;
int *bytes, *type, *root ;
{
  bcast1(buf, *bytes, *type, *root);
}

gand0_(buf, items, datatype, msgtype, root)
char *buf;
int *items, *datatype, *msgtype, *root;
{
  gand0(buf, *items, *datatype, *msgtype, *root);
}

gather0_(vec, n, type, root)
float *vec;
int *n, *type, *root;
{
  gather0(vec, *n, *type, *root);
}

gcomb0_(buf, items, datatype, msgtype, root, comb)
char *buf;
int *items, *datatype, *msgtype, *root, (*comb)();
{
  gcomb0(buf, *items, *datatype, *msgtype, *root, comb);
}

getarc0_(nprocs, top, ord, dir)
int *nprocs, *top, *ord, *dir;
{
  getarc0(nprocs, top, ord, dir);
}   

int ginv0_(i)
int *i;
{
  int ginv0();
  return(ginv0(*i));
}

gmax0_(buf, items, datatype, msgtype, root)
char *buf;
int *items, *datatype, *msgtype, *root;
{
  gmax0(buf, *items, *datatype, *msgtype, *root);
}

gmin0_(buf, items, datatype, msgtype, root)
char *buf;
int *items, *datatype, *msgtype, *root;
{
  gmin0(buf, *items, *datatype, *msgtype, *root);
}

gor0_(buf, items, datatype, msgtype, root)
char *buf;
int *items, *datatype, *msgtype, *root;
{
  gor0(buf, *items, *datatype, *msgtype, *root);
}

gprod0_(buf, items, datatype, msgtype, root)
char *buf;
int *items, *datatype, *msgtype, *root;
{
  gprod0(buf, *items, *datatype, *msgtype, *root);
}

int gray0_(i)
int *i;
{
  int gray0();
  return(gray0(*i));
}

gsum0_(buf, items, datatype, msgtype, root)
char *buf;
int *items, *datatype, *msgtype, *root;
{
  gsum0(buf, *items, *datatype, *msgtype, *root);
}

gxor0_(buf, items, datatype, msgtype, root)
char *buf;
int *items, *datatype, *msgtype, *root;
{
  gxor0(buf, *items, *datatype, *msgtype, *root);
}

int neighbor0_(i, node)
int *i, *node;
{
  int neighbor0();
  return(neighbor0(*i, *node));
}

setarc0_(nprocs, top, ord, dir)
int *nprocs, *top, *ord, *dir;
{
  setarc0(nprocs, top, ord, dir);
}

