/*      PICL to Paragon interface library       */

/* $Log$
 * Revision 1.1  1993/12/29 12:53:47  etseidl
 * Initial revision
 *
 *  */
static char rcsid[]="$Id$";

#include <stdio.h>
#include <comm/picl/picl.h>
#include <nx.h>

/* First let's define some external static variables */
static int checking_flag=1;    /* specifies if parameter checking is on */

/* These are set by the setarc0 and getarc0 call.
 * These are ignored in this version of PICL
 */
static int arc_nproc;
static int arc_top=-1;
static int arc_ord=0;
static int arc_dir=1;

int host0() {return 0;}

void check(checking)
int checking;
{
  checking_flag=checking;
}

double clock0()
{
  static double first;
  static initp = 0;
  double t;
  if (!initp) {
    first = dclock();
    initp = 1;
  }
  t = dclock();
  return (t - first);
}


void close0(release)
int release;
{}

void message0(message)
char *message;
{printf ("%s\n", message);}

void open0(numproc, me, host)
int *numproc, *me, *host;
{
  *numproc = numnodes();
  *me = mynode();
  *host = -1;
}

void setarc0(nprocs,top,ord,dir)
int *nprocs;
int *top;
int *ord;
int *dir;
{
  *nprocs = arc_nproc = numnodes();
  *top = arc_top;
  *ord = arc_ord;
  *dir = arc_dir;
  }

void getarc0(nprocs,top,ord,dir)
int *nprocs;
int *top;
int *ord;
int *dir;
{
  *nprocs = arc_nproc;
  *top = arc_top;
  *ord = arc_ord;
  *dir = arc_dir;
  }

int probe0(msgtype)
int msgtype;
{
  return iprobe(msgtype);
}

void recv0(buf, bytes, type)
int type, bytes;
void *buf;
{
  crecv(type,buf,bytes);
}

picl_check()
{
#if 0
  int source, type, i, j, length, totallength;
  int numproc,me,host;
  who0(&numproc,&me,&host);
  totallength = 0;
  for (i=0; i<32500; i+=500) {
    nrange(i,i+500 -1);
    source = -1;
    type = -1;
    if (ntest(&source,&type) >= 0) {
      for (j=0; j<500; j++) {
        source = -1;
        type = i + j;
        if ((length = ntest(&source,&type)) > 0) {
          if (!totallength) printf("picl_check:\n");
          totallength += length;
          printf(" on %d found type %d from %d of length %d\n",
                 me,type,source,length);
          }
        }
      }
    }
  nrange(0,32767);
  if (totallength > 0) {
    printf("totallength was %d\n",totallength);
    tim_traceback();
    }
#endif
}

void recvinfo0(bytes, msgtype, source)
int *bytes, *msgtype, *source;
{
  *bytes=infocount();
  *source=infonode();
  *msgtype=infotype();
}
 
void send0(buf, bytes, msgtype, dest)
int bytes, msgtype, dest;
void *buf;
{
  csend(msgtype,buf,bytes,dest,0);
}

void sync0()
{
  gsync();
}

void who0(numproc, me, host)
int *numproc, *me, *host;
{
  *numproc = numnodes();
  *me = mynode();
  *host = -1;
}

/*********************************************************************
 * Fortran stubs                                                     *
 *********************************************************************/

void message0_(message) 
char *message; 
{printf ("%s\n", message);}

int host0_() {return host0();}

void check_(checking)
int *checking;
{
  check(*checking);
}

double clock0_()
{
  return(clock0());
}

void close0_(release)
int *release;
{
  close0(*release);
}

void open0_(numproc,me,host)
int *numproc,*me,*host;
{
  open0(numproc,me,host);
}

int probe0_(msgtype)
int *msgtype;
{
  return(probe0(*msgtype));
}

void recv0_(buf,bytes,type)
int *type, *bytes;
char *buf;
{
  recv0(buf,*bytes,*type);
}

void recvinfo0_(bytes,msgtype,source)
int *bytes,*msgtype,*source;
{
  recvinfo0(bytes,msgtype,source);
}

void send0_(buf,bytes,msgtype,dest)
int *bytes,*msgtype,*dest;
char *buf;
{
  send0(buf,*bytes,*msgtype,*dest);
}

void sync0_()
{
  sync0();
}

void who0_(numproc,me,host)
int *numproc,*me,*host;
{
  who0(numproc,me,host);
}

