
/* Simulate PICL in a single unix process */
/* This provides no host.  Has only node routines. */
/* For now only one node is supported.  In the future perhaps more
 * can be supported with some nasty unix programming. */

/* Here's the RCS junk */
/* $Log$
 * Revision 1.2  1994/08/27 00:02:57  etseidl
 * fix some warnings
 *
 * Revision 1.1.1.1  1993/12/29  12:53:49  etseidl
 * SC source tree 0.1
 *
 * Revision 1.6  1992/06/17  22:37:37  jannsen
 * clean up for saber-c
 *
 * Revision 1.5  1992/05/13  18:16:26  jannsen
 * added sync0, which does nothing
 *
 * Revision 1.4  1992/04/16  18:23:38  jannsen
 * eliminated some debugging output
 *
 * Revision 1.3  1992/04/13  11:09:43  seidl
 * *** empty log message ***
 *
 * Revision 1.2  1992/04/06  12:15:50  seidl
 * add definition of CLOCKS_PER_SEC
 * add sandia changes
 *
 * Revision 1.1.1.1  1992/03/31  02:15:54  jannsen
 * The CCE picl routines
 *
 * Revision 1.1  1992/03/31  02:15:53  jannsen
 * Initial revision
 *
 * Revision 1.2  1991/03/02  01:37:30  colvin
 * Removed unused TEST_CASE cpp directives.
 *
 * Revision 1.1  91/02/27  22:08:59  colvin
 * Initial revision
 *  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <util/misc/libmisc.h>
#include <comm/picl/picl.h>

#ifndef OLDCLOCK
#  include <sys/types.h>
#  include <sys/time.h>
#  include <sys/resource.h>
#  if defined(sun) || defined(AIX)
#  include <unistd.h>
#  endif
#endif

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC ((int) 1000000)
#endif

/* First let's define some static variables */
static int last_count;   /* holds the size of the last msg sent or recvd */
static int last_node;    /* holds the id of where last msg came from */
static int last_type;    /* holds the type of where last msg recvd */
static int check_option=1;

/* These are set by the setarc0 and getarc0 call.
 * Since this version of PICL only supports node programs they aren't
 * really used. */
static int arc_top=1;
static int arc_ord=0;
static int arc_dir=1;

struct message_struct {
  void *buf;
  int size;
  int type;
  int source;
  struct message_struct *p;
  };
typedef struct message_struct message_t;

static message_t *messages=NULL;

/* Here are the PICL low-level primitives, in the order they appear in 
tabel 1 of the PICL documentation */

/*******************************************************/

int host0() {return 0;}
int host0_() {return 0;}

void sync0() {}
void sync0_() {}

#ifndef OLDCLOCK
double clock0()
{
  double res;
  struct rusage r;

  getrusage(RUSAGE_SELF,&r);

  res = r.ru_utime.tv_sec + r.ru_stime.tv_sec;
  res += 0.000001 * ( r.ru_utime.tv_usec + r.ru_stime.tv_usec );
  return res;
}
#else
double clock0()
{
  return(((double)clock())/CLOCKS_PER_SEC);
}
#endif

void check0(checking)
int checking;
{
  check_option = checking;
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
  *numproc = 1;
  *me = 0;
  *host = 0xffff;
}

void setarc0(nprocs,top,ord,dir)
int *nprocs;
int *top;
int *ord;
int *dir;
{
  *nprocs = 1;
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
  *nprocs = 1;
  *top = arc_top;
  *ord = arc_ord;
  *dir = arc_dir;
  }

int probe0(msgtype)
int msgtype;
{
  message_t *i;

  for (i=messages; i!=NULL; i = i->p) {
    if (i->type == msgtype || msgtype == -1) {
      last_node=i->source;
      last_type=i->type;
      return 1;
      }
    }
  return 0;
}

/* Print out the queued messages. */
void piclproc_queued()
{
  message_t *I;

  for (I=messages; I!=NULL; I=I->p) {
    printf(" (%d,%d,%d)",I->size,I->type,I->source);
    }
  }

void recv0(buf, bytes, type)
void *buf;
int type, bytes;
{
  message_t *i;
  message_t *last;

  last = NULL;
  for (i=messages; i!=NULL; i = i->p) {
    if (i->type == type || type == -1) {
      if (i->size > bytes) {
        fprintf(stderr,"recv0: message buffer isn't big enough\n");
        exit(1);
        }
      memcpy(buf,i->buf,i->size);
      last_node=i->source;
      last_type=i->type;
      last_count=i->size;

#if 0
      if (type>4000) {
        printf("recv0: %d bytes of type %d, ",i->size,i->type);
        printf("queued:"); piclproc_queued(); printf("\n");
        }
#endif

      /* Remove the message from the list. */
      if (last) {
        last->p = i->p;
        }
      else {
        messages = messages->p;
        }
      free(i->buf);
      free(i);

#if 0
      if (messages) {
        printf("recv0: remaining:"); piclproc_queued(); printf("\n");
        }
#endif

      return;
      }
    last = i;
    }

  fprintf(stderr,"recv0: tried to receive something that isn't there\n");
  fprintf(stderr,"recv0: tried %d bytes of type %d, ",bytes,type);
  printf("queued:"); piclproc_queued(); printf("\n");
#if !defined(SUN)
  debug_start("error in recv0");
#endif
  exit(1);
  }


void recvinfo0(bytes, msgtype, source)
int *bytes, *msgtype, *source;
{
  *bytes=last_count;
  *source=last_node;
  *msgtype=last_type;
}
 

void send0(buf, bytes, msgtype, dest)
void *buf;
int bytes, msgtype, dest;
{
  message_t *msg;
  message_t *I;

#if 0
  if (msgtype>4000) {
    printf("send0: %d bytes of type %d to %d, ",bytes,msgtype,dest);
    printf("queued:"); piclproc_queued(); printf("\n");
    }
#endif

  if (dest != 0) {
    fprintf(stderr,"send0: can only send to 0\n");
    exit(1);
    }

  msg = (message_t *) malloc(sizeof(message_t));
  if (msg) msg->buf = (char *) malloc(bytes);
  if (!msg || !msg->buf) {
    fprintf(stderr,"send0: allocation failed\n");
    exit(1);
    }

  /* Put msg at the end of the linked list, because of some bad
   * assumptions made by the mpscf program and libraries. */
  msg->p = NULL;
  if (!messages) {
    messages = msg;
    }
  else {
    for (I=messages; I->p != NULL; I=I->p);
    I->p = msg;
    }

  memcpy(msg->buf,buf,bytes);
  msg->type = msgtype;
  msg->size = bytes;
  msg->source = 0;
}

void sync()
{
  /* sync is automatic for 1 node. */
}

void who0(numproc, me, host)
int *numproc, *me, *host;
{
  *me = 0;
  *numproc = 1;
}

/* FORTRAN stub routines */

void message0_(message) 
char *message; 
{printf ("%s\n", message);}

void check0_(checking)
int *checking;
{
  check0(*checking);
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

void sync_()
{
  sync();
}

void who0_(numproc,me,host)
int *numproc,*me,*host;
{
  who0(numproc,me,host);
}

