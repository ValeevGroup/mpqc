
/* Simulate PICL using unix messages */
/* This provides no host.  It has only node routines. */

/* Here's the RCS junk */
/* $Log$
 * Revision 1.1  1993/12/29 12:53:47  etseidl
 * Initial revision
 *
 * Revision 1.4  1992/04/16  19:46:53  jannsen
 * debugging output removed
 *
 * Revision 1.3  1992/04/13  11:09:10  seidl
 * *** empty log message ***
 *
 * Revision 1.2  1992/04/06  12:13:48  seidl
 * define clocks_per_sec for sun
 *
 * Revision 1.1.1.1  1992/03/31  02:15:40  jannsen
 * The CCE picl routines
 *
 * Revision 1.1  1992/03/31  02:15:39  jannsen
 * Initial revision
 *
 * Revision 1.1.1.1  1992/03/31  02:11:54  jannsen
 * The CCE picl routines
 *
 * Revision 1.1  1992/03/31  02:11:53  jannsen
 * Initial revision
 *
 * Revision 1.2  1991/03/02  01:37:30  colvin
 * Removed unused TEST_CASE cpp directives.
 *
 * Revision 1.1  91/02/27  22:08:59  colvin
 * Initial revision
 *  */
static char rcsid[]="$Id$";

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/msg.h>
#include <unistd.h>
/* #define _BSD_SIGNALS */
#include <signal.h>

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC ((int) 1000000)
#endif

#define MAX_PROC 8
#define NODE_BITS 3 /* The number of bits needed to encode nodeid's. */
#define PICL_TO_MSG(t) (t + 1)
#define MSG_TO_PICL(t) (t - 1)

/* First let's define some static variables */
static int last_count;   /* holds the size of the last msg sent or recvd */
static int last_node;    /* hold the last node we got a msg from */
static int last_type;    /* holds the type of where last msg recvd */
static int check_option=1;
static int mynodeid;
static int n_process;

static int data_q[MAX_PROC];
static int stat_q[MAX_PROC];

struct piclstat_struct {
  int source;
  int size;
  int type;
  };
typedef struct piclstat_struct piclstat_t;

struct statmsgbuf_struct {
  long mtype;
  piclstat_t piclstat;
  };
typedef struct statmsgbuf_struct statmsgbuf_t;

#define TXTSIZ 4096
struct msgbuf_struct {
  long mtype;
  char mtext[TXTSIZ];
  };
typedef struct msgbuf_struct msgbuf_t;

/* These are set by the setarc0 and getarc0 call.
 * Since this version of PICL only supports node programs they aren't
 * really used. */
static int arc_nproc;
static int arc_top=1;
static int arc_ord=0;
static int arc_dir=1;

/* Here are the PICL low-level primitives, in the order they appear in 
tabel 1 of the PICL documentation */

/*******************************************************/

int host0() {return 0;}
int host0_() {return 0;}

int mynode0() {return mynodeid;}
int mynode0_() {return mynodeid;}

/* Handles interrupts by cleaning up and exiting. */
static void handler(sig,code,scp,addr)
int sig, code;
struct sigcontext *scp;
char *addr;
{
  (void) signal (SIGINT, SIG_IGN);
  printf("cleaning up because of interrupt\n");
  piclcleanup();
  }

double clock0()
{
  return(((double)clock())/CLOCKS_PER_SEC);
}

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
  int i;
  int stat;
  char *numprocch;
  int cleaner = 1;

  *host = 0xffffffff;
  numprocch = getenv("NUMPROC");
  if (numprocch) *numproc = atoi(numprocch);
  else *numproc = 1;
  *me = 0;
  if (! *numproc) *numproc = 1;

  if (*numproc > MAX_PROC) {
    fprintf(stderr,"NUMPROC=%s is too big, %d is max\n",numprocch,MAX_PROC);
    exit(1);
    }

  setpgid(getpid(),getpid());

  /* Set up the message queues. */
  for (i=0; i<*numproc; i++) {
    if ((data_q[i] = msgget(IPC_PRIVATE,0700|IPC_CREAT))==-1) {
      perror("msgget data");
      exit(1);
      }
    if ((stat_q[i] = msgget(IPC_PRIVATE,0700|IPC_CREAT))==-1) {
      perror("msgget stat");
      exit(1);
      }
    }

  arc_nproc = *numproc;
  n_process = *numproc; /* Needed by the cleanup routines. */

  /* Begin forking off the other processes. */
  for (i=0; i<*numproc; i++) {
    int pid;
    if (!(pid = fork())) {
       *me = i;
       cleaner = 0;
       break;
       }
     }
  mynodeid = *me;
#if 0
  printf("mynodeid = %d, cleaner = %d\n",mynodeid,cleaner);
#endif

  /* If this is the original parent process, then clean up when
   * a kid has died. */
  if (cleaner) {
    int pid;
    (void) signal (SIGINT, handler);
    pid = wait(&stat);
    printf("cleaning up because pid %d has died\n",pid);
sleep(1000);
    piclcleanup();
    }
}

piclcleanup()
{
  int i;
  for (i=0; i<n_process; i++) {
    msgctl(stat_q[i],IPC_RMID);
    msgctl(data_q[i],IPC_RMID);
    }
  killpg(getpid(),SIGKILL);
  }

void setarc0(nprocs,top,ord,dir)
int *nprocs;
int *top;
int *ord;
int *dir;
{
  arc_nproc = *nprocs;
  arc_top = *top;
  arc_ord = *ord;
  arc_dir = *dir;
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
  return 0;
}

void recv0(buf, bytes, type)
int type, bytes;
char *buf;
{
  int datatype;
  int i;
  int retcode;
  statmsgbuf_t statmsg;
  msgbuf_t msg;
  int received_type;
  int bytes_remaining;
  int offset;
  int number_to_recv;

#if 0
  printf("recv0: %d bytes of type %d on %d\n",
         bytes,type,mynodeid);
#endif

  statmsg.mtype = PICL_TO_MSG(type);
  if ((retcode = msgrcv(stat_q[mynodeid],(void *)&statmsg,sizeof(piclstat_t),PICL_TO_MSG(type),0)) == -1) {
    perror("msgrcv stat");
    exit(1);
    }
  last_type = MSG_TO_PICL(statmsg.mtype);

#if 0
  printf("got a piclstat on %d:\n type=%d\n source=%d\n size=%d\n retcode=%d\n",
         mynodeid,statmsg.piclstat.type,statmsg.piclstat.source,statmsg.piclstat.size,retcode);
#endif

  if (statmsg.piclstat.type <= 0) {
    printf("got garbage for statmsg.piclstat.type -- aborting\n");
    exit(1);
    }

  msg.mtype = statmsg.piclstat.type;

#if 0
  if (statmsg.piclstat.size > bytes) {
    fprintf(stderr,"recv0: buffer is too small\n");
    exit(1);
    }
#endif

  bytes_remaining = statmsg.piclstat.size;
  offset = 0;
  datatype = statmsg.piclstat.type >> NODE_BITS;
  while(bytes_remaining > 0) {
    number_to_recv = (bytes_remaining > TXTSIZ) ? TXTSIZ: bytes_remaining;

    msg.mtype = ((datatype++)<<NODE_BITS)|statmsg.piclstat.source;
    if (msgrcv(data_q[mynodeid],(void*)&msg,number_to_recv,msg.mtype,0)==-1) {
      perror("msgrcv data");
      exit(1);
      }

    for (i=0; i<number_to_recv; i++) {
      buf[i+offset] = msg.mtext[i];
      }
    bytes_remaining -= number_to_recv;
    offset += number_to_recv;
    }

  last_node = statmsg.piclstat.source;
  last_count= statmsg.piclstat.size;
  }


void recvinfo0(bytes, msgtype, source)
int *bytes, *msgtype, *source;
{
  *bytes=last_count;
  *source=last_node;
  *msgtype=last_type;
}
 

void send0(buf, bytes, msgtype, dest)
int bytes, msgtype, dest;
char *buf;
{
  int i;
  int retcode;
  statmsgbuf_t statmsg;
  msgbuf_t msg;
  static datatype = 100;
  int bytes_remaining;
  int offset;
  int number_to_send;

  if (datatype > 9000000) datatype = 1;

  statmsg.mtype = PICL_TO_MSG(msgtype);
  statmsg.piclstat.type = ((datatype)<<NODE_BITS)|mynodeid;
  statmsg.piclstat.source = mynodeid;
  statmsg.piclstat.size = bytes;
  if ((retcode = msgsnd(stat_q[dest],(void*)&statmsg,sizeof(piclstat_t),0))==-1) {
    perror("msgsnd stat");
    exit(1);
    }

#if 0
  printf("sent a piclstat on %d:\n type=%d\n source=%d\n size=%d retcode=%d\n",
         mynodeid,statmsg.piclstat.type,statmsg.piclstat.source,statmsg.piclstat.size,retcode);
#endif

#if 0
  if (bytes > TXTSIZ) {
    printf("send0: too many bytes\n");
    exit(1);
    }
#endif

  bytes_remaining = bytes;
  offset = 0;
  while(bytes_remaining > 0) {
    number_to_send = (bytes_remaining > TXTSIZ) ? TXTSIZ: bytes_remaining;
    for (i=0; i<number_to_send; i++) {
      msg.mtext[i] = buf[i + offset];
      }

    msg.mtype = ((datatype++)<<NODE_BITS)|mynodeid;
    if (msgsnd(data_q[dest],(void*)&msg,number_to_send,0)==-1) {
      perror("msgsnd data");
      exit(1);
      }
    bytes_remaining -= number_to_send;
    offset += number_to_send;
    }

#if 0
  printf("send0: %d bytes of type %d to %d from %d\n",
         bytes,msgtype,dest,mynodeid);
#endif

}

void sync0()
{
  int sync;
  if (arc_nproc == 1) return;
  bcast0(&sync,sizeof(sync),8882,1);
  bcast0(&sync,sizeof(sync),8881,0);
}

void who0(numproc, me, host)
int *numproc, *me, *host;
{
  *me = mynodeid;
  *numproc = arc_nproc;
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

void sync0_()
{
  sync0();
}

void who0_(numproc,me,host)
int *numproc,*me,*host;
{
  who0(numproc,me,host);
}

