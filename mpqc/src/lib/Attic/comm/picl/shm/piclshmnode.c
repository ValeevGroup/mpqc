/*      PICL to shared memory interface library      
  All questions and criticisms should be directed to 
               Curtis Janssen
     Center for Computational Engineering
	   Sandia National Laboratory
    		Livermore CA 94551                         */



#include <stdio.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <sys/shm.h>

#ifndef OLDCLOCK
#  include <sys/types.h>
#  include <sys/time.h>
#  include <sys/resource.h>
#  if defined(sun) || defined(AIX)
#  include <unistd.h>
#  endif
#endif

/* Set the maximum number of processors (including the host). */
#define MAXPROCS 17

/* NALIGN is the byte boundary that we align data on. */
#define NALIGN 8
#define ROUNDUPTOALIGN(n) (((n) + (NALIGN-1)) & ~(NALIGN-1))

#define SHMCOMMBUFSIZE 1000000

struct commbuf_struct {
  int nmsg;
  char buf[SHMCOMMBUFSIZE];
  };
typedef struct commbuf_struct commbuf_t;

struct msgbuf_struct {
  int type;
  int from;
  int size;
  };
typedef struct msgbuf_struct msgbuf_t;

static commbuf_t *commbuf[MAXPROCS];
static int shmid;
static int sync_semid;
static int sync2_semid;
static int semid;
static int send_semid;
static int recv_semid;
static void* sharedmem;

/* 
  Here we specify a HOSTID.  Its id must be different than all of
  the node ids.
 */

#define HOSTID (MAXPROCS-1)


/* Define some externals that will store info between calls */
static int last_count;   /* holds the size of the last msg sent or recvd */
static int last_node;    /* holds the id of where last msg came from */
static int last_type;    /* holds the type of where last msg recvd */
static int nprocs;       /* holds the number of processors allocated */
static int mynodeid;

/* These are set by the setarc0 and getarc0 call. */
static int arc_nproc;
static int arc_top=-1;
static int arc_ord=0;
static int arc_dir=1;

/* There is never a host program in this version. */
static int hostavail = 0;

static struct sembuf semdec;
static struct sembuf seminc;
static struct sembuf semread;

/* This is wrong: #define NEXT_MESSAGE(m) ((msgbuf_t*)(((char*)m) + sizeof(msgbuf_t) + m->size + m->size%8)) */
static msgbuf_t *
NEXT_MESSAGE(m)
msgbuf_t *m;
{
  msgbuf_t *r;
  if (m->size < 0) {
    printf("NEXT_MESSAGE: m->size = %d (real %d)\n",m->size,sizeof(msgbuf_t) + m->size + m->size%8);
    debug_start("m->size < 0");
    }
  r = ((msgbuf_t*)(((char*)m) + ROUNDUPTOALIGN(sizeof(msgbuf_t) + m->size)));
#if 0
  printf("NEXT_MESSAGE: m = 0x%x, size = %d, r = 0x%x\n",m,m->size,r);
#endif
  return r;
  }

int host0() {return hostavail;}

/* sync_semid must have semval = 0 on entry. */
void sync0()
{
  static struct sembuf semndec;
  if (nprocs == 1) return;
  semndec.sem_num = 0;
  semndec.sem_op = -nprocs + 1;
  semndec.sem_flg = 0;
  if (mynodeid==0) {
    semop(sync_semid,&semndec,1);
    semop(sync2_semid,&semndec,1);
    }
  else {
    semop(sync_semid,&seminc,1);
    semop(sync_semid,&semread,1);
    semop(sync2_semid,&seminc,1);
    semop(sync2_semid,&semread,1);
    }
}

void setarc0(np,top,ord,dir)
int *np;
int *top;
int *ord;
int *dir;
{
  arc_nproc = nprocs;
  *np = arc_nproc;
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

void message0(message)
char *message;
{
  printf ("%s\n", message);
}

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC ((int) 1000000)
#endif
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

void close0(release)
     int release;
{
  static struct sembuf semndec;

  shmdt(sharedmem);

  /* sync the nodes */
  sync0();

  /* Make sure node zero is last to touch the semaphores. */
  if (nprocs > 1) {
    semndec.sem_num = 0;
    semndec.sem_op = -nprocs + 1;
    semndec.sem_flg = 0;
    if (mynodeid==0) {
      semop(sync_semid,&semndec,1);
      }
    else {
      semop(sync_semid,&seminc,1);
      }
    }

  /* Release resources. */
  if (mynodeid == 0) {
    shmctl(shmid,IPC_RMID);
    semctl(sync_semid,0,IPC_RMID);
    semctl(sync2_semid,0,IPC_RMID);
    }
  semctl(semid,mynodeid,IPC_RMID);
  semctl(send_semid,mynodeid,IPC_RMID);
  semctl(recv_semid,mynodeid,IPC_RMID);
}

void open0(numproc, me, host)
     int *numproc, *me, *host;
{
  int i;
  void* nextbuf;

  semdec.sem_num = 0;
  semdec.sem_op = -1;
  semdec.sem_flg = 0;
  seminc.sem_num = 0;
  seminc.sem_op =  1;
  seminc.sem_flg = 0;
  semread.sem_num = 0;
  semread.sem_op =  0;
  semread.sem_flg = 0;

  nprocs = atoi(getenv("NUMPROC"));
  if (!nprocs) nprocs = 1;

  /* Each node gets a buffer for incoming data. */
  shmid = shmget(IPC_PRIVATE,
                 nprocs * sizeof(commbuf_t),
                 IPC_CREAT | SHM_R | SHM_W);

  /* Attach the shared segment. */
  nextbuf = sharedmem = shmat(shmid,NULL,0);

  /* This is used for node synchronization. */
  sync_semid = semget(IPC_PRIVATE,1,IPC_CREAT | SEM_R | SEM_A );
  if (sync_semid == -1) {
    perror("semget");
    exit(-1);
    }
  if (semctl(sync_semid,0,SETVAL,0) == -1) {
    perror("semctl");
    exit(-1);
    }
  sync2_semid = semget(IPC_PRIVATE,1,IPC_CREAT | SEM_R | SEM_A );
  if (sync_semid == -1) {
    perror("semget");
    exit(-1);
    }
  if (semctl(sync2_semid,0,SETVAL,0) == -1) {
    perror("semctl");
    exit(-1);
    }

  /* Each shared memory segment gets a semaphore to synchronize access. */
  semid = semget(IPC_PRIVATE,nprocs,IPC_CREAT | SEM_R | SEM_A );
  if (semid == -1) {
    perror("semget");
    exit(-1);
    }

  send_semid = semget(IPC_PRIVATE,nprocs,IPC_CREAT | SEM_R | SEM_A );
  if (send_semid == -1) {
    perror("semget");
    exit(-1);
    }

  recv_semid = semget(IPC_PRIVATE,nprocs,IPC_CREAT | SEM_R | SEM_A );
  if (recv_semid == -1) {
    perror("semget");
    exit(-1);
    }

  for (i=0; i<nprocs; i++) {

    /* Mark all of the segments as available for writing. */
    if (semctl(semid,i,SETVAL,1) == -1) {
      perror("semctl");
      exit(-1);
      }

    if (semctl(recv_semid,i,SETVAL,0) == -1) {
      perror("semctl");
      exit(-1);
      }

    if (semctl(send_semid,i,SETVAL,0) == -1) {
      perror("semctl");
      exit(-1);
      }

    /* Alloc shm for each node's commbuf. */
    commbuf[i] = nextbuf;
    /* Mark the commbuf as having no messages. */
    commbuf[i]->nmsg = 0;
    nextbuf = (void*)(((char*)nextbuf) + sizeof(commbuf_t));
    }

  /* Create the remaining nodes. */
  mynodeid = 0;
  if (!hostavail) {
    for (i=1; i<nprocs; i++) {
      int pid = fork();
      if (!pid) {
        /* printf("got forked node %d\n",i); */
        mynodeid = i;
        break;
        }
      }
    }

  *me      = mynodeid;
  *numproc = nprocs;
  *host    = HOSTID;

}

static void reset_recv(node)
int node;
{
  semctl(recv_semid,node,SETVAL,0);
  }

static void get_recv(node)
int node;
{
  semdec.sem_num = node;
  semop(recv_semid,&semdec,1);
  semdec.sem_num = 0;
  }

static void put_recv(node)
int node;
{
  seminc.sem_num = node;
  semop(recv_semid,&seminc,1);
  seminc.sem_num = 0;
  }

static void reset_send(node)
int node;
{
  semctl(send_semid,node,SETVAL,0);
  }

static void get_send(node)
int node;
{
  semdec.sem_num = node;
  semop(send_semid,&semdec,1);
  semdec.sem_num = 0;
  }

static void put_send(node)
int node;
{
  seminc.sem_num = node;
  semop(send_semid,&seminc,1);
  seminc.sem_num = 0;
  }

/* Obtain a lock for writing to the node's buffer. */
static void wait_for_write(node)
int node;
{
  void *mem;
  semdec.sem_num = node;
  semop(semid,&semdec,1);
  semdec.sem_num = 0;
  }

/* Release lock for writing to node's buffer. */
static void release_write(node)
int node;
{
  seminc.sem_num = node;
  semop(semid,&seminc,1);
  seminc.sem_num = 0;
  }


static void print_buffer(node)
int node;
{
  int i;
  msgbuf_t *message;
  message = (msgbuf_t*)commbuf[node]->buf;

  printf("Printing buffer on node %d for node %d\n",mynodeid,node);
  for (i=0; i<commbuf[node]->nmsg; i++) {
    printf(" bytes,type,from = %6d, %6d, %3d\n",
           message->size,
           message->type,
           message->from);
    message = NEXT_MESSAGE(message);
    }

  }
 

int probe0(msgtype)
int msgtype;
{
  int i;
  msgbuf_t *message;

  wait_for_write(mynodeid);

  message = (msgbuf_t*)commbuf[mynodeid]->buf;
  for (i=0; i<commbuf[mynodeid]->nmsg; i++) {
    if ((msgtype == -1) || (msgtype == message->type)) {
      last_count = message->size;
      last_type = message->type;
      last_node = message->from;
      release_write(mynodeid);
      return 1;
      }
    message = NEXT_MESSAGE(message);
    }

  release_write(mynodeid);

  return 0;
}

recv0(buf,bytes,type)
     char *buf;
     int bytes, type;
{
  int size;
  int i;
  msgbuf_t *message,*lastmessage;

#if 0
  printf("on %d: receiving %d, size is %d\n",mynodeid,type,bytes);
#endif

  reset_send(mynodeid);

  look_for_message:

  wait_for_write(mynodeid);

  message = (msgbuf_t*)commbuf[mynodeid]->buf;
  for (i=0; i<commbuf[mynodeid]->nmsg; i++) {
    if (message->type == type) break;
    message = NEXT_MESSAGE(message);
    }
  if (i==commbuf[mynodeid]->nmsg) {
#if 0
    printf("on %d: type %d not found (waiting)\n",mynodeid,type);
#endif
    release_write(mynodeid);
    get_send(mynodeid);
#if 0
    printf("on %d: look for message for type %d\n",mynodeid,type);
#endif
    goto look_for_message;
    }

  if (bytes < message->size) {
    print_buffer(mynodeid);
    debug_start("recv0 buffer too small");
    }
  if (bytes < message->size) size = bytes;
  else size = message->size;

  /* Copy the data. */
  /* bcopy(((char*)message) + sizeof(msgbuf_t),buf,size); */
  memcpy(buf,((char*)message) + sizeof(msgbuf_t),size);

  /* Update the locals for recinfo0. */
  last_count = size;
  last_type = message->type;
  last_node = message->from;

  /* Find the lastmessage. */
  lastmessage = (msgbuf_t*)commbuf[mynodeid]->buf;
  for (i=0; i<commbuf[mynodeid]->nmsg; i++) {
    lastmessage = NEXT_MESSAGE(lastmessage);
    }

  /* Repack the message buffer. */
  /* bcopy(NEXT_MESSAGE(message),message,
   *    (int)((char*)lastmessage)-(int)((char*)NEXT_MESSAGE(message)));
   */
  memmove(message,NEXT_MESSAGE(message),
        (int)((char*)lastmessage)-(int)((char*)NEXT_MESSAGE(message)));

  commbuf[mynodeid]->nmsg--;

  release_write(mynodeid);
  put_recv(mynodeid);

#if 0
  printf("on %d: received %d from %d\n",mynodeid,type,last_node);
#endif
}

void recvinfo0(bytes, msgtype, source)
int *bytes, *msgtype, *source;
{
  *bytes = last_count;
  *msgtype = last_type;
  *source = last_node;
}

void send0(buf, bytes, type, dest)
     char *buf;
     int bytes, type, dest;
{
  int i;
  msgbuf_t *availmsg;

#if 0
  printf("on %d: sending %d to %d, size is %d\n",mynodeid,type,dest,bytes);
#endif
  if (dest>=nprocs) {
    debug_start("bad destination");
    }

  reset_recv(dest);

  try_send_again:

  /* Obtain write access to the dest's incoming buffer. */
  wait_for_write(dest);

#if 0
  print_buffer(dest);
#endif
  availmsg = (msgbuf_t*)commbuf[dest]->buf;
#if 0
  printf("on %d: availmsg = 0x%x, sharedmem = 0x%x\n",
         mynodeid,availmsg,sharedmem);
#endif
  for (i=0; i<commbuf[dest]->nmsg; i++) {
    availmsg = NEXT_MESSAGE(availmsg);
#if 0
    printf("on %d: availmsg = 0x%x, sharedmem = 0x%x\n",
           mynodeid,availmsg,sharedmem);
#endif
    }
  if (  (((char*)availmsg) + ROUNDUPTOALIGN(sizeof(msgbuf_t) + bytes))
      > (((char*)commbuf[dest]) + sizeof(commbuf_t))) {
    if (mynodeid == dest) {
      /* sending a message to myself and the buffer is full
       * --cannot recover */
      printf("commbuf size exceeded on %d\n",mynodeid);
      printf(" availmsg = 0x%x\n",availmsg);
      printf(" commbuf[%d] + sizeof(commbuf_t) = 0x%x\n",
             dest,((char*)commbuf[dest]) + sizeof(commbuf_t));
      printf(" size = %d\n",bytes);
      exit(1);
      }
    else {
      /* try to recover from a full buffer by waiting for the dest
       * to read some data. */
      release_write(dest);
      get_recv(dest);
      goto try_send_again;
      }
    }
  availmsg->from = mynodeid;
  availmsg->type = type;
  availmsg->size = bytes;
  /* bcopy(buf,((char*)availmsg) + sizeof(msgbuf_t),bytes); */
  memcpy(((char*)availmsg) + sizeof(msgbuf_t),buf,bytes);
  commbuf[dest]->nmsg++;

  /* print_buffer(dest); */

  /* Release write access to the dest's buffer. */
  release_write(dest);
  put_send(dest);

#if 0
  printf("on %d: sent %d to %d\n",mynodeid,type,dest);
#endif
}

void who0(numproc,me,host)
     int *numproc,*me,*host;
{

  *numproc = nprocs;
  *me      = mynodeid;
  *host    = HOSTID;
}

