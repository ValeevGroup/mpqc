
#include <stdio.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <sys/shm.h>

#include <util/group/messshm.h>

#if defined(OSF) || defined(SUNMOS) || defined(AIX)
union semun {
  int val;
  struct semid_ds *buf;
  u_short *array;
};
#endif

//#define DEBUG

#ifndef SEM_A
#  define SEM_A 0200
#endif

#ifndef SEM_R
#  define SEM_R 0400
#endif

#ifdef L486
#  define SEMCTL_REQUIRES_SEMUN
#endif

#if defined(L486) || defined(PARAGON)
#ifndef SHMCTL_REQUIRES_SHMID
#  define SHMCTL_REQUIRES_SHMID
#endif
#endif

#if defined(L486) || defined(PARAGON)
#ifndef SHMDT_CHAR
#  define SHMDT_CHAR
#endif
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

static struct sembuf semdec;
static struct sembuf seminc;
static struct sembuf semread;

static msgbuf_t *
NEXT_MESSAGE(msgbuf_t *m)
{
  msgbuf_t *r;
  if (m->size < 0) {
      printf("NEXT_MESSAGE: m->size = %d (real %d)\n",
             m->size,sizeof(msgbuf_t) + m->size + m->size%8);
      //debug_start("m->size < 0");
      fprintf(stderr,"messshm.cc: m->size < 0\n");
      abort();
    }
  r = ((msgbuf_t*)(((char*)m) + ROUNDUPTOALIGN(sizeof(msgbuf_t) + m->size)));
  return r;
}

#define CLASSNAME ShmMessageGrp
#define PARENTS public intMessageGrp
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
ShmMessageGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  intMessageGrp::_castdown(cd);
  return do_castdowns(casts,cd);
}

ShmMessageGrp::ShmMessageGrp()
{
  initialize();
}

ShmMessageGrp::ShmMessageGrp(int nprocs)
{
  initialize(nprocs);
}

ShmMessageGrp::ShmMessageGrp(const RefKeyVal& keyval):
  intMessageGrp(keyval)
{
  int nprocs = keyval->intvalue("n");
  if (keyval->error() != KeyVal::OK) initialize();
  else initialize(nprocs);
}

// sync_semid must have semval = 0 on entry.
void
ShmMessageGrp::sync()
{
  static struct sembuf semndec;
  if (n() == 1) return;
  semndec.sem_num = 0;
  semndec.sem_op = -n() + 1;
  semndec.sem_flg = 0;
  if (me()==0) {
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

ShmMessageGrp::~ShmMessageGrp()
{
  static struct sembuf semndec;

#ifdef SHMDT_CHAR
  shmdt((char*)sharedmem);
#else
  shmdt(sharedmem);
#endif

  // sync the nodes
  sync();

  // Make sure node zero is last to touch the semaphores.
  if (n() > 1) {
      semndec.sem_num = 0;
      semndec.sem_op = -n() + 1;
      semndec.sem_flg = 0;
      if (me()==0) {
          semop(sync_semid,&semndec,1);
        }
      else {
          semop(sync_semid,&seminc,1);
        }
    }

  // Release resources.
  if (me() == 0) {
#ifdef SHMCTL_REQUIRES_SHMID
      shmctl(shmid,IPC_RMID,0);
#else
      shmctl(shmid,IPC_RMID);
#endif
#ifdef SEMCTL_REQUIRES_SEMUN
      semun junk;
      junk.val = 0;
      semctl(sync_semid,0,IPC_RMID,junk);
      semctl(sync2_semid,0,IPC_RMID,junk);
#else
      semctl(sync_semid,0,IPC_RMID);
      semctl(sync2_semid,0,IPC_RMID);
#endif
    }
#ifdef SEMCTL_REQUIRES_SEMUN
  semun junk;
  junk.val = 0;
  semctl(semid,me(),IPC_RMID,junk);
  semctl(send_semid,me(),IPC_RMID,junk);
  semctl(recv_semid,me(),IPC_RMID,junk);
#else
  semctl(semid,me(),IPC_RMID);
  semctl(send_semid,me(),IPC_RMID);
  semctl(recv_semid,me(),IPC_RMID);
#endif
}

void
ShmMessageGrp::initialize()
{
  int nprocs = atoi(getenv("NUMPROC"));
  if (!nprocs) nprocs = 1;
  initialize(nprocs);
}

void
ShmMessageGrp::initialize(int nprocs)
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

  // Each node gets a buffer for incoming data.
  shmid = shmget(IPC_PRIVATE,
                 nprocs * sizeof(commbuf_t),
                 IPC_CREAT | SHM_R | SHM_W);

  // Attach the shared segment.
  nextbuf = sharedmem = shmat(shmid,NULL,0);

#ifdef SEMCTL_REQUIRES_SEMUN
  semun semzero;
  semzero.val = 0;
  semun semone;
  semone.val = 1;
#else
  int semzero = 0;
  int semone = 1;
#endif

  // This is used for node synchronization.
  sync_semid = semget(IPC_PRIVATE,1,IPC_CREAT | SEM_R | SEM_A );
  if (sync_semid == -1) {
      perror("semget");
      exit(-1);
    }
  if (semctl(sync_semid,0,SETVAL,semzero) == -1) {
      perror("semctl");
      exit(-1);
    }
  sync2_semid = semget(IPC_PRIVATE,1,IPC_CREAT | SEM_R | SEM_A );
  if (sync_semid == -1) {
      perror("semget");
      exit(-1);
    }
  if (semctl(sync2_semid,0,SETVAL,semzero) == -1) {
      perror("semctl");
      exit(-1);
    }

  // Each shared memory segment gets a semaphore to synchronize access.
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

      // Mark all of the segments as available for writing.
      if (semctl(semid,i,SETVAL,semone) == -1) {
          perror("semctl");
          exit(-1);
        }

      if (semctl(recv_semid,i,SETVAL,semzero) == -1) {
          perror("semctl");
          exit(-1);
        }

      if (semctl(send_semid,i,SETVAL,semzero) == -1) {
          perror("semctl");
          exit(-1);
        }

      // Alloc shm for each node's commbuf.
      commbuf[i] = (commbuf_t*) nextbuf;
      // Mark the commbuf as having no messages.
      commbuf[i]->nmsg = 0;
      nextbuf = (void*)(((char*)nextbuf) + sizeof(commbuf_t));
    }

  // Create the remaining nodes.
  int mynodeid = 0;
  for (i=1; i<nprocs; i++) {
      int pid = fork();
      if (!pid) {
          mynodeid = i;
          break;
        }
    }

  // Initialize the base class.
  intMessageGrp::initialize(mynodeid, nprocs, 30);
}

static void reset_recv(int node)
{
#ifdef SEMCTL_REQUIRES_SEMUN
  semun semzero;
  semzero.val = 0;
#else
  int semzero = 0;
#endif
  semctl(recv_semid,node,SETVAL,semzero);
}

static void get_recv(int node)
{
  semdec.sem_num = node;
  semop(recv_semid,&semdec,1);
  semdec.sem_num = 0;
}

static void put_recv(int node)
{
  seminc.sem_num = node;
  semop(recv_semid,&seminc,1);
  seminc.sem_num = 0;
}

static void reset_send(int node)
{
#ifdef SEMCTL_REQUIRES_SEMUN
  semun semzero;
  semzero.val = 0;
#else
  int semzero = 0;
#endif
  semctl(send_semid,node,SETVAL,semzero);
}

static void get_send(int node)
{
  semdec.sem_num = node;
  semop(send_semid,&semdec,1);
  semdec.sem_num = 0;
}

static void put_send(int node)
{
  seminc.sem_num = node;
  semop(send_semid,&seminc,1);
  seminc.sem_num = 0;
}

// Obtain a lock for writing to the node's buffer.
static void wait_for_write(int node)
{
  semdec.sem_num = node;
  semop(semid,&semdec,1);
  semdec.sem_num = 0;
}

// Release lock for writing to node's buffer.
static void release_write(int node)
{
  seminc.sem_num = node;
  semop(semid,&seminc,1);
  seminc.sem_num = 0;
}

#ifdef DEBUG
static void print_buffer(int node, int me)
{
  int i;
  msgbuf_t *message;
  message = (msgbuf_t*)commbuf[node]->buf;

  printf("Printing buffer for node %d on node %d\n",node,me);
  for (i=0; i<commbuf[node]->nmsg; i++) {
      printf(" on node %2d: to=%2d, bytes=%6d, type=%10d, from=%2d\n",
             me,
             node,
             message->size,
             message->type,
             message->from);
      fflush(stdout);
      message = NEXT_MESSAGE(message);
    }

}
#endif

int
ShmMessageGrp::basic_probe(int msgtype)
{
  int i;
  msgbuf_t *message;

  wait_for_write(me());

  message = (msgbuf_t*)commbuf[me()]->buf;
  for (i=0; i<commbuf[me()]->nmsg; i++) {
      if ((msgtype == -1) || (msgtype == message->type)) {
          set_last_size(message->size);
          set_last_source(message->from);
          set_last_type(message->type);
          release_write(me());
          return 1;
        }
      message = NEXT_MESSAGE(message);
    }

  release_write(me());

  return 0;
}

void
ShmMessageGrp::basic_recv(int type, void* buf, int bytes)
{
  int size;
  int i;
  msgbuf_t *message,*lastmessage;

#ifdef DEBUG
  printf("node %2d recv type 0x%08x length %6d\n",
         me(), type, bytes);
  print_buffer(me(),me());
#endif

  reset_send(me());

  look_for_message:

  wait_for_write(me());

  message = (msgbuf_t*)commbuf[me()]->buf;
  for (i=0; i<commbuf[me()]->nmsg; i++) {
      if (message->type == type) break;
      message = NEXT_MESSAGE(message);
    }
  if (i==commbuf[me()]->nmsg) {
      release_write(me());
      get_send(me());
      goto look_for_message;
    }

  if (bytes < message->size) {
      fprintf(stderr,"messshm.cc: recv buffer too small\n");
      abort();
    }
  if (bytes < message->size) size = bytes;
  else size = message->size;

  // Copy the data.
  memcpy(buf,((char*)message) + sizeof(msgbuf_t),size);

  set_last_size(size);
  set_last_source(message->from);
  set_last_type(message->type);

  // Find the lastmessage.
  lastmessage = (msgbuf_t*)commbuf[me()]->buf;
  for (i=0; i<commbuf[me()]->nmsg; i++) {
      lastmessage = NEXT_MESSAGE(lastmessage);
    }

  // Repack the message buffer.
  memmove(message,NEXT_MESSAGE(message),
          (size_t)((char*)lastmessage)-(size_t)((char*)NEXT_MESSAGE(message)));

  commbuf[me()]->nmsg--;

  release_write(me());
  put_recv(me());
}

void
ShmMessageGrp::basic_send(int dest, int type, void* buf, int bytes)
{
  int i;
  msgbuf_t *availmsg;

#ifdef DEBUG
  printf("node %2d send to %2d type 0x%08x length %6d\n",
         me(), dest, type, bytes);
#endif

  if (dest>=n()) {
      //debug_start("bad destination");
      fprintf(stderr,"ShmMessageGrp::basic_send: bad destination\n");
      abort();
    }

  reset_recv(dest);

  try_send_again:

  // Obtain write access to the dest's incoming buffer.
  wait_for_write(dest);

  availmsg = (msgbuf_t*)commbuf[dest]->buf;
  for (i=0; i<commbuf[dest]->nmsg; i++) {
      availmsg = NEXT_MESSAGE(availmsg);
    }
  if (  (((char*)availmsg) + ROUNDUPTOALIGN(sizeof(msgbuf_t) + bytes))
        > (((char*)commbuf[dest]) + sizeof(commbuf_t))) {
      if (me() == dest) {
          // sending a message to myself and the buffer is full
          // --cannot recover
          printf("commbuf size exceeded on %d\n",me());
          printf(" availmsg = 0x%x\n",availmsg);
          printf(" commbuf[%d] + sizeof(commbuf_t) = 0x%x\n",
                 dest,((char*)commbuf[dest]) + sizeof(commbuf_t));
          printf(" size = %d\n",bytes);
          exit(1);
        }
      else {
          // try to recover from a full buffer by waiting for the dest
          // to read some data.
          release_write(dest);
          get_recv(dest);
          goto try_send_again;
        }
    }
  availmsg->from = me();
  availmsg->type = type;
  availmsg->size = bytes;
  memcpy(((char*)availmsg) + sizeof(msgbuf_t),buf,bytes);
  commbuf[dest]->nmsg++;

  // Release write access to the dest's buffer.
  release_write(dest);
  put_send(dest);
}

int
ShmMessageGrp::last_source()
{
  return last_source_;
}

int
ShmMessageGrp::last_size()
{
  return last_size_;
}

int
ShmMessageGrp::last_type()
{
  return msgtype_typ(last_type_);
}
