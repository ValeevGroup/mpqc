
#include <stdio.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <sys/shm.h>

#include <util/group/message.h>

#ifdef L486
#  define SEM_A 0200
#  define SEM_R 0400
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

#ifdef L486
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
#ifdef L486
      semun junk;
      junk.val = 0;
      shmctl(shmid,IPC_RMID,0);
      semctl(sync_semid,0,IPC_RMID,junk);
      semctl(sync2_semid,0,IPC_RMID,junk);
#else
      shmctl(shmid,IPC_RMID);
      semctl(sync_semid,0,IPC_RMID);
      semctl(sync2_semid,0,IPC_RMID);
#endif
    }
#ifdef L486
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

ShmMessageGrp::ShmMessageGrp()
{
  int nprocs = atoi(getenv("NUMPROC"));
  if (!nprocs) nprocs = 1;
  initialize(nprocs);
}

ShmMessageGrp::ShmMessageGrp(int nprocs)
{
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

  semun semzero;
  semzero.val = 0;
  semun semone;
  semone.val = 1;

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
  semun semzero;
  semzero.val = 0;
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
  semun semzero;
  semzero.val = 0;
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
  void *mem;
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


static void print_buffer(int node)
{
  int i;
  msgbuf_t *message;
  message = (msgbuf_t*)commbuf[node]->buf;

  printf("Printing buffer for node %d\n",node);
  for (i=0; i<commbuf[node]->nmsg; i++) {
      printf(" bytes,type,from = %6d, %6d, %3d\n",
             message->size,
             message->type,
             message->from);
      message = NEXT_MESSAGE(message);
    }

}
 

int
ShmMessageGrp::basic_probe(int msgtype)
{
  int i;
  msgbuf_t *message;

  wait_for_write(me());

  message = (msgbuf_t*)commbuf[me()]->buf;
  for (i=0; i<commbuf[me()]->nmsg; i++) {
      if ((msgtype == -1) || (msgtype == message->type)) {
          last_size_ = message->size;
          last_source_ = message->from;
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

  print_buffer(0);
  print_buffer(1);

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
      //print_buffer(me());
      //debug_start("recv buffer too small");
      fprintf(stderr,"messshm.cc: recv buffer too small\n");
      abort();
    }
  if (bytes < message->size) size = bytes;
  else size = message->size;

  // Copy the data.
  memcpy(buf,((char*)message) + sizeof(msgbuf_t),size);

  last_size_ = size;
  last_source_ = message->from;

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

