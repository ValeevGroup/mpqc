//
// messshm.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//


#include <unistd.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <sys/shm.h>


#include <util/misc/bug.h>
#include <util/misc/formio.h>
#include <util/group/messshm.h>

using namespace std;
using namespace sc;

//#define DEBUG

#ifndef SEM_A
#  define SEM_A 0200
#endif

#ifndef SEM_R
#  define SEM_R 0400
#endif

/* NALIGN is the byte boundary that we align data on. */
#define NALIGN 8
#define ROUNDUPTOALIGN(n) (((n) + (NALIGN-1)) & ~(NALIGN-1))

static ClassDesc ShmMessageGrp_cd(
  typeid(ShmMessageGrp),"ShmMessageGrp",1,"public intMessageGrp",
  0, create<ShmMessageGrp>, 0);

ShmMessageGrp::ShmMessageGrp()
{
  initialize();
}

ShmMessageGrp::ShmMessageGrp(int nprocs)
{
  initialize(nprocs);
}

ShmMessageGrp::ShmMessageGrp(const Ref<KeyVal>& keyval):
  intMessageGrp(keyval)
{
  int nprocs = keyval->intvalue("n");
  if (keyval->error() != KeyVal::OK) initialize();
  else initialize(nprocs);
}

void ShmMessageGrp::sync()
{
  int i;
  for (i=0; i<n(); i++) {
      if (me() == i) continue;
      wait_for_write(i);
      commbuf[i]->n_sync++;
      if (commbuf[i]->n_sync >= n()-1) {
          while(commbuf[i]->n_wait_for_change) {
              put_change(i);
              commbuf[i]->n_wait_for_change--;
            }
        }
      release_write(i);
    }
  wait_for_write(me());
  while (commbuf[me()]->n_sync < n()-1) {
      commbuf[me()]->n_wait_for_change++;
      release_write(me());
      get_change(me());
      wait_for_write(me());
    }
  commbuf[me()]->n_sync -= n()-1;
  while(commbuf[me()]->n_wait_for_change) {
      put_change(me());
      commbuf[me()]->n_wait_for_change--;
    }
  release_write(me());
}

ShmMessageGrp::~ShmMessageGrp()
{
  // sync the nodes
  sync();

  // make sure node zero is las to touch the shared memory
  if (me() == 0) {
      wait_for_write(0);
      while (commbuf[0]->n_sync < n()-1) {
          commbuf[0]->n_wait_for_change++;
          release_write(0);
          get_change(0);
          wait_for_write(0);
        }
      release_write(0);
      shmdt((SHMTYPE)sharedmem);
      // release the memory
      shmctl(shmid,IPC_RMID,0);

      for (int i=0; i<n(); i++) {
#ifdef SEMCTL_REQUIRES_SEMUN
          semun junk;
          junk.val = 0;
#else
          int junk = 0;
#endif
          semctl(semid,i,IPC_RMID,junk);
          semctl(change_semid,i,IPC_RMID,junk);
        }
    }
  else {
      wait_for_write(0);
      commbuf[0]->n_sync++;
      while(commbuf[0]->n_wait_for_change) {
          put_change(0);
          commbuf[0]->n_wait_for_change--;
        }
      shmdt((SHMTYPE)sharedmem);
      release_write(0);
    }
}

void ShmMessageGrp::initialize()
{
  int nprocs = atoi(getenv("NUMPROC"));
  if (!nprocs) nprocs = 1;
  initialize(nprocs);
}

void ShmMessageGrp::initialize(int nprocs)
{
  int i;
  void* nextbuf;

  semdec.sem_num = 0;
  semdec.sem_op = -1;
  semdec.sem_flg = 0;
  seminc.sem_num = 0;
  seminc.sem_op =  1;
  seminc.sem_flg = 0;

  // Each node gets a buffer for incoming data.
  shmid = shmget(IPC_PRIVATE,
                 nprocs * sizeof(commbuf_t),
                 IPC_CREAT | SHM_R | SHM_W);

  // Attach the shared segment.
  nextbuf = sharedmem = shmat(shmid,0,0);

#ifdef SEMCTL_REQUIRES_SEMUN
  semun semzero;
  semzero.val = 0;
  semun semone;
  semone.val = 1;
#else
  int semzero = 0;
  int semone = 1;
#endif

  // Each shared memory segment gets a semaphore to synchronize access.
  semid = semget(IPC_PRIVATE,nprocs,IPC_CREAT | SEM_R | SEM_A );
  if (semid == -1) {
      perror("semget");
      exit(-1);
    }
  
  change_semid = semget(IPC_PRIVATE,nprocs,IPC_CREAT | SEM_R | SEM_A );
  if (change_semid == -1) {
      perror("semget");
      exit(-1);
    }

  for (i=0; i<nprocs; i++) {

      // Mark all of the segments as available for writing.
      if (semctl(semid,i,SETVAL,semone) == -1) {
          perror("semctl");
          exit(-1);
        }

      if (semctl(change_semid,i,SETVAL,semzero) == -1) {
          perror("semctl");
          exit(-1);
        }

      // Alloc shm for each node's commbuf.
      commbuf[i] = (commbuf_t*) nextbuf;
      // Mark the commbuf as having no messages.
      commbuf[i]->nmsg = 0;
      // and no outstanding waits
      commbuf[i]->n_wait_for_change = 0;
      commbuf[i]->n_sync = 0;
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

Ref<MessageGrp>
ShmMessageGrp::clone(void)
{
  Ref<MessageGrp> smgrp = new ShmMessageGrp(n());
  return smgrp;
}

int ShmMessageGrp::basic_probe(int msgtype)
{
  int i;
  msgbuf_t *message;

  wait_for_write(me());

  message = (msgbuf_t*)commbuf[me()]->buf;
  for (i=0; i<commbuf[me()]->nmsg; i++) {
      if ((msgtype == -1) || (msgtype == message->type)) {
          release_write(me());
          return 1;
        }
      message = NEXT_MESSAGE(message);
    }

  release_write(me());

  return 0;
}

void ShmMessageGrp::basic_recv(int type, void* buf, int bytes)
{
  int size;
  int i;
  msgbuf_t *message,*lastmessage;

#ifdef DEBUG
  ExEnv::outn() << "ShmGrp: basic_recv: "
       << "type = " << type << ' '
       << "buf = " << buf << ' '
       << "bytes = " << bytes << ' '
       << "me = " << me() << endl;
  print_buffer(me(),me());
#endif

  look_for_message:

  wait_for_write(me());

  message = (msgbuf_t*)commbuf[me()]->buf;
  for (i=0; i<commbuf[me()]->nmsg; i++) {
      if (message->type == type) break;
      message = NEXT_MESSAGE(message);
    }
  if (i==commbuf[me()]->nmsg) {
      commbuf[me()]->n_wait_for_change++;
      release_write(me());
      get_change(me());
      goto look_for_message;
    }

  if (bytes < message->size) {
      ExEnv::errn() << scprintf("messshm.cc: recv buffer too small\n");
      abort();
    }
  if (bytes < message->size) size = bytes;
  else size = message->size;

  // Copy the data.
  memcpy(buf,((char*)message) + sizeof(msgbuf_t),size);

  // Find the lastmessage.
  lastmessage = (msgbuf_t*)commbuf[me()]->buf;
  for (i=0; i<commbuf[me()]->nmsg; i++) {
      lastmessage = NEXT_MESSAGE(lastmessage);
    }

  // Repack the message buffer.
  memmove(message,NEXT_MESSAGE(message),
          (size_t)((char*)lastmessage)-(size_t)((char*)NEXT_MESSAGE(message)));

  commbuf[me()]->nmsg--;

  while(commbuf[me()]->n_wait_for_change) {
      put_change(me());
      commbuf[me()]->n_wait_for_change--;
    }

  release_write(me());
}

void ShmMessageGrp::basic_send(int dest, int type, const void* buf, int bytes)
{
  int i;
  msgbuf_t *availmsg;

#ifdef DEBUG
  ExEnv::outn() << "ShmGrp: basic_send: "
       << "dest = " << dest << ' '
       << "type = " << type << ' '
       << "buf = " << buf << ' '
       << "bytes = " << bytes << ' '
       << "me = " << me() << endl;
#endif

  if (dest>=n()) {
      //debug_start("bad destination");
      ExEnv::errn() << scprintf("ShmMessageGrp::basic_send: bad destination\n");
      abort();
    }

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
          ExEnv::errn() << scprintf("commbuf size exceeded on %d\n",me());
          ExEnv::errn() << scprintf(" availmsg = 0x%x\n",availmsg);
          ExEnv::errn() << scprintf(" commbuf[%d] + sizeof(commbuf_t) = 0x%x\n",
                           dest,((char*)commbuf[dest]) + sizeof(commbuf_t));
          ExEnv::errn() << scprintf(" size = %d\n",bytes);
          abort();
        }
      else {
          // try to recover from a full buffer by waiting for the dest
          // to read some data.
          commbuf[dest]->n_wait_for_change++;
          release_write(dest);
          get_change(dest);
          goto try_send_again;
        }
    }
  availmsg->from = me();
  availmsg->type = type;
  availmsg->size = bytes;
  memcpy(((char*)availmsg) + sizeof(msgbuf_t),buf,bytes);
  commbuf[dest]->nmsg++;

  // let the dest know that there is more data in the buffer
  while(commbuf[dest]->n_wait_for_change) {
      put_change(dest);
      commbuf[dest]->n_wait_for_change--;
    }

  // Release write access to the dest's buffer.
  release_write(dest);
}

msgbuf_t * ShmMessageGrp::NEXT_MESSAGE(msgbuf_t *m)
{
  msgbuf_t *r;
  if (m->size < 0) {
      ExEnv::errn() << scprintf("NEXT_MESSAGE: m->size = %d (real %d)\n",
                       m->size,sizeof(msgbuf_t) + m->size + m->size%8);
      //debug_start("m->size < 0");
      ExEnv::errn() << scprintf("messshm.cc: m->size < 0\n");
      abort();
    }
  r = ((msgbuf_t*)(((char*)m) + ROUNDUPTOALIGN(sizeof(msgbuf_t) + m->size)));
  return r;
}

void ShmMessageGrp::get_change(int node)
{
  semdec.sem_num = node;
  semop(change_semid,&semdec,1);
  semdec.sem_num = 0;
}

void ShmMessageGrp::put_change(int node)
{
  seminc.sem_num = node;
  semop(change_semid,&seminc,1);
  seminc.sem_num = 0;  

}

// Obtain a lock for writing to the node's buffer.
void ShmMessageGrp::wait_for_write(int node)
{
  semdec.sem_num = node;
  semop(semid,&semdec,1);
  semdec.sem_num = 0;
}

// Release lock for writing to node's buffer.
void ShmMessageGrp::release_write(int node)
{
  seminc.sem_num = node;
  semop(semid,&seminc,1);
  seminc.sem_num = 0;
}

#ifdef DEBUG
void ShmMessageGrp::print_buffer(int node, int me)
{
  int i;
  msgbuf_t *message;
  message = (msgbuf_t*)commbuf[node]->buf;

  ExEnv::outn() << scprintf("Printing buffer for node %d on node %d\n",node,me);
  for (i=0; i<commbuf[node]->nmsg; i++) {
      ExEnv::outn() <<
          scprintf(" on node %2d: to=%2d, bytes=%6d, type=%10d, from=%2d\n",
                   me,
                   node,
                   message->size,
                   message->type,
                   message->from);
      ExEnv::outn().flush();
      message = NEXT_MESSAGE(message);
    }

}
#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
