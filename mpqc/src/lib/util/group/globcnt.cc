
#ifdef __GNUG__
#pragma implementation
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <util/group/globcnt.h>

#if defined(AIX)
union semun {
  int val;
  struct semid_ds *buf;
  u_short *array;
};
#endif

#ifndef SEM_A
#  define SEM_A 0200
#endif

#ifndef SEM_R
#  define SEM_R 0400
#endif

#ifdef L486
#  define SEMCTL_REQUIRES_SEMUN
#endif

GlobalCounter::GlobalCounter()
{
  semid_ = -1;
}

void
GlobalCounter::initialize()
{
  semid_ = semget(IPC_PRIVATE, 1, IPC_CREAT | SEM_R | SEM_A );
  if (semid_ == -1) {
      perror("semget");
      abort();
    }
  controls_release_ = 1;
  operator = (0);
  //printf("got semid = %d\n", semid_);
}

void
GlobalCounter::initialize(const char *stringrep)
{
  semid_ = atoi(stringrep);
  //printf("setting semid (0x%08x) = %d\n", &semid_, semid_);
  controls_release_ = 0;
}

GlobalCounter::~GlobalCounter()
{
  if (semid_ != -1 && controls_release_) {
      //printf("removing semid = %d\n", semid_);
      int ret;
#ifdef SEMCTL_REQUIRES_SEMUN
      semun junk;
      junk.val = 0;
      ret = semctl(semid_, 0, IPC_RMID, junk);
#else
      ret = semctl(semid_, 0, IPC_RMID);
#endif
      if (ret == -1) {
          perror("semctl (IPC_RMID)");
          abort();
        }
    }
}

void
GlobalCounter::operator = (int i)
{
  //printf("GlobalCounter: assigning %5d (0x%08x) to %3d\n",
  //       semid_, &semid_, i);
#ifdef SEMCTL_REQUIRES_SEMUN
  semun val;
  val.val = i;
#else
  int val = i;
#endif
  if (semctl(semid_, 0, SETVAL, val) == -1) {
      perror("semctl (SETVAL)");
      abort();
    }
}

int
GlobalCounter::val()
{
#ifdef SEMCTL_REQUIRES_SEMUN
  semun val;
#else
  int val;
#endif
  int ret;
  if (ret = semctl(semid_, 0, GETVAL, val) == -1) {
      perror("semctl (GETVAL)");
      abort();
    }
  return ret;
}

void
GlobalCounter::wait_for_zero()
{
  operator += (0);
}

void
GlobalCounter::operator+=(int i)
{
  //printf("GlobalCounter: incrementing %5d by %3d\n", semid_, i);
  struct sembuf s;
  s.sem_num = 0;
  s.sem_op = i;
  s.sem_flg = 0;
  if (semop(semid_, &s, 1) == -1) {
      perror("semop");
      abort();
    }
}

void
GlobalCounter::operator--()
{
  operator += (-1);
}

void
GlobalCounter::operator++()
{
  operator += (1);
}

char *
GlobalCounter::stringrep()
{
  char tmp[80];
  sprintf(tmp, "%d", semid_);
  return strcpy(new char[strlen(tmp)+1], tmp);
}

