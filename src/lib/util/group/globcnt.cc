//
// globcnt.cc
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

#ifdef HAVE_CONFIG_H
#include <mpqc_config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <util/group/globcnt.h>

using namespace sc;

#ifndef SEM_A
#  define SEM_A 0200
#endif

#ifndef SEM_R
#  define SEM_R 0400
#endif

GlobalCounter::GlobalCounter()
{
  semid_ = -1;
  controls_release_ = 0;
}

void
GlobalCounter::cleanup()
{
  if (semid_ != -1 && controls_release_) {
      int ret;
#ifdef SEMCTL_REQUIRES_SEMUN
      semun junk;
      junk.val = 0;
#else
      int junk = 0;
#endif
      ret = semctl(semid_, 0, IPC_RMID, junk);
      if (ret == -1) {
          perror("semctl (IPC_RMID)");
          abort();
        }

      semid_ = -1;
    }
}

void
GlobalCounter::initialize()
{
  cleanup();
  semid_ = semget(IPC_PRIVATE, 1, IPC_CREAT | SEM_R | SEM_A );
  if (semid_ == -1) {
      perror("semget");
      abort();
    }
  controls_release_ = 1;
  operator = (0);
}

void
GlobalCounter::initialize(const char *stringrep)
{
  semid_ = atoi(stringrep);
  controls_release_ = 0;
}

GlobalCounter::~GlobalCounter()
{
  cleanup();
}

void
GlobalCounter::operator = (int i)
{
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
  val.val = 0;
#else
  int val = 0;
#endif
  int ret;
  if ((ret = semctl(semid_, 0, GETVAL, val)) == -1) {
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

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
