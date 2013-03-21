//
// memshm.cc
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

#ifndef _util_group_memshm_cc
#define _util_group_memshm_cc

#include <util/misc/formio.h>
#include <util/group/pool.h>
#include <util/group/memshm.h>
#include <limits.h>
#include <errno.h>

using namespace std;
using namespace sc;

#ifndef SHMMAX
// glibc 2.0.3 isn't defining SHMMAX so make set it here
#  ifdef __linux__
#    define SHMMAX 0x1000000
#  else
#    define SHMMAX INT_MAX
#  endif
#endif

#ifndef SHMMIN
#define SHMMIN 1
#endif

#ifndef SIMPLE_LOCK
#define SIMPLE_LOCK 1
#endif

#undef DEBUG

static ClassDesc ShmMemoryGrp_cd(
  typeid(ShmMemoryGrp),"ShmMemoryGrp",1,"public MsgMemoryGrp",
  0, create<ShmMemoryGrp>, 0);

ShmMemoryGrp::ShmMemoryGrp(const Ref<MessageGrp>& msg):
  MsgMemoryGrp(msg)
{
  update_ = 0;
  data_ = 0;
  memory_ = 0;
  pool_ = 0;
  rangelock_ = 0;
  nregion_ = 0;
  shmid_ = 0;
  attach_address_ = 0;
}

ShmMemoryGrp::ShmMemoryGrp(const Ref<KeyVal>& keyval):
  MsgMemoryGrp(keyval)
{
  update_ = 0;
  data_ = 0;
  memory_ = 0;
  pool_ = 0;
  rangelock_ = 0;
  nregion_ = 0;
  shmid_ = 0;
  attach_address_ = 0;
}

int
ShmMemoryGrp::attach_memory(void *ataddress, int size)
{
  int i;
  int fail = 0;
  int isize;
  int rsize = size;
  for (i=0; i<nregion_; i++) {
      attach_address_[i] = 0;
    }
  for (i=0; rsize>0; i++,rsize-=isize) {
      isize = rsize;
      if (isize > SHMMAX) isize = SHMMAX;
      else if (isize < SHMMIN) isize = SHMMIN;
      if (debug_) {
          ExEnv::outn() << me() << ": ";
          ExEnv::outn() << "ShmMemoryGrp: attaching segment with "
               << isize << " bytes at address " << (void*)ataddress
               << " on node " << me()
               << endl;
        }
      attach_address_[i] = shmat(shmid_[i],(SHMTYPE)ataddress,0);
      if (debug_) {
          ExEnv::outn() << me() << ": ";
          ExEnv::outn() << "ShmMemoryGrp: got address "
               << (void*)attach_address_[i]
               << " on node " << me()
               << endl;
        }
      if ((attach_address_[i] == 0)
          || (attach_address_[i] == ((void*) -1))
          || (ataddress && (attach_address_[i] != ataddress))) {
          //ExEnv::outn() << "ShmMemoryGrp: shmat: problem attaching using address: "
          //     << " " << (void*) ataddress
          //     << ": got address: "
          //     << (void*) attach_address_[i]
          //     << " on node " << me()
          //     << endl;
          fail = 1;
        }
      ataddress = (void*)((char*)(attach_address_[i]) + isize);
    }

  memory_ = (void*) attach_address_[0];

  if (fail) detach_memory();

  return fail;
}

void
ShmMemoryGrp::detach_memory()
{
  int i;
  for (i=0; i<nregion_; i++) {
      if (attach_address_[i] != 0 && attach_address_[i] != (void*) -1) {
          if (debug_) {
              ExEnv::outn() << "detaching " << (void*)attach_address_[i]
                   << " on node " << me() << endl;
            }
          shmdt((SHMTYPE)attach_address_[i]);
        }
      attach_address_[i] = 0;
    }
}

void
ShmMemoryGrp::set_localsize(size_t localsize)
{
  int i;

  cleanup();

  MsgMemoryGrp::set_localsize(localsize);

  const int poolallocation = 160000;

  update_ = new GlobalCounter[n()];

  // allocate memory both the data and the Pool
  int size = poolallocation + distsize_to_size(totalsize());
  // compute the number of shared memory regions that will be needed
  nregion_ = size/SHMMAX;
  if (size%SHMMAX) nregion_++;
  shmid_ = new int[nregion_];
  attach_address_ = new void*[nregion_];

  // get the shared memory segments
  if (me() == 0) {
      int rsize = size;
      int isize;
      for (i=0; rsize>0; i++,rsize-=isize) {
          isize = rsize;
          if (isize > SHMMAX) isize = SHMMAX;
          else if (isize < SHMMIN) isize = SHMMIN;
          if (debug_) {
              ExEnv::outn() << me() << ": ";
              ExEnv::outn() << "ShmMemoryGrp: getting segment with " << isize
                   << " bytes" << endl;
            }
          shmid_[i] = shmget(IPC_PRIVATE, isize, IPC_CREAT | SHM_R | SHM_W);
          if (shmid_[i] == -1) {
              ExEnv::outn() << me() << ": ";
              ExEnv::outn() << "ShmMemoryGrp: shmget failed for "
                   << isize << " bytes: "
                   << strerror(errno) << endl;
              abort();
            }
        }
    }

  if (me() == 0) {
      lock_.initialize();
      lock_ = 1;
      char * stringrep = lock_.stringrep();
      int length = strlen(stringrep) + 1;
      msg_->bcast(&length, 1);
      msg_->bcast(stringrep, length);
#ifdef DEBUG
      ExEnv::outn() << scprintf("%d: sent initialize string of \"%s\" (%d)\n",
                       me(), stringrep, length);
      ExEnv::outn().flush();
#endif // DEBUG
      delete[] stringrep;
      for (i=0; i<n(); i++) {
          update_[i].initialize();
          char * stringrep = update_[i].stringrep();
          int length = strlen(stringrep) + 1;
          msg_->bcast(&length, 1);
          msg_->bcast(stringrep, length);
#ifdef DEBUG
          ExEnv::outn() << scprintf("%d: sent initialize string of \"%s\" (%d) for %d\n",
                           me(), stringrep, length, i);
          ExEnv::outn().flush();
#endif // DEBUG
          delete[] stringrep;
        }
    }
  else {
      int length;
      msg_->bcast(&length, 1);
      char * stringrep = new char[length];
      msg_->bcast(stringrep, length);
#ifdef DEBUG
      ExEnv::outn() << scprintf("%d: got initialize string of \"%s\" (%d)\n",
                       me(), stringrep, length);
      ExEnv::outn().flush();
#endif // DEBUG
      lock_.initialize(stringrep);
      delete[] stringrep;
      for (i=0; i<n(); i++) {
          msg_->bcast(&length, 1);
          stringrep = new char[length];
          msg_->bcast(stringrep, length);
#ifdef DEBUG
          ExEnv::outn() << scprintf("%d: got initialize string of \"%s\" (%d) for %d\n",
                           me(), stringrep, length, i);
          ExEnv::outn().flush();
#endif // DEBUG
          update_[i].initialize(stringrep);
          delete[] stringrep;
        }
    }

  // get an initial starting address on node 0
  int ntry = 20;
  int itry;
  void *address;
  if (me() == 0) {
      address = 0;
      itry = 0;
      while (attach_memory(address,size) && itry < ntry) {
          if (address == 0) address = memory_;
          else address = (void*) &((char*)address)[0x1000000];
          itry++;
        }
      if (itry == ntry) {
          ExEnv::errn() << "ShmMemoryGrp: ntry exhausted on node 0" << endl;
          abort();
        }
      // detach again, since we all try together below
      detach_memory();
    }

  msg_->bcast(shmid_, nregion_);
  msg_->raw_bcast((void*)&memory_, sizeof(void*));

  address = memory_;
  itry = 0;
  int fail;
  do {
      fail = attach_memory(address,size);
      msg_->max(fail);
      if (fail) {
          detach_memory();
        }
      address = (void*) &((char*)address)[0x1000000];
      itry++;
    } while(fail && itry < ntry);
  if (itry == ntry) {
      ExEnv::errn() << "ShmMemoryGrp: ntry exhausted on node " << me()
           << " on joint attach phase" << endl;
      abort();
    }

  if (me() == 0) {
      // initialize the pool
      pool_ = new(memory_) Pool(poolallocation);
      rangelock_ = new(pool_->allocate(sizeof(RangeLock))) RangeLock(pool_);
    }

  msg_->raw_bcast((void*)&pool_, sizeof(void*));
  msg_->raw_bcast((void*)&rangelock_, sizeof(void*));

  if (debug_) {
      ExEnv::outn() << scprintf("%d: memory_ = 0x%x shmid_[0] = %d\n",
                       me(), memory_, shmid_[0]);
    }

  data_ = (void *) &((char*)memory_)[poolallocation];
}

void *
ShmMemoryGrp::localdata()
{
  return &((char*)data_)[distsize_to_size(localoffset())];
}

void
ShmMemoryGrp::cleanup()
{

  if (memory_) {
      for (int i=0; i<nregion_; i++) {
          shmdt((SHMTYPE)attach_address_[i]);
        }

      msg_->sync();

      if (me() == 0) {
          for (int i=0; i<nregion_; i++) {
              shmctl(shmid_[i],IPC_RMID,0);
            }
        }
      memory_ = 0;
    }

  delete[] update_;
  update_ = 0;

  nregion_ = 0;
  delete[] shmid_;
  delete[] attach_address_;
  shmid_ = 0;
  attach_address_ = 0;
}

ShmMemoryGrp::~ShmMemoryGrp()
{
  cleanup();

#ifdef DEBUG
  ExEnv::outn() << scprintf("msg_->nreference() = %d\n", msg_->nreference());
  ExEnv::outn().flush();
#endif // DEBUG
  msg_ = 0;
}

void *
ShmMemoryGrp::obtain_readwrite(distsize_t offset, int size)
{
  if (offset + size > totalsize()) {
      ExEnv::errn() << scprintf("ShmMemoryGrp::obtain_readwrite: arg out of range\n");
      abort();
    }

#if SIMPLE_LOCK
  obtain_lock();
#else // SIMPLE_LOCK
#ifdef DEBUG
  ExEnv::outn() << scprintf("%d: clear_release_count\n", me());
  ExEnv::outn().flush();
#endif // DEBUG
  clear_release_count();
#ifdef DEBUG
  ExEnv::outn() << scprintf("%d: obtain_lock\n", me());
  ExEnv::outn().flush();
#endif // DEBUG
  obtain_lock();
#ifdef DEBUG
  ExEnv::outn() << scprintf("%d: checkeq\n", me());
  ExEnv::outn().flush();
#endif
  while (!rangelock_->checkeq(offset, offset + size, 0)) {
#ifdef DEBUG
      ExEnv::outn() << scprintf("%d: range not zero -- waiting for release\n", me());
      ExEnv::outn().flush();
#endif // DEBUG
      //rangelock_->print();
      release_lock();
      wait_for_release();
      obtain_lock();
    }
  rangelock_->decrement(offset, offset + size);
#ifdef DEBUG
  ExEnv::outn() << scprintf("%d: after obtain write\n", me());
  ExEnv::outn().flush();
  //rangelock_->print();
#endif // DEBUG
  release_lock();
#endif // SIMPLE_LOCK

  return &((char*)data_)[distsize_to_size(offset)];
}

void *
ShmMemoryGrp::obtain_readonly(distsize_t offset, int size)
{
  if (offset + size > totalsize()) {
      ExEnv::errn() << scprintf("ShmMemoryGrp::obtain_readonly: arg out of range\n");
      abort();
    }

  return &((char*)data_)[distsize_to_size(offset)];
}

void *
ShmMemoryGrp::obtain_writeonly(distsize_t offset, int size)
{
  if (offset + size > totalsize()) {
      ExEnv::errn() << scprintf("ShmMemoryGrp::obtain_writeonly: arg out of range\n");
      abort();
    }

  return &((char*)data_)[distsize_to_size(offset)];
}

void
ShmMemoryGrp::release_readonly(void *data, distsize_t offset, int size)
{
}

void
ShmMemoryGrp::release_writeonly(void *data, distsize_t offset, int size)
{
}

void
ShmMemoryGrp::release_readwrite(void *data, distsize_t offset, int size)
{
#if SIMPLE_LOCK
  release_lock();
#else // SIMPLE_LOCK
  obtain_lock();
  rangelock_->increment(offset, offset + size);
  note_release();
#ifdef DEBUG
  ExEnv::outn() << scprintf("%d: after release write\n", me());
  //rangelock_->print();
  ExEnv::outn().flush();
#endif // DEBUG
  release_lock();
#endif // SIMPLE_LOCK
}

void
ShmMemoryGrp::obtain_lock()
{
#ifdef DEBUG
  ExEnv::outn() << scprintf("%d: lock val = %d\n", me(), lock_.val());
  ExEnv::outn().flush();
#endif // DEBUG
  lock_--;
#ifdef DEBUG
  ExEnv::outn() << scprintf("%d: lock decremented\n", me());
  ExEnv::outn().flush();
#endif // DEBUG
}

void
ShmMemoryGrp::release_lock()
{
  lock_++;
#ifdef DEBUG
  ExEnv::outn() << scprintf("%d: incremented lock\n", me());
  ExEnv::outn().flush();
#endif // DEBUG
}

void
ShmMemoryGrp::note_release()
{
  for (int i=0; i<n(); i++) {
      update_[i]++;
    }
#ifdef DEBUG
  ExEnv::outn() << scprintf("%d: incremented release flags\n", me());
  ExEnv::outn().flush();
#endif // DEBUG
}

void
ShmMemoryGrp::wait_for_release()
{
#ifdef DEBUG
  ExEnv::outn() << scprintf("%d: decrementing release flag\n", me());
  ExEnv::outn().flush();
#endif // DEBUG
  update_[me()]--;
#ifdef DEBUG
  ExEnv::outn() << scprintf("%d: decremented release flag\n", me());
  ExEnv::outn().flush();
#endif // DEBUG
}

void
ShmMemoryGrp::clear_release_count()
{
  update_[me()] = 0;
#ifdef DEBUG
  ExEnv::outn() << scprintf("%d: clearing release count\n", me());
  ExEnv::outn().flush();
#endif // DEBUG
}

void
ShmMemoryGrp::print(ostream &o) const
{
  MemoryGrp::print(o);
}

#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
