//
// memammap.cc --- implementation of memory group which uses DEC mmap
//
// Copyright (C) 1998 Limit Point Systems, Inc.
//
// Author: Edward T. Seidl <seidl@janed.com>
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <util/misc/formio.h>
#include <util/group/pool.h>
#include <util/group/memammap.h>
#include <limits.h>
#include <errno.h>
#include <sys/param.h>
#include <unistd.h>
#include <sys/mman.h>

////////////////////////////////////////////////////////////////////////////

#define MMAPSIZE 1024

class MmapHeap {
  public:
    void *data;
    int msize;

    MmapHeap() {
      msize=MMAPSIZE;
      char *size = getenv("SCHEAPSIZE");
      if (size) 
        msize = atoi(size);
      
      data = mmap(0, msize, PROT_READ|PROT_WRITE, MAP_SHARED|MAP_ANONYMOUS,
                  -1, 0);
      if (data == (void*) -1) {
        cerr << "MmapHeap: mmap failed" << endl;
        abort();
      }
    }

    ~MmapHeap() { munmap(data, msize); }
};

static MmapHeap mmap_heap;

////////////////////////////////////////////////////////////////////////////

#ifndef SIMPLE_LOCK
#define SIMPLE_LOCK 1
#endif

#undef DEBUG

#define CLASSNAME AlphaMMapMemoryGrp
#define HAVE_KEYVAL_CTOR
#define PARENTS public MsgMemoryGrp
#include <util/class/classi.h>
void *
AlphaMMapMemoryGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  MsgMemoryGrp::_castdown(cd);
  return do_castdowns(casts,cd);
}

AlphaMMapMemoryGrp::AlphaMMapMemoryGrp(const RefMessageGrp& msg):
  MsgMemoryGrp(msg)
{
  update_ = 0;
  data_ = 0;
  memory_ = 0;
  pool_ = 0;
  rangelock_ = 0;
}

AlphaMMapMemoryGrp::AlphaMMapMemoryGrp(const RefKeyVal& keyval):
  MsgMemoryGrp(keyval)
{
  update_ = 0;
  data_ = 0;
  memory_ = 0;
  pool_ = 0;
  rangelock_ = 0;
}

void
AlphaMMapMemoryGrp::set_localsize(int localsize)
{
  int i;

  cleanup();

  MsgMemoryGrp::set_localsize(localsize);

  const int poolallocation = 160000;

  update_ = new GlobalCounter[n()];

  // allocate memory both the data and the Pool
  int size = poolallocation + totalsize();

  // get the shared memory segments
  if (me() == 0) {
    if (size > mmap_heap.msize) {
      cerr << "AlphaMMapMemoryGrp: request too large "
           << size << " > " << mmap_heap.msize
           << endl;
      abort();
    }
    if (mmap_heap.data == (void*) -1) {
      cerr << "AlphaMMapMemoryGrp: mmap_heap not initialized" << endl;
      abort();
    }
    memory_ = mmap_heap.data;
    // initialize the pool
    pool_ = new(memory_) Pool(poolallocation);
    rangelock_ = new(pool_->allocate(sizeof(RangeLock))) RangeLock(pool_);
  }

  if (me() == 0) {
    lock_.initialize();
    lock_ = 1;
    char * stringrep = lock_.stringrep();
    int length = strlen(stringrep) + 1;
    msg_->bcast(&length, 1);
    msg_->bcast(stringrep, length);
#ifdef DEBUG
    cout << scprintf("%d: sent initialize string of \"%s\" (%d)\n",
                     me(), stringrep, length);
    cout.flush();
#endif // DEBUG
    delete[] stringrep;
    for (i=0; i<n(); i++) {
      update_[i].initialize();
      char * stringrep = update_[i].stringrep();
      int length = strlen(stringrep) + 1;
      msg_->bcast(&length, 1);
      msg_->bcast(stringrep, length);
#ifdef DEBUG
      cout << scprintf("%d: sent initialize string of \"%s\" (%d) for %d\n",
                       me(), stringrep, length, i);
      cout.flush();
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
    cout << scprintf("%d: got initialize string of \"%s\" (%d)\n",
                     me(), stringrep, length);
    cout.flush();
#endif // DEBUG
    lock_.initialize(stringrep);
    delete[] stringrep;
    for (i=0; i<n(); i++) {
      msg_->bcast(&length, 1);
      stringrep = new char[length];
      msg_->bcast(stringrep, length);
#ifdef DEBUG
      cout << scprintf("%d: got initialize string of \"%s\" (%d) for %d\n",
                       me(), stringrep, length, i);
      cout.flush();
#endif // DEBUG
      update_[i].initialize(stringrep);
      delete[] stringrep;
    }
  }

  msg_->raw_bcast((void*)&memory_, sizeof(void*));
  msg_->raw_bcast((void*)&pool_, sizeof(void*));
  msg_->raw_bcast((void*)&rangelock_, sizeof(void*));

  if (me() != 0) {
    if (mmap_heap.data != memory_) {
      cerr << scprintf("%d: memory_ = 0x%x 0x%x\n",
                       me(), memory_, mmap_heap.data);
      abort();
    }
    memory_ = mmap_heap.data;
  }

  msg_->sync();

  if (debug_) {
    cout << scprintf("%d: memory_ = 0x%x %d\n",
                     me(), memory_, ((int*)memory_)[0]);
  }

  data_ = (void *) &((char*)memory_)[poolallocation];
}

void
AlphaMMapMemoryGrp::cleanup()
{
  if (memory_) {
    memory_ = 0;
  }

  delete[] update_;
  update_ = 0;
}

AlphaMMapMemoryGrp::~AlphaMMapMemoryGrp()
{
  cleanup();

#ifdef DEBUG
  cout << scprintf("msg_->nreference() = %d\n", msg_->nreference());
  cout.flush();
#endif // DEBUG
  msg_ = 0;
}

void *
AlphaMMapMemoryGrp::obtain_readwrite(int offset, int size)
{
  if (offset + size > totalsize()) {
    cerr << scprintf("AlphaMMapMemoryGrp::obtain_readwrite: arg out of range\n");
    abort();
  }

  if (use_locks_) {
#if SIMPLE_LOCK
    obtain_lock();
#else // SIMPLE_LOCK
#ifdef DEBUG
    cout << scprintf("%d: clear_release_count\n", me());
    cout.flush();
#endif // DEBUG
    clear_release_count();
#ifdef DEBUG
    cout << scprintf("%d: obtain_lock\n", me());
    cout.flush();
#endif // DEBUG
    obtain_lock();
#ifdef DEBUG
    cout << scprintf("%d: checkeq\n", me());
    cout.flush();
#endif
    while (!rangelock_->checkeq(offset, offset + size, 0)) {
#ifdef DEBUG
      cout << scprintf("%d: range not zero -- waiting for release\n", me());
      cout.flush();
#endif // DEBUG
      //rangelock_->print();
      release_lock();
      wait_for_release();
      obtain_lock();
    }
    rangelock_->decrement(offset, offset + size);
#ifdef DEBUG
    cout << scprintf("%d: after obtain write\n", me());
    cout.flush();
    //rangelock_->print();
#endif // DEBUG
    release_lock();
#endif // SIMPLE_LOCK
  }

  return &((char*)data_)[offset];
}

void *
AlphaMMapMemoryGrp::obtain_readonly(int offset, int size)
{
  if (offset + size > totalsize()) {
    cerr << "AlphaMMapMemoryGrp::obtain_readonly: arg out of range" << endl;
    abort();
  }

  if (use_locks_) {
#if SIMPLE_LOCK
    obtain_lock();
#else // SIMPLE_LOCK
    clear_release_count();
    obtain_lock();
    while (!rangelock_->checkgr(offset, offset + size, -1)) {
#ifdef DEBUG
      cout << scprintf("%d: range is -1 -- waiting for release\n", me());
      cout.flush();
      //rangelock_->print();
#endif // DEBUG
      release_lock();
      wait_for_release();
      obtain_lock();
    }
    rangelock_->increment(offset, offset + size);
#ifdef DEBUG
    cout << scprintf("%d: after obtain read\n", me());
    cout.flush();
    //rangelock_->print();
#endif // DEBUG
    release_lock();
#endif // SIMPLE_LOCK
  }

  return &((char*)data_)[offset];
}

void
AlphaMMapMemoryGrp::release_read(void *data, int offset, int size)
{
  if (use_locks_) {
#if SIMPLE_LOCK
    release_lock();
#else // SIMPLE_LOCK
    obtain_lock();
    rangelock_->decrement(offset, offset + size);
    note_release();
#ifdef DEBUG
    cout << scprintf("%d: after release read\n", me());
    //rangelock_->print();
    cout.flush();
#endif // DEBUG
    release_lock();
#endif // SIMPLE_LOCK
  }
}

void
AlphaMMapMemoryGrp::release_write(void *data, int offset, int size)
{
  if (use_locks_) {
#if SIMPLE_LOCK
    release_lock();
#else // SIMPLE_LOCK
    obtain_lock();
    rangelock_->increment(offset, offset + size);
    note_release();
#ifdef DEBUG
    cout << scprintf("%d: after release write\n", me());
    //rangelock_->print();
    cout.flush();
#endif // DEBUG
    release_lock();
#endif // SIMPLE_LOCK
  }
}

void
AlphaMMapMemoryGrp::obtain_lock()
{
#ifdef DEBUG
  cout << scprintf("%d: lock val = %d\n", me(), lock_.val());
  cout.flush();
#endif // DEBUG
  lock_--;
#ifdef DEBUG
  cout << scprintf("%d: lock decremented\n", me());
  cout.flush();
#endif // DEBUG
}

void
AlphaMMapMemoryGrp::release_lock()
{
  lock_++;
#ifdef DEBUG
  cout << scprintf("%d: incremented lock\n", me());
  cout.flush();
#endif // DEBUG
}

void
AlphaMMapMemoryGrp::note_release()
{
  for (int i=0; i<n(); i++) {
    update_[i]++;
  }
#ifdef DEBUG
  cout << scprintf("%d: incremented release flags\n", me());
  cout.flush();
#endif // DEBUG
}

void
AlphaMMapMemoryGrp::wait_for_release()
{
#ifdef DEBUG
  cout << scprintf("%d: decrementing release flag\n", me());
  cout.flush();
#endif // DEBUG
  update_[me()]--;
#ifdef DEBUG
  cout << scprintf("%d: decremented release flag\n", me());
  cout.flush();
#endif // DEBUG
}

void
AlphaMMapMemoryGrp::clear_release_count()
{
  update_[me()] = 0;
#ifdef DEBUG
  cout << scprintf("%d: clearing release count\n", me());
  cout.flush();
#endif // DEBUG
}

void
AlphaMMapMemoryGrp::print(ostream &o)
{
  MemoryGrp::print(o);
  if (me() == 0) {
    if (use_locks_) {
      obtain_lock();
      //rangelock_->print(fp);
      release_lock();
    }
  }
}

void
AlphaMMapMemoryGrp::sum_reduction(double *data, int doffset, int dlength)
{
  int offset = doffset * sizeof(double);
  int length = dlength * sizeof(double);

  if (offset + length > totalsize()) {
    cerr << scprintf("MemoryGrp::sum_reduction: arg out of range\n");
    abort();
  }

  double *source_data = (double*) obtain_readwrite(offset, length);

  for (int i=0; i<dlength; i++) {
    source_data[i] += data[i];
  }

  release_write((void*) source_data, offset, length);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
