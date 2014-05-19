//
// ref.cc --- implementation of the reference counting classes
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

#include <typeinfo>

#include <util/ref/ref.h>
#include <util/misc/exenv.h>

using namespace std;

#if REF_USE_LOCKS

#if HAVE_STHREAD

#include <synch.h>
#include <thread.h>

typedef mutex_t sc_lock_t;

static unsigned int __base_nlock__;
static int __base_locker__;

static int
__init_lock__(sc_lock_t* inLock, int inCount)
{
    for (int i=0; i < inCount; i++)
        mutex_init(&inLock[i], USYNC_THREAD, 0);
    return 1;
}

#define __LOCK(l) mutex_lock(&l)
#define __UNLOCK(l) mutex_unlock(&l)

#elif HAVE_PTHREAD

#include <pthread.h>

typedef pthread_mutex_t sc_lock_t;

static int
__init_lock__(sc_lock_t* inLock, int inCount)
{
    for (int i=0; i < inCount; i++)
        pthread_mutex_init(&inLock[i], 0);
    return 1;
}

#define __LOCK(l) pthread_mutex_lock(&l)
#define __UNLOCK(l) pthread_mutex_unlock(&l)

#elif HAVE_CREATETHREAD

#include <windows.h>

typedef HANDLE sc_lock_t;

HANDLE __base_lock__ = 0;

static int
__init_lock__(sc_lock_t* inLock, int inCount)
{
    for (int i=0; i < inCount; i++)
        inLock[i] = CreateMutex(0, FALSE, 0);
    return 1;
}

// windows threads are recursive, so no fanciness is required
#define __LOCK(l) WaitForSingleObject(l, INFINITE)
#define __UNLOCK(l) ReleaseMutex(l)

#else /* !PTHREAD && !STHREAD  && !CREATETHREAD */

#define __LOCK(l) 0
#define __UNLOCK(l) 0

#endif /* HAVE_STHREAD */

/*
 * this is the number of locks to use in the round-robin.
 * since an unsigned char is used for the lock handle,
 * this cannot be greater than 255.
 */
#define NLOCKS 251  // 251 is the largest prime smaller than 255

static sc_lock_t sRefLocks[NLOCKS];
static int sRefLocksInit = __init_lock__(sRefLocks, NLOCKS);
static unsigned char sRefLock = 0;

#else /* !REF_USE_LOCKS */

#define __LOCK(l) 0
#define __UNLOCK(l) 0

#endif /* !REF_USE_LOCKS */

using namespace sc;

int
RefCount::lock_ptr() const
{
#if REF_USE_LOCKS
    if (ref_lock_ == 0xff)
        return 1;
    return __LOCK(sRefLocks[ref_lock_]);
#else
    return 1;
#endif    
}

int
RefCount::unlock_ptr() const
{
#if REF_USE_LOCKS
    if (ref_lock_ == 0xff)
        return 1;
    return __UNLOCK(sRefLocks[ref_lock_]);
#else
    return 1;
#endif    
}

void
RefCount::use_locks(bool inVal)
{
#if REF_USE_LOCKS
    if (inVal) {
        ref_lock_ = sRefLock;
        unsigned char tmp_sRefLock = sRefLock+1;
        if (tmp_sRefLock >= NLOCKS)
            tmp_sRefLock = 0;
        sRefLock = tmp_sRefLock;
    }
    else
        ref_lock_ = 0xff;
#endif
}


void
RefCount::error(const char * w) const
{
  ExEnv::errn() << "RefCount: ERROR: " << w << endl;
  ExEnv::errn() << "The type name is " << typeid(*this).name()
                << std::endl;
  abort();
}

void
RefCount::too_many_refs() const
{
  error("Too many refs.");
}

void
RefCount::not_enough_refs() const
{
  error("Ref count dropped below zero.");
}

RefCount::~RefCount()
{
#if REF_MANAGE
  if (managed() && nreference()) {
      error("Deleting a referenced object.");
    }
#endif
}

///////////////////////////////////////////////////////////////////////

void
RefBase::warn ( const char * msg) const
{
  ExEnv::errn() << "WARNING: " << msg << endl;
}
void
RefBase::warn_ref_to_stack() const
{
  warn("Ref: creating a reference to stack data");
}
void
RefBase::warn_skip_stack_delete() const
{
  warn("Ref: skipping delete of object on the stack");
}
void
RefBase::warn_bad_ref_count() const
{
  warn("Ref: bad reference count in referenced object\n");
}
void
RefBase::ref_info(RefCount*p, ostream& os) const
{
  if (p)
      os << "nreference() = " << p->nreference() << endl;
  else
      os << "reference is null" << endl;
}

void
RefBase::require_nonnull() const
{
  if (parentpointer() == 0) {
      ExEnv::errn() << "RefBase: needed a nonnull pointer but got null"
           << endl;
      abort();
    }
}

RefBase::~RefBase()
{
}

void
RefBase::check_pointer() const
{
  if (parentpointer() && parentpointer()->nreference() <= 0) {
      warn_bad_ref_count();
    }
}

void
RefBase::ref_info(ostream& os) const
{
  RefBase::ref_info(parentpointer(),os);
}

void
RefBase::reference(RefCount *p)
{
  if (p) {
#if REF_CHECK_STACK
      if (DO_REF_CHECK_STACK(p)) {
          DO_REF_UNMANAGE(p);
          warn_ref_to_stack();
        }
#endif
      p->reference();
    }
}

int
RefBase::dereference(RefCount *p)
{
  if (p) 
      return p->dereference();
  else
      return -1;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
