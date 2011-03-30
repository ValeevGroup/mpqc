//
// pool.cc
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

#include <util/misc/formio.h>
#include <util/group/pool.h>

using namespace std;
using namespace sc;

void
PoolData::check(void* lower_bound, void* upper_bound)
{
  if ((void*)this < lower_bound || (void*)this >= upper_bound) {
      ExEnv::errn() << scprintf("PoolData::check: this out of bounds\n");
      abort();
    }
  if (next_) {
      if ((void*)next_ < lower_bound || (void*)next_ >= upper_bound) {
          ExEnv::errn() << scprintf("PoolData::check: next_ out of bounds\n");
          abort();
        }
      if (next_->prev_ != this) {
          ExEnv::errn() << scprintf("PoolData::check: next pd doesn't point back\n");
          abort();
        }
      if ((char*)next_ != (char*)this + size_ + PoolData_aligned_size) {
          ExEnv::errn() << scprintf("PoolData::check: next_ not consistent with size\n");
          abort();
        }
      if (free_ && next_->free_) {
          ExEnv::errn() << scprintf("PoolData::check: free and next is free\n");
          abort();
        }
    }
  if (prev_) {
      if ((void*)prev_ < lower_bound || (void*)prev_ >= upper_bound) {
          ExEnv::errn() << scprintf("PoolData::check: prev_ out of bounds\n");
          abort();
        }
      if (prev_->next_ != this) {
          ExEnv::errn() << scprintf("PoolData::check: prev pd doesn't point back\n");
          abort();
        }
      if (free_ && prev_->free_) {
          ExEnv::errn() << scprintf("PoolData::check: free and prev is free\n");
          abort();
        }
    }
  if (free_) {
      PoolData* n = f.next_free_;
      PoolData* p = f.prev_free_;
      if (n) {
          if ((void*)n < lower_bound || (void*)n >= upper_bound) {
              ExEnv::errn() << scprintf("PoolData::check: next free out of bounds\n");
              abort();
            }
          if (n->f.prev_free_ != this) {
              ExEnv::errn() << scprintf(
                      "PoolData::check: next free pd doesn't point back\n");
              abort();
            }
        }
      if (p) {
          if ((void*)p < lower_bound || (void*)p >= upper_bound) {
              ExEnv::errn() << scprintf("PoolData::check: prev free out of bounds\n");
              abort();
            }
          if (p->f.next_free_ != this) {
              ExEnv::errn() << scprintf(
                      "PoolData::check: prev free pd doesn't point back\n");
              abort();
            }
        }
    }
}

Pool::Pool(size_t size):
  size_(size)
{

  // Initialize the first and last members of the data list.
  firstdatum_ = (PoolData*)align_pool_data((void*)((char*)this+sizeof(Pool)));

  if ((char*)this + size <= (char*) firstdatum_) {
      ExEnv::errn() << scprintf("Pool::Pool: not given enough space\n");
      abort();
    }
  
  size_t firstdatum_size = align_pool_data_downward((size_t)
                                                    (((char*)this+size)
                                                    - (char*)firstdatum_));
  new(firstdatum_) PoolData(firstdatum_size);

  firstdatum_->prev_next(0,0);

  // Initialize the free lists.
  int i;
  for (i=0; i<freelist_size; i++) freelist_[i] = 0;
  freelist_add(firstdatum_);
}

void
Pool::freelist_add(PoolData*d)
{
  int slot = freelist_find_slot(d->size_);
  d->free_ = 1;
  PoolData* tmp = freelist_[slot];
  d->next_free(tmp);
  d->prev_free(0);
  freelist_[slot] = d;
  if (tmp) tmp->prev_free(d);
#ifdef DEBUG_POOL
  d->check();
  if (d->next()) d->next()->check();
  if (d->prev()) d->prev()->check();
#endif
}

void
Pool::freelist_del(PoolData*d)
{
  if (d->next_free()) d->next_free()->prev_free(d->prev_free());
  if (d->prev_free()) d->prev_free()->next_free(d->next_free());
  else {
      int slot = freelist_find_slot(d->size_);
      freelist_[slot] = d->next_free();
    }
  d->free_ = 0;
#ifdef DEBUG_POOL
  d->check();
  if (d->next()) d->next()->check();
  if (d->prev()) d->prev()->check();
#endif
}

int
Pool::freelist_find_slot(size_t size)
{
  int slot = 0;
  size_t mask = ~ (size_t)0;
  while(mask & size) {
      slot++;
      mask <<= 1;
    }
  return slot;
}

void*
Pool::allocate(size_t size)
{
  int slot = freelist_find_slot(size);
  for (int i=slot; i<freelist_size; i++) {
      PoolData* j;
      for (j=freelist_[i]; j; j = j->next_free()) {
          if (j->size_ >= size) {
              freelist_del(j);
              // Maybe need to break this chunk into two pieces.
              if (j->size_ > size + PoolData_aligned_size) {
                  PoolData* freechunk = (PoolData*)((char*)j
                                                    + PoolData_aligned_size
                                                    + size);
                  new(freechunk) PoolData(j->size_ - size);
                  freechunk->prev_next(j,j->next());
                  if (freechunk->next()) freechunk->next()->prev(freechunk);
                  j->size_ = size;
                  j->next(freechunk);
                  freelist_add(freechunk);
                }
#ifdef DEBUG_POOL
              j->check();
              if (j->next()) j->next()->check();
              if (j->prev()) j->prev()->check();
#endif
              return j->data();
            }
        }
    }
  return 0;
}

void
Pool::release(void* v)
{
  PoolData *d = voidptr_to_pd(v);
  if (d->prev() && d->prev()->free_) {
      freelist_del(d->prev());
      d->prev()->size_ += d->size_ + PoolData_aligned_size;
      d->prev()->next(d->next());
      if (d->next()) d->next()->prev(d->prev());
      d->set_magic(0);
      d = d->prev();
    }
  if (d->next() && d->next()->free_) {
      freelist_del(d->next());
      d->next()->set_magic(0);
      d->size_ += d->next()->size_ + PoolData_aligned_size;
      if (d->next()->next()) d->next()->next()->prev(d);
      d->next(d->next()->next());
    }
  freelist_add(d);
}

static void
print_pooldata(ostream&o,PoolData*d,int free)
{
  PoolData *next,*prev;
  if (free) {
      next = d->next_free();
      prev = d->prev_free();
    }
  else {
      next = d->next();
      prev = d->prev();
    }

  o << scprintf("    PoolData: size=%d", d->size_);
  if (d->free_) o << scprintf(", free");
  else o << scprintf(", allocated");
  if (!prev) o << scprintf(", first");
  if (!next) o << scprintf(", last");
  o << endl;
  if (next) print_pooldata(o,next,free);
}

void
Pool::print(ostream&o)
{
  o << scprintf("Memory Pool:\n");
  o << scprintf("  data chain:\n");
  print_pooldata(o,firstdatum_,0);
  for (int i=0; i<freelist_size; i++) {
      if (freelist_[i]) {
          o << scprintf("  freelist[%d]:\n",i);
          print_pooldata(o,freelist_[i],1);
        }
    }
}

void
Pool::check()
{
  // The bit lost at the beginning to Pool and alignment.
  size_t start = (size_t)
                 ((char*)align_pool_data((void*)((char*)this + sizeof(Pool)))
                  - (char*)this);

  // The bit lost at the end to alignment.
  size_t end = (size_t)
               ((char*)this + size_
                - (char*) align_pool_data_downward((size_t)((void*)this)
                                                   +size_));

  size_t size = start + end;

  PoolData *j;
  for (j=firstdatum_; j; j = j->next()) {
      j->check(this,(void*)((char*)this+size_));
      size += PoolData_aligned_size + j->size_;
    }
  
  if (size != size_) {
      ExEnv::errn() << scprintf("Pool::check(): inconsistent sizes\n");
      ExEnv::errn() << scprintf("  computed: %d\n",size);
      ExEnv::errn() << scprintf("  actual:   %d\n",size_);
      abort();
    }

  // make sure that all data is accounted for
  for (j=firstdatum_; j; j = j->next()) {
      j->flags_ = 1;
    }
  for (int i=0; i<freelist_size; i++) {
      for (j=freelist_[i]; j; j = j->next_free()) {
          j->flags_ = 0;
        }
    }
  for (j=firstdatum_; j; j = j->next()) {
      if (j->free_ && j->flags_) {
          ExEnv::errn() << scprintf("Pool::check: free data not in freelist\n");
          abort();
        }
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
