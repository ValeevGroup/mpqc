
#ifdef __GNUC__
#pragma implementation
#endif

#include <stdio.h>
#include "pool.h"

void
PoolData::check(void* lower_bound, void* upper_bound)
{
  PoolData*tmp = this;
  if ((void*)this < lower_bound || (void*)this >= upper_bound) {
      fprintf(stderr,"PoolData::check: this out of bounds\n");
      abort();
    }
  if (next_) {
      if ((void*)next_ < lower_bound || (void*)next_ >= upper_bound) {
          fprintf(stderr,"PoolData::check: next_ out of bounds\n");
          abort();
        }
      if (next_->prev_ != this) {
          fprintf(stderr,"PoolData::check: next pd doesn't point back\n");
          abort();
        }
      if ((void*)next_ != ((void*)this) + size_ + PoolData_aligned_size) {
          fprintf(stderr,"PoolData::check: next_ not consistent with size\n");
          abort();
        }
      if (free_ && next_->free_) {
          fprintf(stderr,"PoolData::check: free and next is free\n");
          abort();
        }
    }
  if (prev_) {
      if ((void*)prev_ < lower_bound || (void*)prev_ >= upper_bound) {
          fprintf(stderr,"PoolData::check: prev_ out of bounds\n");
          abort();
        }
      if (prev_->next_ != this) {
          fprintf(stderr,"PoolData::check: prev pd doesn't point back\n");
          abort();
        }
      if (free_ && prev_->free_) {
          fprintf(stderr,"PoolData::check: free and prev is free\n");
          abort();
        }
    }
  if (free_) {
      PoolData* n = f.next_free_;
      PoolData* p = f.prev_free_;
      if (n) {
          if ((void*)n < lower_bound || (void*)n >= upper_bound) {
              fprintf(stderr,"PoolData::check: next free out of bounds\n");
              abort();
            }
          if (n->f.prev_free_ != this) {
              fprintf(stderr,
                      "PoolData::check: next free pd doesn't point back\n");
              abort();
            }
        }
      if (p) {
          if ((void*)p < lower_bound || (void*)p >= upper_bound) {
              fprintf(stderr,"PoolData::check: prev free out of bounds\n");
              abort();
            }
          if (p->f.next_free_ != this) {
              fprintf(stderr,
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
  firstdatum_ = (PoolData*)align_pool_data(((void*)this) + sizeof(Pool));

  if (((void*)this) + size <= (void*) firstdatum_) {
      fprintf(stderr,"Pool::Pool: not given enough space\n");
      abort();
    }
  
  size_t firstdatum_size = align_pool_data_downward((size_t)
                                                    ((((void*)this)+size)
                                                    - (void*)firstdatum_));
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
                  PoolData* freechunk = (PoolData*)(((void*)j)
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
print_pooldata(FILE*fp,PoolData*d,int free)
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

  fprintf(fp,"    PoolData: size=%d", d->size_);
  if (d->free_) fprintf(fp,", free");
  else fprintf(fp,", allocated");
  if (!prev) fprintf(fp,", first");
  if (!next) fprintf(fp,", last");
  fprintf(fp,"\n");
  if (next) print_pooldata(fp,next,free);
}

void
Pool::print(FILE*fp)
{
  fprintf(fp,"Memory Pool:\n");
  fprintf(fp,"  data chain:\n");
  print_pooldata(fp,firstdatum_,0);
  for (int i=0; i<freelist_size; i++) {
      if (freelist_[i]) {
          fprintf(fp,"  freelist[%d]:\n",i);
          print_pooldata(fp,freelist_[i],1);
        }
    }
}

void
Pool::check()
{
  // The bit lost at the beginning to Pool and alignment.
  size_t start = (size_t) (align_pool_data(((void*)this) + sizeof(Pool))
                           - (void*)this);

  // The bit lost at the end to alignment.
  size_t end =  (size_t)(((void*)this) + size_
                         - align_pool_data_downward((size_t)((void*)this)
                                                    +size_));

  size_t size = start + end;

  for (PoolData*j=firstdatum_; j; j = j->next()) {
      j->check(this,((void*)this)+size_);
      size += PoolData_aligned_size + j->size_;
    }
  
  if (size != size_) {
      fprintf(stderr,"Pool::check(): inconsistent sizes\n");
      fprintf(stderr,"  computed: %d\n",size);
      fprintf(stderr,"  actual:   %d\n",size_);
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
          fprintf(stderr, "Pool::check: free data not in freelist\n");
          abort();
        }
    }
}
