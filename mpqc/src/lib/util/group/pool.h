
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_pool_h
#define _util_group_pool_h

#include <stdlib.h>
#include <new.h>
#include <iostream.h>

#undef DEBUG_POOL

const int pool_data_alignment_bit = 3;
//const int pool_data_alignment_bit = 14;
const size_t pool_data_alignment = 1<<pool_data_alignment_bit;
inline size_t
align_pool_data(size_t size)
{
  return (size + pool_data_alignment - 1)
      & (~ (pool_data_alignment - 1));
}
inline void*
align_pool_data(void* ptr)
{
  return (void*)( (unsigned long) ((char*)ptr + pool_data_alignment - 1)
                 & (~ (pool_data_alignment - 1)));
}
inline size_t
align_pool_data_downward(size_t size)
{
  return size & (~ (pool_data_alignment - 1));
}
inline void*
align_pool_data_downward(void* ptr)
{
  return (void*) ( (unsigned long) ptr & (~ (pool_data_alignment - 1)));
}

//////////////////////////////////////////////////////////////////////////////

class PoolData;
struct FreeData {
    PoolData* next_free_;
    PoolData* prev_free_;
};

//////////////////////////////////////////////////////////////////////////////

struct UsedData {
    unsigned int flags;
    unsigned int held_:16;
    int priority_:15;
    unsigned int fixed_:1;
};

//////////////////////////////////////////////////////////////////////////////

class PoolData {
  public:
    enum {magic = 0x1f1d1e1c};
    int magic_;
    size_t size_;
    unsigned int free_:1;
    unsigned int flags_:15;
  private:
    PoolData* next_;
    PoolData* prev_;
  public:
    union {
        FreeData f;
        UsedData u;
    };

    // Allocates a chunk of free memory, only initializing the size.
    PoolData(size_t size);

    PoolData* next();
    PoolData* prev();

    void next(PoolData*);
    void prev(PoolData*);
    void prev_next(PoolData*,PoolData*);

    PoolData* next_free();
    PoolData* prev_free();

    void next_free(PoolData*);
    void prev_free(PoolData*);
    void prev_next_free(PoolData*,PoolData*);

    void set_magic(int = magic);

    // This new can only be called with aligned memory.
    //void* operator new(size_t size, void* placement);
    void* data();

    void check(void*lower=(void*)0x0,void*upper=(void*)0x7fffffffL);
};

const int PoolData_aligned_size = (sizeof(PoolData) + pool_data_alignment - 1)
    & (~ (pool_data_alignment - 1));
inline void* PoolData::data()
{
  return (void*)(((char*)this) + PoolData_aligned_size);
}

inline PoolData*
PoolData::next()
{
  return next_;
}

inline PoolData*
PoolData::prev()
{
  return prev_;
}

inline void
PoolData::next(PoolData*p)
{
  next_ = p;
#ifdef DEBUG_POOL
  if (next_ && prev_ && (next_ < prev_)) {
      fprintf(stderr,"PoolData::next(PoolData*): next < prev\n");
      abort();
    }
#endif
}

inline void
PoolData::prev(PoolData*p)
{
  prev_ = p;
#ifdef DEBUG_POOL
  if (next_ && prev_ && (next_ < prev_)) {
      fprintf(stderr,"PoolData::prev(PoolData*): next < prev\n");
      abort();
    }
#endif
}

inline void
PoolData::prev_next(PoolData*p,PoolData*n)
{
  prev_ = p;
  next_ = n;
#ifdef DEBUG_POOL
  if (next_ && prev_ && (next_ < prev_)) {
      fprintf(stderr,"PoolData::prev_next: next < prev\n");
      abort();
    }
#endif
}

//////

inline PoolData*
PoolData::next_free()
{
#ifdef DEBUG_POOL
  if (!free_) {
      fprintf(stderr,"PoolData::next_free(): datum is not free\n");
      abort();
    }
#endif
  return f.next_free_;
}

inline PoolData*
PoolData::prev_free()
{
#ifdef DEBUG_POOL
  if (!free_) {
      fprintf(stderr,"PoolData::prev_free(): datum is not free\n");
      abort();
    }
#endif
  return f.prev_free_;
}

inline void
PoolData::next_free(PoolData*p)
{
#ifdef DEBUG_POOL
  if (!free_) {
      fprintf(stderr,"PoolData::next_free(PoolData*): datum is not free\n");
      abort();
    }
#endif
  f.next_free_ = p;
}

inline void
PoolData::prev_free(PoolData*p)
{
#ifdef DEBUG_POOL
  if (!free_) {
      fprintf(stderr,"PoolData::prev_free(PoolData*): datum is not free\n");
      abort();
    }
#endif
  f.prev_free_ = p;
}

inline void
PoolData::prev_next_free(PoolData*p,PoolData*n)
{
#ifdef DEBUG_POOL
  if (!free_) {
      fprintf(stderr,"PoolData::prev_next_free: datum is not free\n");
      abort();
    }
#endif
  f.prev_free_ = p;
  f.next_free_ = n;
}

inline
PoolData::PoolData(size_t size):
  size_(size-PoolData_aligned_size),
  magic_(magic)
{
}

inline void
PoolData::set_magic(int magic)
{
  magic_ = magic;
}

//////////////////////////////////////////////////////////////////////////////

class Pool {
  protected:
    enum { freelist_size = sizeof(size_t)*8 };
    PoolData* freelist_[freelist_size];

    size_t size_;

    PoolData* firstdatum_;
    PoolData* voidptr_to_pd(void*d);

    int freelist_find_slot(size_t);
    void freelist_add(PoolData*);
    void freelist_del(PoolData*);
  public:
    Pool(size_t);
    ~Pool();
    
//     void* operator new(size_t size, void* placement) { return placement; }

//     Handle& allocate_handle(size_t size, int priority = 0);
//     void release(Handle&);

    void* allocate(size_t size);
    void release(void*d);
    double* allocate_double(size_t n);
    void release(double*d);
    int* allocate_int(size_t n);
    void release(int*d);
    void print(ostream&o=cout);
    void check();
};

inline PoolData*
Pool::voidptr_to_pd(void*d)
{
  return (PoolData*)((char*)d - PoolData_aligned_size);
}

inline double*
Pool::allocate_double(size_t n)
{
  return (double*) allocate(n*sizeof(double));
}

inline void
Pool::release(double*d)
{
  release((void*)d);
}
inline int*
Pool::allocate_int(size_t n)
{
  return (int*) allocate(n*sizeof(int));
}
inline void
Pool::release(int*d)
{
  release((void*)d);
}

#endif

