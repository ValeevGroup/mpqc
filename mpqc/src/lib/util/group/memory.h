
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_memory_h
#define _util_group_memory_h

#include <stdio.h>
#include <util/class/class.h>
#include <util/container/ref.h>
#include <util/group/rnglock.h>

class MemoryGrp: public DescribedClass {
#define CLASSNAME MemoryGrp
#include <util/class/classda.h>
  private:
    RangeLock locks_;

  protected:
    // derived classes must fill in all these
    // ~MemoryGrp deletes the arrays
    int me_;
    int n_;
    int *offsets_; // offsets_[n_] is the fence for all data

    // release locks to the local memory
    void release_read_(int offset, int size);
    void release_write_(int offset, int size);

    // obtains locks to the local memory
    // return 1 if the lock was obtained, otherwise 0
    int obtain_read_(int offset, int size);
    int obtain_write_(int offset, int size);

    int use_locks_;

  public:
    MemoryGrp();
    virtual ~MemoryGrp();
    
    int me() { return me_; }
    int n() { return n_; }
    int localsize() { return offsets_[me_+1] - offsets_[me_]; }
    int localoffset() { return offsets_[me_]; }
    int size(int node) { return offsets_[node+1] - offsets_[node]; }
    int offset(int node) { return offsets_[node]; }
    int totalsize() { return offsets_[n_]; }

    // Activate is called before the memory is to be used and
    // deactivate is called afterwards.  (Necessary to avoid
    // problems with broken message passing routines.)
    virtual void activate();
    virtual void deactivate();

    // Locking of memory regions can be turned off if the application
    // can be sure that only one node has write access to a region
    // at time.  This might improve performance a bit.  lock must
    // be called with the same argument on all nodes.
    virtual void lock(int truefalse);

    // these are should only used by MemoryGrpBuf objects
    // since the arguments must be correct for these to work
    virtual void *obtain_writeonly(int offset, int size);
    virtual void *obtain_readwrite(int offset, int size) = 0;
    virtual void *obtain_readonly(int offset, int size) = 0;
    virtual void release_read(void *data, int offset, int size) = 0;
    virtual void release_write(void *data, int offset, int size) = 0;

    virtual void sum_reduction(double *data, int doffset, int dsize);
    virtual void sum_reduction_on_node(double *data, int doffset, int dsize,
                                       int node = -1);

    // Synchronize all the nodes.  Useful if lock(0) is called.
    virtual void sync() = 0;

    virtual void print(FILE *fp = stdout);

    static MemoryGrp *create_memorygrp(int localsize);
};
DescribedClass_REF_dec(MemoryGrp);

template <class data_t>
class MemoryGrpBuf {
    RefMemoryGrp grp_;
    enum LockType { None, Read, Write };
    LockType locktype_;
    data_t *data_;
    int offset_;
    int length_;
  public:
    MemoryGrpBuf(const RefMemoryGrp &);
    data_t *writeonly(int offset, int length);
    data_t *readwrite(int offset, int length);
    const data_t *readonly(int offset, int length);
    data_t *writeonly_on_node(int offset, int length, int node = -1);
    data_t *readwrite_on_node(int offset, int length, int node = -1);
    const data_t *readonly_on_node(int offset, int length, int node = -1);
    void release();
    int length() const { return length_; }
};

#ifndef __GNUG__
#  include <util/group/memory.cc>
#endif

#endif

