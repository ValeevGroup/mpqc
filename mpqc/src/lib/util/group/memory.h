
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_memory_h
#define _util_group_memory_h

#include <stdio.h>
#include <util/class/class.h>
#include <util/group/rnglock.h>

//. The \clsnm{MemoryGrp} abstract class provides the appearance of global
//. shared memory in a parallel machine.
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
    
    //. Returns who I am and how many nodes there are.
    int me() { return me_; }
    int n() { return n_; }

    //. Returns the amount of memory residing locally on \srccd{me()};
    int localsize() { return offsets_[me_+1] - offsets_[me_]; }
    //. Returns the global offset to this node's memory.
    int localoffset() { return offsets_[me_]; }
    //. Returns the amount of memory residing on \vrbl{node}.
    int size(int node) { return offsets_[node+1] - offsets_[node]; }
    //. Returns the global offset to \vrbl{node}'s memory.
    int offset(int node) { return offsets_[node]; }
    //. Returns the sum of all memory allocated on all nodes.
    int totalsize() { return offsets_[n_]; }

    //. Activate is called before the memory is to be used and
    //. deactivate is called afterwards.  (Necessary to avoid
    //. problems with broken message passing routines.)
    virtual void activate();
    virtual void deactivate();

    //. Locking of memory regions can be turned off if the application
    //. can be sure that only one node has write access to a region
    //. at time.  This might improve performance a bit.  lock must
    //. be called with the same argument on all nodes.
    virtual void lock(int truefalse);

    //. These are should only used by MemoryGrpBuf objects.
    virtual void *obtain_writeonly(int offset, int size);
    virtual void *obtain_readwrite(int offset, int size) = 0;
    virtual void *obtain_readonly(int offset, int size) = 0;
    virtual void release_read(void *data, int offset, int size) = 0;
    virtual void release_write(void *data, int offset, int size) = 0;

    virtual void sum_reduction(double *data, int doffset, int dsize);
    virtual void sum_reduction_on_node(double *data, int doffset, int dsize,
                                       int node = -1);

    //. Synchronizes all the nodes.  Consider using this when the way you
    //. you access memory changes.
    virtual void sync() = 0;

    //. Prints out information about the object.
    virtual void print(FILE *fp = stdout);

    //. Create a memory group with \vrbl{localsize} bytes on this
    //. node.  When memory is accessed using a global offset counting
    //. starts at node 0 and proceeds up to node \srccd{n()} - 1.
    //. This routine looks at the default \clsnmref{MessageGrp} object
    //. to decide which specialization of \clsnm{MemoryGrp} would be
    //. appropriate.
    static MemoryGrp *create_memorygrp(int localsize);
};
DescribedClass_REF_dec(MemoryGrp);

//. The \clsnm{MemoryGrpBug} class provides access to pieces of the
//. global shared memory that have been obtained with \clsnmref{MemoryGrp}.
//. \clsnm{MemoryGrpBug} is a template class that is parameterized on
//. \srccd{data\_t}.  All lengths and offsets of given in terms
//. of \srccd{sizeof(data\_t)}.
template <class data_t>
class MemoryGrpBuf {
    RefMemoryGrp grp_;
    enum LockType { None, Read, Write };
    LockType locktype_;
    data_t *data_;
    int offset_;
    int length_;
  public:
    //. Creates a new \clsnm{MemoryGrpBuf} given a \clsnmref{MemoryGrp}
    //. reference.  This is a template class parameterized on
    //. \srccd{data\_t}.
    MemoryGrpBuf(const RefMemoryGrp &);
    //. Request write only access to global memory at the global address
    //. \vrbl{offset} and with size \vrbl{length}.
    data_t *writeonly(int offset, int length);
    //. Request read write access to global memory at the global address
    //. \vrbl{offset} and with size \vrbl{length}.
    data_t *readwrite(int offset, int length);
    //. Request read only access to global memory at the global address
    //. \vrbl{offset} and with size \vrbl{length}.
    const data_t *readonly(int offset, int length);
    //. These behave like \srccd{writeonly}, \srccd{readwrite}, and
    //. \srccd{readonly}, except the \vrbl{offset} is local to the
    //. node specified by \vrbl{node}.  If \vrbl{node} = -1, then
    //. the local node is used.
    data_t *writeonly_on_node(int offset, int length, int node = -1);
    data_t *readwrite_on_node(int offset, int length, int node = -1);
    const data_t *readonly_on_node(int offset, int length, int node = -1);
    //. Release the access to the chunk of global memory that was
    //. obtained with \srccd{writeonly}, \srccd{readwrite},
    //. \srccd{readonly}, \srccd{writeonly\_on\_node},
    //. \srccd{readwrite\_on\_node}, and \srccd{readonly\_on\_node}.
    void release();
    //. The length of the current bit of memory.
    int length() const { return length_; }
};

#ifndef __GNUG__
#  include <util/group/memory.cc>
#endif

#endif

