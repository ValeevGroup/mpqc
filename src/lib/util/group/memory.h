//
// memory.h
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

#ifndef _util_group_memory_h
#define _util_group_memory_h

#include <iostream>

#include <scconfig.h>
#include <util/class/class.h>
#include <util/group/thread.h>

namespace sc {

#if 0 // this can be used to catch accidental conversions to int
class distsize_t {
    friend size_t distsize_to_size(const distsize_t &a);
    friend distsize_t operator *(const int &a,const distsize_t &b);
    friend distsize_t operator +(const int &a,const distsize_t &b);
    friend distsize_t operator -(const int &a,const distsize_t &b);
    friend distsize_t operator /(const int &a,const distsize_t &b);
    friend distsize_t operator %(const int &a,const distsize_t &b);
    friend ostream& operator <<(ostream& o, const distsize_t &s);
  private:
    unsigned long long s;
  public:
    distsize_t(): s(999999999999999LL) {}
    distsize_t(int a): s(a) {}
    distsize_t(unsigned int a): s(a) {}
    distsize_t(unsigned long long a): s(a) {}
    distsize_t &operator =(const distsize_t &a)
        { s=a.s; return *this; }
    distsize_t &operator +=(const distsize_t &a)
        { s+=a.s; return *this; }
    distsize_t operator *(const distsize_t &a) const
        { return s*a.s; }
    distsize_t operator +(const distsize_t &a) const
        { return s+a.s; }
    distsize_t operator -(const distsize_t &a) const
        { return s-a.s; }
    distsize_t operator /(const distsize_t &a) const
        { return s/a.s; }
    distsize_t operator %(const distsize_t &a) const
        { return s%a.s; }
    bool operator <(const distsize_t &a) const
        { return s<a.s; }
    bool operator <=(const distsize_t &a) const
        { return s<=a.s; }
    bool operator >(const distsize_t &a) const
        { return s>a.s; }
    bool operator >=(const distsize_t &a) const
        { return s>=a.s; }
    bool operator ==(const distsize_t &a) const
        { return s==a.s; }
    distsize_t operator *(const int &a) const
        { return s*a; }
    distsize_t operator +(const int &a) const
        { return s+a; }
    distsize_t operator -(const int &a) const
        { return s-a; }
    distsize_t operator /(const int &a) const
        { return s/a; }
    distsize_t operator %(const int &a) const
        { return s%a; }
};
inline distsize_t operator *(const int &a,const distsize_t &b)
{ return a*b.s; }
inline distsize_t operator +(const int &a,const distsize_t &b)
{ return a+b.s; }
inline distsize_t operator -(const int &a,const distsize_t &b)
{ return a-b.s; }
inline distsize_t operator /(const int &a,const distsize_t &b)
{ return a/b.s; }
inline distsize_t operator %(const int &a,const distsize_t &b)
{ return a%b.s; }
inline ostream& operator <<(ostream& o, const distsize_t &s) { return o<<s.s; }
inline size_t distsize_to_size(const distsize_t &a) {return a.s;}
#elif defined(HAVE_LONG_LONG)
typedef unsigned long long distsize_t;
typedef long long distssize_t;
inline size_t distsize_to_size(const distsize_t &a) {return a;}
#else
typedef unsigned long distsize_t;
typedef long distssize_t;
inline size_t distsize_to_size(const distsize_t &a) {return a;}
#endif

/** The MemoryGrp abstract class provides a way of accessing distributed
memory in a parallel machine.  Several specializations are available.  For
one processor, ProcMemoryGrp provides a simple stub implementation.
Parallel specializations include ShmMemoryGrp, MTMPIMemoryGrp, and
ARMCIMemoryGrp.  The particular specializations that work depend highly on
the target hardware and software environment.

*/
class MemoryGrp: public DescribedClass {
  private:
    Ref<ThreadLock> *locks_;
    int nlock_;
 
    void init_locks();


  protected:
    // derived classes must fill in all these
    // ~MemoryGrp deletes the arrays
    int me_;
    int n_;
    distsize_t *offsets_; // offsets_[n_] is the fence for all data

    // set to nonzero for debugging information
    int debug_;

    void obtain_local_lock(size_t start, size_t fence);
    void release_local_lock(size_t start, size_t fence);
  public:
    MemoryGrp();
    MemoryGrp(const Ref<KeyVal>&);
    virtual ~MemoryGrp();
    
    /// Returns who I am.
    int me() const { return me_; }
    /// Returns how many nodes there are.
    int n() const { return n_; }

    /** Set the size of locally held memory.
        When memory is accessed using a global offset counting
        starts at node 0 and proceeds up to node n() - 1. */
    virtual void set_localsize(size_t) = 0;
    /// Returns the amount of memory residing locally on me().
    size_t localsize() { return distsize_to_size(offsets_[me_+1]-offsets_[me_]); }
    /// Returns a pointer to the local data.
    virtual void *localdata() = 0;
    /// Returns the global offset to this node's memory.
    distsize_t localoffset() { return offsets_[me_]; }
    /// Returns the amount of memory residing on node.
    int size(int node)
        { return distsize_to_size(offsets_[node+1] - offsets_[node]); }
    /// Returns the global offset to node's memory.
    distsize_t offset(int node) { return offsets_[node]; }
    /// Returns the sum of all memory allocated on all nodes.
    distsize_t totalsize() { return offsets_[n_]; }

    /// Activate is called before the memory is to be used.
    virtual void activate();
    /// Deactivate is called after the memory has been used.
    virtual void deactivate();

    /// This gives write access to the memory location.  No locking is done.
    virtual void *obtain_writeonly(distsize_t offset, int size) = 0;
    /** Only one thread can have an unreleased obtain_readwrite at a time.
        The actual memory region locked can be larger than that requested.
        If the memory region is already locked this will block.  For this
        reason, data should be held as read/write for as short a time as
        possible. */
    virtual void *obtain_readwrite(distsize_t offset, int size) = 0;
    /// This gives read access to the memory location.  No locking is done.
    virtual void *obtain_readonly(distsize_t offset, int size) = 0;
    /// This is called when read access is no longer needed.
    virtual void release_readonly(void *data, distsize_t offset, int size) = 0;
    /// This is called when write access is no longer needed.
    virtual void release_writeonly(void *data, distsize_t offset, int size)=0;
    /** This is called when read/write access is no longer needed.
        The memory will be unlocked. */
    virtual void release_readwrite(void *data, distsize_t offset, int size)=0;

    /** Perform a sum reduction on double data.
        @param data the contribution to sum into the global array.
        @param doffset the global offset in terms of doubles.
        @param dsize the size in terms of doubles.
    */
    virtual void sum_reduction(double *data, distsize_t doffset, int dsize);
    /** Perform a sum reduction on double data localized to a single node.
        @param data the contribution to sum into the global array.
        @param doffset the local offset on the node in terms of doubles.
        @param dsize the size in terms of doubles.
        @param node the node holding the data. The default is the local
        node.
    */
    virtual void sum_reduction_on_node(double *data, size_t doffset, int dsize,
                                       int node = -1);

    /** Synchronizes all the nodes.  This is useful after remote memory
        writes to be certain that all of the writes have completed and the
        data can be accessed locally, for example. */
    virtual void sync() = 0;

    /** Allocate data that will be accessed locally only.  Using this
        for data that will be used for global operations can improve
        efficiency.  Data allocated in this way must be freed with
        free_local_double.  */
    virtual void* malloc_local(size_t nbyte);
    /** Allocate double data that will be accessed locally only.
        \sa malloc_local
    */
    double* malloc_local_double(size_t ndouble);

    /** Free data that was allocated with malloc_local, and set the pointer to 0 */
    virtual void free_local(void * & data);
    /** Free data that was allocated with malloc_local_double, and set the pointer to 0  */
    void free_local_double(double * & data);

    /** Processes outstanding requests. Some memory group implementations
        don't have access to real shared memory or even active messages.
        Instead, requests are processed whenever certain memory group
        routines are called.  This can cause large latencies and buffer
        overflows.  If this is a problem, then the catchup member can be
        called to process all outstanding requests. */
    virtual void catchup();

    /** Returns a copy of this MemoryGrp specialization that provides an
        independent communication context. This will throw for
        specializations that do not reimplement it. */
    virtual Ref<MemoryGrp> clone(void);

    /// Prints out information about the object.
    virtual void print(std::ostream &o = ExEnv::out0()) const;

    /** Create a memory group.  This routine looks for a -memorygrp
        argument, and then the environmental variable MEMORYGRP to decide
        which specialization of MemoryGrp would be appropriate.  The
        argument to -memorygrp or the value of the environmental variable
        should be either string for a ParsedKeyVal constructor or a
        classname.  The default ThreadGrp and MessageGrp objects should be
        initialized before this is called. */
    static MemoryGrp* initial_memorygrp(int &argc, char** argv);
    /// Create a memory group. \sa initial_memorygrp(int&,char**)
    static MemoryGrp* initial_memorygrp();
    /** The default memory group contains the primary memory group to
        be used by an application. */
    static void set_default_memorygrp(const Ref<MemoryGrp>&);
    /** Returns the default memory group.  If the default memory
        group has not yet been set, then one is created.
        The particular specialization used is determined by
        configuration options and which specializations are being
        used for MessageGrp and ThreadGrp. */
    static MemoryGrp* get_default_memorygrp();
};


/** The MemoryGrpBuf class provides access to pieces of the
    global shared memory that have been obtained with MemoryGrp.
    MemoryGrpBuf is a template class that is parameterized on
    data_t.  All lengths and offsets of given in terms
    of sizeof(data_t). */
template <class data_t>
class MemoryGrpBuf {
    Ref<MemoryGrp> grp_;
    enum AccessType { None, Read, Write, ReadWrite };
    AccessType accesstype_;
    data_t *data_;
    distsize_t offset_;
    int length_;
  public:
    /** Creates a new MemoryGrpBuf given a MemoryGrp
        reference.  This is a template class parameterized on
        data_t. */
    MemoryGrpBuf(const Ref<MemoryGrp> &);
    /** Request write only access to global memory at the global address
        offset and with size length.  Writing the same bit of memory twice
        without an intervening sync of the MemoryGrp will have undefined
        results. */
    data_t *writeonly(distsize_t offset, int length);
    /** Request read write access to global memory at the global address
        offset and with size length.  This will lock the memory it uses
        until release is called unless locking has been turned off in the
        MemoryGrp object. */
    data_t *readwrite(distsize_t offset, int length);
    /** Request read only access to global memory at the global address
        offset and with size length.  Writing to the
        specified region without an intervening sync of the MemoryGrp
        will have undefined results. */
    const data_t *readonly(distsize_t offset, int length);
    /** This behaves like writeonly, except the
        offset is local to the node specified by node.  If node = -1, then
        the local node is used. */
    data_t *writeonly_on_node(size_t offset, int length, int node = -1);
    /** This behaves like readwrite, except the
        offset is local to the node specified by node.  If node = -1, then
        the local node is used. */
    data_t *readwrite_on_node(size_t offset, int length, int node = -1);
    /** This behaves like readonly, except the
        offset is local to the node specified by node.  If node = -1, then
        the local node is used. */
    const data_t *readonly_on_node(size_t offset, int length, int node = -1);
    /** Release the access to the chunk of global memory that was obtained
        with writeonly, readwrite, readonly, writeonly_on_node,
        readwrite_on_node, and readonly_on_node. */
    void release();
    /// The length of the current bit of memory.
    int length() const { return length_; }
};

//////////////////////////////////////////////////////////////////////
// MemoryGrpBuf members

template <class data_t>
MemoryGrpBuf<data_t>::MemoryGrpBuf(const Ref<MemoryGrp> & grp)
{
  grp_ = grp;
  accesstype_ = None;
}

template <class data_t>
data_t *
MemoryGrpBuf<data_t>::writeonly(distsize_t offset, int length)
{
  if (accesstype_ != None) release();
  data_ = (data_t *) grp_->obtain_writeonly(sizeof(data_t)*offset,
                                            sizeof(data_t)*length);
  offset_ = offset;
  length_ = length;
  accesstype_ = Write;
  return data_;
}

template <class data_t>
data_t *
MemoryGrpBuf<data_t>::readwrite(distsize_t offset, int length)
{
  if (accesstype_ != None) release();
  data_ = (data_t *) grp_->obtain_readwrite(sizeof(data_t)*offset,
                                            sizeof(data_t)*length);
  offset_ = offset;
  length_ = length;
  accesstype_ = ReadWrite;
  return data_;
}

template <class data_t>
const data_t *
MemoryGrpBuf<data_t>::readonly(distsize_t offset, int length)
{
  if (accesstype_ != None) release();
  data_ = (data_t *) grp_->obtain_readonly(sizeof(data_t)*offset,
                                           sizeof(data_t)*length);
  offset_ = offset;
  length_ = length;
  accesstype_ = Read;
  return data_;
}

template <class data_t>
data_t *
MemoryGrpBuf<data_t>::writeonly_on_node(size_t offset, int length, int node)
{
  if (node == -1) node = grp_->me();
  return writeonly(offset + grp_->offset(node)/sizeof(data_t), length);
}

template <class data_t>
data_t *
MemoryGrpBuf<data_t>::readwrite_on_node(size_t offset, int length, int node)
{
  if (node == -1) node = grp_->me();
  return readwrite(offset + grp_->offset(node)/sizeof(data_t), length);
}

template <class data_t>
const data_t *
MemoryGrpBuf<data_t>::readonly_on_node(size_t offset, int length, int node)
{
  if (node == -1) node = grp_->me();
  return readonly(offset + grp_->offset(node)/sizeof(data_t), length);
}

template <class data_t>
void
MemoryGrpBuf<data_t>::release()
{
  if (accesstype_ == Write)
      grp_->release_writeonly((data_t *)data_,
                              sizeof(data_t)*offset_, sizeof(data_t)*length_);
  if (accesstype_ == Read)
      grp_->release_readonly(data_, sizeof(data_t)*offset_,
                             sizeof(data_t)*length_);
  if (accesstype_ == ReadWrite)
      grp_->release_readwrite(data_, sizeof(data_t)*offset_,
                              sizeof(data_t)*length_);

  accesstype_ = None;
}

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
