//
// file.h
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

#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_file_h
#define _util_group_file_h

#include <iostream>

#include <scconfig.h>
#include <util/class/class.h>
#include <util/group/thread.h>
#include <util/group/memory.h>

namespace sc {

/** The FileGrp abstract class provides a way of accessing distributed
file in a parallel machine.  Several specializations are available.  For
one processor, ProcFileGrp provides a simple stub implementation.
Otherwise, the specializations that should work are MPIIOIFileGrp and MTMPIFileGrp.

If a FileGrp is not given to the program, then one will be automatically
chosen depending on which MessageGrp is used by default, the type of
machine on which the code was compiled, and what options were given at
configuration time.
*/

class FileGrp: public DescribedClass {
  private:
    int datafile_;
    char *filename_;

    Ref<ThreadLock> *locks_;
    int nlock_;
 
    void init_locks();


  protected:
    
    // derived classes must fill in all these
    // ~FileGrp deletes the arrays
    int me_;
    int n_;
    distsize_t *offsets_; // offsets_[n_] is the fence for all data

    // set to nonzero for debugging information
    int debug_;

    void obtain_local_lock(size_t start, size_t fence);
    void release_local_lock(size_t start, size_t fence);
  public:
    FileGrp();
    FileGrp(const Ref<KeyVal>&);
    virtual ~FileGrp();

    /// Opens the files
    void open();
    /// Closes the files
    void close();
    /// Sets the filename for the FileGrp
    void set_filename(char *name);
    /// Returns the filename for the FileGrp
    const char* get_filename() const { return datafile_; };
    
    /// Returns who I am.
    int me() const { return me_; }
    /// Returns how many nodes there are.
    int n() const { return n_; }

    /** Set the size of locally held data.
        When data is accessed using a global offset counting
        starts at node 0 and proceeds up to node n() - 1. */
    virtual void set_localsize(size_t) = 0;
    /// Returns the amount of data residing locally on me().
    size_t localsize() { return distsize_to_size(offsets_[me_+1]-offsets_[me_]); }
    /// Returns the global offset to this node's data.
    distsize_t localoffset() { return offsets_[me_]; }
    /// Returns the amount of data residing on node.
    int size(int node)
        { return distsize_to_size(offsets_[node+1] - offsets_[node]); }
    /// Returns the global offset to node's data.
    distsize_t offset(int node) { return offsets_[node]; }
    /// Returns the sum of all data allocated on all nodes.
    distsize_t totalsize() { return offsets_[n_]; }

    /// Activate is called before the data is to be used.
    virtual void activate();
    /// Deactivate is called after the data has been used.
    virtual void deactivate();

    /// This gives write access to the data location.  No locking is done.
    virtual void *obtain_writeonly(distsize_t offset, int size) = 0;
    /** Only one thread can have an unreleased obtain_readwrite at a time.
        The actual file region locked can be larger than that requested.
        If the file region is already locked this will block.  For this
        reason, data should be held as read/write for as short a time as
        possible. */
    virtual void *obtain_readwrite(distsize_t offset, int size) = 0;
    /// This gives read access to the file location.  No locking is done.
    virtual void *obtain_readonly(distsize_t offset, int size) = 0;
    /// This is called when read access is no longer needed.
    virtual void release_readonly(void *data, distsize_t offset, int size) = 0;
    /// This is called when write access is no longer needed.
    virtual void release_writeonly(void *data, distsize_t offset, int size)=0;
    /** This is called when read/write access is no longer needed.
        The data will be unlocked. */
    virtual void release_readwrite(void *data, distsize_t offset, int size)=0;

    virtual void sum_reduction(double *data, distsize_t doffset, int dsize);
    virtual void sum_reduction_on_node(double *data, size_t doffset, int dsize,
                                       int node = -1);

    /** Synchronizes all the nodes.  Consider using this when the way you
        you access data changes. */
    virtual void sync() = 0;

    /** Processes outstanding requests. Some file group implementations
        don't have access to real shared memory or even active messages.
        Instead, requests are processed whenever certain file group
        routines are called.  This can cause large latencies and buffer
        overflows.  If this is a problem, then the catchup member can be
        called to process all outstanding requests. */
    virtual void catchup();

    /// Prints out information about the object.
    virtual void print(std::ostream &o = ExEnv::out0()) const;

    /** Create a file group.  This routine looks for a -filegrp
        argument, then the environmental variable FILEGRP, and, finally,
        the default MessageGrp object to decide which specialization of
        FileGrp would be appropriate.  The argument to -integralgrp should
        be either string for a ParsedKeyVal constructor or a classname.
        The default ThreadGrp and MessageGrp objects should be initialized
        before this is called. */
    static FileGrp* initial_filegrp(int &argc, char** argv);
    static FileGrp* initial_filegrp();
    /** The default file group contains the primary file group to
        be used by an application. */
    static void set_default_filegrp(const Ref<FileGrp>&);
    /// Returns the default file group.
    static FileGrp* get_default_filegrp();
    /// Clones the given FileGrp. The new FileGrp may need to be initialized additionally.
    virtual FileGrp* clone() =0;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
