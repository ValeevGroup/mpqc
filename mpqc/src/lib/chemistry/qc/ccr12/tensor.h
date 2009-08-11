//
// tensor.h --- tensor class that contains an interface to MemoryGrp etc
//
// Copyright (C) 2008 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@qtp.ufl.edu>
// Maintainer: TS
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

#ifndef _chemistry_qc_ccr12_tensor_h
#define _chemistry_qc_ccr12_tensor_h

#include <string>
#include <vector>
#include <util/misc/compute.h>
#include <util/group/memory.h>
#include <util/group/memregion.h>
#include <util/group/message.h>
#include <util/group/thread.h>
#include <chemistry/qc/mbptr12/distarray4.h>

namespace sc {

class Tensor : virtual public RefCount {
  protected:
    const Ref<MemoryGrp>& mem_;

    /// for use in disk-based algorithm
    std::string filename_;

    /// hash table: hash_table_.insert(key,offset)
    std::map<long,long> hash_table_;

    /// tensor size in double
    long filesize_;

    /// data area
    MemoryGrpRegion* file_;
    bool file_allocated_;

    /// determines the distribution of blocks to nodes
    std::vector<long> determine_filesizes();

  public:
    Tensor(std::string filename,const Ref<MemoryGrp>& mem);
    ~Tensor();

    /// returns filename_
    std::string filename() const {return filename_;};

    /// returns MemoryGrpRegion for this tensor
    MemoryGrpRegion* file() const {return file_;};

    /// set/get the filesize of the tensor
    void set_filesize(long i);
    long get_filesize() const {return filesize_;};

    /// create/delete the distributed file for the tensor
    void createfile();
    void deletefile();

    /// input for the hash table
    void input_offset(long tag, long offset);

    /// get a block from the distributed file (non-blocking)
    void get_block(long tag, double* data);
    /// add a block to the distributed file (non-blocking); double* data will be destroyed
    void add_block(long tag, double* data);
    /// put a block to the distributed file (non-blocking)
    void put_block(long tag, double* data);
    /// does this block exist?
    bool exists(long tag) const;
    /// returns if a block associated with long tag resides in a local memory or not
    bool is_this_local(long tag);

    /// initialize/clear tensors to zero
    void zero();

    /// assigns all values to a
    void assign(double a);

    /// obtain the norm of a tensor
    double norm() const;

    /// obtain the ddot of two tensors
    double ddot(Ref<Tensor>&) const;

    /// return a copy of self
    Ref<Tensor> copy() const;

    /// return a tensor that have the same hashtable, and is intilaized to zero
    Ref<Tensor> clone() const;

    /// perform daxpy of tensors (self+=a*other)
    void daxpy(const Ref<Tensor>&, double);

    /// scale self by a
    void scale(double a);

    /// sync
    void sync() const { const_cast<MemoryGrpRegion*>(file_)->sync();};

    /// print
    void print(const std::string& label, std::ostream& os = ExEnv::out0()) const;
};

/** Computes the ``RMS norm'' of the tensor, defined as tensor->norm() divided by the size of the tensor.
 *  The advantage of RMS norm over norm is that it's size-independent.
 */
double RMS(const Tensor& t);

}

#endif

