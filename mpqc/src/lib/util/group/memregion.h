//
// memregion.h
//
// Copyright (C) 2008 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
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

#ifndef _util_group_memregion_h
#define _util_group_memregion_h

#include <list>
#include <vector>
#include <memory>
#include <sys/types.h>

#include <util/group/memory.h>

namespace sc {

/** The MemoryGrpRegion is a MemoryGrp proxy to a region of a MemoryGrp. It is assumed that the host MemoryGrp
    is not resized during lifetime of all MemoryGrpRegion objects bound to it.

    Each region has a static local offset and size. set_localsize() then cannot ask for more than the maximum size,
    but can be called multiple times.
  */
class MemoryGrpRegion: public MemoryGrp {
  public:
    /// MemoryGrpRegion is defined by the host MemoryGrp and position and size of the local segment within the local memory of the host
    MemoryGrpRegion(const Ref<MemoryGrp>& mem, size_t mem_start, size_t max_size);
    /// same as above, except mem_start is determined automatically
    MemoryGrpRegion(const Ref<MemoryGrp>& mem, size_t max_size);
    ~MemoryGrpRegion();

    void set_localsize(size_t);
    void* localdata();

    void* obtain_readwrite(distsize_t offset, int size);
    void* obtain_readonly(distsize_t offset, int size);
    void* obtain_writeonly(distsize_t offset, int size);
    void release_readonly(void* data, distsize_t offset, int size);
    void release_writeonly(void* data, distsize_t offset, int size);
    void release_readwrite(void* data, distsize_t offset, int size);

    void activate();
    void deactivate();
    void sync();
    void catchup();
    void* malloc_local(size_t nbytes);
    void free_local(void* data);

    /// reimplementation of MemoryGrp::clone(). always throws ProgrammingError -- this type of MemoryGrp is not clonable
    Ref<MemoryGrp> clone();
    void print(std::ostream &o=ExEnv::out0()) const;

  private:

    // init() should be called in every constructor to finish up construction duties
    void init();

    /// LocalRegion is an interval [start, start+size)
    class LocalRegion {
      public:
        LocalRegion(size_t start, size_t size) : start_(start), size_(size) {}
        size_t start() const { return start_; }
        size_t size() const { return size_; }
        bool operator==(const LocalRegion& other) const {
          return start() == other.start() && size() == other.size();
        }
        bool operator>(const LocalRegion& other) const {
          return start() > other.start();
        }
      private:
        size_t start_;
        size_t size_;
    };
    typedef std::list<LocalRegion> LocalRegions;

    Ref<MemoryGrp> host_;
    LocalRegion reserve_;              // reserved space
    std::vector<distsize_t> host_offsets_;  // start of regions on each node using the host MemoryGrp coordinates
    // converts offset from this object's coordinates to host coordinates
    distsize_t offset_to_host_offset(const distsize_t& offset) const;

    /// this maps MemoryGrp to the ordered list of its allocated regions, in order of increasing starts
    class HostToRegionsMap {
      public:
        HostToRegionsMap();
        ~HostToRegionsMap();
        void insert(const MemoryGrp* host, const LocalRegion& region);
        void erase(const MemoryGrp* host, const LocalRegion& region);
        /// suggests where to put a block of size sz. Does not check whether host has enough memory.
        size_t find_free(const MemoryGrp* host, size_t sz) const;
      private:
        const LocalRegions* regions(const MemoryGrp* host) const;
        typedef std::map<const MemoryGrp*,LocalRegions*> ImplType;
        ImplType impl_;
        Ref<ThreadLock> lock_;
    };
    static HostToRegionsMap map_;

};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
