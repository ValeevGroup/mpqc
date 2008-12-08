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
#include <sys/types.h>

#include <util/group/memory.h>

namespace sc {

/** The MemoryGrpRegion is a MemoryGrp proxy to a region of a MemoryGrp. It is assumed that the host MemoryGrp
    is not resized.
    
    Each region has a static start and size.
    set_localsize() then cannot ask for more than the maximum size.
  */
class MemoryGrpRegion: public MemoryGrp {
  public:
    /// MemoryGrpRegion is defined by the host MemoryGrp and position and size of the local segment within the local memory of the host
    MemoryGrpRegion(const Ref<MemoryGrp>& mem, size_t mem_start, size_t max_size);
    ~MemoryGrpRegion();

    void set_localsize(size_t);
    void* localdata();

    void* obtain_readwrite(distsize_t offset, int size);
    void* obtain_readonly(distsize_t offset, int size);
    void* obtain_writeonly(distsize_t offset, int size);
    void release_readonly(void* data, distsize_t offset, int size);
    void release_writeonly(void* data, distsize_t offset, int size);
    void release_readwrite(void* data, distsize_t offset, int size);

    void sync();
    
  private:
    class Region {
      public:
        Region(size_t start, size_t size) : start_(start), size_(size) {}
        size_t start() const { return start_; }
        size_t size() const { return size_; }
      private:
        size_t start_;
        size_t size_;
    };
    typedef std::list<Region> Regions;
    
    Ref<MemoryGrp> host_;
    Region reserve_;              // reserved space
    distsize_t* memgrp_offsets_;  // offsets within the host MemoryGrp
    
    // this maps MemoryGrp to the ordered list of its allocated regions, in order of increasing starts
    class HostToRegionsMap {
      public:
        const Regions* regions(const MemoryGrp* host) const;
        void insert(const MemoryGrp* host, const Region& region);
      private:
        std::map<MemoryGrp*,Regions*> impl_;
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
