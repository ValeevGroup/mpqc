//
// dim.h
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

#ifndef _math_scmat_dim_h
#define _math_scmat_dim_h

#include <util/keyval/keyval.h>
#include <util/state/state.h>

namespace sc {

class RefSCDimension;
/** SCBlockInfo contains blocking information for the SCDimension class.
    There are really two ways that it can contain blocking information.  In
    the first way, a vector of block offsets and block sizes is stored.
    The second method is only used by those specializations created by the
    BlockedSCMatrixKit class.  In this method the blocking information is
    stored as subdimensions of type SCDimension.  If both methods are used,
    they must be used consistently.  That is, the number, sizes, and order
    of the blocks must match the number, sizes, and order of the
    SCDimension objects.  */
class SCBlockInfo: public SavableState {
  protected:
    int n_;
    int nblocks_;
    int *start_;
    int *size_;
    RefSCDimension *subdims_;
    void init_start();
  public:
    /// Create a SCBlockInfo object.
    SCBlockInfo(int n, int nblocks = 0, const int *blocksizes = 0);
    SCBlockInfo(StateIn&);
    /** The KeyVal constructor.
        <dl>

        <dt><tt>sizes</tt><dd> This is a vector giving the size of each
        subblock.  There is no default.

        <dt><tt>subdims</tt><dd> If this vector is given there is must be
        entry for each entry in the sizes vector.  Each entry is an
        SCDimension object.  The default is to not store subdimension
        information.

        </dl> */
    SCBlockInfo(const Ref<KeyVal>& keyval);

    ~SCBlockInfo();
    void save_data_state(StateOut&);

    /// Return nonzero if this is equivalent to bi.
    int equiv(SCBlockInfo *bi);
    /// Return the total number of elements.
    int nelem() const { return n_; }
    /// Return the number of blocks.
    int nblock() const { return nblocks_; }
    /// Return the starting index for block i.
    int start(int i) const { return start_[i]; }
    /// Return the size of block i.
    int size(int i) const { return size_[i]; }
    ///  Return the last index + 1 for block i.
    int fence(int i) const { return start_[i] + size_[i]; }

    void elem_to_block(int i, int &block, int &offset);

    /// Retreive subdimension information.
    RefSCDimension subdim(int i);
    /** Set subdimension information.  The dimension dim and index i must
        be consistent with the nblocks and blocksizes information given to
        the constructor. */
    void set_subdim(int i, const RefSCDimension &dim);

    /// Print the object to the stream o.
    void print(std::ostream&o=ExEnv::out0()) const;
};


/** The SCDimension class is used to determine the size and blocking of
    matrices.  The blocking information is stored by an object of class
    SCBlockInfo.  */
class SCDimension: public SavableState {
  protected:
    std::string name_;
    int n_;
    Ref<SCBlockInfo> blocks_;
    SCDimension(const char* name = 0);
  public:
    /** Create a dimension with an optional name.  The name is a copy of
        the '0' terminated string name. */
    SCDimension(int n, const char* name = 0);
    SCDimension(const Ref<SCBlockInfo>&, const char *name = 0);
    SCDimension(int n, int nblocks, const int *blocksizes = 0,
                const char* name = 0);
    /** The KeyVal constructor.
        <dl>

        <dt><tt>n</tt><dd> This gives size of the dimension.  One of n or
        blocks is required.

        <dt><tt>blocks</tt><dd> The block information for the dimension can
        be given as a SCBlockInfo object.  One of n or blocks is required.

        </dl> */
    SCDimension(const Ref<KeyVal>&);
    SCDimension(StateIn&s);

    ~SCDimension();
    void save_data_state(StateOut&);

    /// Test to see if two dimensions are equivalent.
    int equiv(const SCDimension*) const;
    
    /// Return the dimension.
    int n() const { return n_; }
    /** Return the name of the dimension.  If no name was given to the
        constructor, then return 0. */
    const char* name() const { return name_.c_str(); }

    /// Return the blocking information for this dimension.
    const Ref<SCBlockInfo>& blocks() const { return blocks_; }

    /// Print information about this dimension to o.
    void print(std::ostream&o=ExEnv::out0()) const;
};

/** The RefSCDimension class is a smart pointer to an SCDimension
    specialization. */
class RefSCDimension: public Ref<SCDimension> {
    // standard overrides
  public:
    /** Initializes the dimension pointer to 0.  The
        reference must be initialized before it is used. */
    RefSCDimension();
    /// Make this and d refer to the same SCDimension.
    RefSCDimension(const RefSCDimension& d);
    /// Make this refer to d.
    RefSCDimension(SCDimension *d);

    ~RefSCDimension();
    /// Make this refer to d.
    RefSCDimension& operator=(SCDimension* d);

    RefSCDimension& operator<<(RefCount*);
    RefSCDimension& operator<<(const RefBase &);
    /// Make this and d refer to the same SCDimension.
    RefSCDimension& operator=(const RefSCDimension & d);

    // dimension specific functions
  public:
    /// Return the dimension.
    operator int() const;
    int n() const;

    void print(std::ostream&o=ExEnv::out0()) const;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
