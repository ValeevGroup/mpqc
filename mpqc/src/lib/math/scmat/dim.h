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

#ifdef __GNUC__
#pragma interface
#endif

#ifndef _math_scmat_dim_h
#define _math_scmat_dim_h

#include <util/keyval/keyval.h>
#include <util/state/state.h>

class RefSCDimension;
//. \clsnm{SCBlockInfo} contains blocking information for the
//\clsnmref{SCDimension} class.  There are really two ways that it can
//contain blocking information.  In the first way, a vector of block
//offsets and block sizes is stored.  The second method is only used by
//those specializations created by the \clsnm{BlockedSCMatrixKit} class.
//In this method the blocking information is stored as subdimensions of
//type \clsnmref{SCDimension}.  If both methods are used, they must be used
//consistently.  That is, the number, sizes, and order of the blocks must
//match the number, sizes, and order of the \clsnm{SCDimension} objects.
class SCBlockInfo: public SavableState {
#   define CLASSNAME SCBlockInfo
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    int n_;
    int nblocks_;
    int *start_;
    int *size_;
    RefSCDimension *subdims_;
    void init_start();
  public:
    //. Create a \clsnm{SCBlockInfo} object.
    SCBlockInfo(int n, int nblocks = 0, const int *blocksizes = 0);
    SCBlockInfo(StateIn&);
    SCBlockInfo(const RefKeyVal& keyval);

    ~SCBlockInfo();
    void save_data_state(StateOut&);

    //. Return nonzero if \vrbl{this} is equivalent to \vrbl{bi}.
    int equiv(SCBlockInfo *bi);
    //. Return the total number of elements.
    int nelem() const { return n_; }
    //. Return the number of blocks.
    int nblock() const { return nblocks_; }
    //. Return the starting index for block \vrbl{i}.
    int start(int i) const { return start_[i]; }
    //. Return the size of block \vrbl{i}.
    int size(int i) const { return size_[i]; }
    //. Return the last index $+ 1$ for block \vrbl{i}.
    int fence(int i) const { return start_[i] + size_[i]; }

    void elem_to_block(int i, int &block, int &offset);

    //. Retreive subdimension information.
    RefSCDimension subdim(int i);
    //. Set subdimension information.  The dimension \vrbl{dim} and
    // index \vrbl{i} must be consistent with the \vrbl{nblocks} and
    // \vrbl{blocksizes} information given to the constructor.
    void set_subdim(int i, const RefSCDimension &dim);

    //. Print the object to the stream \vrbl{o}.
    void print(ostream&o=cout);
};
SavableState_REF_dec(SCBlockInfo);

//. The \clsnm{SCDimension} class is used to determine the size and
// blocking of matrices.  The blocking information is stored by
// an object of class \clsnmref{SCBlockInfo}.
class SCDimension: public SavableState {
#   define CLASSNAME SCDimension
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    char *name_;
    int n_;
    RefSCBlockInfo blocks_;
    SCDimension(const char* name = 0);
  public:
    //. Create a dimension with an optional name.  The
    //. name is a copy of the \srccd{'0'} terminated string \vrbl{name}.
    SCDimension(int n, const char* name = 0);
    SCDimension(const RefSCBlockInfo&, const char *name = 0);
    SCDimension(int n, int nblocks, const int *blocksizes = 0,
                const char* name = 0);
    SCDimension(const RefKeyVal&);
    SCDimension(StateIn&s);

    ~SCDimension();
    void save_data_state(StateOut&);

    //. Test to see if two dimensions are equivalent.
    int equiv(const SCDimension*) const;
    
    //. Return the dimension.
    int n() const { return n_; }
    //. Return the name of the dimension.  If no name was given
    //. to the constructor, then return \srccd{0}.
    const char* name() const { return name_; }

    //. Return the blocking information for this dimension.
    RefSCBlockInfo blocks() { return blocks_; }

    //. Print information about this dimension to \vrbl{o}.
    void print(ostream&o=cout);
};

DCRef_declare(SCDimension);
SSRef_declare(SCDimension);

//. The \clsnm{RefSCDimension} class is a smart pointer to an
//. \clsnm{SCDimension} specialization.
class RefSCDimension: public SSRefSCDimension {
    // standard overrides
  public:
    //. Initializes the dimension pointer to \srccd{0}.  The
    //. reference must be initialized before it is used.
    RefSCDimension();
    //. Make this and \vrbl{d} refer to the same \clsnmref{SCDimension}.
    RefSCDimension(const RefSCDimension& d);
    //. Make this refer to \vrbl{d}.
    RefSCDimension(SCDimension *d);

    RefSCDimension(const DCRefBase&);
    ~RefSCDimension();
    //. Make this refer to \vrbl{d}.
    RefSCDimension& operator=(SCDimension* d);

    RefSCDimension& operator=(const DCRefBase & c);
    //. Make this and \vrbl{d} refer to the same \clsnmref{SCDimension}.
    RefSCDimension& operator=(const RefSCDimension & d);

    // dimension specific functions
  public:
    //. Return the dimension.
    operator int() const;
    int n() const;

    void print(ostream&o=cout);
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
