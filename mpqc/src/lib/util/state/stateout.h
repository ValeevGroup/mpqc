//
// stateout.h
//
// Copyright (C) 1998 Limit Point Systems, Inc.
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

#ifndef _util_state_stateout_h
#define _util_state_stateout_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/class/class.h>
#include <util/state/state.h>
#include <util/container/avlmap.h>

class StateOutData {
  public:
    int num;
    int size;
    int type;
    int offset;

    StateOutData(): num(0), size(0), type(0), offset(0) {}
};

//. The \clsnm{StateOut} class serializes objects of types
//. that derive from \clsnmref{SavableState}.  It keeps track
//. of pointers to data so that two references to the same
//. piece of data do not result in that data being sent to the
//. output device two times.
class StateOut: public DescribedClass {
#   define CLASSNAME StateOut
#   include <util/class/classda.h>
    friend class SavableState;
    friend class TranslateDataOut;
  private:
    // do not allow copy constructor or assignment
    StateOut(const StateOut&);
    void operator=(const StateOut&);
    int have_cd_;
  protected:
    int dir_loc_loc_;
    TranslateDataOut *translate_;
    int copy_references_;
    int next_object_number_;
    AVLMap<RefSavableState,StateOutData> ps_;
    AVLMap<ClassDescP,int> classidmap_;
    int nextclassid_;
    int node_to_node_;
    virtual int put_array_void(const void*,int);
    virtual int putparents(const ClassDesc*);

    void put_directory();

    // The following members are called by friend SavableState

    void have_classdesc() { have_cd_ = 1; }
    int need_classdesc() { int tmp = have_cd_; have_cd_ = 0; return !tmp; }

    //. This will prepare \clsnm{StateOut} to output a pointer to
    //. data.  It first checks to see if the data has already been
    //. saved.  If it has, then a reference to this data is saved.
    //. Otherwise the object is written out.
    virtual int putobject(const RefSavableState &);

    //. Write out information about the given \clsnmref{ClassDesc}.
    virtual int put(const ClassDesc*);
  public:
    StateOut();
    virtual ~StateOut();

    //. Write out header information.
    virtual void put_header();

    //. This is like \srccd{put} except the length of the \srccd{char}
    //. array is determined by interpreting the character array as
    //. a character string.
    virtual int putstring(const char*);

    //. Write the given datum.
    virtual int put(char r);
    virtual int put(unsigned int r);
    virtual int put(int r);
    virtual int put(float r);
    virtual int put(double r);
    //. Write the given array data.  Size information is also saved.  The
    //data is allocated and read by the \srccd{get(T*&)} routines.
    virtual int put(const char*,int);
    virtual int put(const unsigned int*,int);
    virtual int put(const int*,int);
    virtual int put(const float*,int);
    virtual int put(const double*,int);
    //. Put arrays of data.  No size information is stored.  This
    //data is read by the \srccd{get\_array\_T} routines.
    virtual int put_array_char(const char*p,int size);
    virtual int put_array_uint(const unsigned int*p,int size);
    virtual int put_array_int(const int*p,int size);
    virtual int put_array_float(const float*p,int size);
    virtual int put_array_double(const double*p,int size);

    //. Don't keep track of pointers to objects.  Calling this
    //. causes duplicated references to objects to be copied.
    //. The directory will not contain the forgotten objects.
    void forget_references();
    //. If a reference to an object that has already been written
    //. is encountered, copy it instead of generating a reference
    //. to the first object.
    //. The directory will not be updated with new objects.
    void copy_references();

    //. Returns true if this object uses a directory.
    virtual int use_directory();

    //. Flush out any remaining data.
    virtual void flush();

    //. True if this is a node to node save/restore.  This is
    //necessary for classes that try to avoid saving databases
    //to files that can otherwise be read in, but want to avoid
    //reading the database from disk on all nodes.
    int node_to_node() const { return node_to_node_; }

    //. Returns the current position in the file.  The default
    //implementation returns 0.
    virtual int tell();
    //. Set the current position in the file.  The default implementation
    //does nothing.
    virtual void seek(int loc);
    //. Return non-zero if tell and seek do anything sensible.  The
    //default implementation returns 0.
    virtual int seekable();
  };
DescribedClass_REF_dec(StateOut);

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
