//
// state.h
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

#ifndef _util_state_state_h
#define _util_state_state_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/class/class.h>
#include <util/container/array.h>

#define SavableState_REF_dec(T) SavableState_named_REF_dec(Ref ## T,T)
#define SavableState_REF_def(T) SavableState_named_REF_def(Ref ## T,T)

#ifdef USE_REF_MACROS
#define SavableState_named_REF_dec(refname,T)				      \
   DCRef_declare(T); SSRef_declare(T); typedef class SSRef ## T refname;
#else
#define SSRef_declare(T) typedef class SSRef<T> SSRef ## T;
#define SavableState_named_REF_dec(refname,T) typedef class SSRef<T> refname;
#endif
#define SavableState_named_REF_def(refname,T)

// This does forward declarations of REF classes.
#ifdef USE_REF_MACROS
#define SavableState_REF_fwddec(T) class SSRef ## T; \
                                   typedef class SSRef ## T Ref ## T;
#else
#define SavableState_REF_fwddec(T) class T; \
                                   typedef class SSRef<T> Ref ## T;
#endif

class StateIn;
class StateOut;
class TranslateDataIn;
class TranslateDataOut;

// If multiple inheritance is used, SavableState should be a virtual
// parent.  This breaks some compilers so the virtual_base macro should
// be used instead of the virtual keyword (see identity.h).  Since for StateIn
// CTORs the SavableState base class must be initialized with the
// StateIn object, the following macro should be used to determine if
// nondirect descendants of SavableState need to call the SavableState
// CTOR.  s is the StateIn
#ifdef NO_VIRTUAL_BASES
#  define maybe_SavableState(s)
#else
#  define maybe_SavableState(s) ,SavableState(s)
#endif

//. \clsnm{SavableState} give objects of derivative classes the
//. ability to save and restore their state or to send their state.
class SavableState: public DescribedClass {
#   define CLASSNAME SavableState
#   include <util/class/classda.h>
  protected:
    SavableState();
    SavableState(const SavableState&);
    // helper for save_object_state overrides
    void save_object_state_(StateOut&, const ClassDesc *);
#ifndef __GNUC__
  public:
#endif
    SavableState& operator=(const SavableState&);
  public:
    virtual ~SavableState();

    //// save functions

    //. Save the state of the object as specified by the
    //. \clsnmref{StateOut} object.  This routine saves the state
    //. of the object (which includes the nonvirtual bases),
    //. the virtual bases, and type information.  The default
    //. implementation should be adequate.
    void save_state(StateOut&);

    // Like save_state(StateOut&), but will handle null pointers correctly.
    static void save_state(SavableState*s, StateOut&);

    //. This can be used for saving state when the exact type of
    //. the object is known for both the save and the restore.  To
    //. restore objects saved in this way the user must directly
    //. invoke the object's \srccd{\clsnmref{StateIn}\&} constructor.
    void save_object_state(StateOut&);

    //. Save the virtual bases for the object.
    //. This must be done in the same order that the ctor
    //. initializes the virtual bases.  This does not include
    //. the DescribedClass and SavableState virtual base classes.
    //. This must be implemented by the user if the class has other
    //. virtual bases.  (These virtual bases must come after
    //. \clsnm{SavableState}, if \clsnm{SavableState} is virtual.)
    virtual void save_vbase_state(StateOut&);

    //. Save the base classes (with \srccd{save\_data\_state})
    //. and the members
    //. in the same order that the \clsnmref{StateIn} CTOR
    //. initializes them.  This must be implemented by the derived
    //. class if the class has data.
    virtual void save_data_state(StateOut&);

    //// restore functions

    //. This restores objects saved with \srccd{save\_state}.  The
    //. exact type of the next object in \srccd{si} can be any
    //. type publically derived from the \clsnm{SavableState}.
    //. Derived classes implement a similar static function that
    //. returns a pointer to the derived class.  If the name is
    //. given the directory will be consulted to find and restore
    //. that object.
    static SavableState* restore_state(StateIn& si, const char *name = 0);

  protected:

    //. Each derived class \clsnmref{StateIn} CTOR handles the
    //. restore corresponding to calling \srccd{save\_object\_state},
    //. \srccd{save\_vbase\_state}, and \srccd{save\_data\_state}
    //. listed above.  All derived class
    //. \srccd{\clsnmref{StateIn}\&} constructors must invoke the
    //. \srccd{\clsnm{SavableState}(\srccd{StateIn}\&)}
    //. constructor.
    SavableState(StateIn&);
  };

// just do a fwddec here, since StateIn and StateOut are not
// yet declared.
SavableState_REF_fwddec(SavableState);

////////////////////////////////////////////////////////////////////

//. \clsnm{SSRefBase} provides a few utility routines common to all
//. \clsnmref{SSRef} template instantiations.
class SSRefBase {
  protected:
    void check_castdown_result(void*, SavableState *, const ClassDesc *);
  public:
    virtual SavableState *sspointer() = 0;
    virtual void restore_state(StateIn&, const char *name = 0) = 0;
    void save_data_state(StateOut&);
    //. Save the state of the reference.
    void save_state(StateOut&);
};

// Include the smart pointer to SavableState templates and macros.
#include <util/state/stattmpl.h>
#ifdef USE_REF_MACROS
#include <util/state/statmacr.h>
#endif

////////////////////////////////////////////////////////////////////

class StateDataPtrSet;
class StateDataNumSet;
class ClassDescPintMap;

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
    StateDataPtrSet* ps_;
    ClassDescPintMap* classidmap_;
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
    virtual int put(int r);
    virtual int put(float r);
    virtual int put(double r);
    //. Write the given array data.  Size information is also saved.  The
    //data is allocated and read by the \srccd{get(T*&)} routines.
    virtual int put(const char*,int);
    virtual int put(const int*,int);
    virtual int put(const float*,int);
    virtual int put(const double*,int);
    //. Put arrays of data.  No size information is stored.  This
    //data is read by the \srccd{get\_array\_T} routines.
    virtual int put_array_char(const char*p,int size);
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

class StateClassData {
  public:
    int version;
    char *name;
    const ClassDesc *classdesc;
    int ninstance;
  public:
    StateClassData(int v=-1, const ClassDesc *c=0, char *name=0):
      version(v), classdesc(c), name(name), ninstance(0) {}
    StateClassData(const StateClassData &d) { operator=(d); }
    ~StateClassData();
    StateClassData &operator=(const StateClassData &d);
};
class intStateClassDataMap;

//. Objects of a type derived from \clsnmref{SavableState} can be
//. restored from a \clsnm{StateIn} object.
class StateIn: public DescribedClass {
#   define CLASSNAME StateIn
#   include <util/class/classda.h>
    friend class SavableState;
    friend class TranslateDataIn;
  private:
    // do not allow copy constructor or assignment
    StateIn(const StateIn&);
    void operator=(const StateIn&);
    int have_cd_;
    int dir_loc_;
  protected:
    TranslateDataIn *translate_;
    StateDataNumSet* ps_;
    int expected_object_num_;
    ClassDescPintMap* classidmap_;
    intStateClassDataMap* classdatamap_;
    int nextclassid_;
    int node_to_node_;
    int version_;
    int date_;
    char userid_[9];
    char format_;
    virtual int get_array_void(void*,int);

    void get_directory();
    int directory_location() const { return dir_loc_; }
    void find_and_get_directory();

    // The following members are called by friend SavableState

    //. This is used to restore an object.  It is called with the
    //. reference to the reference being restored.  If the data being
    //. restored has previously been restored, then the pointer
    //. being restored is set to a reference to the previously
    //. restored object.
    virtual int getobject(RefSavableState &);

    //. This restores objects that are listed in the directory.
    virtual int dir_getobject(RefSavableState &, const char *name);

    //. When storage has been allocated during object restoration,
    //. this routine is called with the object reference number
    //. and the pointer to the new storage so \srccd{getpointer}
    //. can find the data if it is referenced again.
    virtual void haveobject(int,const RefSavableState &);

    //. A call to nextobject followed by havepointer(int) is equiv
    //. to havepointer(int,void**);
    virtual void nextobject(int);
    virtual void haveobject(const RefSavableState &);

    void have_classdesc() { have_cd_ = 1; }
    int need_classdesc() { int tmp = have_cd_; have_cd_ = 0; return !tmp; }

    //. This restores \srccd{ClassDesc}'s.  It will set the
    //. pointer to the address of the static \srccd{ClassDesc} for
    //. the class which has the same name as the class that had
    //. the \clsnm{ClassDesc} that was saved by \srccd{put(const
    //. \clsnmref{ClassDesc}*)}.
    virtual int get(const ClassDesc**);
  public:
    StateIn();
    virtual ~StateIn();

    //. Read in the header information.  Changes the translation
    //scheme if necessary.
    virtual void get_header();

    //. Returns the version of the ClassDesc in the persistent object
    //. or -1 if info on the ClassDesc doesn't exist
    virtual int version(const ClassDesc*);
    
    //. This restores strings saved with
    //. \srccd{\clsnmref{StateOut}::putstring}.
    virtual int getstring(char*&);

    //. These restore data saved with \clsnmref{StateOut}'s \srccd{put}.
    //members.
    virtual int get(char&r);
    virtual int get(int&r);
    virtual int get(float&r);
    virtual int get(double&r);
    //. These restore data saved with \clsnmref{StateOut}'s \srccd{put}.
    //members.  The data is allocated by StateIn.
    virtual int get(char*&);
    virtual int get(int*&);
    virtual int get(float*&);
    virtual int get(double*&);
    //. These restore data saved with \clsnmref{StateOut}'s \srccd{put}.
    //members.  The data must be preallocated by the user.
    virtual int get_array_char(char*p,int size);
    virtual int get_array_int(int*p,int size);
    virtual int get_array_float(float*p,int size);
    virtual int get_array_double(double*p,int size);

    //. True if this is a node to node save/restore.  This is
    //necessary for classes that try to avoid saving databases
    //to files that can otherwise be read in, but want to avoid
    //reading the database from disk on all nodes.
    int node_to_node() const { return node_to_node_; }

    //. Returns true of this object uses a directory.
    virtual int use_directory();

    //. Return the current position in the file.
    virtual int tell();
    //. Set the current position in the file.  The default implementation
    //does nothing.
    virtual void seek(int);
    //. Return non-zero if seek does anything sensible.  The
    //default implementation returns 0.
    virtual int seekable();

    virtual void list_objects(ostream& =cout);
  };
DescribedClass_REF_dec(StateIn);

SavableState_REF_dec(SavableState);

////////////////////////////////////////////////////////////////////

#ifndef __GNUC__
static SavableState * att_hack_job(StateIn&si)
{
  return SavableState::restore_state(si);
}
#endif

// Now that StateIn and StateOut are complete, some savable arrays of
// basic types can be declared.

SSB_ARRAY_dec(int);
SSB_ARRAY2_dec(int);
SSB_ARRAY_dec(double);
SSB_ARRAY2_dec(double);

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
