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

#ifndef _libqc_state_h
#define _libqc_state_h

#ifdef __GNUC__
#pragma interface
#endif

#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>

#include <util/class/class.h>
#include <util/state/qc_xdr.h>
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

// If multiple inheritance is used, SavableState should be a virtual
// parent.  This breaks some compilers so the virtual_base macro should
// be used instead of the virtual keyword (see identity.h).  Since for StateIn
// CTORs the SavableState base class must be initialized with the
// StateIn object, the following macro should be used to determine if
// nondirect descendants of SavableState need to call the SavableState
// CTOR.  b might be a : or , and s is the StateIn
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
    //. returns a pointer to the derived class.
    static SavableState* restore_state(StateIn& si);

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

////////////////////////////////////////////////////////////////////

//. \clsnm{SSRefBase} provides a few utility routines common to all
//. \clsnmref{SSRef} template instantiations.
class SSRefBase {
  protected:
    void check_castdown_result(void*, SavableState *);
  public:
    virtual SavableState *sspointer() = 0;
    virtual void restore_state(StateIn&) = 0;
    void save_data_state(StateOut&);
    SavableState *restore_ss(StateIn&);
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
  private:
    // do not allow copy constructor or assignment
    StateOut(const StateOut&);
    void operator=(const StateOut&);
    int have_cd_;
  protected:
    int next_pointer_number;
    StateDataPtrSet* ps_;
    ClassDescPintMap* _classidmap;
    int _nextclassid;
    virtual int put_array_void(const void*,int);
    virtual void putparents(const ClassDesc*);
  public:
    StateOut();
    virtual ~StateOut();

    void have_classdesc() { have_cd_ = 1; }
    int need_classdesc() { int tmp = have_cd_; have_cd_ = 0; return tmp; }

    //. This will prepare \clsnm{StateOut} to output a pointer to
    //. data.  It first checks to see if the data has already been
    //. saved.  If it has, then a reference to this data is saved
    //. and 1 is returned.  Otherwise, 0 is returned and the class
    //. must continue with the save of the data referenced by this
    //. pointer.
    virtual int putpointer(void*);

    //. This is like \srccd{put} except the length of the \srccd{char}
    //. array is determined by interpreting the character array as
    //. a character string.
    virtual int putstring(char*);

    //. Write the version of the given \clsnmref{ClassDesc}, if
    //. it has not already been written.
    //. It shouldn't be necessary to call this member.
    virtual int put_version(const ClassDesc*);

    //. Write out information about the given \clsnmref{ClassDesc}.
    //. It shouldn't be necessary to call this member.
    virtual int put(const ClassDesc*);

    //. Write the given data.
    //. The member functions taking both a pointer and integer
    //. argument save a vector of the specified integer length of
    //. the appropiate type.  Additionally, the address of the
    //. data is kept in a table, so additional references to the
    //. data will not result in duplicated data being saved.
    //. Also, the restored object will contain the duplicated
    //. references as they appeared in the original.
    virtual int put(char r);
    virtual int put(int r);
    virtual int put(float r);
    virtual int put(double r);
    virtual int put(char*,int);
    virtual int put(int*,int);
    virtual int put(float*,int);
    virtual int put(double*,int);
    virtual int put_array_char(const char*p,int size);
    virtual int put_array_int(const int*p,int size);
    virtual int put_array_float(const float*p,int size);
    virtual int put_array_double(const double*p,int size);

    //. Don't keep track of pointers to objects.  Calling this
    //. causes duplicated references to objects to be copied.
    void forget_references();
    //. If a reference to an object that has already been written
    //. is encountered, copy it instead of generating a reference
    //. to the first object.
    void copy_references();

    //. Flush out any remaining data.
    virtual void flush();
  };
DescribedClass_REF_dec(StateOut);

//. Objects of a type derived from \clsnmref{SavableState} can be
//. restored from a \clsnm{StateIn} object.
class StateIn: public DescribedClass {
#   define CLASSNAME StateIn
#   include <util/class/classda.h>
  private:
    // do not allow copy constructor or assignment
    StateIn(const StateIn&);
    void operator=(const StateIn&);
    int have_cd_;
  protected:
    StateDataNumSet* ps_;
    int _nextobject;
    ArraysetCClassDescP _cd;
    Arrayint _version;
    virtual int get_array_void(void*,int);
  public:
    StateIn();
    virtual ~StateIn();

    //. This is used to restore a pointer.  It is called with the
    //. pointer to the pointer being restored.  If the data being
    //. restored has previously been restored, then the pointer
    //. being restored is set to a reference to the previously
    //. restored object and 0 is returned.  If the return value is
    //. nonzero then storage must be allocated for the data and
    //. the pointer to the new data along with the return value
    //. from \srccd{getpointer} must be given as arguments to the
    //. \srccd{havepointer} routine.  Note that \srccd{getpointer}
    //. and \srccd{havepointer} are automatically called for all
    //. of the above data types and need not be used, except for
    //. two dimensional arrays and other special situations.
    virtual int  getpointer(void**);

    //. When storage has been allocated during object restoration,
    //. this routine is called with the object reference number
    //. and the pointer to the new storage so \srccd{getpointer}
    //. can find the data if it is referenced again.
    virtual void havepointer(int,void*);

    //. A call to nextobject followed by havepointer(int) is equiv
    //. to havepointer(int,void**);
    virtual void nextobject(int);
    virtual void havepointer(void*);

    void have_classdesc() { have_cd_ = 1; }
    int need_classdesc() { int tmp = have_cd_; have_cd_ = 0; return tmp; }

    //. Returns the version of the ClassDesc in the persistent object
    //. or -1 if info on the ClassDesc doesn't exist
    int version(const ClassDesc*);
    
    //. This restores strings saved with
    //. \srccd{\clsnmref{StateOut}::putstring}.
    virtual int getstring(char*&);

    //. If the version of the ClassDesc in the persistent object
    //. has been read in yet, read it in.
    //. It shouldn't be necessary to call this member.
    virtual int get_version(const ClassDesc*);

    //. This restores \srccd{ClassDesc}'s.  It will set the
    //. pointer to the address of the static \srccd{ClassDesc} for
    //. the class which has the same name as the class that had
    //. the \clsnm{ClassDesc} that was saved by \srccd{put(const
    //. \clsnmref{ClassDesc}*)}.  It is not necessary for the user
    //. call this member.
    virtual int get(const ClassDesc**);

    //. These restore data saved with \srccd{\clsnmref{StateOut}::put}.
    virtual int get(char&r);
    virtual int get(int&r);
    virtual int get(float&r);
    virtual int get(double&r);
    virtual int get(char*&);
    virtual int get(int*&);
    virtual int get(float*&);
    virtual int get(double*&);
    virtual int get_array_char(char*p,int size);
    virtual int get_array_int(int*p,int size);
    virtual int get_array_float(float*p,int size);
    virtual int get_array_double(double*p,int size);

    //. Don't keep track of pointers to objects.  Calling this
    //. causes duplicated references to objects to be copied.
    void forget_references();
    //. If a reference to an object that has already been written
    //. is encountered, copy it instead of generating a reference
    //. to the first object.
    void copy_references();
  };
DescribedClass_REF_dec(StateIn);

SavableState_REF_dec(SavableState);

////////////////////////////////////////////////////////////////////

//. The \clsnmref{StateOutFile} provides a \clsnmref{StateOut}
//. which writes to files.  It is still abstract---one of its
//. derived classes, \clsnmref{StateOutFileText} or
//. \clsnmref{StateOutFileBin}, must be used to obtain a
//. \clsnmref{StateOut} object.  The
//. \clsnmref{StateOutFileText} class writes in a text format
//. and the \clsnmref{StateOutFileBin} writes in a binary
//. format.
class StateOutFile: public StateOut {
  private:
    // do not allow copy constructor or assignment
    StateOutFile(const StateOutFile&);
    operator=(const StateOutFile&);
  protected:
    int opened_;
    FILE* fp_;
  public:
    //. State information will be written to \srccd{stdout}.
    StateOutFile();
    //. State information will be written to \vrbl{fp}.
    StateOutFile(FILE*fp);
    //. State information will be written to \filnm{name}.
    StateOutFile(const char *name, const char *access = "w");

    ~StateOutFile();

    //. State information will be written to \filnm{name}.
    virtual int open(const char *name, const char *access = "w");
    //. Miscellaneous file operations.
    virtual void flush();
    virtual void close();
    virtual void rewind();
  };

//. The \clsnm{StateInFile} provides a \clsnmref{StateIn} which
//. reads from files.  It is still abstract---one of its
//. derived classes, \clsnmref{StateInFileText} or
//. \clsnmref{StateInFileBin}, must be used to obtain a
//. \clsnmref{StateIn} object.  The \clsnmref{StateInFileText} class
//. reads with a text format and the \clsnmref{StateInFileBin}
//. reads with a binary format.
class StateInFile: public StateIn {
  private:
    // do not allow copy constructor or assignment
    StateInFile(const StateInFile&);
    operator=(const StateInFile&);
  protected:
    int opened_;
    FILE* fp_;
  public:
    //. State information will be obtained from \srccd{stdin}.
    StateInFile();
    //. State information will be obtained from \vrbl{fp}.
    StateInFile(FILE*fp);
    //. State information will be obtained from \filnm{name}.
    StateInFile(const char *name, const char *access = "r");

    ~StateInFile();

    //. State information will be obtained from \filnm{name}.
    virtual int open(const char *name, const char *access = "r");
    //. Miscellaneous file operations.
    virtual void flush();
    virtual void close();
    virtual void rewind();
  };

////////////////////////////////////////////////////////////////////

//. \clsnm{StateOutText} writes out state information in an
//. almost human readable format.  It is intended for debugging
//. only.  The state information can read in again with
//. \clsnmref{StateInText}.
class StateOutText: public StateOutFile {
  private:
    // do not allow copy constructor or assignment
    StateOutText(const StateOutText&);
    operator=(const StateOutText&);
  protected:
    void newline();
    void comment(const char*,...);
    void start_array();
    void end_array();
    int putpointer(void*);
    void putparents(const ClassDesc*);
  public:
    StateOutText();
    StateOutText(FILE*);
    StateOutText(const char *, const char * = "w");
    ~StateOutText();
    int putstring(char*);
    int put_array_char(const char*,int);
    int put_array_int(const int*,int);
    int put_array_float(const float*,int);
    int put_array_double(const double*,int);
    int put(const ClassDesc*);
    int put(char r);
    int put(int r);
    int put(float r);
    int put(double r);
    int put(char*,int);
    int put(int*,int);
    int put(float*,int);
    int put(double*,int);
  };

//. \clsnm{StateInText} reads state information written
//. with \clsnm{StateOutText}.
class StateInText: public StateInFile {
  private:
    // do not allow copy constructor or assignment
    StateInText(const StateInText&);
    operator=(const StateInText&);
  protected:
    int _newlines;
    
    int read(char*);
    int read(int&);
    int read(float&);
    int read(double&);
    void newline();
    void comment();
    void start_array();
    void end_array();
    int  getpointer(void**);

    void abort();
  public:
    StateInText();
    StateInText(FILE*);
    StateInText(const char *, const char * = "r");
    ~StateInText();
    int getstring(char*&);
    int get_array_char(char*,int);
    int get_array_int(int*,int);
    int get_array_float(float*,int);
    int get_array_double(double*,int);
    int get(const ClassDesc**);
    int get(char&r);
    int get(int&r);
    int get(float&r);
    int get(double&r);
    int get(char*&);
    int get(int*&);
    int get(float*&);
    int get(double*&);
  };

////////////////////////////////////////////////////////////////////

//. \clsnm{StateOutBin} is used to write binary files.
class StateOutBin: public StateOutFile {
  private:
    // do not allow copy constructor or assignment
    StateOutBin(const StateOutBin&);
    operator=(const StateOutBin&);
  protected:
    int put_array_void(const void*,int);
  public:
    StateOutBin();
    StateOutBin(FILE*);
    StateOutBin(const char *, const char * = "w");
    ~StateOutBin();
  };

//. \clsnm{StateInBin} is used to read objects written with
//. \clsnm{StateOutBin}.
class StateInBin: public StateInFile {
  private:
    // do not allow copy constructor or assignment
    StateInBin(const StateInBin&);
    operator=(const StateInBin&);
  protected:
    int get_array_void(void*,int);
  public:
    StateInBin();
    StateInBin(FILE*);
    StateInBin(const char *, const char * = "r");
    ~StateInBin();
  };

////////////////////////////////////////////////////////////////////

//. \clsnm{StateOutBinXDR} is used to write binary files that
//. can be moved between machines with different endianess.
class StateOutBinXDR : public StateOutBin, public QCXDR
{
  private:
    // do not allow copy constructor or assignment
    StateOutBinXDR(const StateOutBinXDR&);
    operator=(const StateOutBinXDR&);
  protected:
    //this is needed for a mips-sgi-irix4 gcc 2.5.2 bug
    int put_array_void(const void*v,int i);
  public:
    StateOutBinXDR();
    StateOutBinXDR(FILE*);
    StateOutBinXDR(const char *, const char * = "w");
    ~StateOutBinXDR();
    int put_array_char(const char*,int);
    int put_array_int(const int*,int);
    int put_array_float(const float*,int);
    int put_array_double(const double*,int);
};

//. \clsnm{StateInBinXDR} is used to read objects written with
//. \clsnm{StateOutBinXDR}.
class StateInBinXDR : public StateInBin, public QCXDR
{
  private:
    // do not allow copy constructor or assignment
    StateInBinXDR(const StateInBinXDR&);
    operator=(const StateInBinXDR&);
  protected:
    //this is needed for a mips-sgi-irix4 gcc 2.5.2 bug
    int get_array_void(void*v,int i);
  public:
    StateInBinXDR();
    StateInBinXDR(FILE*);
    StateInBinXDR(const char *, const char * = "r");
    ~StateInBinXDR();
    int get_array_char(char*,int);
    int get_array_int(int*,int);
    int get_array_float(float*,int);
    int get_array_double(double*,int);
  };

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
