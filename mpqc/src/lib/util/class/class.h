//
// class.h
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

#ifdef __GNUG__
#pragma interface
#endif

#ifndef _util_class_class_h
#define _util_class_class_h

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <iostream.h>
#include <iomanip.h>
#include <util/ref/ref.h>
#include <util/container/array.h>
#include <util/container/set.h>

template <class T, class C>
class DescribedMemberDatum {
  private:
    T C::*member_;
  public:
    DescribedMemberDatum(T C::*member): member_(member) {}
    //T &member(C *c) { return c->*member_; }
};

class ClassKeyClassDescPMap;
class ClassKeySet;
class DescribedClass;
class ClassDesc;
typedef ClassDesc* ClassDescP;
typedef const ClassDesc* CClassDescP;

class ClassKey {
  private:
    char* classname_;
  public:
    ClassKey();
    ClassKey(const char* name);
    ClassKey(const ClassKey&);
    ~ClassKey();
    ClassKey& operator=(const ClassKey&);
    int operator==(ClassKey& ck);
    int hash() const;
    int cmp(ClassKey&ck) const;
    char* name() const;
  };

class ClassDesc;

class ParentClass
{
  public:
    enum Access { Private, Protected, Public };
  private:
    Access _access;
    int _is_virtual;
    ClassDesc* _classdesc;
  public:
    ParentClass(ClassDesc*,Access access = Private,int is_virtual = 0);
    ParentClass(const ParentClass&);
    ~ParentClass();
    int is_virtual() const;
    Access access() const { return _access; }
    const ClassDesc* classdesc() const;
    void change_classdesc(ClassDesc*n);
};

class ParentClasses
{
  private:
    int _n;
    ParentClass** _classes;
    void add(ParentClass*);
    // do not allow copy constructor or assignment
    ParentClasses(const ParentClasses&);
    void operator=(const ParentClasses&);
  public:
    ParentClasses();
    void init(const char*);
    ~ParentClasses();
    ParentClass& parent(int i) { return *_classes[i]; }
    const ParentClass& parent(int i) const { return *_classes[i]; }
    ParentClass& operator[](int i) { return *_classes[i]; }
    const ParentClass& operator[](int i) const { return *_classes[i]; }
    int n() const { return _n; }
    void change_parent(ClassDesc*oldcd,ClassDesc*newcd);
};
    

REF_fwddec(KeyVal);
class StateIn;

//. This class is used to contain information about classes.
//. Each \clsnmref{DescribedClass} type has a static \clsnm{ClassDesc}
//. member.  This member has lists of the parents, children
//. and virtual parents for each class.  The
//. \clsnm{ClassDesc} class also has a static member that is
//. a list of all described classes in the system.  These
//. lists are constructed as the constructors for the static
//. \clsnm{ClassDesc} members for each class are called and
//. are completed before \clsnm{main} is entered.
class ClassDesc: public Identity {
    friend class ParentClasses;
  private:
    static ClassKeyClassDescPMap* all_;
    static char * classlib_search_path_;
    static ClassKeySet* unresolved_parents_;

    char* classname_;
    int version_;
    ParentClasses parents_;
    ClassKeySet* children_;
    DescribedClass* (*ctor_)();
    DescribedClass* (*keyvalctor_)(const RefKeyVal&);
    DescribedClass* (*stateinctor_)(StateIn&);

    void change_parent(ClassDesc*oldcd,ClassDesc*newcd);

    // do not allow copy constructor or assignment
    ClassDesc(const ClassDesc&);
    void operator=(const ClassDesc&);
  public:
    ClassDesc(char*,int=1,char* p=0,
              DescribedClass* (*ctor)()=0,
              DescribedClass* (*keyvalctor)(const RefKeyVal&)=0,
              DescribedClass* (*stateinctor)(StateIn&)=0);
    ~ClassDesc();

    static ClassKeyClassDescPMap& all();
    const ParentClasses& parents() const { return parents_; }

    //. Writes a list of all of the classes to \srccd{stdout}.
    static void list_all_classes();
    //. Given the name of the class, return a pointer to the
    //. class descriptor.
    static ClassDesc* name_to_class_desc(const char*);
    //. Returns the name of the class.
    const char* name() const { return classname_; }
    //. Returns the version number of the class.
    int version() const { return version_; }
    //. This member has been replaced by \srccd{create()}.
    DescribedClass* create_described_class() const;
    //. Create an instance of \clsnmref{DescribedClass} with
    //. exact type equal to the class to which this class
    //. descriptor belongs.  The constructor which takes no
    //. arguments is used.  If this constructor doesn't exist or
    //. a static function that calls it with \srccd{new} wasn't
    //. given to this \clsnm{ClassDesc} when it was created, then
    //. 0 will be returned.
    virtual DescribedClass* create() const;
    //. Create an instance of \clsnmref{DescribedClass} with
    //. exact type equal to the class to which this class
    //. descriptor belongs.  The \clsnmref{KeyVal}\srccd{\&}
    //. constructor is used.  If this constructor doesn't exist
    //. or a static function that calls it with \srccd{new}
    //. wasn't passed to this \clsnm{ClassDesc}, then 0 will be
    //. returned.
    virtual DescribedClass* create(const RefKeyVal&) const;
    //. Create an instance of \clsnmref{DescribedClass} with exact
    //. type equal to the class to which this class descriptor
    //. belongs.  The \clsnmref{StateIn}\srccd{\&} constructor is used.
    //. If this constructor doesn't
    //. exist or a static function that calls it with \srccd{new}
    //. wasn't passed to this \clsnm{ClassDesc}, then 0
    //. will be returned.
    virtual DescribedClass* create(StateIn&) const;

    //. Attempt to dynamically load the shared object file for
    //. \vrbl{classname}
    static int load_class(const char* classname);
};

//. Classes which need runtime information about themselves and
//. their relationship to other classes can virtually inherit
//. from \clsnm{DescribedClass}.  This will provide the class
//. with the ability to query its name, query its version, and
//. perform safe castdown operations.  Furthermore, the class's
//. static \clsnmref{ClassDesc} can be obtained which permits
//. several other operations.
class DescribedClass : public VRefCount {
  private:
    static ClassDesc class_desc_;
  public:
    DescribedClass();
    DescribedClass(const DescribedClass&);
    DescribedClass& operator=(const DescribedClass&);
    virtual ~DescribedClass();
    //. Returns the argument.  This member is more interesting
    //. for types that derive from \clsnm{DescribedClass}.
    //. In this case the \srccd{castdown} member converts a
    //. \clsnm{DescribedClass} pointer to the derived type.
    //. If the type of the pointer is not the same as or
    //. (perhaps indirectly) derived from the type whose
    //. \srccd{castdown} member is being called, then 0 is return.
    //. The return type for derived class is a pointer to the
    //. derived class.  This member is implemented by the
    //. include files for derived types.
    static DescribedClass* castdown(DescribedClass*);
    //. Returns a pointer to the \clsnmref{ClassDesc} for the
    //. \clsnm{DescribedClass}.  Similar static members are
    //. implemented for derivatives of \clsnm{DescribedClass}
    //. by the provided include files.
    static const ClassDesc* static_class_desc();
    //. This returns the unique pointer to the
    //. \clsnmref{ClassDesc} for the object.  For derived types,
    //. this is declared and implemented for the programmer
    //. by the provided include files.
    virtual const ClassDesc* class_desc() const;
    //. Return the name of the object's exact type.
    const char* class_name() const;
    //. Return the version of the class.
    int class_version() const;
    //. This is a helper function that the programmer must override
    //. for derived types.
    virtual void* _castdown(const ClassDesc*);
    //. These print the object.
    virtual void print(ostream& =cout) const;
    virtual void print(ostream& =cout);
  };

//. \clsnm{DCRefBase} provides a few utility routines common to all
//. \clsnmref{DCRef} template instantiations.
class DCRefBase: private RefBase {
  protected:
    void reference(VRefCount *p);
    void dereference(VRefCount *p);
  public:
    DCRefBase() {}
    virtual ~DCRefBase();
    //. Returns the \clsnmref{DescribedClass} pointer for the
    //. contained object.
    virtual DescribedClass* parentpointer() const = 0;
    //. Requires that a nonnull reference is held.  If not,
    //. the program will abort.
    void require_nonnull() const;
    //. Ordering relations are provided using the \clsnmref{Identity}
    //. class.
    int operator==(const DescribedClass*a) const;
    int operator!=(const DescribedClass*a) const;
    int operator>=(const DescribedClass*a) const;
    int operator<=(const DescribedClass*a) const;
    int operator> (const DescribedClass*a) const;
    int operator< (const DescribedClass*a) const;
    int operator==(const DCRefBase &a) const;
    int operator!=(const DCRefBase &a) const;
    int operator>=(const DCRefBase &a) const;
    int operator<=(const DCRefBase &a) const;
    int operator> (const DCRefBase &a) const;
    int operator< (const DCRefBase &a) const;
    //. Miscellaneous utility functions.
    void warn(const char * msg) const;
    void warn_ref_to_stack() const;
    void warn_skip_stack_delete() const;
    void warn_bad_ref_count() const;
    void ref_info(VRefCount*p, ostream& os) const;
    void ref_info(ostream& os) const;
    void check_pointer() const;
};
ostream &operator<<(ostream&,const DCRefBase&);

inline void
DCRefBase::reference(VRefCount *p)
{
  if (p) {
#if REF_CHECK_STACK
      if (DO_REF_CHECK_STACK(p)) {
          DO_REF_UNMANAGE(p);
          warn_ref_to_stack();
        }
#endif
      p->reference();
    }
}

inline void
DCRefBase::dereference(VRefCount *p)
{
  if (p && p->dereference()<=0) {
      delete p;
    }
}

// These files declare template and macro smart pointer classes for
// DescribedClass objects.  They use macros from util/ref/ref.h.
#include <util/class/clastmpl.h>
#ifdef USE_REF_MACROS
#include <util/class/clasmacr.h>
#endif

#ifdef USE_REF_MACROS
#define DescribedClass_named_REF_dec(name,T) DCRef_declare(T); \
                                             typedef class DCRef ## T name;
#else
#define DCRef_declare(T) typedef class DCRef<T> DCRef ## T;
#define DescribedClass_named_REF_dec(name,T) typedef class DCRef<T> name;
#endif
#define DescribedClass_named_REF_def(name,T)

// These macros choose a default name for the reference class formed from
// "Ref" followed by the type name.
#define DescribedClass_REF_dec(T) DescribedClass_named_REF_dec(Ref ## T,T)
#define DescribedClass_REF_def(T) DescribedClass_named_REF_def(Ref ## T,T)

// This does forward declarations of REF classes.
#ifdef USE_REF_MACROS
#define DescribedClass_REF_fwddec(T) class DCRef ## T; \
                                     typedef class DCRef ## T Ref ## T;
#else
#define DescribedClass_REF_fwddec(T) class T; typedef class DCRef<T> Ref ## T;
#endif

DescribedClass_REF_dec(DescribedClass);

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
