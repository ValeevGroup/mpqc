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

#ifndef _util_class_class_h
#define _util_class_class_h

#include <map>
#include <set>
#include <string>

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <util/ref/ref.h>
#include <util/misc/exenv.h>
#include <util/misc/exception.h>

namespace sc {

class DescribedClass;
class ClassDesc;
typedef ClassDesc* ClassDescP;
typedef const ClassDesc* CClassDescP;

class ClassDesc;

/// Gives one parent class of a class.
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

/// Gives a list of parent classes of a class.
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
    

class KeyVal;
class StateIn;

/** This is used to pass a function that make void constructor calls to the
    ClassDesc constructor. */
template <class T>
DescribedClass* create()
{
  return new T;
}

/** This is used to pass a function that make KeyVal constructor calls to
    the ClassDesc constructor. */
template <class T>
DescribedClass* create(const Ref<KeyVal>& keyval)
{
  return new T(keyval);
}

/** This is used to pass a function that make StateIn constructor calls to
    the ClassDesc constructor. */
template <class T>
DescribedClass* create(StateIn& statein)
{
  return new T(statein);
}

class type_info_key {
  private:
    const std::type_info *ti_;
  public:
    type_info_key(): ti_(0) {}
    type_info_key(const std::type_info *ti): ti_(ti) {}
    type_info_key& operator=(const type_info_key&);
    int operator==(const type_info_key&) const;
    int operator<(const type_info_key&) const;
    int cmp(const type_info_key&) const;
};

/** This class is used to contain information about classes.
 Each DescribedClass type has a static ClassDesc
 member.  This member has lists of the parents, children
 and virtual parents for each class.  The
 ClassDesc class also has a static member that is
 a list of all described classes in the system.  These
 lists are constructed as the constructors for the static
 ClassDesc members for each class are called and
 are completed before main is entered.  See \ref class for
 more information.
*/
class ClassDesc {
    friend class ParentClasses;
  private:
    static std::map<std::string,ClassDescP> *all_;
    static std::map<type_info_key,ClassDescP> *type_info_all_;
    static char * classlib_search_path_;
    static std::set<std::string> *unresolved_parents_;

    char* classname_;
    int version_;
    ParentClasses parents_;
    std::set<std::string> *children_;
    DescribedClass* (*ctor_)();
    DescribedClass* (*keyvalctor_)(const Ref<KeyVal>&);
    DescribedClass* (*stateinctor_)(StateIn&);
    const std::type_info *ti_;

    void change_parent(ClassDesc*oldcd,ClassDesc*newcd);

    // do not allow copy constructor or assignment
    ClassDesc(const ClassDesc&);
    void operator=(const ClassDesc&);

    // this is used for temporary parent class descriptors
    ClassDesc(const char*);
    void init(const char*,int=1,const char* p=0,
              const std::type_info *ti=0,
              DescribedClass* (*ctor)()=0,
              DescribedClass* (*keyvalctor)(const Ref<KeyVal>&)=0,
              DescribedClass* (*stateinctor)(StateIn&)=0);
  public:
    ClassDesc(const std::type_info&, const char*,int=1,const char* p=0,
              DescribedClass* (*ctor)()=0,
              DescribedClass* (*keyvalctor)(const Ref<KeyVal>&)=0,
              DescribedClass* (*stateinctor)(StateIn&)=0);
    ~ClassDesc();

    static std::map<std::string,ClassDescP>& all();
    const ParentClasses& parents() const { return parents_; }

    /// Writes a list of all of the classes to ExEnv::out0().
    static void list_all_classes();
    /** Given the name of the class, return a pointer to the
        class descriptor. */
    static ClassDesc* name_to_class_desc(const char*);
    /** Given a type_info object return a pointer to the ClassDesc. */
    static ClassDesc *class_desc(const std::type_info &);
    /// Returns the name of the class.
    const char* name() const { return classname_; }
    /// Returns the version number of the class.
    int version() const { return version_; }
    /// This member has been replaced by create().
    DescribedClass* create_described_class() const;
    /** Create an instance of DescribedClass with
        exact type equal to the class to which this class
        descriptor belongs.  The constructor which takes no
        arguments is used.  If this constructor doesn't exist or
        a static function that calls it with new wasn't
        given to this ClassDesc when it was created, then
        0 will be returned. */
    virtual DescribedClass* create() const;
    /** Create an instance of DescribedClass with exact type equal to the
        class to which this class descriptor belongs.  The KeyVal&
        constructor is used.  If this constructor doesn't exist or a static
        function that calls it with new wasn't passed to this ClassDesc,
        then 0 will be returned. */
    virtual DescribedClass* create(const Ref<KeyVal>&) const;
    /** Create an instance of DescribedClass with exact type equal to the
        class to which this class descriptor belongs.  The StateIn&
        constructor is used.  If this constructor doesn't exist or a static
        function that calls it with new wasn't passed to this ClassDesc,
        then 0 will be returned. */
    virtual DescribedClass* create(StateIn&) const;

    /** Attempt to dynamically load the shared object file for
        classname. */
    static int load_class(const char* classname);
};

/** Classes which need runtime information about themselves and their
    relationship to other classes can virtually inherit from
    DescribedClass.  This will provide the class with the ability to query
    its name and its version.
    Furthermore, the class's static ClassDesc can be obtained
    which permits several other operations.  See \ref class for
    more information. */
class DescribedClass : virtual public RefCount {
  public:
    DescribedClass();
    DescribedClass(const DescribedClass&);
    DescribedClass& operator=(const DescribedClass&);
    virtual ~DescribedClass();
    /** This returns the unique pointer to the ClassDesc corresponding
        to the given type_info object.  Null is returned if it fails. */
    ClassDesc* class_desc() const MPQC__NOEXCEPT;
    /// Return the name of the object's exact type.
    const char* class_name() const;
    /// Return the version of the class.
    int class_version() const;
    /// Print the object.
    virtual void print(std::ostream& = ExEnv::out0()) const;
    /// Return this object wrapped up in a Ref smart pointer.  This
    /// member is mainly a convenience function for the Python MPQC
    /// interface.
    Ref<DescribedClass> ref() { return Ref<DescribedClass>(this); }
  };

/** Return the ClassDesc corresponding to template argument. */
template <class T>
inline ClassDesc *
class_desc()
{
  return ClassDesc::class_desc(typeid(T));
}

/** Return the ClassDesc corresponding to the exact type for the
    argument. */
inline ClassDesc *
class_desc(DescribedClass *d)
{
  return ClassDesc::class_desc(typeid(*d));
}

/** Attempt to cast a DescribedClass pointer to a DescribedClass
    descendent.  It is an error for the result to be a null pointer. */
template<class T>
inline T
require_dynamic_cast(DescribedClass*p,const char * errmsg,...)
{
  T t = dynamic_cast<T>(p);
  if (p && !t) {
      va_list args;
      va_start(args,errmsg);
      fprintf(stderr,"A required dynamic_cast failed in: ");
      vfprintf(stderr,errmsg,args);
      fprintf(stderr,"\nwanted type \"%s\" but got \"%s\"\n",
              typeid(T).name(),p->class_desc()->name());
      fflush(stderr);
      va_end(args);
      abort();
  }
  return t;
}

/** Attempt to cast a const DescribedClass pointer to a DescribedClass
    descendent.  It is an error for the result to be a null pointer. */
template<class T>
inline T
require_dynamic_cast(const DescribedClass*p,const char * errmsg,...)
{
  T t = dynamic_cast<T>(p);
  if (p && !t) {
      va_list args;
      va_start(args,errmsg);
      fprintf(stderr,"A required dynamic_cast failed in: ");
      vfprintf(stderr,errmsg,args);
      fprintf(stderr,"\nwanted type \"%s\" but got \"%s\"\n",
              typeid(T).name(),p->class_desc()->name());
      fflush(stderr);
      va_end(args);
      abort();
  }
  return t;
}

/** This, together with ForceLink, is used to force code for particular
    classes to be linked into executables. */
template <class A>
class ForceLinkBase {
  public:
    ForceLinkBase() {};
    virtual ~ForceLinkBase() {};
    virtual DescribedClass *create(A) = 0;
};

/** This, together with ForceLinkBase, is used to force code for particular
classes to be linked into executables.  Objects are created from input and
checkpoint files by using class name lookup to find that class's ClassDesc
object.  The ClassDesc object has members that can create the class.
Unfortunately, linking in a library doesn't cause code for the
ClassDesc, and thus the class itself, to be linked.  ForceLink objects are
created in linkage.h files for each library.  The code containing the main
routine for an executable can include these linkage files to force code for
that library's classes to be linked. */
template <class T, class A = const Ref<KeyVal> &>
class ForceLink: public ForceLinkBase<A> {
  public:
    ForceLink() {};
    virtual ~ForceLink() {};
    DescribedClass *create(A a) { return new T(a); }
};


/**
 * computes Damerau–Levenshtein distance between two strings
 * @param str1 a string
 * @param str2 another string
 * @return the distance between the strings
 */
std::string::size_type
string_distance(const std::string& str1,
                const std::string& str2);

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
