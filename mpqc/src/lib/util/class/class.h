
#ifdef __GNUG__
#pragma interface
#endif

#ifndef _libQC_class_h
#define _libQC_class_h

#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <util/misc/identity.h>
#include <util/container/ref.h>
#include <util/container/array.h>
#include <util/container/set.h>

class ClassKeyClassDescPMap;
class ClassKeySet;
class DescribedClass;
class ClassDesc;
typedef ClassDesc* ClassDescP;

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
    operator=(const ParentClasses&);
  public:
    ParentClasses(const char*);
    ~ParentClasses();
    ParentClass& parent(int i) { return *_classes[i]; }
    const ParentClass& parent(int i) const { return *_classes[i]; }
    ParentClass& operator[](int i) { return *_classes[i]; }
    const ParentClass& operator[](int i) const { return *_classes[i]; }
    int n() const { return _n; }
    void change_parent(ClassDesc*oldcd,ClassDesc*newcd);
};
    

class RefKeyVal;
class StateIn;

class ClassDesc: public Identity {
    friend class ParentClasses; // ParentClasses initializes all_
  private:
    static ClassKeyClassDescPMap* all_;

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
    operator=(const ClassDesc&);
  public:
    ClassDesc(char*,int=1,char* p=0,
              DescribedClass* (*ctor)()=0,
              DescribedClass* (*keyvalctor)(const RefKeyVal&)=0,
              DescribedClass* (*stateinctor)(StateIn&)=0);
    ~ClassDesc();
    static void list_all_classes();
    static ClassKeyClassDescPMap& all();
    static ClassDesc* name_to_class_desc(const char*);
    const ParentClasses& parents() const { return parents_; }
    const char* name() const { return classname_; }
    int version() const { return version_; }
    DescribedClass* create_described_class() const;

    // create an object using the default constructor
    DescribedClass* create() const;

    // create an object using the keyval constructor
    DescribedClass* create(const RefKeyVal&) const;

    // create an object using the statein constructor
    DescribedClass* create(StateIn&) const;

};

ARRAY_dec(ClassDescP);
SET_dec(ClassDescP);
ARRAYSET_dec(ClassDescP);

// This makes info about the class available.
class DescribedClass : public VRefCount {
  private:
    static ClassDesc class_desc_;
  public:
    DescribedClass();
    DescribedClass(const DescribedClass&);
    DescribedClass& operator=(const DescribedClass&);
    static DescribedClass* castdown(DescribedClass*);
    static const ClassDesc* static_class_desc();
    virtual ~DescribedClass();
    virtual const ClassDesc* class_desc() const;
    virtual void* _castdown(const ClassDesc*);
    const char* class_name() const;
    int class_version() const;
  };

class  RefDescribedClassBase: private RefBase {
  public:
    RefDescribedClassBase() {}
    virtual DescribedClass* parentpointer() const = 0;
    virtual ~RefDescribedClassBase () {}
    void warn ( const char * msg) const { RefBase::warn(msg); }
    void warn_ref_to_stack() const { RefBase::warn_ref_to_stack(); }
    void warn_skip_stack_delete() const { RefBase::warn_skip_stack_delete(); }
    void warn_bad_ref_count() const { RefBase::warn_bad_ref_count(); }
    void ref_info(VRefCount*p,FILE*fp) const { RefBase::ref_info(p,fp); }
    void require_nonnull() const;
};

// These files declare template and macro smart pointer classes for
// DescribedClass objects.  They use macros from util/container/ref.h.
#include <util/class/clastmpl.h>
#include <util/class/clasmacr.h>

#define DescribedClass_named_REF_dec(name,T) DCRef_declare(T); \
                                             typedef class DCRef ## T name;
#define DescribedClass_named_REF_def(name,T)

// These macros choose a default name for the reference class formed from
// "Ref" followed by the type name.
#define DescribedClass_REF_dec(T) DescribedClass_named_REF_dec(Ref ## T,T)
#define DescribedClass_REF_def(T) DescribedClass_named_REF_def(Ref ## T,T)

DescribedClass_REF_dec(DescribedClass);
ARRAY_dec(RefDescribedClass);
SET_dec(RefDescribedClass);
ARRAYSET_dec(RefDescribedClass);

#endif

