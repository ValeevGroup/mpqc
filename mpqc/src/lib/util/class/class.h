
#ifdef __GNUG__
#pragma interface
#endif

#ifndef _util_class_class_h
#define _util_class_class_h

#include <string.h>
#include <stdio.h>
#include <stdarg.h>
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
    operator=(const ParentClasses&);
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
    

class RefKeyVal;
class StateIn;

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
    virtual DescribedClass* create() const;

    // create an object using the keyval constructor
    virtual DescribedClass* create(const RefKeyVal&) const;

    // create an object using the statein constructor
    virtual DescribedClass* create(StateIn&) const;

    static int load_class(const char* classname);
};

ARRAY_dec(ClassDescP);
SET_dec(ClassDescP);
ARRAYSET_dec(ClassDescP);

ARRAY_dec(CClassDescP);
SET_dec(CClassDescP);
ARRAYSET_dec(CClassDescP);

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

class DCRefBase: private RefBase {
  public:
    DCRefBase() {}
    virtual DescribedClass* parentpointer() const = 0;
    virtual ~DCRefBase();
    void warn(const char * msg) const;
    void warn_ref_to_stack() const;
    void warn_skip_stack_delete() const;
    void warn_bad_ref_count() const;
    void ref_info(VRefCount*p,FILE*fp) const;
    void require_nonnull() const;
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
    void check_pointer() const;
    void ref_info(FILE*fp) const;
    void reference(VRefCount *p);
    void dereference(VRefCount *p);
};

// These files declare template and macro smart pointer classes for
// DescribedClass objects.  They use macros from util/ref/ref.h.
#include <util/class/clastmpl.h>
#include <util/class/clasmacr.h>

#define DescribedClass_named_REF_dec(name,T) DCRef_declare(T); \
                                             typedef class DCRef ## T name;
#define DescribedClass_named_REF_def(name,T)

// These macros choose a default name for the reference class formed from
// "Ref" followed by the type name.
#define DescribedClass_REF_dec(T) DescribedClass_named_REF_dec(Ref ## T,T)
#define DescribedClass_REF_def(T) DescribedClass_named_REF_def(Ref ## T,T)

// This does forward declarations of REF classes.
#define DescribedClass_REF_fwddec(T) class DCRef ## T; \
                                     typedef class DCRef ## T Ref ## T;

DescribedClass_REF_dec(DescribedClass);
ARRAY_dec(RefDescribedClass);
SET_dec(RefDescribedClass);
ARRAYSET_dec(RefDescribedClass);

#endif

