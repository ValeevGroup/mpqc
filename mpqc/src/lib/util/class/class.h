
#ifndef _libQC_class_h
#define _libQC_class_h

#include <string.h>
#include <stdio.h>
#include <stdarg.h>
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
    Access access() const;
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
    ParentClass& parent(int i);
    const ParentClass& parent(int i) const;
    ParentClass& operator[](int i);
    const ParentClass& operator[](int i) const;
    int n() const;
    void change_parent(ClassDesc*oldcd,ClassDesc*newcd);
};
    

class KeyVal;
class StateIn;

class ClassDesc {
  private:
    char* classname_;
    int version_;
    ParentClasses parents_;
    ClassKeySet* children_;
    DescribedClass* (*ctor_)();
    DescribedClass* (*keyvalctor_)(KeyVal&);
    DescribedClass* (*stateinctor_)(StateIn&);

    void change_parent(ClassDesc*oldcd,ClassDesc*newcd);

    // do not allow copy constructor or assignment
    ClassDesc(const ClassDesc&);
    operator=(const ClassDesc&);
  public:
    ClassDesc(char*,int=1,char* p=0,
              DescribedClass* (*ctor)()=0,
              DescribedClass* (*keyvalctor)(KeyVal&)=0,
              DescribedClass* (*stateinctor)(StateIn&)=0);
    ~ClassDesc();
    static void list_all_classes();
    static ClassDesc* name_to_class_desc(const char*);
    const ParentClasses& parents() const;
    const char* name() const;
    int version() const;
    DescribedClass* create_described_class() const;

    // create an object using the default constructor
    DescribedClass* create() const;

    // create an object using the keyval constructor
    DescribedClass* create(KeyVal&) const;

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

class  RefDescribedClassBase {
  public:
    RefDescribedClassBase();
    RefDescribedClassBase(const RefDescribedClassBase&);
    RefDescribedClassBase& operator=(const RefDescribedClassBase&);
    virtual DescribedClass* parentpointer() = 0;
    virtual ~RefDescribedClassBase ();
    int operator==( RefDescribedClassBase &a);
    int operator>=( RefDescribedClassBase &a);
    int operator<=( RefDescribedClassBase &a);
    int operator>( RefDescribedClassBase &a);
    int operator<( RefDescribedClassBase &a);
};

// this uses macros from util/container/ref.h
#define DescribedClass_named_REF_dec(refname,T)				      \
class  refname : public RefDescribedClassBase  {			      \
  private:								      \
    T* p;								      \
  public:								      \
    DescribedClass* parentpointer();					      \
    T* operator->();							      \
    const T* operator->() const;					      \
    T* pointer();							      \
    const T* pointer() const;						      \
    operator T*();							      \
    operator const T*() const;						      \
    T& operator *();							      \
    const T& operator *() const;					      \
    refname ();								      \
    refname (T*a);							      \
    refname ( refname &a);						      \
    refname ( RefDescribedClassBase &);					      \
    ~refname ();							      \
    int null();								      \
    int nonnull();							      \
    void require_nonnull();						      \
    refname& operator=(T* cr);						      \
    refname& operator=( RefDescribedClassBase & c);			      \
    refname& operator=( refname & c);					      \
    void assign_pointer(T* cr);						      \
    void  ref_info(FILE*fp=stdout);					      \
    void warn(const char *);						      \
    void clear();							      \
    void check_pointer();						      \
}
#define DescribedClass_named_REF_def(refname,T)				      \
T* refname :: operator->() { return p; };				      \
const T* refname :: operator->() const { return p; };			      \
T* refname :: pointer() { return p; };					      \
const T* refname :: pointer() const { return p; };			      \
refname :: operator T*() { return p; };					      \
refname :: operator const T*() const { return p; };			      \
T& refname :: operator *() { return *p; };				      \
const T& refname :: operator *() const { return *p; };			      \
int refname :: null() { return p == 0; };				      \
int refname :: nonnull() { return p != 0; };				      \
void refname :: require_nonnull()					      \
{									      \
  if (p == 0) {								      \
      fprintf(stderr,#refname": have null reference where nonnull needed\n"); \
      abort();								      \
    }									      \
};									      \
DescribedClass* refname :: parentpointer() { return p; }		      \
refname :: refname (): p(0) {}						      \
refname :: refname (T*a): p(a)						      \
{									      \
  if (DO_REF_CHECK_STACK(p)) {						      \
      warn("Ref" # T ": creating a reference to stack data");		      \
    }									      \
  if (p) p->reference();						      \
  if (REF_CHECK_POINTER) check_pointer();				      \
}									      \
refname :: refname ( refname &a): p(a.p)				      \
{									      \
  if (p) p->reference();						      \
  if (REF_CHECK_POINTER) check_pointer();				      \
}									      \
refname :: refname ( RefDescribedClassBase &a)				      \
{									      \
  p = T::castdown(a.parentpointer());					      \
  if (p) p->reference();						      \
  if (REF_CHECK_POINTER) check_pointer();				      \
}									      \
refname :: ~refname ()							      \
{									      \
  clear();								      \
}									      \
void									      \
refname :: clear()							      \
{									      \
  if (REF_CHECK_POINTER) check_pointer();				      \
  if (p && p->dereference()<=0) {					      \
      if (DO_REF_CHECK_STACK(p)) {					      \
          warn("Ref" # T ": skipping delete of object on the stack");	      \
        }								      \
      else {								      \
           delete p;							      \
         }								      \
    }									      \
  p = 0;								      \
}									      \
void									      \
refname :: warn ( const char * msg)					      \
{									      \
  fprintf(stderr,"WARNING: %s\n",msg);					      \
}									      \
refname& refname :: operator=( refname & c)				      \
{									      \
  if (c.p) c.p->reference();						      \
  clear();								      \
  p=c.p;								      \
  if (REF_CHECK_POINTER) check_pointer();				      \
  return *this;								      \
}									      \
refname& refname :: operator=(T* cr)					      \
{									      \
  if (cr) cr->reference();						      \
  clear();								      \
  p = cr;								      \
  if (REF_CHECK_POINTER) check_pointer();				      \
  return *this;								      \
}									      \
refname& refname :: operator=( RefDescribedClassBase & c)		      \
{									      \
  T* cr = T::castdown(c.parentpointer());				      \
  if (cr) cr->reference();						      \
  clear();								      \
  p = cr;								      \
  if (REF_CHECK_POINTER) check_pointer();				      \
  return *this;								      \
}									      \
void refname :: assign_pointer(T* cr)					      \
{									      \
  if (cr) cr->reference();						      \
  clear();								      \
  p = cr;								      \
  if (REF_CHECK_POINTER) check_pointer();				      \
}									      \
void refname :: check_pointer()						      \
{									      \
  if (p && p->nreference() <= 0) {					      \
      warn("Ref" # T ": bad reference count in referenced object\n");	      \
    }									      \
}									      \
void refname :: ref_info(FILE*fp)					      \
{									      \
  if (nonnull()) fprintf(fp,"nreference() = %d\n",p->nreference());	      \
  else fprintf(fp,"reference is null\n");				      \
}

// These macros choose a default name for the reference class formed from
// "Ref" followed by the type name.
#define DescribedClass_REF_dec(T) DescribedClass_named_REF_dec(Ref ## T, T)
#define DescribedClass_REF_def(T) DescribedClass_named_REF_def(Ref ## T, T)


DescribedClass_REF_dec(DescribedClass);
ARRAY_dec(RefDescribedClass);
SET_dec(RefDescribedClass);
ARRAYSET_dec(RefDescribedClass);

#endif

