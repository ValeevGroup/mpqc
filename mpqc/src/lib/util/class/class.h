
#ifdef __GNUG__
#pragma interface
#endif

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
    virtual DescribedClass* parentpointer() const = 0;
    virtual ~RefDescribedClassBase ();
    int operator==( const RefDescribedClassBase &a) const;
    int operator!=( const RefDescribedClassBase &a) const;
    int operator>=( const RefDescribedClassBase &a) const;
    int operator<=( const RefDescribedClassBase &a) const;
    int operator>( const RefDescribedClassBase &a) const;
    int operator<( const RefDescribedClassBase &a) const;
};

#ifdef TYPE_CONV_BUG
#  define DCREF_TYPE_CAST_DEC(refname,T)
#  define DCREF_TYPE_CAST_DEF(refname,T)
#else
#  define DCREF_TYPE_CAST_DEC(refname,T) operator T*() const
#  define DCREF_TYPE_CAST_DEF(refname,T) \
    refname :: operator T*() const { return p; }
#endif


// this uses macros from util/container/ref.h
#define DescribedClass_named_REF_dec(refname,T)				      \
class  refname : public RefDescribedClassBase  {			      \
  private:								      \
    T* p;								      \
  public:								      \
    DescribedClass* parentpointer() const;			      \
    T* operator->() const;					      \
    T* pointer() const;						      \
    DCREF_TYPE_CAST_DEC(refname,T);			\
    T& operator *() const;					      \
    refname ();								      \
    refname (T*a);							      \
    refname (const refname &a);						      \
    refname (const RefDescribedClassBase &);				      \
    ~refname ();							      \
    int null() const;							      \
    int nonnull() const;						      \
    void require_nonnull() const;					      \
    refname& operator=(T* cr);						      \
    refname& operator=(const RefDescribedClassBase & c);		      \
    refname& operator=(const refname & c);				      \
    void assign_pointer(T* cr);						      \
    void ref_info(FILE*fp=stdout) const;				      \
    void warn(const char *) const;					      \
    void clear();							      \
    void check_pointer() const;						      \
}
#define DescribedClass_named_REF_def(refname,T)				      \
T* refname :: operator->() const { return p; };			      \
T* refname :: pointer() const { return p; };			      \
DCREF_TYPE_CAST_DEF(refname,T);				\
T& refname :: operator *() const { return *p; };			      \
int refname :: null() const { return p == 0; };				      \
int refname :: nonnull() const { return p != 0; };			      \
void refname :: require_nonnull() const					      \
{									      \
  if (p == 0) {								      \
      fprintf(stderr,#refname": have null reference where nonnull needed\n"); \
      abort();								      \
    }									      \
};									      \
DescribedClass* refname :: parentpointer() const { return p; }	      \
refname :: refname (): p(0) {}						      \
refname :: refname (T*a): p(a)						      \
{									      \
  if (DO_REF_CHECK_STACK(p)) {						      \
      warn("Ref" # T ": creating a reference to stack data");		      \
    }									      \
  if (p) p->reference();						      \
  if (REF_CHECK_POINTER) check_pointer();				      \
}									      \
refname :: refname (const refname &a): p(a.p)				      \
{									      \
  if (p) p->reference();						      \
  if (REF_CHECK_POINTER) check_pointer();				      \
}									      \
refname :: refname (const RefDescribedClassBase &a)			      \
{									      \
  p = T::castdown((DescribedClass*)a.parentpointer());		              \
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
refname :: warn ( const char * msg) const				      \
{									      \
  fprintf(stderr,"WARNING: %s\n",msg);					      \
}									      \
refname& refname :: operator=(const refname & c)			      \
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
refname& refname :: operator=(const RefDescribedClassBase & c)		      \
{									      \
  T* cr = T::castdown((DescribedClass*)c.parentpointer());		      \
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
void refname :: check_pointer() const					      \
{									      \
  if (p && p->nreference() <= 0) {					      \
      warn("Ref" # T ": bad reference count in referenced object\n");	      \
    }									      \
}									      \
void refname :: ref_info(FILE*fp) const					      \
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

