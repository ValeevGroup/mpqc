
#ifndef _libQC_class_h
#define _libQC_class_h

#include <string.h>
#include <stdio.h>
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
    inline int cmp(ClassKey&ck) const { return strcmp(classname_,ck.classname_); }
    inline char* name() const {return classname_;}
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
    ~ParentClass();
    inline int is_virtual() const { return _is_virtual; };
    inline Access access() const { return _access; };
    inline const ClassDesc* classdesc() const { return _classdesc; };
    inline void change_classdesc(ClassDesc*n) { _classdesc = n; };
};

class ParentClasses
{
  private:
    int _n;
    ParentClass** _classes;
    void add(ParentClass*);
  public:
    ParentClasses(const char*);
    ~ParentClasses();
    inline const ParentClass& parent(int i) const { return *_classes[i]; };
    inline const ParentClass& operator[](int i) const { return *_classes[i]; };
    inline int n() const { return _n; };
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
  public:
    ClassDesc(char*,int=1,char* p=0,
              DescribedClass* (*ctor)()=0,
              DescribedClass* (*keyvalctor)(KeyVal&)=0,
              DescribedClass* (*stateinctor)(StateIn&)=0);
    ~ClassDesc();
    inline const ParentClasses& parents() const { return parents_; };
    inline const char* name() const { return classname_; }
    inline int version() const { return version_; }
    static void list_all_classes();
    static const ClassDesc* name_to_class_desc(const char*);
    inline DescribedClass* create_described_class() const { return create(); };

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
    static DescribedClass* castdown(DescribedClass*);
    static const ClassDesc* static_class_desc();
    virtual ~DescribedClass();
    virtual const ClassDesc* class_desc() const;
    virtual void* _castdown(const ClassDesc*);
    inline const char* class_name() const { return class_desc()->name(); }
    inline int class_version() const { return class_desc()->version(); }
  };

class  RefDescribedClassBase {
  public:
    inline RefDescribedClassBase() {};
    virtual DescribedClass* parentpointer() = 0;
    virtual ~RefDescribedClassBase ();
    inline int operator==( RefDescribedClassBase &a) {
        return parentpointer() == a.parentpointer();
      };
    inline int operator>=( RefDescribedClassBase &a) {
        return parentpointer() >= a.parentpointer();
      };
    inline int operator<=( RefDescribedClassBase &a) {
        return parentpointer() <= a.parentpointer();
      };
    inline int operator>( RefDescribedClassBase &a) {
        return parentpointer() > a.parentpointer();
      };
    inline int operator<( RefDescribedClassBase &a) {
        return parentpointer() < a.parentpointer();
      };
};
#define DescribedClass_REF_dec(T)					      \
class  Ref ## T : public RefDescribedClassBase  {			      \
  private:								      \
    T* p;								      \
  public:								      \
    DescribedClass* parentpointer();					      \
    inline T* operator->() { return p; };				      \
    inline const T* operator->() const { return p; };			      \
    inline T* pointer() { return p; };					      \
    inline const T* pointer() const { return p; };			      \
    inline operator T*() { return p; };					      \
    inline const operator T*() const { return p; };			      \
    inline T& operator *() { return *p; };				      \
    inline const T& operator *() const { return *p; };			      \
    Ref ## T ();							      \
    Ref ## T (T*a);							      \
    Ref ## T ( Ref ## T &a);						      \
    Ref ## T ( RefDescribedClassBase &);				      \
    ~Ref ## T ();							      \
    inline int null() { return p == 0; };				      \
    inline int nonnull() { return p != 0; };				      \
    Ref ## T  operator=(T* cr);						      \
    Ref ## T  operator=( RefDescribedClassBase & c);			      \
    Ref ## T  operator=( Ref ## T & c);					      \
    void assign_pointer(T* cr);						      \
    void  ref_info(FILE*fp=stdout);					      \
}
#define DescribedClass_REF_def(T)					      \
DescribedClass* Ref ## T :: parentpointer() { return p; }		      \
Ref ## T :: Ref ## T (): p(0) {}					      \
Ref ## T :: Ref ## T (T*a): p(a) { if (p) p->count++; }			      \
Ref ## T :: Ref ## T ( Ref ## T &a): p(a.p) { if (p) p->count++; }	      \
Ref ## T :: Ref ## T ( RefDescribedClassBase &a)			      \
{									      \
  p = T::castdown(a.parentpointer());					      \
  if (p) p->count++;							      \
}									      \
Ref ## T :: ~Ref ## T () { if (p && --p->count<=0) delete p; }		      \
Ref ## T  Ref ## T :: operator=( Ref ## T & c)				      \
{									      \
  if (   p								      \
         && --p->count <= 0						      \
         && p != c.p)							      \
    delete p;								      \
  p=c.p;								      \
  if (p) p->count++;							      \
  return *this;								      \
}									      \
Ref ## T  Ref ## T :: operator=( RefDescribedClassBase & c)		      \
{									      \
  T* newp = T::castdown(c.parentpointer());				      \
  if (   p								      \
         && --p->count <= 0						      \
         && p != newp)							      \
    delete p;								      \
  p=newp;								      \
  if (p) p->count++;							      \
  return *this;								      \
}									      \
Ref ## T  Ref ## T :: operator=(T* cr)					      \
{									      \
  if (p && --p->count <= 0 && p != cr) delete p;			      \
  p=cr;									      \
  if (p) p->count++;							      \
  return *this;								      \
}									      \
void Ref ## T :: assign_pointer(T* cr)					      \
{									      \
  if (p && --p->count <= 0 && p != cr) delete p;			      \
  p=cr;									      \
  if (p) p->count++;							      \
}									      \
void Ref ## T :: ref_info(FILE*fp=stdout)				      \
{									      \
  if (nonnull()) fprintf(fp,"count = %d\n",p->count);			      \
  else fprintf(fp,"reference is null\n");				      \
}

DescribedClass_REF_dec(DescribedClass);
ARRAY_dec(RefDescribedClass);
SET_dec(RefDescribedClass);
ARRAYSET_dec(RefDescribedClass);

#endif

