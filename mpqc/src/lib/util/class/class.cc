
#ifdef __GNUG__
#pragma implementation
#endif

#include <stdio.h>
#include <stdlib.h>

#include "class.h"

#include "classMap.h"
#include "classImplMap.h"

#include "classkeySet.h"
#include "classkeyImplSet.h"

static ClassKeyClassDescPMap* all_;

/////////////////////////////////////////////////////////////////

ClassKey::ClassKey() :
  classname_(0)
{
}

ClassKey::ClassKey(const char* name):
  classname_(::strcpy(new char[strlen(name)+1],name))
{
}

ClassKey::ClassKey(const ClassKey& key)
{
  if (key.classname_) {
      classname_ = ::strcpy(new char[strlen(key.classname_)+1],key.classname_);
    }
  else {
      classname_ = 0;
    }
}

ClassKey::~ClassKey()
{
  delete[] classname_;
}

ClassKey& ClassKey::operator=(const ClassKey& key)
{
  if (classname_ && classname_ != key.classname_) {
      delete[] classname_;
      classname_ = ::strcpy(new char[strlen(key.classname_)+1],key.classname_);
    }
  return *this;
}

int ClassKey::operator==(ClassKey& ck)
{
  if (!classname_) {
      if (!ck.classname_) return 1;
      else return 0;
    }
  else if (!ck.classname_) return 0;
  
  return !::strcmp(classname_,ck.classname_);
}

int
ClassKey::hash() const
{
  int r=0;
  int i;

  // Even numbered bytes make up the lower part of the hash index
  for (i=0; i<::strlen(classname_); i+=2) {
      r ^= classname_[i];
    }

  // Odd numbered bytes make up the upper part of the hash index
  for (i=1; i<::strlen(classname_); i+=2) {
      r ^= classname_[i]<<8;
    }

  return r;
}

int ClassKey::cmp(ClassKey&ck) const
{
  return strcmp(classname_,ck.classname_);
}

char* ClassKey::name() const
{
  return classname_;
}

/////////////////////////////////////////////////////////////////

ParentClass::ParentClass(ClassDesc*classdesc,Access access,int is_virtual):
  _access(access),
  _is_virtual(is_virtual),
  _classdesc(classdesc)
{
}

ParentClass::ParentClass(const ParentClass&p):
  _access(p._access),
  _is_virtual(p._is_virtual),
  _classdesc(p._classdesc)
{
}

ParentClass::~ParentClass()
{
}

int ParentClass::is_virtual() const
{
  return _is_virtual;
}

ParentClass::Access ParentClass::access() const
{
  return _access;
}

const ClassDesc* ParentClass::classdesc() const
{
  return _classdesc;
}

void ParentClass::change_classdesc(ClassDesc*n)
{
  _classdesc = n;
}

/////////////////////////////////////////////////////////////////

ParentClasses::ParentClasses(const char* parents):
  _n(0),
  _classes(0)
{

  // if this the first ParentClasses to be initialized all_ will
  // be zero and I must initialize it.
  if (!all_) {
      all_ = new MAPCTOR;
    }

  // if parents is empty then we are done
  if (!parents || strlen(parents) == 0) return;

  char* tokens = ::strcpy(new char[strlen(parents)+1],parents);
  const char* whitesp = "\t\n,() ";
  char* token;
  int is_virtual = 0;
  ParentClass::Access access = ParentClass::Private;
  for (token = ::strtok(tokens,whitesp);
       token;
       token = ::strtok(0,whitesp)) {
      if (!strcmp(token,"virtual")) {
          is_virtual = 1;
        }
      else if (!strcmp(token,"public")) {
          access = ParentClass::Public;
        }
      else if (!strcmp(token,"protected")) {
          access = ParentClass::Protected;
        }
      else if (!strcmp(token,"private")) {
          access = ParentClass::Private;
        }
      else {
          ClassKey parentkey(token);
          // if the parents class desc does not exist create a temporary
          // the temporary will be incorrect,because it does not have the
          // parent's parents   ' for braindead compiler
          if (!(*all_)[parentkey]) {
              (*all_)[parentkey] = new ClassDesc(token);
            }
          ParentClass* p = new ParentClass((*all_)[parentkey],
                                           access,
                                           is_virtual);
          add(p);
          access = ParentClass::Private;
          is_virtual = 0;
        }
    }
  delete[] tokens;
  
}

ParentClasses::~ParentClasses()
{
  for (int i=0; i<_n; i++) delete _classes[i];
  
  if (_classes) delete[] _classes;
  _classes = 0;
  _n = 0;
}

void
ParentClasses::add(ParentClass*p)
{
  ParentClass** newpp = new ParentClass*[_n+1];
  for (int i=0; i<_n; i++) newpp[i] = _classes[i];
  newpp[_n] = p;
  _n++;
  delete[] _classes;
  _classes = newpp;
}

void
ParentClasses::change_parent(ClassDesc*oldcd,ClassDesc*newcd)
{
  for (int i=0; i<_n; i++) {
      if (parent(i).classdesc() == oldcd) parent(i).change_classdesc(newcd);
    }
}

ParentClass& ParentClasses::parent(int i)
{
  return *_classes[i];
}

const ParentClass& ParentClasses::parent(int i) const
{
  return *_classes[i];
}

ParentClass& ParentClasses::operator[](int i)
{
  return *_classes[i];
}

const ParentClass& ParentClasses::operator[](int i) const
{
  return *_classes[i];
}

int ParentClasses::n() const
{
  return _n;
}

////////////////////////////////////////////////////////////////////////

ARRAY_def(ClassDescP);
SET_def(ClassDescP);
ARRAYSET_def(ClassDescP);

ClassDesc::ClassDesc(char* name, int version,
                     char* parents,
                     DescribedClass* (*ctor)(),
                     DescribedClass* (*keyvalctor)(KeyVal&),
                     DescribedClass* (*stateinctor)(StateIn&)
                     ):
  classname_(0),
  version_(version),
  parents_(parents),
  children_(0),
  ctor_(ctor),
  keyvalctor_(keyvalctor),
  stateinctor_(stateinctor)
{

  classname_ = ::strcpy(new char[strlen(name)+1],name);

  // test the version number to see if it is valid
  if (version <= 0) {
      fprintf(stderr,"error in ClassDesc ctor: version <= 0\n");
      exit(1);
    }

  ClassKey key(name);

  // let each of the parents know that this is a child
  for (int i=0; i<parents_.n(); i++) {
      ClassKey parentkey(parents_[i].classdesc()->name());
      if (!(*all_)[parentkey]->children_)
        (*all_)[parentkey]->children_ = new SETCTOR;
      // let the parents know about the child
#ifndef GNUBUG
      ((*all_)[parentkey]->children_)->add(key);
#else
      ClassKeySet *ppp=(*all_)[parentkey]->children_;
      ppp->add(key);
#endif
    }

  // if this class is aleady in all_, then it was put there by a child
  // preserve children info, destroy the old entry, and put this there
  if ((*all_)[key]) {
      children_ = (*all_)[key]->children_;
      (*all_)[key]->children_ = 0;

      // go thru the list of children and correct their
      // parent class descriptors
      for (Pix i=children_->first(); i; children_->next(i)) {
          (*all_)[children_->operator()(i)]->change_parent((*all_)[key],this);
        }

      delete (*all_)[key];
    }
  (*all_)[key] = this;
}

ClassDesc::~ClassDesc()
{
  delete[] classname_;
  if (children_) delete children_;
}

ClassDesc*
ClassDesc::name_to_class_desc(const char* name)
{
  ClassKey key(name);
  return (*all_)[key];
}

DescribedClass*
ClassDesc::create() const
{
  if (ctor_) return (*ctor_)();
  return 0;
}

DescribedClass*
ClassDesc::create(KeyVal&keyval) const
{
  DescribedClass* result;
  if (keyvalctor_) {
      result = (*keyvalctor_)(keyval);
    }
  else result = 0;
  return result;
}

DescribedClass*
ClassDesc::create(StateIn&statein) const
{
  if (stateinctor_) return (*stateinctor_)(statein);
  return 0;
}

void
ClassDesc::change_parent(ClassDesc*oldcd,ClassDesc*newcd)
{
  parents_.change_parent(oldcd,newcd);
}

void
ClassDesc::list_all_classes()
{
  printf("Listing all classes:\n");
  for (Pix ind=all_->first(); ind!=0; all_->next(ind)) {
      ClassDesc* classdesc = all_->contents(ind);
      printf("class %s\n",classdesc->name());
      ParentClasses& parents = classdesc->parents_;
      if (parents.n()) {
          printf("  parents:");
          for (int i=0; i<parents.n(); i++) {
              if (parents[i].is_virtual()) {
                  printf(" virtual");
                }
              if (parents[i].access() == ParentClass::Public) {
                  printf(" public");
                }
              else if (parents[i].access() == ParentClass::Protected) {
                  printf(" protected");
                }
              printf(" %s",parents[i].classdesc()->name());
            }
          printf("\n");
        }
      ClassKeySet* children = classdesc->children_;
      if (children) {
          printf("  children:");
          for (Pix pind=children->first(); pind!=0; children->next(pind)) {
              printf(" %s",(*children)(pind).name());
            }
          printf("\n");
        }
    }
}

const ParentClasses& ClassDesc::parents() const
{
  return parents_;
}

const char* ClassDesc::name() const
{
  return classname_;
}

int ClassDesc::version() const
{ return version_;
}

DescribedClass* ClassDesc::create_described_class() const
{
    return create();
}

////////////////////////////////////////////////////

ClassDesc DescribedClass::class_desc_("DescribedClass");

DescribedClass::DescribedClass() {}
DescribedClass::DescribedClass(const DescribedClass&) {}
DescribedClass& DescribedClass::operator=(const DescribedClass&)
{
  return *this;
}

DescribedClass::~DescribedClass()
{
}

const ClassDesc*
DescribedClass::class_desc() const
{
  return &class_desc_;
}

void*
DescribedClass::_castdown(const ClassDesc*cd)
{
  if (cd == &class_desc_) return this;
  return 0;
}

DescribedClass*
DescribedClass::castdown(DescribedClass*p)
{
  return (DescribedClass*) p->_castdown(DescribedClass::static_class_desc());
}

const ClassDesc*
DescribedClass::static_class_desc()
{
  return &class_desc_;
}

const char* DescribedClass::class_name() const
{
    return class_desc()->name();
}

int DescribedClass::class_version() const
{
    return class_desc()->version();
}

RefDescribedClassBase::RefDescribedClassBase()
{
}
RefDescribedClassBase::~RefDescribedClassBase()
{
}
int RefDescribedClassBase::operator==( const RefDescribedClassBase &a) const {
  return parentpointer() == a.parentpointer();
};
int RefDescribedClassBase::operator!=( const RefDescribedClassBase &a) const {
  return parentpointer() != a.parentpointer();
};
int RefDescribedClassBase::operator>=( const RefDescribedClassBase &a) const {
  return parentpointer() >= a.parentpointer();
};
int RefDescribedClassBase::operator<=( const RefDescribedClassBase &a) const {
  return parentpointer() <= a.parentpointer();
};
int RefDescribedClassBase::operator>( const RefDescribedClassBase &a) const {
  return parentpointer() > a.parentpointer();
};
int RefDescribedClassBase::operator<( const RefDescribedClassBase &a) const {
  return parentpointer() < a.parentpointer();
};

RefDescribedClassBase::RefDescribedClassBase(const RefDescribedClassBase&)
{
}
RefDescribedClassBase&
RefDescribedClassBase::operator=(const RefDescribedClassBase&)
{
  return *this;
}

DescribedClass_REF_def(DescribedClass);
ARRAY_def(RefDescribedClass);
SET_def(RefDescribedClass);
ARRAYSET_def(RefDescribedClass);

