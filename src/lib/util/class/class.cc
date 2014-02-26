//
// class.cc
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

#ifdef HAVE_CONFIG_H
#include <mpqc_config.h>
#endif

#include <string>
#include <vector>

#include <stdlib.h>
#include <string.h>

#ifdef HAVE_DLFCN_H
#include <dlfcn.h>
#endif // HAVE_DLFCN_H

#include <util/misc/formio.h>
#include <util/misc/scexception.h>

#include <util/class/class.h>
#include <util/class/proxy.h>

using namespace std;
using namespace sc;

std::map<std::string,ClassDescP>* ClassDesc::all_ = 0;
std::map<type_info_key,ClassDescP>* ClassDesc::type_info_all_ = 0;
char * ClassDesc::classlib_search_path_ = 0;
std::set<std::string>* ClassDesc::unresolved_parents_ = 0;

/////////////////////////////////////////////////////////////////

static sc::ClassDesc DescribedClassProxy_cd(
    typeid(sc::DescribedClassProxy), "DescribedClassProxy",1,"public DescribedClass");

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

const ClassDesc* ParentClass::classdesc() const
{
  return _classdesc;
}

void ParentClass::change_classdesc(ClassDesc*n)
{
  _classdesc = n;
}

/////////////////////////////////////////////////////////////////

ParentClasses::ParentClasses():
  _n(0),
  _classes(0)
{
}

void
ParentClasses::init(const char* parents)
{
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
          std::string parentkey(token);
          // if the parents class desc does not exist create a temporary
          // the temporary will be incorrect,because it does not have the
          // parent's parents
          if (ClassDesc::all().find(parentkey) == ClassDesc::all().end()) {
              ClassDesc *tmp_classdesc = new ClassDesc(token);
              ClassDesc::all()[parentkey] = tmp_classdesc;
              if (ClassDesc::unresolved_parents_ == 0) {
                  ClassDesc::unresolved_parents_ = new std::set<std::string>;
                }
              ClassDesc::unresolved_parents_->insert(token);
            }
          ParentClass* p = new ParentClass(ClassDesc::all()[parentkey],
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

////////////////////////////////////////////////////////////////////////

type_info_key&
type_info_key::operator=(const type_info_key&t)
{
  ti_ = t.ti_;
  return *this;
}

int
type_info_key::operator==(const type_info_key&t) const
{
  if (!ti_ && !t.ti_) return 1;
  if (!ti_ || !t.ti_) return 0;

  return *ti_ == *t.ti_;
}

int
type_info_key::operator<(const type_info_key&t) const
{
  if (!ti_ && !t.ti_) return 0;
  if (!ti_) return 0;
  if (!t.ti_) return 1;

  return ti_->before(*t.ti_);
}

int
type_info_key::cmp(const type_info_key&t) const
{
  if (*this == t) return 0;
  if (*this < t) return -1;
  return 1;
}

////////////////////////////////////////////////////////////////////////

ClassDesc::ClassDesc(const type_info &ti,
                     const char* name, int version,
                     const char* parents,
                     DescribedClass* (*ctor)(),
                     DescribedClass* (*keyvalctor)(const Ref<KeyVal>&),
                     DescribedClass* (*stateinctor)(StateIn&)
                     )
{
  typedef std::map<type_info_key,ClassDescP> class_registry;
  if (!type_info_all_) {
      type_info_all_ = new class_registry;
    }
  // look for duplicates, first by type_info ...
  type_info_key key(&ti);
  if (type_info_all_->find(key) != type_info_all_->end()) {
      std::stringstream oss;
      oss << "duplicate ClassDesc detected for class "
          << name << " type_info name = " << ti.name() << endl;
      throw sc::ProgrammingError(oss.str().c_str(), __FILE__, __LINE__);
  }
  else {
    // look for duplicates, second by name
    for(class_registry::iterator i = type_info_all_->begin();
        i != type_info_all_->end();
        ++i){
      if (strcmp(i->second->name(),name) == 0) {
        std::stringstream oss;
        oss << "duplicate ClassDesc detected for class named "
            << name << endl;
        oss << "ClassDesc being added            : type_info name = " << i->first.type_info()->name() << endl;
        oss << "ClassDesc already in the registry: type_info name = " << ti.name() << endl;
        throw sc::ProgrammingError(oss.str().c_str(), __FILE__, __LINE__);
      }
    }

    (*type_info_all_)[key] = this;
  }

  // test the version number to see if it is valid
  if (version <= 0)
      throw ProgrammingError("ERROR: ClassDesc ctor: version <= 0", __FILE__, __LINE__);

  init(name,version,parents,&ti,ctor,keyvalctor,stateinctor);
}

ClassDesc::ClassDesc(const char* name)
{
  init(name, 0);
}

void
ClassDesc::init(const char* name, int version,
                const char* parents,
                const type_info *ti,
                DescribedClass* (*ctor)(),
                DescribedClass* (*keyvalctor)(const Ref<KeyVal>&),
                DescribedClass* (*stateinctor)(StateIn&))
{
  classname_ = 0;
  version_ = version;
  children_ = 0;
  ctor_ = ctor;
  keyvalctor_ = keyvalctor;
  stateinctor_ = stateinctor;
  ti_ = ti;

  // make sure that the static members have been initialized
  if (!all_) {
      all_ = new std::map<std::string,ClassDescP>;
      const char* tmp = getenv("LD_LIBRARY_PATH");
      if (tmp) {
          // Needed for misbehaving getenv's.
          if (strncmp(tmp, "LD_LIBRARY_PATH=", 16) == 0) {
              tmp = ::strchr(tmp,'=');
              tmp++;
            }
        }
      else tmp = ".";
      classlib_search_path_ = ::strcpy(new char[strlen(tmp)+1],tmp);
    }
  
  // see if I'm already in the list
  ClassDesc *me = name_to_class_desc(name);
  int temp_copy_present = 0;
  if (me && me->version() != 0) {
    std::stringstream oss;
    oss << "duplicate ClassDesc detected for class named "
        << name << endl;
    throw sc::ProgrammingError(oss.str().c_str(), __FILE__, __LINE__);
  }
  else if (me) {
      temp_copy_present = 1;
    }

  parents_.init(parents);

  if (!temp_copy_present && name_to_class_desc(name)) {
      // I wasn't in the list before, but am in it now
    std::stringstream oss;
    oss << "inheritance loop detected for class "
        << name << endl;
    throw sc::ProgrammingError(oss.str().c_str(), __FILE__, __LINE__);
  }

  classname_ = ::strcpy(new char[strlen(name)+1],name);

  std::string key(name);

  // let each of the parents know that this is a child
  for (int i=0; i<parents_.n(); i++) {
      std::string parentkey(parents_[i].classdesc()->name());
      if (!(*all_)[parentkey]->children_)
        (*all_)[parentkey]->children_ = new std::set<std::string>;
      // let the parents know about the child
      ((*all_)[parentkey]->children_)->insert(key);
    }

  // if this class is aleady in all_, then it was put there by a child
  // preserve children info, destroy the old entry, and put this there
  if (all_->find(key) != all_->end()) {
      children_ = (*all_)[key]->children_;
      (*all_)[key]->children_ = 0;

      if (!children_) {
        std::ostringstream oss;
        oss << "ERROR: ClassDesc: inconsistency in initialization for "
            << key
            << "--perhaps a duplicated CTOR call" << endl;
        throw ProgrammingError(oss.str().c_str(), __FILE__, __LINE__);
      }

      // go thru the list of children and correct their
      // parent class descriptors
      for (std::set<std::string>::iterator i=children_->begin();
           i!=children_->end(); i++) {
          (*all_)[*i]->change_parent((*all_)[key],this);
        }

      delete (*all_)[key];
      unresolved_parents_->erase(key);
      if (unresolved_parents_->size() == 0) {
          delete unresolved_parents_;
          unresolved_parents_ = 0;
        }
    }
  (*all_)[key] = this;
}

ClassDesc::~ClassDesc()
{
  // remove references to this class descriptor
  if (children_) {
      for (std::set<std::string>::iterator i=children_->begin();
           i!=children_->end(); i++) {
          if (all_->find(*i) != all_->end()) {
              (*all_)[*i]->change_parent(this,0);
            }
        }
    }
  // delete this ClassDesc from the list of all ClassDesc's
  std::string key(classname_);
  all_->erase(key);

  // if the list of all ClassDesc's is empty, delete it
  if (all_->size() == 0) {
      delete all_;
      all_ = 0;
      delete[] classlib_search_path_;
      classlib_search_path_ = 0;
    }

  // delete this ClassDesc entry from the type_info map
  if (ti_ != 0) {
      type_info_key key(ti_);
      type_info_all_->erase(key);
      if (type_info_all_->size() == 0) {
          delete type_info_all_;
          type_info_all_ = 0;
        }
    }

  // delete local data
  delete[] classname_;
  if (children_) delete children_;
}

ClassDesc*
ClassDesc::class_desc(const type_info &ti)
{
  if (type_info_all_->find(type_info_key(&ti))
      == type_info_all_->end()) return 0;
  return (*type_info_all_)[type_info_key(&ti)];
}

std::map<std::string,ClassDescP>&
ClassDesc::all()
{
  if (!all_) throw sc::ProgrammingError("ClassDesc::all(): all not initialized", __FILE__, __LINE__);
  return *all_;
}

ClassDesc*
ClassDesc::name_to_class_desc(const char* name)
{
  std::string key(name);
  if (all_->find(key) == all_->end()) return 0;
  return (*all_)[key];
}

DescribedClass*
ClassDesc::create() const
{
  if (ctor_) return (*ctor_)();
  return 0;
}

DescribedClass*
ClassDesc::create(const Ref<KeyVal>&keyval) const
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
  ExEnv::out0() << "Listing all classes:" << endl;
  for (std::map<std::string,ClassDescP>::iterator ind=all_->begin();
       ind!=all_->end(); ind++) {
      ClassDesc* classdesc = ind->second;
      ExEnv::out0() << "class " << classdesc->name() << endl;
      ParentClasses& parents = classdesc->parents_;
      if (parents.n()) {
          ExEnv::out0() << "  parents:";
          for (int i=0; i<parents.n(); i++) {
              if (parents[i].is_virtual()) {
                  ExEnv::out0() << " virtual";
                }
              if (parents[i].access() == ParentClass::Public) {
                  ExEnv::out0() << " public";
                }
              else if (parents[i].access() == ParentClass::Protected) {
                  ExEnv::out0() << " protected";
                }
              if (parents[i].classdesc() == 0) {
                std::ostringstream oss;
                oss << "ERROR: parent " << i
                    << " for " << classdesc->name()
                    << " is missing" << endl;
                throw ProgrammingError(oss.str().c_str(), __FILE__, __LINE__);
              }
              const char *n = parents[i].classdesc()->name();
              ExEnv::out0() << " " << parents[i].classdesc()->name();
            }
          ExEnv::out0() << endl;
        }
      std::set<std::string>* children = classdesc->children_;
      if (children) {
          ExEnv::out0() << "  children:";
          for (std::set<std::string>::iterator pind=children->begin();
               pind!=children->end(); pind++) {
              ExEnv::out0() << " " << (*pind);
            }
          ExEnv::out0() << endl;
        }
    }
}

DescribedClass* ClassDesc::create_described_class() const
{
    return create();
}

// Returns 0 for success and -1 for failure.
int
ClassDesc::load_class(const char* classname)
{
  // See if the class has already been loaded.
  if (name_to_class_desc(classname) != 0) {
      return 0;
    }
  
#ifdef HAVE_DLFCN_H
  // make a copy of the library search list
  char* path = new char[strlen(classlib_search_path_) + 1];
  strcpy(path, classlib_search_path_);

  // go through each directory in the library search list
  char* dir = strtok(path,":");
  while (dir) {
      // find the 'classes' files
      char* filename = new char[strlen(dir) + 8 + 1];
      strcpy(filename,dir);
      strcat(filename,"/classes");
      ExEnv::outn() << "ClassDesc::load_class looking for \"" << filename << "\""
           << endl;
      FILE* fp = fopen(filename, "r");
      delete[] filename;

      if (fp) {
          // read the lines in the classes file
          const int bufsize = 10000;
          char buf[bufsize];
          while(fgets(buf, bufsize, fp)) {
              if (buf[0] != '\0' && buf[strlen(buf)-1] == '\n') {
                  buf[strlen(buf)-1] = '\0';
                }
              char* lib = strtok(buf," ");
              char* testclassname = strtok(0," ");
              ExEnv::outn() << "lib = \"" << lib << "\"" << endl;
              while(testclassname) {
                  ExEnv::outn() << "classname = \"" << testclassname << "\"" << endl;
                  if (strcmp(testclassname,classname) == 0) {
                      // found it
                      char* libname = new char[strlen(lib) + strlen(dir) + 2];
                      strcpy(libname, dir);
                      strcat(libname, "/");
                      strcat(libname, lib);
                      // load the libraries this lib depends upon

                      // i should look in the library's .dep file to
                      // get the dependencies, but this makes it a little
                      // difficult to make sure the same library doesn't
                      // get loaded twice (which is important) so for now
                      // i'll just wait until after i load the library and
                      // then look in the unresolved parents set
                      // and load parents until nothing is left

                      // load the library
                      ExEnv::outn() << "loading \"" << libname << "\"" << endl;
                      dlopen(libname, RTLD_LAZY);

                      // load code for parents
                      while (unresolved_parents_
                             && unresolved_parents_->size()) {
                          load_class((*unresolved_parents_->begin()).c_str());
                        }

                      fclose(fp);
                      delete[] path;
                      // make sure it worked.
                      if (name_to_class_desc(classname) == 0) {
                          ExEnv::errn() << "load of \"" << classname << "\" from \""
                               << libname << "\" failed" << endl;
                          delete[] libname;
                          return -1;
                        }
                      ExEnv::outn() << "loaded \"" << classname << "\" from \""
                           << libname << "\"" << endl;
                      delete[] libname;
                      return 0;
                    }
                  testclassname = strtok(0," ");
                }
            }
          fclose(fp);
        }
      
      dir = strtok(0, ":");
    }

  delete[] path;
#endif // HAVE_DLFCN_H

  ExEnv::outn() << "ClassDesc::load_class(\"" << classname << "\"): load failed!"
       << endl
       << "Possible causes:" << std::endl;

  // check if the name was simply misspelled
  int reason = 1;
  for (std::map<std::string,ClassDescP>::iterator ind=all_->begin();
       ind!=all_->end(); ind++) {
      ClassDesc* classdesc = ind->second;
      if (sc::string_distance(classname, classdesc->name()) <= std::size_t(4)) {
        ExEnv::outn() << "(" << reason << ") \"" << classname
            << "\" could be a misspelling of the known class name \""
            << classdesc->name() << "\"" << std::endl;
        ++reason;
      }
  }

  ExEnv::outn() << "(" << reason++ << ") the code for \"" << classname
                << "\" was not linked into the executable" << std::endl;
  ExEnv::outn() << "(" << reason++ << ") \"" << classname << "\" is an invalid class name" << std::endl;

  return -1;
}

////////////////////////////////////////////////////

static ClassDesc DescribedClass_cd(
    typeid(DescribedClass),"DescribedClass");

DescribedClass::DescribedClass()
{
}

DescribedClass::DescribedClass(const DescribedClass&) {}
DescribedClass& DescribedClass::operator=(const DescribedClass&)
{
  return *this;
}

DescribedClass::~DescribedClass()
{
}

ClassDesc*
DescribedClass::class_desc() const MPQC__NOEXCEPT
{
  ClassDesc *cd;
  try {
      cd = ClassDesc::class_desc(typeid(*this));
    }
  catch (...) {
      cd = 0;
    }
  return cd;
}

const char* DescribedClass::class_name() const
{
    return class_desc()->name();
}

int DescribedClass::class_version() const
{
    return class_desc()->version();
}

void
DescribedClass::print(ostream &o) const
{
  o << indent << "Object of type " << class_name() << endl;
}

ostream &
operator <<(ostream&o, const RefBase &ref)
{
  DescribedClass *dc = dynamic_cast<DescribedClass*>(ref.parentpointer());
  if (dc) {
      dc->print(o);
    }
  else {
      o << indent << "reference to null" << endl;
    }

  return o;
}

#ifdef EXPLICIT_TEMPLATE_INSTANTIATION

template class std::map<std::string,ClassDescP>;
               
template class std::map<std::string,int>;
template class std::set<std::string>;

#endif

/////////////////////////////////////////////////////////////////////////////

namespace {
  template <typename T>
  T
  min (const T& t1,
       const T& t2,
       const T& t3) {
    return std::min(std::min(t1,t2),t3);
  }
}

std::string::size_type
sc::string_distance(const std::string& str1,
                    const std::string& str2) {

  typedef std::string::size_type intsize;
  const intsize lenStr1 = str1.length();
  const intsize lenStr2 = str2.length();

  std::vector< std::vector<intsize> > d(lenStr1+1);
  for(intsize i=0; i<=lenStr1; ++i)
    d[i].resize(lenStr2+1);

  for(intsize i=0; i<=lenStr1; ++i)
    d[i][0] = i;
  for(intsize j=1; j<=lenStr2; ++j)
    d[0][j] = j;

  for(intsize i=1; i<=lenStr1; ++i) {
    for(intsize j=1; j<=lenStr2; ++j) {

      intsize cost;
      if (str1[i-1] == str2[j-1])
        cost = 0;
      else
        cost = 1;

      d[i][j] = min(
                   d[i-1][j] + 1,     // deletion
                   d[i][j-1] + 1,     // insertion
                   d[i-1][j-1] + cost   // substitution
                );

      if (i > 1 and j > 1 and str1[i] == str2[j-1] and str1[i-1] == str2[j])
        d[i][j] = std::min(
            d[i][j],
            d[i-2][j-2] + cost   // transposition
        );
    }
  }

  return d[lenStr1][lenStr2];
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
