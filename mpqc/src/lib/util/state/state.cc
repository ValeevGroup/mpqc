//
// state.cc
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>

#include <util/misc/formio.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/state/translate.h>
#include <util/state/state_ptr.h>
#include <util/state/stateptrSet.h>
#include <util/state/statenumSet.h>
#include <util/state/stateptrImplSet.h>
#include <util/state/statenumImplSet.h>
#include <util/state/classdImplMap.h>
#include <util/state/classdatImplMap.h>

/////////////////////////////////////////////////////////////////

DescribedClass_REF_def(SavableState);

#define CLASSNAME SavableState
#define PARENTS public DescribedClass
#include <util/class/classia.h>
void *
SavableState::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

SavableState::SavableState()
{
}

SavableState::SavableState(const SavableState&)
{
}

SavableState& SavableState::operator=(const SavableState&)
{
  return *this;
}

SavableState::SavableState(StateIn&si)
{
  // In case si is looking for the next pointer, let it know i
  // have one.
  reference();
  RefSavableState th(this);
  si.haveobject(th);
  th = 0;
  dereference();

  // The following gets the version of this class and all of the
  // parent classes.  This is only needed for restoring objects
  // that were saved with save_object_state and don't necessarily
  // have all of their version information already restored.
  if (si.need_classdesc()) {
      const ClassDesc* tcd;
      si.get(&tcd);
    }
}

SavableState::~SavableState()
{
}

void
SavableState::save_state(StateOut&so)
{
  save_state(this,so);
}

void
SavableState::save_state(SavableState*th,StateOut&so)
{
  so.putobject(th);
}

SavableState*
SavableState::restore_state(StateIn& si)
{
  return dir_restore_state(si,0,0);
}

SavableState*
SavableState::key_restore_state(StateIn& si,
                                const char *keyword)
{
  return dir_restore_state(si,0,keyword);
}

SavableState*
SavableState::dir_restore_state(StateIn&si, const char *objectname,
                                const char *keyword)
{
  RefKeyVal old_override;
  RefSavableState overriding_value;
  int p = si.push_key(keyword);
  const int can_override_objects = 0;
  if (can_override_objects && keyword && si.override().nonnull()) {
      overriding_value = si.override()->describedclassvalue(si.key());
      old_override = si.override();
      if (overriding_value.nonnull()) {
          si.set_override(0);
        }
    }
  // restore the pointer
  RefSavableState ss;
  if (objectname) si.dir_getobject(ss, objectname);
  else si.getobject(ss);
  if (overriding_value.nonnull()) {
      cout << node0 << indent
           << "overriding \"" << si.key() << "\": object of type ";
      if (ss.null()) cout << node0 << "(null)";
      else cout << node0 << ss->class_name();
      cout << node0 << " -> object of type "
           << overriding_value->class_name()
           << endl;
      ss = overriding_value;
    }
  SavableState *ret = ss.pointer();
  if (ret) {
      ret->reference();
      ss = 0;
      ret->dereference();
    }
  if (old_override.nonnull()) {
      si.set_override(old_override);
    }
  si.pop_key(p);
  return ret;
}

void
SavableState::save_object_state_(StateOut&so, const ClassDesc *cd)
{
  if (class_desc() != cd) {
      cerr <<  "Warning:"
           << cd->name()
           << "::save_object_state: "
           << "exact type not known -- object not saved" << endl;
      return;
    }
  save_vbase_state(so);
  save_data_state(so);
}

void
SavableState::save_object_state(StateOut&)
{
  cerr << "SavableState::save_object_state(StateOut&):" << endl
       << " only can be used when exact type is known" << endl
       << " otherwise use save_state(StateOut&)" << endl
       << " (object not saved)" << endl;
}

void
SavableState::save_vbase_state(StateOut&so)
{
  SavableState::save_data_state(so);
}
void
SavableState::save_data_state(StateOut& so)
{
  if (so.need_classdesc()) so.put(class_desc());
}

/////////////////////////////////////////////////////////////////

DescribedClass_REF_def(StateOut);

#define CLASSNAME StateOut
#define PARENTS public DescribedClass
#include <util/class/classia.h>
void *
StateOut::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

StateOut::StateOut() :
  next_object_number_(1),
  ps_(new StateDataPtr_CTOR),
  copy_references_(0),
  classidmap_(new ClassDescPintMap_CTOR),
  nextclassid_(0),
  have_cd_(0),
  node_to_node_(0),
  translate_(new TranslateDataOut(this, new TranslateDataBigEndian))
{
}

StateOut::StateOut(const StateOut&)
{
  cerr << "StateOut: private copy ctor called???" << endl;
  abort();
}

void
StateOut::operator=(const StateOut&)
{
  cerr << "StateOut: private assignment called???" << endl;
  abort();
}

StateOut::~StateOut()
{
  delete ps_;
  delete classidmap_;
  delete translate_;
}

void
StateOut::flush()
{
}

int
StateOut::tell()
{
  return 0;
}

void
StateOut::seek(int loc)
{
}

int
StateOut::seekable()
{
  return 0;
}

int
StateOut::use_directory()
{
  return 0;
}

int
StateOut::put_array_char(const char*p,int size)
{
  return translate_->put(p,size);
}

int
StateOut::put_array_int(const int*p,int size)
{
  return translate_->put(p,size);
}

int
StateOut::put_array_float(const float*p,int size)
{
  return translate_->put(p,size);
}

int
StateOut::put_array_double(const double*p,int size)
{
  return translate_->put(p,size);
}

int StateOut::put(char r) { return put_array_char((char*)&r,1); }
int StateOut::put(int r) { return put_array_int((int*)&r,1); }
int StateOut::put(float r) { return put_array_float((float*)&r,1); }
int StateOut::put(double r) { return put_array_double((double*)&r,1); }

/////////////////////////////////////////////////////////////////

DescribedClass_REF_def(StateIn);

#define CLASSNAME StateIn
#define PARENTS public DescribedClass
#include <util/class/classia.h>
void *
StateIn::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

StateIn::StateIn(const StateIn&)
{
  cerr << "StateIn: private copy ctor called???" << endl;
  abort();
}

void
StateIn::operator=(const StateIn&)
{
  cerr << "StateIn: private assignment called???" << endl;
  abort();
}

StateIn::StateIn() :
  ps_(new StateDataNum_CTOR),
  classidmap_(new ClassDescPintMap_CTOR),
  classdatamap_(new intStateClassDataMap_CTOR),
  expected_object_num_(0),
  have_cd_(0),
  node_to_node_(0),
  nextclassid_(0),
  translate_(new TranslateDataIn(this, new TranslateDataBigEndian))
{
  key_[0] = '\0';
  keylength_ = 0;
}

StateIn::~StateIn()
{
  delete ps_;
  delete translate_;
  delete classidmap_;
  delete classdatamap_;
}

int
StateIn::push_key(const char *keyword)
{
  if (!keyword || override_.null()) return keylength_;

  int length = strlen(keyword);
  if (keylength_ + length + 1 >= KeyVal::MaxKeywordLength) {
      cerr << "StateIn: KeyVal::MaxKeywordLength exceeded" << endl;
      abort();
    }
  int old_keylength = keylength_;
  if (keylength_) key_[keylength_++] = ':';
  char *tmp = &key_[keylength_];
  for (int i=0; i<length; i++) tmp[i] = keyword[i];
  keylength_ += length;
  key_[keylength_] = '\0';

  return old_keylength;
}

int
StateIn::tell()
{
  return 0;
}

void
StateIn::seek(int loc)
{
}

int
StateIn::seekable()
{
  return 0;
}

int
StateIn::use_directory()
{
  return 0;
}

int
StateIn::get_array_char(char*p,int size)
{
  return translate_->get(p,size);
}

int
StateIn::get_array_int(int*p,int size)
{
  return translate_->get(p,size);
}

int
StateIn::get_array_float(float*p,int size)
{
  return translate_->get(p,size);
}

int
StateIn::get_array_double(double*p,int size)
{
  return translate_->get(p,size);
}

int
StateIn::get(char&r, const char *keyword)
{
  int n = get_array_char(&r,1);
  if (keyword && override().nonnull()) {
      int p = push_key(keyword);
      char roverride = override()->charvalue(key());
      if (override()->error() == KeyVal::OK) {
          cout << node0 << indent << "overriding \"" << key()
               << "\": " << r << " -> " << roverride << endl;
          r = roverride;
        }
      pop_key(p);
    }
  return n;
}

int
StateIn::get(int&r, const char *keyword)
{
  int n = get_array_int(&r,1);
  if (keyword && override().nonnull()) {
      int p = push_key(keyword);
      int roverride = override()->intvalue(key());
      if (override()->error() == KeyVal::OK) {
          cout << node0 << indent << "overriding \"" << key()
               << "\": " << r << " -> " << roverride << endl;
          r = roverride;
        }
      pop_key(p);
    }
  return n;
}

int
StateIn::get(float&r, const char *keyword)
{
  int n = get_array_float(&r,1);
  if (keyword && override().nonnull()) {
      int p = push_key(keyword);
      float roverride = override()->floatvalue(key());
      if (override()->error() == KeyVal::OK) {
          cout << node0 << indent << "overriding \"" << key()
               << "\": " << r << " -> " << roverride << endl;
          r = roverride;
        }
      pop_key(p);
    }
  return n;
}

int
StateIn::get(double&r, const char *keyword)
{
  int n = get_array_double(&r,1);
  if (keyword && override().nonnull()) {
      int p = push_key(keyword);
      double roverride = override()->doublevalue(key());
      if (override()->error() == KeyVal::OK) {
          cout << node0 << indent << "overriding \"" << key()
               << "\": " << r << " -> " << roverride << endl;
          r = roverride;
        }
      pop_key(p);
    }
  return n;
}

/////////////////////////////////////////////////////////////////

// This deletes all references to objects, so if they are output
// again, they will be written in their entirety.
void StateOut::forget_references()
{
  ps_->clear();
}

// This deletes all references to objects, so if they are output
// again, they will be written in their entirety.  These also
// cause all future reference information to be ignored.  All
// referenced objects will be copied.
void StateOut::copy_references()
{
  copy_references_ = 1;
}

/////////////////////////////////////////////////////////////////

int StateOut::put_array_void(const void*p,int s)
{
  cerr << "StateOut::put_array_void(const void*p,int s) "
       << "is a derived class responsiblility" << endl
       << "  exact type is \"" << class_name() << "\"" << endl;
  abort();
  return -1;
}

int StateIn::get_array_void(void*p,int s)
{
  cerr << "StateIn::get_array_void(void*p,int s) "
       << "is a derived class responsiblility" << endl
       << "  exact type is \"" << class_name() << "\"" << endl;
  abort();
  return -1;
}

/////////////////////////////////////////////////////////////////

void
StateOut::put_header()
{
  const char *magic = "\001SCSO\002";
  put_array_char(magic,6);

  // Switch to the native format and get_header will figure it out when read
  delete translate_;
  translate_ = new TranslateDataOut(this,new TranslateData);

  char format = translate_->translator()->format_code();
  put_array_char(&format,1);

  const int version = 1;
  put_array_int(&version,1);

  char userid[L_cuserid+9];
  memset(userid, 0, 9);
  cuserid(userid);
  userid[8] = 0;
  put_array_char(userid,9);

  timeval tv;
  gettimeofday(&tv,0);
  int date = (int) tv.tv_sec;
  put_array_int(&date,1);

  // record the position of the directory locator
  dir_loc_loc_ = tell();

  // the directory location defaults to 0 (no directory)
  int dir_loc = 0;
  // however, if a directory is to be used make dir_loc -1 (invalid)
  if (use_directory()) dir_loc = -1;
  put_array_int(&dir_loc,1);
}

void
StateIn::get_header()
{
  const char *magic = "\001SCSO\002";
  char tmp[7];
  get_array_char(tmp,6);
  tmp[6] = '\0';
  if (strcmp(tmp,magic)) {
      cerr << "StateIn: bad magic number" << endl;
      abort();
    }

  get_array_char(&format_,1);
  // switch to the new format
  if (translate_->translator()->format_code() != format_) {
      delete translate_;
      translate_ = new TranslateDataIn(this,TranslateData::vctor(format_));
    }

  get_array_int(&version_,1);

  get_array_char(userid_,9);

  get_array_int(&date_,1);

  // get the directory location
  get_array_int(&dir_loc_,1);
  if (dir_loc_ == -1) {
      cerr << "ERROR: StateIn: directory corrupted" << endl;
      abort();
    }
}

/////////////////////////////////////////////////////////////////

void
StateOut::put_directory()
{
  Pix i;

  // write the class information
  put(classidmap_->length());
  for (i=classidmap_->first(); i; classidmap_->next(i)) {
      const ClassDesc *cd = classidmap_->key(i);
      int classid = classidmap_->contents(i);
      putstring(cd->name());
      put(cd->version());
      put(classid);
      //cout << "PUT CLASS:"
      //     << " NAME = " << cd->name()
      //     << " VERSION = " << cd->version()
      //     << " ID = " << classid << endl;
    }

  // write the object information
  put(ps_->length());
  for (i=ps_->first(); i; ps_->next(i)) {
      const StateDataPtr& ptr = ps_->operator()(i);
      put(ptr.num);
      put(ptr.type);
      put(ptr.offset);
      put(ptr.size);
      //cout << "PUT OBJECT:"
      //     << " NUM = " << ptr.num
      //     << " TYPE = " << ptr.type
      //     << " OFFSET = " << ptr.offset
      //     << " SIZE = " << ptr.size
      //     << endl;
    }
}

void
StateIn::get_directory()
{
  int i, length;

  // read the type information
  get(length);
  for (i=0; i<length; i++) {
      char *name;
      int version, classid;
      getstring(name);
      get(version);
      get(classid);
      //cout << "GET CLASS:"
      //     << " NAME = " << name
      //     << " VERSION = " << version
      //     << " ID = " << classid << endl;
      ClassDesc* tmp = ClassDesc::name_to_class_desc(name);

      classidmap_->operator[](tmp) = classid;
      StateClassData classdat(version,tmp,name);
      classdatamap_->operator[](classid) = classdat;
    }

  // read the object information
  get(length);
  for (i=0; i<length; i++) {
      int n;
      get(n);
      StateDataNum num(n);
      get(num.type);
      get(num.offset);
      get(num.size);
      //cout << "GET OBJECT:"
      //     << " NUM = " << num.num
      //     << " TYPE = " << num.type
      //     << " OFFSET = " << num.offset
      //     << " SIZE = " << num.size
      //     << endl;
      ps_->add(num);
      classdatamap_->operator[](num.type).ninstance++;
    }
}

void
StateIn::find_and_get_directory()
{
  if (directory_location()) {
      int original_loc = tell();
      seek(directory_location());
      get_directory();
      seek(original_loc);
    }
}

/////////////////////////////////////////////////////////////////

int StateOut::putstring(const char*s)
{
  int r=0;
  if (s) {
      int size = strlen(s)+1;
      r += put(size);
      r += put_array_char(s,size-1);
    }
  else {
      r += put((int)0);
    }
  return r;
}

int StateIn::getstring(char*&s)
{
  int r=0;
  int size;
  r += get(size);
  if (size) {
      s = new char[size];
      r += get_array_char(s,size-1);
      s[size-1] = '\0';
    }
  else {
      s = 0;
    }
  return r;
}

/////////////////////////////////////////////////////////////////

int StateOut::put(const char*s,int size)
{
  int r=0;
  if (s) {
      r += put(size);
      r += put_array_char(s,size);
    }
  else {
      r += put((int)0);
    }
  return r;
}

int StateIn::get(char*&s)
{
  int r=0;
  int size;
  r += get(size);
  if (size) {
      s = new char[size];
      r += get_array_char(s,size);
    }
  else {
      s = 0;
    }
  return r;
}

/////////////////////////////////////////////////////////////////

int StateOut::put(const int*s,int size)
{
  int r=0;
  if (s) {
      r += put(size);
      r += put_array_int(s,size);
    }
  else {
      r += put((int)0);
    }
  return r;
}

int StateIn::get(int*&s)
{
  int r=0;
  int size;
  r += get(size);
  if (size) {
      s = new int[size];
      r += get_array_int(s,size);
    }
  else {
      s = 0;
    }
  return r;
}

/////////////////////////////////////////////////////////////////

int StateOut::put(const float*s,int size)
{
  int r=0;
  if (s) {
      r += put(size);
      r += put_array_float(s,size);
    }
  else {
      r += put((int)0);
    }
  return r;
}

int StateIn::get(float*&s)
{
  int r=0;
  int size;
  r += get(size);
  if (size) {
      s = new float[size];
      r += get_array_float(s,size);
    }
  else {
      s = 0;
    }
  return r;
}

/////////////////////////////////////////////////////////////////

int StateOut::put(const double*s,int size)
{
  int r=0;
  if (s) {
      r += put(size);
      r += put_array_double(s,size);
    }
  else {
      r += put((int)0);
    }
  return r;
}

int StateIn::get(double*&s)
{
  int r=0;
  int size;
  r += get(size);
  if (size) {
      s = new double[size];
      r += get_array_double(s,size);
    }
  else {
      s = 0;
    }
  return r;
}

/////////////////////////////////////////////////////////////////

int StateIn::version(const ClassDesc* cd)
{
  int classid = classidmap_->operator[]((ClassDesc*)cd);
  return classdatamap_->operator[](classid).version;
}

int StateIn::get(const ClassDesc**cd)
{
  int r=0;

  // if a list of class descriptors exist then read it in 
  if (!use_directory()) {
      int size;
      r += get(size);
      while (size) {
          char* name = new char[size+1];
          r += get_array_char(name,size);
          name[size] = '\0';
          int version;
          r += get(version);
          ClassDesc* tmp = ClassDesc::name_to_class_desc(name);
          // save the class descriptor and the version
          int classid = nextclassid_++;
          classidmap_->operator[]((ClassDesc*)tmp) = classid;
          StateClassData classdat(version,tmp,name);
          classdatamap_->operator[](classid) = classdat;
          r += get(size);
        }
    }

  // get the class id for the object
  int classid;
  r += get(classid);

  Pix ind = classdatamap_->seek(classid);
  if (ind == 0) {
      cerr << "ERROR: StateIn: couldn't find class descriptor for classid "
           << classid << endl;
      abort();
    }

  // convert the class id into the class descriptor
  *cd = classdatamap_->contents(ind).classdesc;

  return r;
}

int StateOut::put(const ClassDesc*cd)
{
  int r=0;
  // write out parent info
  if (!classidmap_->contains((ClassDesc*)cd)) {
      r += putparents(cd);
      if (!use_directory()) {
          const char* name = cd->name();
          int size = strlen(name);
          r += put(size);
          r += put_array_char(name,size);
          r += put(cd->version());
        }
      classidmap_->operator[]((ClassDesc*)cd) = nextclassid_++;
    }
  if (!use_directory()) {
      // write out a 0 to indicate the end of the list
      r += put((int)0);
    }
  // the cast is needed to de-const-ify cd
  r += put(classidmap_->operator[]((ClassDesc*)cd));
  return r;
  }

int
StateOut::putparents(const ClassDesc*cd)
{
  int r=0;
  const ParentClasses& parents = cd->parents();

  for (int i=0; i<parents.n(); i++) {
      // the cast is needed to de-const-ify the class descriptor
      ClassDesc*tmp = (ClassDesc*) parents[i].classdesc();
      if (!classidmap_->contains(tmp)) {
          r += putparents(tmp);
          if (!use_directory()) {
              const char* name = tmp->name();
              int size = strlen(name);
              r += put(size);
              r += put_array_char(name,size);
              r += put(tmp->version());
            }
          classidmap_->operator[]((ClassDesc*)tmp) = nextclassid_++;
        }
    }

  return r;
}

////////////////////////////////////////////////////////////////////////

void
StateIn::list_objects(ostream &o)
{
  Pix i;
  if (SCFormIO::getverbose(o)) {
      int ii = 1;
      for (i=ps_->first(); i; ps_->next(i),ii++) {
          StateDataNum &num(ps_->operator()(i));
          Pix dati = classdatamap_->seek(num.type);
          const char *classname = classdatamap_->contents(dati).name;
          o << indent
            << "object " << ii
            << " at offset " << num.offset
            << " is of type " << classname
            << endl;
        }
    }
  else {
      int ntot = 0;
      for (i=classdatamap_->first(); i; classdatamap_->next(i)) {
          StateClassData &dat = classdatamap_->contents(i);
          if (dat.ninstance > 0) {
              o << indent << dat.ninstance
                << " "
                << dat.name
                << endl;
              ntot += dat.ninstance;
            }
        }
      o << indent << "total of " << ntot << endl;
    }
}

////////////////////////////////////////////////////////////////////////
// SSRefBase members

void
SSRefBase::save_data_state(StateOut&s)
{
  save_state(s);
}

void
SSRefBase::save_state(StateOut&so)
{
  SavableState::save_state(sspointer(),so);
}

void
SSRefBase::check_castdown_result(void* t, SavableState *ss,
                                 const ClassDesc *cd)
{
  if (!t && ss) {
      cerr << node0
           << "SSRef::restore_state() got type \"" << ss->class_name()
           << "\""
           << " but expected \""
           << cd->name()
           << "\""
           << endl;
        abort();
    }
}

/////////////////////////////////////////////////////////////////////////////

StateClassData::~StateClassData()
{
  delete[] name;
}

StateClassData &
StateClassData::operator=(const StateClassData &d)
{
  version = d.version;
  classdesc = d.classdesc;
  ninstance = d.ninstance;
  if (d.name) name = strcpy(new char[strlen(d.name)+1], d.name);
  else name = 0;
  return *this;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
