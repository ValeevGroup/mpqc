//
// statein.cc
//
// Copyright (C) 1998 Limit Point Systems, Inc.
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

#include <ctype.h>
#include <util/misc/formio.h>
#include <util/state/translate.h>
#include <util/state/statein.h>

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
  have_cd_(0),
  translate_(new TranslateDataIn(this, new TranslateDataBigEndian)),
  expected_object_num_(0),
  nextclassid_(0),
  node_to_node_(0)
{
  key_[0] = '\0';
  keylength_ = 0;
}

StateIn::~StateIn()
{
  delete translate_;
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
StateIn::get_array_uint(unsigned int*p,int size)
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
StateIn::get(unsigned int&r, const char *keyword)
{
  int n = get_array_uint(&r,1);
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

int
StateIn::get_array_void(void*p,int s)
{
  cerr << "StateIn::get_array_void(void*p,int s) "
       << "is a derived class responsiblility" << endl
       << "  exact type is \"" << class_name() << "\"" << endl;
  abort();
  return -1;
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
#if DEBUG
      cout << "GET CLASS:"
           << " NAME = " << name
           << " VERSION = " << version
           << " ID = " << classid << endl;
#endif
      ClassDesc* tmp = ClassDesc::name_to_class_desc(name);

      classidmap_[tmp] = classid;
      StateClassData classdat(version,tmp,name);
      classdatamap_[classid] = classdat;
    }

  // read the object information
  get(length);
  for (i=0; i<length; i++) {
      int n;
      get(n);
      StateInData num;
      get(num.type);
      get(num.offset);
      get(num.size);
#if DEBUG
      cout << "GET OBJECT:"
           << " NUM=" << setw(2) << n
           << " OFF=" << setw(5) << num.offset
           << " SZ=" << setw(4) << num.size
           << " ID=" << setw(2) << num.type
           << " (" << classdatamap_[num.type].name << ")"
           << endl;
#endif
      ps_[n] = num;
      classdatamap_[num.type].ninstance++;
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

int
StateIn::getstring(char*&s)
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

int
StateIn::get(char*&s)
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

int
StateIn::get(unsigned int*&s)
{
  int r=0;
  int size;
  r += get(size);
  if (size) {
      s = new unsigned int[size];
      r += get_array_uint(s,size);
    }
  else {
      s = 0;
    }
  return r;
}

int
StateIn::get(int*&s)
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

int
StateIn::get(float*&s)
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

int
StateIn::get(double*&s)
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

int
StateIn::version(const ClassDesc* cd)
{
  int classid = classidmap_[(ClassDesc*)cd];
  return classdatamap_[classid].version;
}

int
StateIn::get(const ClassDesc**cd)
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
          classidmap_[tmp] = classid;
          StateClassData classdat(version,tmp,name);
          classdatamap_[classid] = classdat;
          r += get(size);
        }
    }

  // get the class id for the object
  int classid;
  r += get(classid);

  if (!classdatamap_.contains(classid)) {
      cerr << "ERROR: StateIn: couldn't find class descriptor for classid "
           << classid << endl;
      abort();
    }

  // convert the class id into the class descriptor
  *cd = classdatamap_[classid].classdesc;

  return r;
}

void
StateIn::list_objects(ostream &o)
{
  if (SCFormIO::getverbose(o)) {
      int ii = 1;
      for (AVLMap<int,StateInData>::iterator i=ps_.begin(); i!=ps_.end();
           i++,ii++) {
          const StateInData &num(i.data());
          const char *classname = classdatamap_[num.type].name;
          o << indent
            << "object " << ii
            << " at offset " << num.offset
            << " is of type " << classname
            << endl;
        }
    }
  else {
      int ntot = 0;
      for (AVLMap<int,StateClassData>::iterator i=classdatamap_.begin();
           i!=classdatamap_.end(); i++) {
          StateClassData &dat = i.data();
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

int
StateIn::dir_getobject(RefSavableState &p, const char *name)
{
  int r=0;

  p = 0;

  if (!use_directory()) {
      cerr << "ERROR: StateIn: no directory to get object from" << endl;
      abort();
    }

  // find the class name and/or object number
  const char *colon = strrchr(name,':');
  int number = 1;
  char *classname = 0;
  if (colon == 0) {
      if (isdigit(*name)) number = atoi(name);
      else classname = strcpy(new char[strlen(name)+1], name);
    }
  else {
      number = atoi(&colon[1]);
      classname = strcpy(new char[strlen(name)+1], name);
      *strrchr(classname,':') = '\0';
    }

  const ClassDesc *cd = 0;
  if (classname) {
      cd = ClassDesc::name_to_class_desc(classname);
      if (!cd) {
          cerr << "ERROR: StateIn: class " << classname << " unknown" << endl;
          abort();
        }
      delete[] classname;
    }

  int classid;
  if (cd) classid = classidmap_[(ClassDesc*)cd];
  else classid = -1;

  AVLMap<int,StateInData>::iterator i;
  int nfound = 0;
  for (i=ps_.begin(); i!=ps_.end(); i++) {
      if (classid == -1 || i.data().type == classid) nfound++;
      if (nfound == number) {
          if (i.data().ptr.nonnull()) {
              p = i.data().ptr;
            }
          else {
              seek(i.data().offset);
              r += getobject(p);
            }
          return r;
        }
    }

  return r;
}

int
StateIn::getobject(RefSavableState &p)
{
  int use_dir = use_directory();
  int r=0;
  int refnum;
  int original_loc=0;
  if (use_dir) original_loc = tell();
  int size_refnum;
  r += (size_refnum = get(refnum));
  if (refnum == 0) {
      // reference to null
#if DEBUG
      cout << indent << "getting null object" << endl;
#endif
      p = 0;
    }
  else {
#if DEBUG
      cout << indent << "getting object number " << setw(2) << refnum << endl;
      cout << incindent;
#endif
      AVLMap<int,StateInData>::iterator ind = ps_.find(refnum);
      if (ind == ps_.end() && use_dir) {
          cerr << "ERROR: StateIn: directory missing object number "
               << refnum << endl;
          abort();
        }
      if (ind == ps_.end() || ind.data().ptr.null()) {
#if DEBUG
          cout << indent << "reading object" << endl;
#endif
          // object has not yet been read in
          int need_seek = 0;
          if (use_dir) {
              if (original_loc != ind.data().offset) {
                  need_seek = 1;
                  original_loc = tell();
#if DEBUG
                  cout << indent << "seeking to"
                       << setw(5) << ind.data().offset << endl;
#endif
                  seek(ind.data().offset);
                  int trefnum;
                  get(trefnum);
                  if (trefnum != refnum) {
                      cerr << "StateIn: didn't find expected reference"<<endl;
                      abort();
                    }
                }
            }
          const ClassDesc *cd;
          r += get(&cd);
          have_classdesc();
          nextobject(refnum);
          DescribedClass *dc = cd->create(*this);
          p = SavableState::castdown(dc);
          if (use_dir) {
              ind.data().ptr = p;
              if (need_seek) seek(original_loc);
            }
#if DEBUG
          cout << indent << "got object with type = "
               << p->class_name() << endl;
#endif
        }
      else {
          // object already exists
          p = ind.data().ptr;
#if DEBUG
          cout << indent << "object already existed, type = "
               << p->class_name() << endl;
          cout << indent
               << "  use_dir = " << use_dir
               << " tell() = " << setw(5) << tell()
               << " offset = " << setw(5) << ind.data().offset
               << " size_refnum = " << setw(1) << size_refnum
               << endl;
#endif
          if (use_dir && tell() - size_refnum == ind.data().offset) {
              seek(tell() - size_refnum + ind.data().size);
#if DEBUG
              cout << indent << "  seeking to "
                   << tell() - size_refnum + ind.data().offset
                   << endl;
#endif
            }
        }
#if DEBUG
      cout << decindent;
#endif
    }
  return r;
}

void
StateIn::nextobject(int objnum)
{
  expected_object_num_ = objnum;
}

void
StateIn::haveobject(const RefSavableState &p)
{
  if (expected_object_num_) {
      haveobject(expected_object_num_,p);
      expected_object_num_ = 0;
    }
}

void
StateIn::haveobject(int objnum,const RefSavableState &p)
{
  AVLMap<int,StateInData>::iterator ind = ps_.find(objnum);
  if (ind == ps_.end()) {
      ps_[objnum].ptr = p;
#if DEBUG
      cout << indent << "have object adding number " << objnum << endl;
#endif
    }
  else {
      ind.data().ptr = p;
#if DEBUG
      cout << indent << "have object updating number " << objnum << endl;
#endif
    }
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
