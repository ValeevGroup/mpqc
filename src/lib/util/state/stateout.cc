//
// stateout.cc
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

#include <limits.h>

#include <mpqc_config.h>

#include <unistd.h>
#include <string.h>
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#ifdef HAVE_PWD_H
#include <pwd.h>
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif

#include <util/misc/scexception.h>
#include <util/state/translate.h>
#include <util/state/stateout.h>

using namespace std;
using namespace sc;

#define DEBUG 0

static ClassDesc StateOut_cd(
    typeid(StateOut),"StateOut",1,"public DescribedClass");

StateOut::StateOut() :
  have_cd_(0),
  translate_(new TranslateDataOut(this, new TranslateDataBigEndian)),
  copy_references_(0),
  next_object_number_(1),
  nextclassid_(0),
  node_to_node_(0)
{
}

StateOut::StateOut(const StateOut&)
{
  ExEnv::errn() << "StateOut: private copy ctor called???" << endl;
  abort();
}

void
StateOut::operator=(const StateOut&)
{
  ExEnv::errn() << "StateOut: private assignment called???" << endl;
  abort();
}

StateOut::~StateOut()
{
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
StateOut::put_array_uint(const unsigned int*p,int size)
{
  return translate_->put(p,size);
}

int
StateOut::put_array_int(const int*p,int size)
{
  return translate_->put(p,size);
}

int
StateOut::put_array_ulong(const unsigned long*p,int size)
{
  return translate_->put(p,size);
}

int
StateOut::put_array_long(const long*p,int size)
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

int StateOut::put(char r) { return put_array_char(&r,1); }
int StateOut::put(unsigned int r) { return put_array_uint(&r,1); }
int StateOut::put(bool r) { return put(int(r)); }
int StateOut::put(int r) { return put_array_int(&r,1); }
int StateOut::put(unsigned long r) { return put_array_ulong(&r,1); }
int StateOut::put(long r) { return put_array_long(&r,1); }
int StateOut::put(float r) { return put_array_float(&r,1); }
int StateOut::put(double r) { return put_array_double(&r,1); }

// This deletes all references to objects, so if they are output
// again, they will be written in their entirety.
void
StateOut::forget_references()
{
  ps_.clear();
}

// This deletes all references to objects, so if they are output
// again, they will be written in their entirety.  These also
// cause all future reference information to be ignored.  All
// referenced objects will be copied.
void
StateOut::copy_references()
{
  copy_references_ = 1;
}

int
StateOut::put_array_void(const void*p,int s)
{
  ExEnv::errn() << "StateOut::put_array_void(const void*p,int s) "
       << "is a derived class responsiblility" << endl
       << "  exact type is \"" << class_name() << "\"" << endl;
  abort();
  return -1;
}

void
StateOut::put_header()
{
  putstring("\001MPQCSO\002");

  // Switch to the native format and get_header will figure it out when read
  delete translate_;
  translate_ = new TranslateDataOut(this,new TranslateData);

  char format = translate_->translator()->format_code();
  put_array_char(&format,1);

  const int version = 1;
  put_array_int(&version,1);

  char userid[9];
  memset(userid,0,9);
#if defined(HAVE_GETPWUID) && defined(HAVE_GETEUID)
  const char *pw_name = getpwuid(geteuid())->pw_name;
  if (pw_name) {
      strncpy(userid, pw_name, 9);
      userid[8] = 0;
    }
#else
  strcpy(userid,"UNKNOWN");
#endif
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
StateOut::put_directory()
{
  std::map<ClassDescP,int>::iterator iid;
  std::map<Ref<SavableState>,StateOutData>::iterator isd;

  // write the type information
#if DEBUG
  ExEnv::outn() << "StateOut::put_directory(): directory length location = " << tell() << endl;
  ExEnv::outn() << "StateOut::put_directory(): directory length = " << classidmap_.size() << endl;
#endif
  put(classidmap_.size());
#if DEBUG
  ExEnv::outn() << "StateOut::put_directory(): directory entries location = " << tell() << endl;
#endif
  for (iid=classidmap_.begin(); iid!=classidmap_.end(); iid++) {
      const ClassDesc *cd = iid->first;
      int classid = iid->second;
      putstring(cd->name());
      put(cd->version());
      put(classid);
#if DEBUG
      ExEnv::outn() << "PUT CLASS:"
                   << " NAME = " << cd->name()
                   << " VERSION = " << cd->version()
                   << " ID = " << classid << endl;
#endif
    }

  // write the object information
  put(ps_.size());
  for (isd=ps_.begin(); isd!=ps_.end(); isd++) {
      const StateOutData& ptr = isd->second;
      put(ptr.num);
      put(ptr.type);
      put(ptr.offset);
      put(ptr.size);
#if DEBUG
      ExEnv::outn() << "PUT OBJECT:"
                   << " NUM = " << ptr.num
                   << " TYPE = " << ptr.type
                   << " OFFSET = " << ptr.offset
                   << " SIZE = " << ptr.size
                   << endl;
#endif
    }
}

int
StateOut::putstring(const char*s)
{
  int r=0;
  if (s) {
      int size = strlen(s)+1;
#if DEBUG
  ExEnv::outn() << "StateOut::putstring: string length location = " << tell() << endl;
  ExEnv::outn() << "StateOut::putstring: string length = " << size << endl;
#endif
      r += put(size);
#if DEBUG
      ExEnv::outn() << "StateOut::putstring: string location = " << tell() << endl;
#endif
      r += put_array_char(s,size-1);
    }
  else {
      r += put((int)0);
    }
  return r;
}

int
StateOut::put(const std::string &s)
{
  int r = putstring(s.c_str());
  return r;
}

int
StateOut::put(const char*s,int size)
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

int
StateOut::put(const unsigned int*s,int size)
{
  int r=0;
  if (s) {
      r += put(size);
      r += put_array_uint(s,size);
    }
  else {
      r += put((int)0);
    }
  return r;
}

int
StateOut::put(const int*s,int size)
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

int
StateOut::put(const unsigned long*s,int size)
{
  int r=0;
  if (s) {
      r += put(size);
      r += put_array_ulong(s,size);
    }
  else {
      r += put((int)0);
    }
  return r;
}

int
StateOut::put(const long*s,int size)
{
  int r=0;
  if (s) {
      r += put(size);
      r += put_array_long(s,size);
    }
  else {
      r += put((int)0);
    }
  return r;
}

int
StateOut::put(const float*s,int size)
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

int
StateOut::put(const double*s,int size)
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

int
StateOut::put(const ClassDesc*cd)
{
  int r=0;
  // write out parent info
  if (classidmap_.find((ClassDesc*)cd) == classidmap_.end()) {
      r += putparents(cd);
      if (!use_directory()) {
          const char* name = cd->name();
          int size = strlen(name);
          r += put(size);
          r += put_array_char(name,size);
          r += put(cd->version());
        }
      classidmap_[(ClassDesc*)cd] = nextclassid_++;
    }
  if (!use_directory()) {
      // write out a 0 to indicate the end of the list
      r += put((int)0);
    }
  // the cast is needed to de-const-ify cd
  r += put(classidmap_[(ClassDesc*)cd]);
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
      if (classidmap_.find(tmp) == classidmap_.end()) {
          r += putparents(tmp);
          if (!use_directory()) {
              const char* name = tmp->name();
              int size = strlen(name);
              r += put(size);
              r += put_array_char(name,size);
              r += put(tmp->version());
            }
          classidmap_[(ClassDesc*)tmp] = nextclassid_++;
        }
    }

  return r;
}

int
StateOut::putobject(const Ref<SavableState> &p)
{
  int r=0;
  if (p == 0) {
      // reference to null
      r += put(0);
    }
  else {
      std::map<Ref<SavableState>,StateOutData>::iterator ind = ps_.find(p);
      if (ind == ps_.end() || copy_references_) {
          StateOutData dp;
          // object has not been written yet
          dp.num = next_object_number_++;
          dp.offset = tell();
          r += put(dp.num);
          const ClassDesc *cd = p->class_desc();
          r += put(cd);
          dp.type = classidmap_[(ClassDesc*)cd];
          if (!copy_references_) ps_[p] = dp;
          have_classdesc();
          p->save_vbase_state(*this);
          p->save_data_state(*this);
          if (!copy_references_) {
              ind = ps_.find(p);
              ind->second.size = tell() - ind->second.offset;
            }
        }
      else {
          // object has already been written
          r += put(ind->second.num);
        }
    }
  return r;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
