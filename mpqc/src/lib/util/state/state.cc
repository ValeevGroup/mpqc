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

#include <util/class/class.h>
#include <util/state/state.h>
#include <util/state/translate.h>
#include <util/state/state_ptr.h>
#include <util/state/stateptrSet.h>
#include <util/state/statenumSet.h>
#include <util/state/stateptrImplSet.h>
#include <util/state/statenumImplSet.h>

#include <util/state/classdImplMap.h>

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
  si.havepointer(this);

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
  // save the pointer
  if (so.putpointer(th)) {
      // The pointer hasn't been written yet, so write it.
      
      // Save the class descriptor for the exact type and all base classes
      // (base classes are needed to get the version information correct).
      so.put(th->class_desc());

      so.have_classdesc();
      
      // save the object
      th->save_vbase_state(so);
      th->save_data_state(so);
    }
}

SavableState*
SavableState::restore_state(StateIn&si)
{
  // restore the pointer
  SavableState* ss;
  int objnum = si.getpointer((void**)&ss);
  if (objnum) {
      // The object doesn't yet exist, so create it

      // Get the class descriptor (and all parent class descriptors are
      // retrieved too).
      const ClassDesc* cd;
      si.get(&cd);

      // Tell si that the next pointer corresponds to objnum.
      // (The SavableState(StateIn&) CTOR will call havepointer
      // with the this pointer.)
      si.nextobject(objnum);

      si.have_classdesc();
      DescribedClass* dc = cd->create(si);
      //cout << "dc = 0x" << setbase(16) << dc << endl;
      ss = SavableState::castdown(dc);
    }
  //cout << "ss = 0x" << setbase(16) << ss << endl;
  return ss;
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
  next_pointer_number(1),
  ps_(new StateDataPtr_CTOR),
  copy_references_(0),
  _classidmap(new ClassDescPintMap_CTOR),
  _nextclassid(0),
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
  delete _classidmap;
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
  _nextobject(0),
  have_cd_(0),
  node_to_node_(0),
  translate_(new TranslateDataIn(this, new TranslateDataBigEndian))
{
}

StateIn::~StateIn()
{
  delete ps_;
  delete translate_;
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

int StateIn::get(char&r) { return get_array_char(&r,1); }
int StateIn::get(int&r) { return get_array_int(&r,1); }
int StateIn::get(float&r) { return get_array_float(&r,1); }
int StateIn::get(double&r) { return get_array_double(&r,1); }

/////////////////////////////////////////////////////////////////

// This deletes all references to objects, so if they are output
// again, they will be written in their entirety.
void StateOut::forget_references()
{
  for (Pix i=ps_->first(); i; ps_->next(i)) {
      ps_->operator()(i).can_refer = 0;
    }
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

  char format;
  get_array_char(&format,1);
  // switch to the new format
  if (translate_->translator()->format_code() != format) {
      delete translate_;
      translate_ = new TranslateDataIn(this,TranslateData::vctor(format));
    }

  int version;
  get_array_int(&version,1);

  char userid[9];
  get_array_char(userid,9);

  int date;
  get_array_int(&date,1);
}

/////////////////////////////////////////////////////////////////

int StateOut::putstring(char*s)
{
  if (putpointer((void*)s)) {
      if (s) {
	  int size = strlen(s)+1;
	  put(size);
	  return put_array_char(s,size);
	}
      else {
	  put((int)0);
	}
    }
  return 0;
}

int StateIn::getstring(char*&s)
{
  int objnum;
  if (objnum = getpointer((void**)&s)) {
      int size;
      get(size);
      if (size) {
	  s = new char[size];
	  havepointer(objnum,s);
	  return get_array_char(s,size);
	}
      else {
	  s = 0;
	}
    }
  return 0;
}

/////////////////////////////////////////////////////////////////

int StateOut::put(char*s,int size)
{
  if (putpointer((void*)s)) {
    if (s) {
      put(size);
      return put_array_char(s,size);
      }
    else {
      put((int)0);
      }
    }
  return 0;
  }

int StateIn::get(char*&s)
{
  int objnum;
  if (objnum = getpointer((void**)&s)) {
    int size;
    get(size);
    if (size) {
      s = new char[size];
      havepointer(objnum,s);
      return get_array_char(s,size);
      }
    else {
      s = 0;
      }
    }
  return 0;
  }

/////////////////////////////////////////////////////////////////

int StateOut::put(int*s,int size)
{
  if (putpointer((void*)s)) {
    if (s) { put(size); return put_array_int(s,size); }
    else put((int)0);
    }
  return 0;
  }

int StateIn::get(int*&s)
{
  int objnum;
  if (objnum = getpointer((void**)&s)) {
    int size; get(size);
    if (size) {
      s = new int[size];
      havepointer(objnum,(void*)s);
      return get_array_int(s,size);
      }
    else s = 0;
    }
  return 0;
  }

/////////////////////////////////////////////////////////////////

int StateOut::put(float*s,int size)
{
  if (putpointer((void*)s)) {
    if (s) { put(size); return put_array_float(s,size); }
    else put((int)0);
    }
  return 0;
  }

int StateIn::get(float*&s)
{
  int objnum;
  if (objnum = getpointer((void**)&s)) {
    int size; get(size);
    if (size) {
      s = new float[size];
      havepointer(objnum,(void*)s);
      return get_array_float(s,size);
      }
    else s = 0;
    }
  return 0;
  }

/////////////////////////////////////////////////////////////////

int StateOut::put(double*s,int size)
{
  if (putpointer((void*)s)) {
    if (s) { put(size); return put_array_double(s,size); }
    else put((int)0);
    }
  return 0;
  }

int StateIn::get(double*&s)
{
  int objnum;
  if (objnum = getpointer((void**)&s)) {
    int size; get(size);
    if (size) {
      s = new double[size];
      havepointer(objnum,(void*)s);
      return get_array_double(s,size);
      }
    else s = 0;
    }
  return 0;
  }

int StateIn::version(const ClassDesc* cd)
{
  int position = _cd.iseek(cd);
  if (position >= 0) {
      return _version[position];
    }
  return -1;
}

int StateIn::get_version(const ClassDesc*cd)
{
  if (!_cd.contains((const ClassDescP&)cd)) {
      const ClassDesc *tmp;
      get(&tmp);
    }
  return 0;
}
int StateIn::get(const ClassDesc**cd)
{

  // if a list of class descriptors exist then read it in 
  int size;
  get(size);
  while (size) {
      char* name = new char[size+1];
      get_array_char(name,size);
      name[size] = '\0';
      int version;
      get(version);
      //cout << "just got \"" << name << "\" " << version << endl;
      ClassDesc* tmp = ClassDesc::name_to_class_desc(name);
      // save the class descriptor and the version
      _cd.add(tmp);
      int position = _cd.iseek(tmp);
      if (_version.length() <= position) {
          _version.reset_length(position + 10);
        }
      _version[position] = version;
      delete[] name;
      get(size);
    }

  // get the class id for the object
  int classid;
  get(classid);

  // convert the class id into the class descriptor
  *cd = _cd[classid];
  
  return 0;
  }

int StateOut::put_version(const ClassDesc*cd)
{
  if (!_classidmap->contains((ClassDesc*)cd)) {
      put(cd);
    }
  return 0;
}
int StateOut::put(const ClassDesc*cd)
{
  // write out parent info
  if (!_classidmap->contains((ClassDesc*)cd)) {
      putparents(cd);
      const char* name = cd->name();
      int size = strlen(name);
      put(size);
      put_array_char(name,size);
      put(cd->version());
      _classidmap->operator[]((ClassDesc*)cd) = _nextclassid++;
    }
  // write out a 0 to indicate the end of the list
  put((int)0);
  // the cast is needed to de-const-ify cd
  put(_classidmap->operator[]((ClassDesc*)cd));
  return 0;
  }

void
StateOut::putparents(const ClassDesc*cd)
{
  const ParentClasses& parents = cd->parents();

  for (int i=0; i<parents.n(); i++) {
      // the cast is needed to de-const-ify the class descriptor
      ClassDesc*tmp = (ClassDesc*) parents[i].classdesc();
      if (!_classidmap->contains(tmp)) {
          putparents(tmp);
          const char* name = tmp->name();
          int size = strlen(name);
          put(size);
          put_array_char(name,size);
          put(tmp->version());
          _classidmap->operator[](tmp) = _nextclassid++;
        }
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
SSRefBase::check_castdown_result(void* t, SavableState *ss)
{
  if (!t && ss) {
      cerr << "SSRef::restore_state() got type \"" << ss->class_name()
           << "\"" << endl;
        abort();
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
