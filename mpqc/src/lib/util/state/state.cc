
#include <stdlib.h>

#include <util/class/class.h>
#include "state.h"
#include "state_ptr.h"
#include "stateptrSet.h"
#include "statenumSet.h"
#include "stateptrImplSet.h"
#include "statenumImplSet.h"

#include "classdImplMap.h"

/////////////////////////////////////////////////////////////////

DescribedClass_REF_def(SavableState);

#define CLASSNAME SavableState
#define PARENTS virtual public DescribedClass
#include <util/class/classia.h>
void *
SavableState::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { DescribedClass::_castdown(cd) };
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

SavableState::SavableState(StateIn&si,const ClassDesc& cd)
{
  // In case si is looking for the next pointer, let it know i
  // have one.
  si.havepointer(this);
  si.get_version(&cd);
  // the order of the above two doesn't matter since when i are looking
  // for a pointer i already have the version info and the get_version
  // is ignored
}

SavableState::~SavableState()
{
}

void
SavableState::save_state(StateOut&so)
{
  // save the pointer
  if (so.putpointer(this)) {
      // The pointer hasn't been written yet, so write it.
      
      // Save the class descriptor for the exact type and all base classes
      // (base classes are needed to get the version information correct).
      so.put(class_desc());

      // save the object
      save_vbase_state(so);
      save_data_state(so);
    }
}

SavableState*
SavableState::restore_state(StateIn&si)
{
  // restore the pointer
  SavableState* ss;
  int objnum = si.getpointer((void**)&ss);
  if (objnum) {
      // The object doesn't yet exist.

      // Get the class descriptor (and all parent class descriptors are
      // retrieved too).
      const ClassDesc* cd;
      si.get(&cd);

      // Tell si that the next pointer corresponds to objnum.
      // (The SavableState(StateIn&) CTOR will call havepointer
      // with the this pointer.)
      si.nextobject(objnum);

      DescribedClass* dc = cd->create(si);
      //printf("dc = 0x%x\n",dc);
      ss = SavableState::castdown(dc);
    }
  //printf("ss = 0x%x\n",ss);
  return ss;
}

void
SavableState::save_object_state(StateOut&)
{
  fprintf(stderr,"SavableState::save_object_state(StateOut&):\n");
  fprintf(stderr," only can be used when exact type is known\n");
  fprintf(stderr," otherwise use save_state(StateOut&)\n");
  fprintf(stderr," (object not saved)\n");
}

void
SavableState::save_vbase_state(StateOut&)
{
}
void
SavableState::save_data_state(StateOut&)
{
}

/////////////////////////////////////////////////////////////////

DescribedClass_REF_def(StateOut);

#define CLASSNAME StateOut
#define PARENTS virtual public DescribedClass
#include <util/class/classia.h>
void *
StateOut::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { DescribedClass::_castdown(cd) };
  return do_castdowns(casts,cd);
}

StateOut::StateOut() :
  ps_(new StateDataPtr_CTOR),
  _classidmap(new ClassDescPintMap_CTOR),
  _nextclassid(0)
{
}

StateOut::StateOut(const StateOut&) {
    fprintf(stderr,"StateOut: private copy ctor called???\n");
    abort();
}
StateOut::operator=(const StateOut&) {
    fprintf(stderr,"StateOut: private assignment called???\n");
    abort();
}

StateOut::~StateOut()
{
  if(ps_) delete ps_; ps_=0;
}

int StateOut::put_array_char(const char*p,int size)
    { return put_array_void((const void*)p,size*sizeof(char)); }
int StateOut::put_array_int(const int*p,int size)
    { return put_array_void((const void*)p,size*sizeof(int)); }
int StateOut::put_array_float(const float*p,int size)
    { return put_array_void((const void*)p,size*sizeof(float)); }
int StateOut::put_array_double(const double*p,int size)
    { return put_array_void((const void*)p,size*sizeof(double)); }
int StateOut::put(char r) { return put_array_char((char*)&r,1); }
int StateOut::put(int r) { return put_array_int((int*)&r,1); }
int StateOut::put(float r) { return put_array_float((float*)&r,1); }
int StateOut::put(double r) { return put_array_double((double*)&r,1); }

/////////////////////////////////////////////////////////////////

DescribedClass_REF_def(StateIn);

#define CLASSNAME StateIn
#define PARENTS virtual public DescribedClass
#include <util/class/classia.h>
void *
StateIn::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { DescribedClass::_castdown(cd) };
  return do_castdowns(casts,cd);
}

StateIn::StateIn(const StateIn&) {
    fprintf(stderr,"StateIn: private copy ctor called???\n");
    abort();
}
StateIn::operator=(const StateIn&) {
    fprintf(stderr,"StateIn: private assignment called???\n");
    abort();
}

StateIn::StateIn() :
  ps_(new StateDataNum_CTOR),
  _nextobject(0)
{
}

StateIn::~StateIn()
{
  if(ps_) delete ps_; ps_=0;
}

int StateIn::get_array_char(char*p,int size)
{ return get_array_void((void*)p,size*sizeof(char)); }
int StateIn::get_array_int(int*p,int size)
{ return get_array_void((void*)p,size*sizeof(int)); }
int StateIn::get_array_float(float*p,int size)
{ return get_array_void((void*)p,size*sizeof(float)); }
int StateIn::get_array_double(double*p,int size)
{ return get_array_void((void*)p,size*sizeof(double)); }
int StateIn::get(char&r) { return get_array_char(&r,1); }
int StateIn::get(int&r) { return get_array_int(&r,1); }
int StateIn::get(float&r) { return get_array_float(&r,1); }
int StateIn::get(double&r) { return get_array_double(&r,1); }

/////////////////////////////////////////////////////////////////

// This deletes all references to objects, so if they are output
// again, they will be written in their entirety.
void StateOut::forget()
{
  ps_->clear();
}

void StateIn::forget()
{
  ps_->clear();
}

/////////////////////////////////////////////////////////////////

int StateOut::put_array_void(const void*p,int s)
{
  fprintf(stderr,"StateOut::put_array_void(const void*p,int s)"
    " is a derived class responsiblility\n");
  fprintf(stderr,"  exact type is \"%s\"\n",class_name());
  abort();
  return -1;
}

int StateIn::get_array_void(void*p,int s)
{
  fprintf(stderr,"StateIn::get_array_void(void*p,int s)"
    " is a derived class responsiblility\n");
  fprintf(stderr,"  exact type is \"%s\"\n",class_name());
  abort();
  return -1;
}

/////////////////////////////////////////////////////////////////

int StateOut::putstring(char*s)
{
  if (putpointer((void*)s)) {
      if (s) {
	  int size = strlen(s)+1;
	  put(size);
	  return put_array_void((void*)s,size);
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

/////////////////////////////////////////////////////////////////
// 
// int StateOut::put(SavableState*ss)
// {
//   if (ss) ss->save_state(*this);
//   else putpointer((void*)(DescribedClass*)ss);
//   return 0;
// }
// 
/////////////////////////////////////////////////////////////////
// 
// int StateOut::put(SavableState&ss)
// {
//   ss.save_member_state(*this);
//   return 0;
// }
// 
// int StateIn::get(SavableState&ss)
// {
//   ss.restore_member_state(*this);
//   return 0;
// }
// 
/////////////////////////////////////////////////////////////////

int StateIn::get_version(const ClassDesc*cd)
{
  if (!_cd.contains(cd)) {
      ClassDesc *tmp;
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
      //printf("just got \"%s\" %d\n",name,version);
      const ClassDesc* tmp = ClassDesc::name_to_class_desc(name);
      // save the class descriptor and the version
      _cd.add(tmp);
      _version.add(version);
      delete[] name;
      get(size);
    };

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
  ParentClasses& parents = cd->parents();

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

  
