
#include <builtin.h>

#include "symm.h"
#include "symmQCList.h"

// DescribedClass_IMPL(SymmCoPtr,1,"","",SymmCoPtr::create,0,0,0,0)
// create_IMPL(SymmCoPtr)
// SavableState_IMPL(SymmCoPtr)
// void * SymmCoPtr::_castdown(const ClassDesc *cd)
// {
//   if(&class_desc_ == cd) return this;
//   return 0;
//   }
// 
// SymmCoPtr::SymmCoPtr() : p(0) {}
// 
// SymmCoPtr::SymmCoPtr(SymmCo *sc) : p(sc) { if(p) p->count++; }
// 
// SymmCoPtr::SymmCoPtr(SymmCoPtr& sc) : p(sc.p) { if(p) p->count++; }
// 
// SymmCoPtr::~SymmCoPtr() { if (p && --(p->count) <= 0) delete p; }
// 
// SymmCoPtr& SymmCoPtr::operator=(SymmCo* sc)
// {
//   if( p && --(p->count) <= 0 && p != sc) delete p;
// 
//   p=sc;
//   if(p) p->count++;
//   return *this;
//   }
// 
// SymmCoPtr& SymmCoPtr::operator=(SymmCoPtr& sc)
// {
//   if( p && --(p->count) <= 0 && p != sc.p) delete p;
// 
//   p=sc.p;
//   if(p) p->count++;
//   return *this;
//   }
// 
// void SymmCoPtr::save_data_state(StateOut& so)
// {
//   so.put(p);
//   }
// 
// void SymmCoPtr::restore_data_state(int v, StateIn& si)
// {
//   p = SymmCo::restore_state(si);
//   if(p) p->count++;
//   }

///////////////////////////////////////////////////////////////////////

#define CLASSNAME SymmCoListLink
#define PARENTS virtual public SavableState
#define HAVE_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
SymmCoListLink::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}

SymmCoListLink::SymmCoListLink() : next(0), prev(0) {}

SymmCoListLink::SymmCoListLink(SymmCo *o) : obj(o), next(0), prev(0) {}

SymmCoListLink::SymmCoListLink(RefSymmCo& o) : obj(o), next(0), prev(0) {}

SymmCoListLink::SymmCoListLink(const SymmCoListLink& bll)
  : next(0), prev(0)
{
  *this=bll;
  }

SymmCoListLink::~SymmCoListLink()
{
  if(next) delete next; next=0; prev=0;
  }

int SymmCoListLink::operator==(SymmCoListLink& bl)
{
  if(obj==bl.obj) return 1;
  else if((void*)obj && !bl.obj || !obj && (void*)bl.obj) return 0;
  else if(obj) return(*obj == *bl.obj);
  return 1;
  }

void SymmCoListLink::save_data_state(StateOut&so)
{
  obj->save_state(so);
  next->save_state(so);
  prev->save_state(so);
}

SymmCoListLink::SymmCoListLink(StateIn&si):
  SavableState(si,class_desc_)
{
  obj = SymmCo::restore_state(si);
  next = SymmCoListLink::restore_state(si);
  prev = SymmCoListLink::restore_state(si);
  }

void SymmCoListLink::print(ostream& os, const char *pad) const
{
  char *mypad = new char[strlen(pad)+3]; sprintf(mypad,"%s  ",pad);
  if(obj.pointer()) obj->print(os,mypad);
  if(next) next->print(os,pad);
  delete[] mypad;
  os.flush();
  }

void SymmCoListLink::print(FILE *of, const char *pad) const
{
  char *mypad = new char[strlen(pad)+3]; sprintf(mypad,"%s  ",pad);
  if(obj.pointer()) obj->print(of,mypad);
  if(next) next->print(of,pad);
  delete[] mypad;
  fflush(of);
  }

///////////////////////////////////////////////////////////////////////

DescribedClass_REF_def(SymmCoList);

#define CLASSNAME SymmCoList
#define PARENTS virtual public SavableState
#define HAVE_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
SymmCoList::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}

SymmCoList::SymmCoList()
  : current(0), head(0) {}

SymmCoList::SymmCoList(const SymmCoList& bl)
  : current(0), head(0)
{
  *this = bl;
  }

SymmCoList::SymmCoList(SymmCo *bo)
  : current(0), head(0)
{
  head = new SymmCoListLink(bo);
  current=head;
  }

SymmCoList::SymmCoList(RefSymmCo&bo)
  : current(0), head(0)
{
  head = new SymmCoListLink(bo);
  current=head;
  }

SymmCoList::~SymmCoList()
{
  current=0;
  if(head) delete head; head=0;
  }

// add object to tail end of list
void SymmCoList::add(SymmCo *bo)
{
  SymmCoListLink *bll = new SymmCoListLink(bo);

  if(!head)
    head=current=bll;
  else {
   // set current to end of list
    for(; current->next; current=current->next) ;
   // now tack bll onto current
    current->next = bll;
    bll->prev = current;
    }
  }

void SymmCoList::add(RefSymmCo& bo)
{
  SymmCoListLink *bll = new SymmCoListLink(bo);

  if(!head)
    head=current=bll;
  else {
   // set current to end of list
    for(; current->next; current=current->next) ;
   // now tack bll onto current
    current->next = bll;
    bll->prev = current;
    }
  }

// remove object from list, only finds first occurrence
void SymmCoList::remove_first(RefSymmCo &bo)
{
 // search through list for bo
  for(current=head; current; ) {
    if (*bo == *current->obj) {
      SymmCoListLink *thislink=current;

      if(!current->prev) { // this is head of list
        if(current->next) current->next->prev = 0;
        head = current->next;
        current = current->next;
        }
      else if(!current->next) { // this is tail of list
        if(current->prev) current->prev->next = current->next;
        current = 0;
        }
      else { // somewhere inbetween
        current->next->prev = current->prev;
        current->prev->next = current->next;
        current = current->next;
        }

      thislink->next = thislink->prev = 0;
      delete thislink;
      return;
      }
    else
      current=current->next;
    }
  }

// remove all occurences of object from list
void SymmCoList::remove_all(RefSymmCo &bo)
{
 // search through list for bo
  for(current=head; current; ) {
    if (*bo == *current->obj) {
      SymmCoListLink *thislink=current;

      if(!current->prev) { // this is head of list
        if(current->next) current->next->prev = 0;
        head = current->next;
        current = current->next;
        }
      else if(!current->next) { // this is tail of list
        if(current->prev) current->prev->next = current->next;
        current = 0;
        }
      else { // somewhere inbetween
        current->next->prev = current->prev;
        current->prev->next = current->next;
        current = current->next;
        }

      thislink->next = thislink->prev = 0;
      delete thislink;
      }
    else
      current=current->next;
    }
  }

void SymmCoList::save_data_state(StateOut& so)
{
  head->save_state(so);
  current->save_state(so);
}

SymmCoList::SymmCoList(StateIn& si):
  SavableState(si,class_desc_)
{
  head = SymmCoListLink::restore_state(si);
  current = SymmCoListLink::restore_state(si);
}

RefSymmCo& SymmCoList::operator[](int i)
{
  int j=0;
  for(current=head; current; j++,current=current->next ) {
    if(i==j) return current->obj;
    }

  static RefSymmCo ret;
  return ret;
  }

void SymmCoList::print(ostream& os, const char *pad) const
{
  os << pad << "SymmCoList:\n";
  if(head) head->print(os,pad);
  os.flush();
  }

void SymmCoList::print(FILE *of, const char *pad) const
{
  fprintf(of,"%sSymmCoList:\n",pad);
  if(head) head->print(of,pad);
  fflush(of);
  }

///////////////////////////////////////////////////////////////////////

SymmCoListIter::SymmCoListIter() : p(0) {}

SymmCoListIter::SymmCoListIter(SymmCoList* l)
  : list(l), p(0)
{
  p=l->head;
  }

SymmCoListIter& SymmCoListIter::operator=(int i)
{
  if(!list) return *this;

  int j=0;
  for(p=list->head; j!=i && p ; p=p->next)
    j++;
  return *this;
  }

SymmCoListIter& SymmCoListIter::operator=(SymmCoList& b)
{
  list = &b; p = list->head;
  return *this;
  }

SymmCoListIter& SymmCoListIter::operator=(SymmCoList *b)
{
  list = b; p = list->head;
  return *this;
  }

RefSymmCo& SymmCoListIter::this_object()
{
  if(!p)
    err_quit("SymmCoListIter::this_object(): p is null");

  return p->obj;
  }

int SymmCoListIter::operator==(SymmCo& b)
{
  if(p && p->obj.pointer()) return (*p->obj.pointer()==b);
  else return 0;
  }
