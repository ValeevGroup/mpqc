
#include <builtin.h>

#include "simple.h"
#include "simpleQCList.h"

DescribedClass_IMPL(SimpleCoPtr,1,"","")
SavableState_IMPL(SimpleCoPtr)
void * SimpleCoPtr::_castdown(const ClassDesc *cd)
{
  if(&class_desc_ == cd) return this;
  return 0;
  }

SimpleCoPtr::SimpleCoPtr() : p(0) {}

SimpleCoPtr::SimpleCoPtr(SimpleCo *sc) : p(sc) { if(p) p->count++; }

SimpleCoPtr::SimpleCoPtr(SimpleCoPtr& sc) : p(sc.p) { if(p) p->count++; }

SimpleCoPtr::~SimpleCoPtr() { if (p && --(p->count) <= 0) delete p; }

SimpleCoPtr& SimpleCoPtr::operator=(SimpleCo* sc)
{
  if( p && --(p->count) <= 0 && p != sc) delete p;

  p=sc;
  if(p) p->count++;
  return *this;
  }

SimpleCoPtr& SimpleCoPtr::operator=(SimpleCoPtr& sc)
{
  if( p && --(p->count) <= 0 && p != sc.p) delete p;

  p=sc.p;
  if(p) p->count++;
  return *this;
  }

void SimpleCoPtr::save_data_state(StateOut& so)
{
  so.put(p);
  }

void SimpleCoPtr::restore_data_state(int v, StateIn& si)
{
  p = SimpleCo::restore_state(si);
  if(p) p->count++;
  }

///////////////////////////////////////////////////////////////////////

DescribedClass_IMPL(SimpleCoListLink,1,"","")
SavableState_IMPL(SimpleCoListLink)
void * SimpleCoListLink::_castdown(const ClassDesc *cd)
{
  if(&class_desc_ == cd) return this;
  return 0;
  }

SimpleCoListLink::SimpleCoListLink() : next(0), prev(0) {}

SimpleCoListLink::SimpleCoListLink(SimpleCo *o) : obj(o), next(0), prev(0) {}

SimpleCoListLink::SimpleCoListLink(SimpleCoPtr& o) : obj(o), next(0), prev(0) {}

SimpleCoListLink::SimpleCoListLink(const SimpleCoListLink& bll)
  : next(0), prev(0)
{
  *this=bll;
  }

SimpleCoListLink::~SimpleCoListLink()
{
  if(next) delete next; next=0; prev=0;
  }

int SimpleCoListLink::operator==(SimpleCoListLink& bl)
{
  if(obj==bl.obj) return 1;
  else if((void*)obj && !bl.obj || !obj && (void*)bl.obj) return 0;
  else if(obj) return(*obj == *bl.obj);
  return 1;
  }

void SimpleCoListLink::save_data_state(StateOut&so)
{
  so.put(obj); so.put(next); so.put(prev);
  }

void SimpleCoListLink::restore_data_state(int v,StateIn&si)
{
  si.get(obj);
  next = SimpleCoListLink::restore_state(si);
  prev = SimpleCoListLink::restore_state(si);
  }

void SimpleCoListLink::print(ostream& os, const char *pad) const
{
  char *mypad = new char[strlen(pad)+3]; sprintf(mypad,"%s  ",pad);
  if(obj.const_pointer()) obj->print(os,mypad);
  if(next) next->print(os,pad);
  delete[] mypad;
  os.flush();
  }

void SimpleCoListLink::print(FILE *of, const char *pad) const
{
  char *mypad = new char[strlen(pad)+3]; sprintf(mypad,"%s  ",pad);
  if(obj.const_pointer()) obj->print(of,mypad);
  if(next) next->print(of,pad);
  delete[] mypad;
  fflush(of);
  }

///////////////////////////////////////////////////////////////////////

DescribedClass_IMPL(SimpleCoList,1,"","")
SavableState_IMPL(SimpleCoList)
void * SimpleCoList::_castdown(const ClassDesc *cd)
{
  if(&class_desc_ == cd) return this;
  return 0;
  }

SimpleCoList::SimpleCoList()
  : current(0), head(0) {}

SimpleCoList::SimpleCoList(const SimpleCoList& bl)
  : current(0), head(0)
{
  *this = bl;
  }

SimpleCoList::SimpleCoList(SimpleCo *bo)
  : current(0), head(0)
{
  head = new SimpleCoListLink(bo);
  current=head;
  }

SimpleCoList::SimpleCoList(SimpleCoPtr&bo)
  : current(0), head(0)
{
  head = new SimpleCoListLink(bo);
  current=head;
  }

SimpleCoList::~SimpleCoList()
{
  current=0;
  if(head) delete head; head=0;
  }

// add object to tail end of list
void SimpleCoList::add(SimpleCo *bo)
{
  SimpleCoListLink *bll = new SimpleCoListLink(bo);

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

void SimpleCoList::add(SimpleCoPtr& bo)
{
  SimpleCoListLink *bll = new SimpleCoListLink(bo);

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
void SimpleCoList::remove_first(SimpleCoPtr &bo)
{
 // search through list for bo
  for(current=head; current; ) {
    if (*bo == *current->obj) {
      SimpleCoListLink *thislink=current;

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
void SimpleCoList::remove_all(SimpleCoPtr &bo)
{
 // search through list for bo
  for(current=head; current; ) {
    if (*bo == *current->obj) {
      SimpleCoListLink *thislink=current;

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

void SimpleCoList::save_data_state(StateOut& so)
{
  so.put(head);
  so.put(current);
  }

void SimpleCoList::restore_data_state(int v, StateIn& si)
{
  head = SimpleCoListLink::restore_state(si);
  current = SimpleCoListLink::restore_state(si);
  }

SimpleCoPtr& SimpleCoList::operator[](int i)
{
  int j=0;
  for(current=head; current; j++,current=current->next ) {
    if(i==j) return current->obj;
    }

  static SimpleCoPtr ret;
  return ret;
  }

void SimpleCoList::print(ostream& os, const char *pad) const
{
  os << pad << "SimpleCoList:\n";
  if(head) head->print(os,pad);
  os.flush();
  }

void SimpleCoList::print(FILE *of, const char *pad) const
{
  fprintf(of,"%sSimpleCoList:\n",pad);
  if(head) head->print(of,pad);
  fflush(of);
  }

///////////////////////////////////////////////////////////////////////

SimpleCoListIter::SimpleCoListIter() : p(0) {}

SimpleCoListIter::SimpleCoListIter(SimpleCoList* l)
  : list(l), p(0)
{
  p=l->head;
  }

SimpleCoListIter& SimpleCoListIter::operator=(int i)
{
  if(!list) return *this;

  int j=0;
  for(p=list->head; j!=i && p ; p=p->next)
    j++;
  return *this;
  }

SimpleCoListIter& SimpleCoListIter::operator=(SimpleCoList& b)
{
  list = &b; p = list->head;
  return *this;
  }

SimpleCoListIter& SimpleCoListIter::operator=(SimpleCoList *b)
{
  list = b; p = list->head;
  return *this;
  }

SimpleCoPtr& SimpleCoListIter::this_object()
{
  if(!p)
    err_quit("SimpleCoListIter::this_object(): p is null");

  return p->obj;
  }

int SimpleCoListIter::operator==(SimpleCo& b)
{
  if(p && p->obj.pointer()) return (*p->obj.pointer()==b);
  else return 0;
  }
