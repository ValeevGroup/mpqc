
#include <builtin.h>

#include "symm.h"
#include "symmQCList.h"

///////////////////////////////////////////////////////////////////////

#define CLASSNAME SymmCoListLink
#define PARENTS virtual public SavableState
//#define HAVE_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
SymmCoListLink::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}

//SymmCoListLink::SymmCoListLink() : next(0), prev(0) {}

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
#define HAVE_KEYVAL_CTOR
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

SymmCoList::SymmCoList(KeyVal& kv):
  head(0),
  current(0)
{
  int nco = kv.count();
  for (int i = 0; i<nco; i++) {
      RefDescribedClass val = kv.describedclassvalue(i);
      RefSymmCo sc = val;
      if (val.nonnull() && sc.null()) {
          fprintf(stderr,"could not convert type %s to SymmCoList\n",
                  val->class_name());
          abort();
        }
      val = 0;
      add(sc.pointer());
    }
}

int SymmCoList::length()
{
  int j=0;
  for(current=head; current; j++,current=current->next );
  return j;
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
  if (l) p=l->head;
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
  list = b;
  if (list) p = list->head;
  else p = 0;
  return *this;
  }

RefSymmCo& SymmCoListIter::this_object()
{
  if(!p)
    err_quit("SymmCoListIter::this_object(): p is null");

  return p->obj;
  }

SymmCo * SymmCoListIter::operator->()
{
  if (!p || p->obj.null())
    err_quit("SymmCoListIter::operator->(): p or p->obj is null\n");

  return p->obj.pointer();
}

int SymmCoListIter::operator==(SymmCo& b)
{
  if(p && p->obj.pointer()) return (*p->obj.pointer()==b);
  else return 0;
  }
