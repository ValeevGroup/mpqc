

#ifndef _SimpleCoQCList_h
#define _SimpleCoQCList_h

#include <stdio.h>
#include <iostream.h>

#include <util/unix/cct_cprot.h>
#include <util/class/class.h>
#include <util/state/state.h>

class SimpleCo;

// class SimpleCoPtr
//   : virtual public SavableState {
// #   define CLASSNAME A
// #   define HAVE_CTOR
// #   define HAVE_KEYVAL_CTOR
// #   define HAVE_STATEIN_CTOR
// #   include <util/state/stated.h>
// #   include <util/class/classd.h>
//   protected:
//     SimpleCo *p;
// 
//   public:
//     SimpleCoPtr();
//     SimpleCoPtr(SimpleCo*);
//     SimpleCoPtr(SimpleCoPtr&);
//     virtual ~SimpleCoPtr();
// 
//     SimpleCoPtr& operator=(SimpleCo*);
//     SimpleCoPtr& operator=(SimpleCoPtr&);
// 
//     inline int operator==(const SimpleCoPtr&a) const { return p == a.p; }
//     inline int operator!=(const SimpleCoPtr&a) const { return p != a.p; }
//     inline int operator>=(const SimpleCoPtr&a) const { return p >= a.p; }
//     inline int operator<=(const SimpleCoPtr&a) const { return p <= a.p; }
//     inline int operator>(const SimpleCoPtr&a) const { return p > a.p; }
//     inline int operator<(const SimpleCoPtr&a) const { return p < a.p; }
// 
//     inline int nrefs() const { if(p) return p->count; else return 0; }
// 
//     inline SimpleCo* pointer() { return p; }
//     inline SimpleCo* const_pointer() const { return p; }
// 
//     inline int null() { return p == 0; }
//     inline int nonnull() { return p != 0; }
// 
//     inline operator void*() { return p; }
//     inline operator SavableState*() { return p; }
//     inline operator DescribedClass*() { return p; }
//     inline operator SimpleCo*() { return p; }
//     inline SimpleCo* operator->() const { return p; }
// 
//     inline SimpleCo& operator*() { return *p; }
// 
//     void save_data_state(StateOut&);
//     void restore_data_state(int,StateIn&);
//   };

class SimpleCoList;
class SimpleCoListIter;

class SimpleCoListLink
  : virtual public SavableState {
#   define CLASSNAME SimpleCoListLink
#   define HAVE_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
    friend class SimpleCoList;
    friend class SimpleCoListIter;
  protected:
    RefSimpleCo obj;
    SimpleCoListLink *next;
    SimpleCoListLink *prev;

    SimpleCoListLink();
    SimpleCoListLink(SimpleCo*);
    SimpleCoListLink(RefSimpleCo&);
    SimpleCoListLink(const SimpleCoListLink&);
    ~SimpleCoListLink();
  public:
    void save_data_state(StateOut&);
    SimpleCoListLink(StateIn&);

    int operator==(SimpleCoListLink&);
    inline int operator!=(SimpleCoListLink&bl) { return !(*this == bl); }

    virtual void print(ostream&,const char* = " ") const;
    virtual void print(FILE* =stdout,const char* = " ") const;
  };

class SimpleCoList
  : virtual public SavableState {
#   define CLASSNAME SimpleCoList
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
    friend class SimpleCoListIter;
  protected:
    SimpleCoListLink *head;
    SimpleCoListLink *current;
  public:
    SimpleCoList();
    SimpleCoList(SimpleCo*);
    SimpleCoList(RefSimpleCo&);
    SimpleCoList(KeyVal&);
    SimpleCoList(const SimpleCoList&);
    virtual ~SimpleCoList();

    void add(SimpleCo*);
    void add(RefSimpleCo&);
    void remove_first(RefSimpleCo&);
    void remove_all(RefSimpleCo&);

    void save_data_state(StateOut&so);
    SimpleCoList(StateIn&);

    RefSimpleCo& operator[](int);

    virtual void print(ostream&,const char* = " ") const;
    virtual void print(FILE* =stdout,const char* = " ") const;
  };

class SimpleCoListIter {
    friend class SimpleCoList;
  protected:
    SimpleCoList *list;
    SimpleCoListLink *p;
  public:
    SimpleCoListIter();
    SimpleCoListIter(SimpleCoList*);

    inline operator void*() { return p; }
    inline operator DescribedClass*()
                        { if(p) return p->obj.pointer(); else return 0; }

    SimpleCoListIter& operator=(int);
    SimpleCoListIter& operator=(SimpleCoList&);
    SimpleCoListIter& operator=(SimpleCoList*);

    inline SimpleCoListIter& operator++() { if(p) p=p->next; return *this; }
    inline SimpleCoListIter& operator--() { if(p) p=p->prev; return *this; }

    inline SimpleCo * operator->()
                        { if(p) return p->obj.pointer(); else return 0; }

    RefSimpleCo& this_object();

    int operator==(SimpleCo&);
    inline int operator!=(SimpleCo&b) { return !(*this == b); }
  };

#endif
