

#ifndef _SymmCoQCList_h
#define _SymmCoQCList_h

#include <stdio.h>
#include <iostream.h>

#include <util/unix/cct_cprot.h>
#include <util/class/class.h>
#include <util/state/state.h>

class SymmCo;

// class SymmCoPtr
//   : virtual public SavableState {
// #   define CLASSNAME A
// #   define HAVE_CTOR
// #   define HAVE_KEYVAL_CTOR
// #   define HAVE_STATEIN_CTOR
// #   include <util/state/stated.h>
// #   include <util/class/classd.h>
//   protected:
//     SymmCo *p;
// 
//   public:
//     SymmCoPtr();
//     SymmCoPtr(SymmCo*);
//     SymmCoPtr(SymmCoPtr&);
//     virtual ~SymmCoPtr();
// 
//     SymmCoPtr& operator=(SymmCo*);
//     SymmCoPtr& operator=(SymmCoPtr&);
// 
//     inline int operator==(const SymmCoPtr&a) const { return p == a.p; }
//     inline int operator!=(const SymmCoPtr&a) const { return p != a.p; }
//     inline int operator>=(const SymmCoPtr&a) const { return p >= a.p; }
//     inline int operator<=(const SymmCoPtr&a) const { return p <= a.p; }
//     inline int operator>(const SymmCoPtr&a) const { return p > a.p; }
//     inline int operator<(const SymmCoPtr&a) const { return p < a.p; }
// 
//     inline int nrefs() const { if(p) return p->count; else return 0; }
// 
//     inline SymmCo* pointer() { return p; }
//     inline SymmCo* const_pointer() const { return p; }
// 
//     inline int null() { return p == 0; }
//     inline int nonnull() { return p != 0; }
// 
//     inline operator void*() { return p; }
//     inline operator SavableState*() { return p; }
//     inline operator DescribedClass*() { return p; }
//     inline operator SymmCo*() { return p; }
//     inline SymmCo* operator->() const { return p; }
// 
//     inline SymmCo& operator*() { return *p; }
// 
//     void save_data_state(StateOut&);
//     void restore_data_state(int,StateIn&);
//   };

class SymmCoList;
class SymmCoListIter;

class SymmCoListLink
  : virtual public SavableState {
#   define CLASSNAME SymmCoListLink
//#   define HAVE_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
    friend class SymmCoList;
    friend class SymmCoListIter;
  protected:
    RefSymmCo obj;
    SymmCoListLink *next;
    SymmCoListLink *prev;

    //SymmCoListLink();
    //SymmCoListLink(SymmCo*);
    SymmCoListLink(RefSymmCo&);
    SymmCoListLink(const SymmCoListLink&);
    ~SymmCoListLink();
  public:
    void save_data_state(StateOut&);
    SymmCoListLink(StateIn&);

    int operator==(SymmCoListLink&);
    inline int operator!=(SymmCoListLink&bl) { return !(*this == bl); }

    virtual void print(ostream&,const char* = " ") const;
    virtual void print(FILE* =stdout,const char* = " ") const;
  };

class SymmCoList
  : virtual public SavableState {
#   define CLASSNAME SymmCoList
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
    friend class SymmCoListIter;
  protected:
    SymmCoListLink *head;
    SymmCoListLink *current;
  public:
    SymmCoList();
    SymmCoList(SymmCo*);
    SymmCoList(RefSymmCo&);
    SymmCoList(const SymmCoList&);
    virtual ~SymmCoList();

    //void add(SymmCo*);
    void add(RefSymmCo&);
    void remove_first(RefSymmCo&);
    void remove_all(RefSymmCo&);

    void save_data_state(StateOut&so);
    SymmCoList(StateIn&);
    SymmCoList(KeyVal&);

    int length();

    RefSymmCo& operator[](int);

    virtual void print(ostream&,const char* = " ") const;
    virtual void print(FILE* =stdout,const char* = " ") const;
  };
DescribedClass_REF_dec(SymmCoList);

class SymmCoListIter {
    friend class SymmCoList;
  protected:
    RefSymmCoList list;
    SymmCoListLink *p;
  public:
    SymmCoListIter();
    SymmCoListIter(SymmCoList*);

    inline operator void*() { return p; }
    inline operator DescribedClass*()
                        { if(p) return p->obj.pointer(); else return 0; }

    SymmCoListIter& operator=(int);
    SymmCoListIter& operator=(SymmCoList&);
    SymmCoListIter& operator=(SymmCoList*);

    inline SymmCoListIter& operator++() { if(p) p=p->next; return *this; }
    inline SymmCoListIter& operator--() { if(p) p=p->prev; return *this; }

    SymmCo * operator->();

    RefSymmCo& this_object();

    int operator==(SymmCo&);
    inline int operator!=(SymmCo&b) { return !(*this == b); }
  };

#endif
