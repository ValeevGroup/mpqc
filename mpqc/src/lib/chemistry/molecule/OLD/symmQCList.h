

#ifndef _SymmCoQCList_h
#define _SymmCoQCList_h

#include <stdio.h>
#include <iostream.h>

#include <util/unix/cct_cprot.h>
#include <util/class/class.h>
#include <util/state/state.h>

class SymmCo;

class SymmCoPtr
  : virtual public DescribedClass, virtual public SavableState {
DescribedClass_DECLARE(SymmCoPtr)
SavableState_DECLARE(SymmCoPtr)
  protected:
    SymmCo *p;

  public:
    SymmCoPtr();
    SymmCoPtr(SymmCo*);
    SymmCoPtr(SymmCoPtr&);
    virtual ~SymmCoPtr();

    SymmCoPtr& operator=(SymmCo*);
    SymmCoPtr& operator=(SymmCoPtr&);

    inline int operator==(const SymmCoPtr&a) const { return p == a.p; }
    inline int operator!=(const SymmCoPtr&a) const { return p != a.p; }
    inline int operator>=(const SymmCoPtr&a) const { return p >= a.p; }
    inline int operator<=(const SymmCoPtr&a) const { return p <= a.p; }
    inline int operator>(const SymmCoPtr&a) const { return p > a.p; }
    inline int operator<(const SymmCoPtr&a) const { return p < a.p; }

    inline int nrefs() const { if(p) return p->count; else return 0; }

    inline SymmCo* pointer() { return p; }
    inline SymmCo* const_pointer() const { return p; }

    inline int null() { return p == 0; }
    inline int nonnull() { return p != 0; }

    inline operator void*() { return p; }
    inline operator SavableState*() { return p; }
    inline operator DescribedClass*() { return p; }
    inline operator SymmCo*() { return p; }
    inline SymmCo* operator->() const { return p; }

    inline SymmCo& operator*() { return *p; }

    void save_data_state(StateOut&);
    void restore_data_state(int,StateIn&);
  };

class SymmCoList;
class SymmCoListIter;

class SymmCoListLink
  : virtual public DescribedClass, virtual public SavableState {
DescribedClass_DECLARE(SymmCoListLink)
SavableState_DECLARE(SymmCoListLink)
    friend class SymmCoList;
    friend class SymmCoListIter;
  protected:
    SymmCoPtr obj;
    SymmCoListLink *next;
    SymmCoListLink *prev;

    SymmCoListLink();
    SymmCoListLink(SymmCo*);
    SymmCoListLink(SymmCoPtr&);
    SymmCoListLink(const SymmCoListLink&);
    ~SymmCoListLink();
  public:
    void save_data_state(StateOut&);
    void restore_data_state(int,StateIn&);

    int operator==(SymmCoListLink&);
    inline int operator!=(SymmCoListLink&bl) { return !(*this == bl); }

    virtual void print(ostream&,const char* = " ") const;
    virtual void print(FILE* =stdout,const char* = " ") const;
  };

class SymmCoList
  : virtual public DescribedClass, virtual public SavableState {
DescribedClass_DECLARE(SymmCoList)
SavableState_DECLARE(SymmCoList)
    friend class SymmCoListIter;
  protected:
    SymmCoListLink *head;
    SymmCoListLink *current;
  public:
    SymmCoList();
    SymmCoList(SymmCo*);
    SymmCoList(SymmCoPtr&);
    SymmCoList(const SymmCoList&);
    virtual ~SymmCoList();

    void add(SymmCo*);
    void add(SymmCoPtr&);
    void remove_first(SymmCoPtr&);
    void remove_all(SymmCoPtr&);

    void save_data_state(StateOut&so);
    void restore_data_state(int,StateIn&);

    SymmCoPtr& operator[](int);

    virtual void print(ostream&,const char* = " ") const;
    virtual void print(FILE* =stdout,const char* = " ") const;
  };

class SymmCoListIter {
    friend class SymmCoList;
  protected:
    SymmCoList *list;
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

    inline SymmCo * operator->()
                        { if(p) return p->obj.pointer(); else return 0; }

    SymmCoPtr& this_object();

    int operator==(SymmCo&);
    inline int operator!=(SymmCo&b) { return !(*this == b); }
  };

#endif
