

#ifndef _SimpleCoQCList_h
#define _SimpleCoQCList_h

#include <stdio.h>
#include <iostream.h>

#include <util/unix/cct_cprot.h>
#include <util/class/class.h>
#include <util/state/state.h>

class SimpleCo;

class SimpleCoPtr
  : virtual public DescribedClass, virtual public SavableState {
DescribedClass_DECLARE(SimpleCoPtr)
SavableState_DECLARE(SimpleCoPtr)
  protected:
    SimpleCo *p;

  public:
    SimpleCoPtr();
    SimpleCoPtr(SimpleCo*);
    SimpleCoPtr(SimpleCoPtr&);
    virtual ~SimpleCoPtr();

    SimpleCoPtr& operator=(SimpleCo*);
    SimpleCoPtr& operator=(SimpleCoPtr&);

    inline int operator==(const SimpleCoPtr&a) const { return p == a.p; }
    inline int operator!=(const SimpleCoPtr&a) const { return p != a.p; }
    inline int operator>=(const SimpleCoPtr&a) const { return p >= a.p; }
    inline int operator<=(const SimpleCoPtr&a) const { return p <= a.p; }
    inline int operator>(const SimpleCoPtr&a) const { return p > a.p; }
    inline int operator<(const SimpleCoPtr&a) const { return p < a.p; }

    inline int nrefs() const { if(p) return p->count; else return 0; }

    inline SimpleCo* pointer() { return p; }
    inline SimpleCo* const_pointer() const { return p; }

    inline int null() { return p == 0; }
    inline int nonnull() { return p != 0; }

    inline operator void*() { return p; }
    inline operator SavableState*() { return p; }
    inline operator DescribedClass*() { return p; }
    inline operator SimpleCo*() { return p; }
    inline SimpleCo* operator->() const { return p; }

    inline SimpleCo& operator*() { return *p; }

    void save_data_state(StateOut&);
    void restore_data_state(int,StateIn&);
  };

class SimpleCoList;
class SimpleCoListIter;

class SimpleCoListLink
  : virtual public DescribedClass, virtual public SavableState {
DescribedClass_DECLARE(SimpleCoListLink)
SavableState_DECLARE(SimpleCoListLink)
    friend class SimpleCoList;
    friend class SimpleCoListIter;
  protected:
    SimpleCoPtr obj;
    SimpleCoListLink *next;
    SimpleCoListLink *prev;

    SimpleCoListLink();
    SimpleCoListLink(SimpleCo*);
    SimpleCoListLink(SimpleCoPtr&);
    SimpleCoListLink(const SimpleCoListLink&);
    ~SimpleCoListLink();
  public:
    void save_data_state(StateOut&);
    void restore_data_state(int,StateIn&);

    int operator==(SimpleCoListLink&);
    inline int operator!=(SimpleCoListLink&bl) { return !(*this == bl); }

    virtual void print(ostream&,const char* = " ") const;
    virtual void print(FILE* =stdout,const char* = " ") const;
  };

class SimpleCoList
  : virtual public DescribedClass, virtual public SavableState {
DescribedClass_DECLARE(SimpleCoList)
SavableState_DECLARE(SimpleCoList)
    friend class SimpleCoListIter;
  protected:
    SimpleCoListLink *head;
    SimpleCoListLink *current;
  public:
    SimpleCoList();
    SimpleCoList(SimpleCo*);
    SimpleCoList(SimpleCoPtr&);
    SimpleCoList(const SimpleCoList&);
    virtual ~SimpleCoList();

    void add(SimpleCo*);
    void add(SimpleCoPtr&);
    void remove_first(SimpleCoPtr&);
    void remove_all(SimpleCoPtr&);

    void save_data_state(StateOut&so);
    void restore_data_state(int,StateIn&);

    SimpleCoPtr& operator[](int);

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

    SimpleCoPtr& this_object();

    int operator==(SimpleCo&);
    inline int operator!=(SimpleCo&b) { return !(*this == b); }
  };

#endif
