// This may look like C code, but it is really -*- C++ -*-
// Modified by C L Janssen
/* 
Copyright (C) 1988 Free Software Foundation
    written by Doug Lea (dl@rocky.oswego.edu)

This file is part of the GNU C++ Library.  This library is free
software; you can redistribute it and/or modify it under the terms of
the GNU Library General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your
option) any later version.  This library is distributed in the hope
that it will be useful, but WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the GNU Library General Public License for more details.
You should have received a copy of the GNU Library General Public
License along with this library; if not, write to the Free Software
Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#ifndef _util_container_gnuavlse_h
#ifdef __GNUG__
#pragma interface
#endif

#include <util/container/gnuset.h>

class AVLNodeBase
{
  protected:
    // cannot permit delete since this is a base class without a virtual DTOR
    ~AVLNodeBase() {};
  public:
    AVLNodeBase*         lt;
    AVLNodeBase*         rt;
    char                stat;

    AVLNodeBase(AVLNodeBase*l, AVLNodeBase*r): lt(l), rt(r), stat(0) {}
};

typedef AVLNodeBase* AVLNodeBasePtr;

template <class T>
class AVLNode: public AVLNodeBase
{
  public:
    T                 item;
    AVLNode(T& h, AVLNode<T>* l=0, AVLNode<T>* r=0):
      item(h), AVLNodeBase(l, r) {}
    ~AVLNode() {}

static //template <class T>
void* get_item(AVLNodeBase*p) {
    return (void*)&((AVLNode<T>*)p)->item;
};

static //template <class T>
int equal(void*a1, void*a2)
{
  //printf("equal(%d,%d)\n",*(T*)a1, *(T*)a2);
  return *(T*)a1 == *(T*)a2;
}

static //template <class T>
int compare_items(void*a1, void*a2)
{
  //printf("compare_items(%d,%d)\n",*(T*)a1, *(T*)a2);
  T& d1 = *((T*)a1);
  T& d2 = *((T*)a2);
  return (d1 < d2)? -1: ((d1 == d2)? 0: 1);
}

static //template <class T>
AVLNodeBase* new_from_item(void*d)
{
  //printf("new_from_item(%d)\n",*(T*)d);
  return new AVLNode<T>(*((T*)d));
}

static //template <class T>
void delete_node(AVLNodeBase*b)
{
  delete (AVLNode<T>*)b;
}

static //template <class T>
void item_assign(void*i1,void*i2)
{
  *(T*)i1 = *(T*)i2;
}
};

typedef AVLNodeBase* (*AVLNodeNew)(void*);
typedef void (*AVLNodeDelete)(AVLNodeBase*);
typedef void* (*AVLNodeBaseGetItem)(AVLNodeBase*);
typedef int (*AVLNodeBaseInequality)(void*,void*);
typedef void (*AVLNodeItemAssign)(void*,void*);

class AVLSetBase
{
  protected:
    AVLNodeItemAssign     itemassign;
    AVLNodeNew            newnode;
    AVLNodeDelete         deletenode;
    AVLNodeBaseGetItem    item;
    AVLNodeBaseInequality eq;
    AVLNodeBaseInequality cmpfunc;
    
    AVLNodeBase* root;

    int count;

    AVLNodeBase* leftmost();
    AVLNodeBase* rightmost();
    AVLNodeBase* pred(AVLNodeBase* t);
    AVLNodeBase* succ(AVLNodeBase* t);
    void _kill(AVLNodeBase* t);
    void _add(AVLNodeBase*& t);
    void _del(AVLNodeBase* p, AVLNodeBase*& t);

    AVLSetBase(AVLNodeBaseInequality equal,
               AVLNodeBaseInequality compare,
               AVLNodeNew n,
               AVLNodeDelete d,
               AVLNodeBaseGetItem g,
               AVLNodeItemAssign ia):
      eq(equal),
      cmpfunc(compare),
      newnode(n),
      deletenode(d),
      item(g),
      itemassign(ia),
      root(0) {};
  public:
    virtual ~AVLSetBase();
    void clear() {
        _kill(root);
        root = 0;
        count = 0;
      }
    Pix last() {
        return Pix(rightmost());
      }

    void prev(Pix& i) {
        if (i != 0) i = Pix(pred((AVLNodeBase*)i));
      }
    Pix first() {
        return Pix(leftmost());
      }
    void next(Pix& i) {
        if (i != 0) i = Pix(succ((AVLNodeBase*)i));
      }
    int owns(Pix i);
    int OK();
    Pix seek(void*);
    Pix add(void*);
    void del(void*);
    AVLSetBase(AVLSetBase&);

    void operator |= (AVLSetBase& b);
    void operator -= (AVLSetBase& b);
    void operator &= (AVLSetBase& b);
         
    int operator == (AVLSetBase& b);
    int operator <= (AVLSetBase& b);
    void error(const char* msg) {
        (*lib_error_handler)("AVLSetBase", msg);
      }
};

template <class T>
class AVLSet: public Set<T>, private AVLSetBase
{

  public:
    AVLSet(): AVLSetBase(AVLNode<T>::equal,
                         AVLNode<T>::compare_items,
                         AVLNode<T>::new_from_item,
                         AVLNode<T>::delete_node,
                         AVLNode<T>::get_item,
                         AVLNode<T>::item_assign) {}
    AVLSet(AVLSet<T>& a): AVLSetBase(a) {}
    ~AVLSet() {}

    int length() { return count; }

    Pix add(T& item) { return AVLSetBase::add((void*)&item); }
    void del(T& item) { AVLSetBase::del((void*)&item); }
    int contains(T& item) {
        return seek(item) != 0;
      }
         
    void clear() {
        AVLSetBase::clear();
      }
         
    T& operator () (Pix i) {
        if (i == 0) AVLSetBase::error("null Pix");
        return *(T*)(*item)((AVLNodeBase*)i);
      }
    Pix seek(T& item) { return AVLSetBase::seek((void*)&item); }
         
    void operator |= (AVLSet<T>& b) { AVLSetBase::operator |= (b); }
    void operator -= (AVLSet<T>& b) { AVLSetBase::operator -= (b); }
    void operator &= (AVLSet<T>& b) { AVLSetBase::operator &= (b); }
         
    int operator == (AVLSet<T>& b) { AVLSetBase::operator == (b); }
    int operator != (AVLSet<T>& b) {
        return ! ((*this) == b);
      }
    int operator <= (AVLSet<T>& b) { AVLSetBase::operator <= (b); }

    Pix first() { return AVLSetBase::first(); }
    void next(Pix& i) { AVLSetBase::next(i); }
    int OK() { return AVLSetBase::OK(); }
    
};

#endif
