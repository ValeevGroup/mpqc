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

#ifdef __GNUG__
#pragma implementation
#endif

#include <stream.h>
#include <stdlib.h>
#include "gnuavlse.h"

AVLSetBase::~AVLSetBase()
{
  _kill(root);
}

/*
 constants & inlines for maintaining balance & thread status in tree nodes
*/

#define AVLBALANCEMASK    3
#define AVLBALANCED       0
#define AVLLEFTHEAVY      1
#define AVLRIGHTHEAVY     2

#define LTHREADBIT        4
#define RTHREADBIT        8


static inline int bf(AVLNodeBase* t)
{
  return t->stat & AVLBALANCEMASK;
}

static inline void set_bf(AVLNodeBase* t, int b)
{
  t->stat = (t->stat & ~AVLBALANCEMASK) | (b & AVLBALANCEMASK);
}


static inline int rthread(AVLNodeBase* t)
{
  return t->stat & RTHREADBIT;
}

static inline void set_rthread(AVLNodeBase* t, int b)
{
  if (b)
    t->stat |= RTHREADBIT;
  else
    t->stat &= ~RTHREADBIT;
}

static inline int lthread(AVLNodeBase* t)
{
  return t->stat & LTHREADBIT;
}

static inline void set_lthread(AVLNodeBase* t, int b)
{
  if (b)
    t->stat |= LTHREADBIT;
  else
    t->stat &= ~LTHREADBIT;
}

/*
 traversal primitives
*/


AVLNodeBase* AVLSetBase::leftmost()
{
  AVLNodeBase* t = root;
  if (t != 0) while (t->lt != 0) t = t->lt;
  return t;
}

AVLNodeBase* AVLSetBase::rightmost()
{
  AVLNodeBase* t = root;
  if (t != 0) while (t->rt != 0) t = t->rt;
  return t;
}

AVLNodeBase* AVLSetBase::succ(AVLNodeBase* t)
{
  AVLNodeBase* r = t->rt;
  if (!rthread(t)) while (!lthread(r)) r = r->lt;
  return r;
}

AVLNodeBase* AVLSetBase::pred(AVLNodeBase* t)
{
  AVLNodeBase* l = t->lt;
  if (!lthread(t)) while (!rthread(l)) l = l->rt;
  return l;
}


Pix AVLSetBase::seek(void*key)
{
  AVLNodeBase* t = root;
  if (t == 0)
    return 0;
  for (;;)
  {
    int cmp = (*cmpfunc)(key, (*item)(t));
    if (cmp == 0)
      return Pix(t);
    else if (cmp < 0)
    {
      if (lthread(t))
        return 0;
      else
        t = t->lt;
    }
    else if (rthread(t))
      return 0;
    else
      t = t->rt;
  }
}


/*
 The combination of threads and AVL bits make adding & deleting
 interesting, but very awkward.

 We use the following statics to avoid passing them around recursively
*/

static int _need_rebalancing;   // to send back balance info from rec. calls
static void*   _target_item;     // add/del_item target
static AVLNodeBase* _found_node; // returned added/deleted node
static int    _already_found;   // for deletion subcases

static AVLNodeBase** _hold_nodes;       // used for rebuilding trees
static int  _max_hold_index;              // # elements-1 in _hold_nodes


void AVLSetBase:: _add(AVLNodeBase*& t)
{
  int cmp = (*cmpfunc)(_target_item, (*item)(t));
  if (cmp == 0)
  {
    _found_node = t;
    return;
  }
  else if (cmp < 0)
  {
    if (lthread(t))
    {
      ++count;
      _found_node = (*newnode)(_target_item);
      set_lthread(_found_node, 1);
      set_rthread(_found_node, 1);
      _found_node->lt = t->lt;
      _found_node->rt = t;
      t->lt = _found_node;
      set_lthread(t, 0);
      _need_rebalancing = 1;
    }
    else
      _add(t->lt);
    if (_need_rebalancing)
    {
      switch(bf(t))
      {
      case AVLRIGHTHEAVY:
        set_bf(t, AVLBALANCED);
        _need_rebalancing = 0;
        return;
      case AVLBALANCED:
        set_bf(t, AVLLEFTHEAVY);
        return;
      case AVLLEFTHEAVY:
	{
        AVLNodeBase* l = t->lt;
        if (bf(l) == AVLLEFTHEAVY)
        {
          if (rthread(l))
            t->lt = l;
          else
            t->lt = l->rt;
          set_lthread(t, rthread(l));
          l->rt = t;
          set_rthread(l, 0);
          set_bf(t, AVLBALANCED);
          set_bf(l, AVLBALANCED);
          t = l;
          _need_rebalancing = 0;
        }
        else
        {
          AVLNodeBase* r = l->rt;
          set_rthread(l, lthread(r));
          if (lthread(r))
            l->rt = r;
          else
            l->rt = r->lt;
          r->lt = l;
          set_lthread(r, 0);
          set_lthread(t, rthread(r));
          if (rthread(r))
            t->lt = r;
          else
            t->lt = r->rt;
          r->rt = t;
          set_rthread(r, 0);
          if (bf(r) == AVLLEFTHEAVY)
            set_bf(t, AVLRIGHTHEAVY);
          else
            set_bf(t, AVLBALANCED);
          if (bf(r) == AVLRIGHTHEAVY)
            set_bf(l, AVLLEFTHEAVY);
          else
            set_bf(l, AVLBALANCED);
          set_bf(r, AVLBALANCED);
          t = r;
          _need_rebalancing = 0;
          return;
        }
	}
      }
    }
  }
  else
  {
    if (rthread(t))
    {
      ++count;
      _found_node = (*newnode)(_target_item);
      set_rthread(t, 0);
      set_lthread(_found_node, 1);
      set_rthread(_found_node, 1);
      _found_node->lt = t;
      _found_node->rt = t->rt;
      t->rt = _found_node;
      _need_rebalancing = 1;
    }
    else
      _add(t->rt);
    if (_need_rebalancing)
    {
      switch(bf(t))
      {
      case AVLLEFTHEAVY:
        set_bf(t, AVLBALANCED);
        _need_rebalancing = 0;
        return;
      case AVLBALANCED:
        set_bf(t, AVLRIGHTHEAVY);
        return;
      case AVLRIGHTHEAVY:
	{
        AVLNodeBase* r = t->rt;
        if (bf(r) == AVLRIGHTHEAVY)
        {
          if (lthread(r))
            t->rt = r;
          else
            t->rt = r->lt;
          set_rthread(t, lthread(r));
          r->lt = t;
          set_lthread(r, 0);
          set_bf(t, AVLBALANCED);
          set_bf(r, AVLBALANCED);
          t = r;
          _need_rebalancing = 0;
        }
        else
        {
          AVLNodeBase* l = r->lt;
          set_lthread(r, rthread(l));
          if (rthread(l))
            r->lt = l;
          else
            r->lt = l->rt;
          l->rt = r;
          set_rthread(l, 0);
          set_rthread(t, lthread(l));
          if (lthread(l))
            t->rt = l;
          else
            t->rt = l->lt;
          l->lt = t;
          set_lthread(l, 0);
          if (bf(l) == AVLRIGHTHEAVY)
            set_bf(t, AVLLEFTHEAVY);
          else
            set_bf(t, AVLBALANCED);
          if (bf(l) == AVLLEFTHEAVY)
            set_bf(r, AVLRIGHTHEAVY);
          else
            set_bf(r, AVLBALANCED);
          set_bf(l, AVLBALANCED);
          t = l;
          _need_rebalancing = 0;
          return;
        }
	}
      }
    }
  }
}

    
Pix AVLSetBase::add(void* item)
{
  if (root == 0)
  {
    ++count;
    root = (*newnode)(item);
    set_rthread(root, 1);
    set_lthread(root, 1);
    return Pix(root);
  }
  else
  {
    _target_item = item;
    _need_rebalancing = 0;
    _add(root);
    return Pix(_found_node);
  }
}


void AVLSetBase::_del(AVLNodeBase* par, AVLNodeBase*& t)
{
  int comp;
  if (_already_found)
  {
    if (rthread(t))
      comp = 0;
    else
      comp = 1;
  }
  else 
    comp = (*cmpfunc)(_target_item, (*item)(t));
  if (comp == 0)
  {
    if (lthread(t) && rthread(t))
    {
      _found_node = t;
      if (t == par->lt)
      {
        set_lthread(par, 1);
        par->lt = t->lt;
      }
      else
      {
        set_rthread(par, 1);
        par->rt = t->rt;
      }
      _need_rebalancing = 1;
      return;
    }
    else if (lthread(t))
    {
      _found_node = t;
      AVLNodeBase* s = succ(t);
      if (s != 0 && lthread(s))
        s->lt = t->lt;
      t = t->rt;
      _need_rebalancing = 1;
      return;
    }
    else if (rthread(t))
    {
      _found_node = t;
      AVLNodeBase* p = pred(t);
      if (p != 0 && rthread(p))
        p->rt = t->rt;
      t = t->lt;
      _need_rebalancing = 1;
      return;
    }
    else                        // replace item & find someone deletable
    {
      AVLNodeBase* p = pred(t);
      (*itemassign)((*item)(t), (*item)(p));
      _already_found = 1;
      comp = -1;                // fall through below to left
    }
  }

  if (comp < 0)
  {
    if (lthread(t))
      return;
    _del(t, t->lt);
    if (!_need_rebalancing)
      return;
    switch (bf(t))
    {
    case AVLLEFTHEAVY:
      set_bf(t, AVLBALANCED);
      return;
    case AVLBALANCED:
      set_bf(t, AVLRIGHTHEAVY);
      _need_rebalancing = 0;
      return;
    case AVLRIGHTHEAVY:
      {
      AVLNodeBase* r = t->rt;
      switch (bf(r))
      {
      case AVLBALANCED:
        if (lthread(r))
          t->rt = r;
        else
          t->rt = r->lt;
        set_rthread(t, lthread(r));
        r->lt = t;
        set_lthread(r, 0);
        set_bf(t, AVLRIGHTHEAVY);
        set_bf(r, AVLLEFTHEAVY);
        _need_rebalancing = 0;
        t = r;
        return;
      case AVLRIGHTHEAVY:
        if (lthread(r))
          t->rt = r;
        else
          t->rt = r->lt;
        set_rthread(t, lthread(r));
        r->lt = t;
        set_lthread(r, 0);
        set_bf(t, AVLBALANCED);
        set_bf(r, AVLBALANCED);
        t = r;
        return;
      case AVLLEFTHEAVY:
	{
        AVLNodeBase* l = r->lt;
        set_lthread(r, rthread(l));
        if (rthread(l))
          r->lt = l;
        else
          r->lt = l->rt;
        l->rt = r;
        set_rthread(l, 0);
        set_rthread(t, lthread(l));
        if (lthread(l))
          t->rt = l;
        else
          t->rt = l->lt;
        l->lt = t;
        set_lthread(l, 0);
        if (bf(l) == AVLRIGHTHEAVY)
          set_bf(t, AVLLEFTHEAVY);
        else
          set_bf(t, AVLBALANCED);
        if (bf(l) == AVLLEFTHEAVY)
          set_bf(r, AVLRIGHTHEAVY);
        else
          set_bf(r, AVLBALANCED);
        set_bf(l, AVLBALANCED);
        t = l;
        return;
	}
      }
    }
    }
  }
  else
  {
    if (rthread(t))
      return;
    _del(t, t->rt);
    if (!_need_rebalancing)
      return;
    switch (bf(t))
    {
    case AVLRIGHTHEAVY:
      set_bf(t, AVLBALANCED);
      return;
    case AVLBALANCED:
      set_bf(t, AVLLEFTHEAVY);
      _need_rebalancing = 0;
      return;
    case AVLLEFTHEAVY:
      {
      AVLNodeBase* l = t->lt;
      switch (bf(l))
      {
      case AVLBALANCED:
        if (rthread(l))
          t->lt = l;
        else
          t->lt = l->rt;
        set_lthread(t, rthread(l));
        l->rt = t;
        set_rthread(l, 0);
        set_bf(t, AVLLEFTHEAVY);
        set_bf(l, AVLRIGHTHEAVY);
        _need_rebalancing = 0;
        t = l;
        return;
      case AVLLEFTHEAVY:
        if (rthread(l))
          t->lt = l;
        else
          t->lt = l->rt;
        set_lthread(t, rthread(l));
        l->rt = t;
        set_rthread(l, 0);
        set_bf(t, AVLBALANCED);
        set_bf(l, AVLBALANCED);
        t = l;
        return;
      case AVLRIGHTHEAVY:
	{
        AVLNodeBase* r = l->rt;
        set_rthread(l, lthread(r));
        if (lthread(r))
          l->rt = r;
        else
          l->rt = r->lt;
        r->lt = l;
        set_lthread(r, 0);
        set_lthread(t, rthread(r));
        if (rthread(r))
          t->lt = r;
        else
          t->lt = r->rt;
        r->rt = t;
        set_rthread(r, 0);
        if (bf(r) == AVLLEFTHEAVY)
          set_bf(t, AVLRIGHTHEAVY);
        else
          set_bf(t, AVLBALANCED);
        if (bf(r) == AVLRIGHTHEAVY)
          set_bf(l, AVLLEFTHEAVY);
        else
          set_bf(l, AVLBALANCED);
        set_bf(r, AVLBALANCED);
        t = r;
        return;
	}
      }
      }
    }
  }
}

        

void AVLSetBase::del(void* item)
{
  if (root == 0) return;
  _need_rebalancing = 0;
  _already_found = 0;
  _found_node = 0;
  _target_item = item;
  _del(root, root);
  if (_found_node)
  {
    (*deletenode)(_found_node);
    if (--count == 0)
      root = 0;
  }
}

// build an ordered array of pointers to tree nodes back into a tree
// we know that at least one element exists

static AVLNodeBase* _do_treeify(int lo, int hi, int& h)
{
  int lh, rh;
  int mid = (lo + hi) / 2;
  AVLNodeBase* t = _hold_nodes[mid];
  if (lo > mid - 1)
  {
    set_lthread(t, 1);
    if (mid == 0)
      t->lt = 0;
    else
      t->lt = _hold_nodes[mid-1];
    lh = 0;
  }
  else
  {
    set_lthread(t, 0);
    t->lt = _do_treeify(lo, mid-1, lh);
  }
  if (hi < mid + 1)
  {
    set_rthread(t, 1);
    if (mid == _max_hold_index)
      t->rt = 0;
    else
      t->rt = _hold_nodes[mid+1];
    rh = 0;
  }
  else 
  {
    set_rthread(t, 0);
    t->rt = _do_treeify(mid+1, hi, rh);
  }
  if (lh == rh)
  {
    set_bf(t, AVLBALANCED);
    h = lh + 1;
  }
  else if (lh == rh - 1)
  {
    set_bf(t, AVLRIGHTHEAVY);
    h = rh + 1;
  }
  else if (rh == lh - 1)
  {
    set_bf(t, AVLLEFTHEAVY);
    h = lh + 1;
  }
  else                          // can't happen
    abort();

  return t;
}

static AVLNodeBase* _treeify(int n)
{
  AVLNodeBase* t;
  if (n == 0)
    t = 0;
  else
  {
    int b;
    _max_hold_index = n-1;
    t = _do_treeify(0, _max_hold_index, b);
  }
  delete _hold_nodes;
  return t;
}


void AVLSetBase::_kill(AVLNodeBase* t)
{
  if (t != 0)
  {
    if (!lthread(t)) _kill(t->lt);
    if (!rthread(t)) _kill(t->rt);
    (*deletenode)(t);
  }
}


AVLSetBase::AVLSetBase(AVLSetBase& b)
{
  if ((count = b.count) == 0)
  {
    root = 0;
  }
  else
  {
    _hold_nodes = new AVLNodeBasePtr [count];
    AVLNodeBase* t = b.leftmost();
    int i = 0;
    while (t != 0)
    {
      _hold_nodes[i++] = (*newnode)((*item)(t));
      t = b.succ(t);
    }
    root = _treeify(count);
  }
}


int AVLSetBase::operator == (AVLSetBase& y)
{
  if (count != y.count)
    return 0;
  else
  {
    AVLNodeBase* t = leftmost();
    AVLNodeBase* u = y.leftmost();
    for (;;)
    {
      if (t == 0)
        return 1;
      else if (!((*eq)((*item)(t), (*item)(u))))
        return 0;
      else
      {
        t = succ(t);
        u = y.succ(u);
      }
    }
  }
}

int AVLSetBase::operator <= (AVLSetBase& y)
{
  if (count > y.count)
    return 0;
  else
  {
    AVLNodeBase* t = leftmost();
    AVLNodeBase* u = y.leftmost();
    for (;;)
    {
      if (t == 0)
        return 1;
      else if (u == 0)
        return 0;
      int cmp = (*cmpfunc)((*item)(t), (*item)(u));
      if (cmp == 0)
      {
        t = succ(t);
        u = y.succ(u);
      }
      else if (cmp < 0)
        return 0;
      else
        u = y.succ(u);
    }
  }
}

void AVLSetBase::operator |=(AVLSetBase& y)
{
  AVLNodeBase* t = leftmost();
  AVLNodeBase* u = y.leftmost();
  int rsize = count + y.count;
  _hold_nodes = new AVLNodeBasePtr [rsize];
  int k = 0;
  for (;;)
  {
    if (t == 0)
    {
      while (u != 0)
      {
        _hold_nodes[k++] = (*newnode)((*item)(u));
        u = y.succ(u);
      }
      break;
    }
    else if (u == 0)
    {
      while (t != 0)
      {
        _hold_nodes[k++] = t;
        t = succ(t);
      }
      break;
    }
    int cmp = (*cmpfunc)((*item)(t), (*item)(u));
    if (cmp == 0)
    {
      _hold_nodes[k++] = t;
      t = succ(t);
      u = y.succ(u);
    }
    else if (cmp < 0)
    {
      _hold_nodes[k++] = t;
      t = succ(t);
    }
    else
    {
      _hold_nodes[k++] = (*newnode)((*item)(u));
      u = y.succ(u);
    }
  }
  root = _treeify(k);
  count = k;
}

void AVLSetBase::operator &= (AVLSetBase& y)
{
  AVLNodeBase* t = leftmost();
  AVLNodeBase* u = y.leftmost();
  int rsize = (count < y.count)? count : y.count;
  _hold_nodes = new AVLNodeBasePtr [rsize];
  int k = 0;
  for (;;)
  {
    if (t == 0)
      break;
    if (u == 0)
    {
      while (t != 0)
      {
        AVLNodeBase* tmp = succ(t);
        (*deletenode)(t);
        t = tmp;
      }
      break;
    }
    int cmp = (*cmpfunc)((*item)(t), (*item)(u));
    if (cmp == 0)
    {
      _hold_nodes[k++] = t;
      t = succ(t);
      u = y.succ(u);
    }
    else if (cmp < 0)
    {
      AVLNodeBase* tmp = succ(t);
      (*deletenode)(t);
      t = tmp;
    }
    else
      u = y.succ(u);
  }
  root = _treeify(k);
  count = k;
}


void AVLSetBase::operator -=(AVLSetBase& y)
{
  AVLNodeBase* t = leftmost();
  AVLNodeBase* u = y.leftmost();
  int rsize = count;
  _hold_nodes = new AVLNodeBasePtr [rsize];
  int k = 0;
  for (;;)
  {
    if (t == 0)
      break;
    else if (u == 0)
    {
      while (t != 0)
      {
        _hold_nodes[k++] = t;
        t = succ(t);
      }
      break;
    }
    int cmp = (*cmpfunc)((*item)(t), (*item)(u));
    if (cmp == 0)
    {
      AVLNodeBase* tmp = succ(t);
      (*deletenode)(t);
      t = tmp;
      u = y.succ(u);
    }
    else if (cmp < 0)
    {
      _hold_nodes[k++] = t;
      t = succ(t);
    }
    else
      u = y.succ(u);
  }
  root = _treeify(k);
  count = k;
}

int AVLSetBase::owns(Pix i)
{
  if (i == 0) return 0;
  for (AVLNodeBase* t = leftmost(); t != 0; t = succ(t)) 
    if (Pix(t) == i) return 1;
  return 0;
}

int AVLSetBase::OK()
{
  int v = 1;
  if (root == 0) 
    v = count == 0;
  else
  {
    int n = 1;
    AVLNodeBase* trail = leftmost();
    AVLNodeBase* t = succ(trail);
    while (t != 0)
    {
      ++n;
      v &= (*cmpfunc)((*item)(trail), (*item)(t)) < 0;
      trail = t;
      t = succ(t);
    }
    v &= n == count;
  }
  if (!v) error("invariant failure");
  return v;
}
