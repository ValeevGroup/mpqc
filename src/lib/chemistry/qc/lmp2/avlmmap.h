//
// avlmmap.h --- definition for avl multimap class
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
// Maintainer: SNL
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifndef _util_container_avlmmap_h
#define _util_container_avlmmap_h

#include <iostream>
#include <functional>
#include <vector>
#include <stdexcept>

#include <unistd.h> // for size_t on solaris
#include <stdlib.h>

namespace sc {

template <class K, class T>
class AVLMMapNode {
  public:
    std::pair<K,T> keyval;
    AVLMMapNode<K,T>* lt;
    AVLMMapNode<K,T>* rt;
    AVLMMapNode<K,T>* up;
    int balance;
  public:
    AVLMMapNode(const std::pair<K,T> &kv): keyval(kv) {}
    AVLMMapNode(const K& k, const T&v): keyval(std::make_pair(k,v)) {}
    K& key() { return keyval.first; }
    const K& key() const { return keyval.first; }
    T& val() { return keyval.second; }
    const T& val() const { return keyval.second; }
};

template <class T>
class chunk_allocator {
    enum {default_chunksize = 4096};
    struct chunk {
        void *data;
        chunk *next;
        int current;
        int chunksize;
        chunk(int _chunksize, int datasize) {
          data = new char[_chunksize*datasize];
          next = 0;
          current = 0;
          chunksize = _chunksize;
        }
        ~chunk() {
          delete[] (char*)data;
        }
    };
    chunk *head;
    long total;
  public:
    chunk_allocator() {
      total = 0;
      head = new chunk(default_chunksize, sizeof(T));
    }
    ~chunk_allocator() {
      if (total != 0) {
          std::cout << "chunks remain: " << total << std::endl;
          abort();
        }
      for (chunk *i = head; i;) {
          chunk *next = i->next;
          delete i;
          i = next;
        }
    }
    void reinit() {
      if (total != 0) {
          std::cout << "reinit: chunks remain: " << total << std::endl;
          abort();
        }
      for (chunk *i = head; i;) {
          chunk *next = i->next;
          delete i;
          i = next;
        }
      head = new chunk(default_chunksize, sizeof(T));
    }
    void destroy(T*d) {
      d->~T();
    }
    void deallocate(T*d,int n) {
      // deallocated memory is not reclaimed until the chunk_allocator dtor
      // is called
      total -= n;
    }
    T *allocate(int n) {
      if (n==1) {
          if (head->current == head->chunksize) {
              chunk *newhead = new chunk(default_chunksize, sizeof(T));
              newhead->next = head;
              head = newhead;
            }
          T *r = &((T*)head->data)[head->current];
          head->current++;
          total += n;
          return r;
        }
      else if (n>0) {
          chunk *newchunk = new chunk(n, sizeof(T));
          newchunk->next = head->next;
          head->next = newchunk;
          total += n;
          return (T*)newchunk->data;
        }
      else {
          return 0;
        }
    }
};

template <class T>
class tristate_less {
  public:
    bool operator()(const T&a, const T&b) const {
      return a < b;
    }
    int compare(const T&a, const T&b) const {
      if (a < b) return -1;
      else if (a > b) return 1;
      return 0;
    }
};

template <class K, class T, class C = tristate_less<K>,
          class A = chunk_allocator<AVLMMapNode<K,T> > >
//          class A = std::allocator<AVLMMapNode<K,T> > >
class AVLMMap {
  private:
    typedef AVLMMapNode<K,T> * link_type;
    C key_comp_;
    A alloc_;
    size_t length_;
    AVLMMapNode<K,T> *root_;
    AVLMMapNode<K,T> *start_;
  public:
    static AVLMMapNode<K,T>*& rlink(AVLMMapNode<K,T>* n) { return n->rt; }
    static AVLMMapNode<K,T>*  rlink(const AVLMMapNode<K,T>* n) { return n->rt; }
    static AVLMMapNode<K,T>*& llink(AVLMMapNode<K,T>* n) { return n->lt; }
    static AVLMMapNode<K,T>*  llink(const AVLMMapNode<K,T>* n) { return n->lt; }
    static AVLMMapNode<K,T>*& uplink(AVLMMapNode<K,T>* n) { return n->up; }
    static AVLMMapNode<K,T>*  uplink(const AVLMMapNode<K,T>* n) { return n->up; }
    static int& balance(AVLMMapNode<K,T>* n) { return n->balance; }
    static int  balance(const AVLMMapNode<K,T>* n) { return n->balance; }
    static const K& key(const AVLMMapNode<K,T>* n) { return n->key(); }
    int compare(const K&n,const K&m) const {
        return key_comp_.compare(n,m);
//       This is for a standard boolean comparision operation
//       if (key_comp_(n, m)) return -1;
//       else if (key_comp_(m, n)) return 1;
//       return 0;
    }
  private:
    void adjust_balance_insert(AVLMMapNode<K,T>* _A, AVLMMapNode<K,T>* child);
    void adjust_balance_remove(AVLMMapNode<K,T>* _A, AVLMMapNode<K,T>* child);
    void clear(AVLMMapNode<K,T>*);
  public:
    typedef C key_compare;
    typedef std::pair<K,T> value_type;
    class iterator {
        friend class AVLMMap<K,T,C,A>;
      public:
        AVLMMapNode<K,T> *node;
      public:
        iterator(AVLMMapNode<K,T> *n = 0):node(n){}
        iterator(const typename AVLMMap<K,T,C>::iterator &i) { node=i.node; }
        void operator++() { AVLMMap<K,T,C>::next(node); }
        void operator++(int) { operator++(); }
        int operator == (const typename AVLMMap<K,T,C>::iterator &i) const
            { return node == i.node; }
        int operator != (const typename AVLMMap<K,T,C>::iterator &i) const
            { return node != i.node; }
        void operator = (const typename AVLMMap<K,T,C>::iterator &i)
            { node = i.node; }
        const K &key() const { return node->key(); }
        std::pair<K,T> & operator *() { return node->keyval; }
        std::pair<K,T> * operator->() { return &node->keyval; }
    };
    class const_iterator {
        friend class AVLMMap<K,T,C,A>;
      private:
        const AVLMMapNode<K,T> *node;
      public:
        const_iterator(const AVLMMapNode<K,T> *n = 0):node(n){}
        const_iterator(const typename AVLMMap<K,T,C>::const_iterator &i) { node=i.node; }
        const_iterator(const typename AVLMMap<K,T,C>::iterator &i) { node=i.node; }
        void operator++() { AVLMMap<K,T,C>::next(node); }
        void operator++(int) { operator++(); }
        int operator == (const typename AVLMMap<K,T,C>::const_iterator &i) const
            { return node == i.node; }
        int operator != (const typename AVLMMap<K,T,C>::const_iterator &i) const
            { return node != i.node; }
        void operator = (const typename AVLMMap<K,T,C>::const_iterator &i)
            { node = i.node; }
        const K &key() const { return node->key(); }
        const std::pair<K,T> & operator *() { return node->keyval; }
        const std::pair<K,T> * operator->() { return &node->keyval; }
    };
  private:
    iterator subtree_insert_new_(AVLMMapNode<K,T> *root, AVLMMapNode<K,T> *node);
    iterator subtree_insert_equal_(AVLMMapNode<K,T> *root, AVLMMapNode<K,T> *node);
    iterator subtree_insert_unique_(AVLMMapNode<K,T> *root,
                                    const std::pair<K,T> &ipair);
    iterator subtree_find_(AVLMMapNode<K,T> *root, const K&);
    const_iterator subtree_find_(const AVLMMapNode<K,T> *root, const K&) const;
  public:
    AVLMMap();
    AVLMMap(const C& c);
    AVLMMap(const AVLMMap<K,T,C,A>&m);
    ~AVLMMap() { clear(); }
    key_compare key_comp() const { return key_comp_; }
    void initialize();
    void clear() {
      clear(root_);
      length_ = 0;
      root_ = 0;
      start_ = 0;
      // destroy the allocator
      // this is needed to force deletion of all data if the chunk_allocator is used
      alloc_.~A();
      // reconstruct the allocator
      new((void*)&alloc_)A;
//       alloc_.reinit();
    }
    // Insert an element. The element must not already exist.
    iterator insert_new(const iterator &hint, const std::pair<K,T> &ipair);
//     // Insert an element. If an element with the same key already exists
//     // then a duplicate element will be added.
//     iterator insert_equal(const iterator &hint, const std::pair<K,T> &ipair);
//     // Insert an element.  If an entry with the key already exists then
//     // this will do nothing and return the iterator for the pre-existing
//     // element.
//     iterator insert_unique(const iterator &hint, const std::pair<K,T> &ipair);
    // slow; do not use:
    iterator insert_equal(const iterator &hint, const std::pair<K,T> &ipair);
    // Insert an element. The element must not already exist.
    iterator insert_new(const std::pair<K,T> &ipair);
    // Insert an element. If an element with the same key already exists
    // then a duplicate element will be added.
    iterator insert_equal(const std::pair<K,T> &ipair);
    // Insert an element.  If an entry with the key already exists then
    // this will do nothing and return the iterator for the pre-existing
    // element.
    iterator insert_unique(const std::pair<K,T> &ipair);
    // This is equivalent to insert_equal.
    iterator insert(const std::pair<K,T> &ipair) { return insert_equal(ipair); }
    // This is equivalent to insert_equal with a hint.
    iterator insert(const iterator &hint, const std::pair<K,T> &ipair) { return insert_equal(hint, ipair); }

    // Insert elements in the given range.
    void insert(const const_iterator &b, const const_iterator &e) {
        for (const_iterator i = b; i != e; i++) insert(*i);
    }

    void remove(AVLMMapNode<K,T>*);
    iterator find(const K&);
    const_iterator find(const K&) const;
    iterator find(const iterator &hint, const K&);
    const_iterator find(const const_iterator &hint, const K&) const;

    void erase(const iterator &i) { remove(i.node); }

    int height(AVLMMapNode<K,T>* node);
    int height() { return height(root_); }
    void check();
    void check_node(AVLMMapNode<K,T>*) const;

    AVLMMapNode<K,T>* start() const { return start_; }
    static void next(const AVLMMapNode<K,T>*&);
    static void next(AVLMMapNode<K,T>*&);

    iterator begin() { return iterator(start()); }
    const_iterator begin() const { return const_iterator(start()); }

    iterator end() { return iterator(0); }
    const_iterator end() const { return const_iterator(0); }

    const_iterator lower_bound(const K& k) const;
    const_iterator upper_bound(const K& k) const;
    std::pair<const_iterator,const_iterator> equal_range(const K& k) const;

    void print(std::ostream &o = std::cout) const;
    size_t size() const { return length_; }
    int depth(AVLMMapNode<K,T>*) const;
};

template <class K, class T, class C, class A>
std::pair<typename AVLMMap<K,T,C,A>::const_iterator,
          typename AVLMMap<K,T,C,A>::const_iterator>
AVLMMap<K,T,C,A>::equal_range(const K& __k) const
{
  link_type __lb = 0; /* Last node which is not less than __k. */
  link_type __ub = 0; /* Last node which is greater than __k. */
  link_type __x = root_; /* Current node. */
      
  while (__x != 0) {
      int cmp = compare(key(__x), __k);
      if (cmp == 1) {
	  __lb = __x;
          __ub = __x;
          __x = llink(__x);
        }
      else if (cmp == -1) {
	  __x = rlink(__x);
        }
      else {
          // behaviours of lower_bound and upper_bound diverge here
          link_type __saved_x = __x;
          // find the lowerbound
          __lb = __x, __x = llink(__x);
          while (__x != 0)
              if (!key_comp_(key(__x), __k))
                  __lb = __x, __x = llink(__x);
              else
                  __x = rlink(__x);
          // find the upperbound
          __x = __saved_x;
          __x = rlink(__x);
          while (__x != 0)
              if (key_comp_(__k, key(__x)))
                  __ub = __x, __x = llink(__x);
              else
                  __x = rlink(__x);
        }
    }

  typedef typename AVLMMap<K,T,C,A>::const_iterator iterator;
  return std::pair<iterator,iterator>(__lb, __ub);
}

template <class K, class T, class C, class A>
typename AVLMMap<K,T,C,A>::const_iterator
AVLMMap<K,T,C,A>::lower_bound(const K& __k) const
{
  // this is the GNU STL RB tree algorithm
  link_type __y = 0; /* Last node which is not less than __k. */
  link_type __x = root_; /* Current node. */
      
  while (__x != 0)
      if (!key_comp_(key(__x), __k))
	  __y = __x, __x = llink(__x);
      else
	  __x = rlink(__x);
      
  return const_iterator(__y);
}

template <class K, class T, class C, class A>
typename AVLMMap<K,T,C,A>::const_iterator
AVLMMap<K,T,C,A>::upper_bound(const K& __k) const
{
  // this the GNU STL RB tree algorithm
  link_type __y = 0; /* Last node which is greater than __k. */
  link_type __x = root_; /* Current node. */
      
  while (__x != 0)
      if (key_comp_(__k, key(__x)))
	  __y = __x, __x = llink(__x);
      else
	  __x = rlink(__x);
      
  return const_iterator(__y);
}

template <class K, class T, class C, class A>
inline typename AVLMMap<K,T,C,A>::iterator
AVLMMap<K,T,C,A>::find(const K& k)
{
  return subtree_find_(root_,k);
}

template <class K, class T, class C, class A>
inline typename AVLMMap<K,T,C,A>::const_iterator
AVLMMap<K,T,C,A>::find(const K& k) const
{
  return subtree_find_(root_,k);
}

template <class K, class T, class C, class A>
inline typename AVLMMap<K,T,C,A>::iterator
AVLMMap<K,T,C,A>::find(const iterator &hint, const K& k)
{
  AVLMMapNode<K,T> *x = hint.node;
  if (x == 0) return subtree_find_(root_,k);
  AVLMMapNode<K,T> *y = uplink(x);
  while (y != 0) {
      int cmp = compare(k, key(x));
      if (cmp == 0) return x;
      int cmp2 = compare(k, key(y));
      if (cmp2 == 0) return y;
      if (cmp == -1 && cmp2 == 1) {
          return subtree_find_(llink(x), k);
        }
      else if (cmp == 1 && cmp2 == -1) {
          return subtree_find_(rlink(x), k);
        }
      x = y;
      cmp = cmp2;
      y = uplink(x);
    }
  return subtree_find_(x,k);
}

template <class K, class T, class C, class A>
inline typename AVLMMap<K,T,C,A>::const_iterator
AVLMMap<K,T,C,A>::find(const const_iterator &hint, const K& k) const
{
  const AVLMMapNode<K,T> *x = hint.node;
  if (x == 0) return subtree_find_(root_,k);
  const AVLMMapNode<K,T> *y = uplink(x);
  while (y != 0) {
      int cmp = compare(k, key(x));
      if (cmp == 0) return x;
      int cmp2 = compare(k, key(y));
      if (cmp2 == 0) return y;
      if (cmp == -1 && cmp2 == 1) {
          return subtree_find_(llink(x), k);
        }
      else if (cmp == 1 && cmp2 == -1) {
          return subtree_find_(rlink(x), k);
        }
      x = y;
      cmp = cmp2;
      y = uplink(x);
    }
  return subtree_find_(x,k);
}

template <class K, class T, class C, class A>
typename AVLMMap<K,T,C,A>::iterator
AVLMMap<K,T,C,A>::subtree_find_(AVLMMapNode<K,T> *root, const K& k)
{
  AVLMMapNode<K,T>* n = root;

  while (n) {
      int cmp = compare(key(n), k);
      if (cmp < 0) n = rlink(n);
      else if (cmp > 0) n = llink(n);
      else return iterator(n);
    }

  return iterator(0);
}

template <class K, class T, class C, class A>
typename AVLMMap<K,T,C,A>::const_iterator
AVLMMap<K,T,C,A>::subtree_find_(const AVLMMapNode<K,T> *root, const K& k) const
{
  const AVLMMapNode<K,T>* n = root;

  while (n) {
      int cmp = compare(key(n), k);
      if (cmp < 0) n = rlink(n);
      else if (cmp > 0) n = llink(n);
      else return const_iterator(n);
    }

  return const_iterator(0);
}

template <class K, class T, class C, class A>
void
AVLMMap<K,T,C,A>::remove(AVLMMapNode<K,T>* node)
{
  if (!node) return;

  length_--;

  if (node == start_) {
      next(start_);
    }

  AVLMMapNode<K,T> *rebalance_point;
  AVLMMapNode<K,T> *q;

  if (llink(node) == 0) {
      q = rlink(node);
      rebalance_point = uplink(node);

      if (q) uplink(q) = rebalance_point;

      if (rebalance_point) {
          if (rlink(rebalance_point) == node) rlink(rebalance_point) = q;
          else llink(rebalance_point) = q;
        }
      else root_ = q;
    }
  else if (rlink(node) == 0) {
      q = llink(node);
      rebalance_point = uplink(node);

      if (q) uplink(q) = rebalance_point;

      if (rebalance_point) {
          if (rlink(rebalance_point) == node) rlink(rebalance_point) = q;
          else llink(rebalance_point) = q;
        }
      else root_ = q;
    }
  else {
      AVLMMapNode<K,T> *r = node;
      next(r);

      if (r == 0 || llink(r) != 0) {
          throw std::runtime_error("AVLMMap::remove: inconsistency");
        }

      if (r == rlink(node)) {
          llink(r) = llink(node);
          if (llink(r)) uplink(llink(r)) = r;
          balance(r) = balance(node);
          rebalance_point = r;
          q = rlink(r);
        }
      else {
          q = rlink(r);

          rebalance_point = uplink(r);

          if (llink(rebalance_point) == r) llink(rebalance_point) = q;
          else rlink(rebalance_point) = q;

          if (q) uplink(q) = rebalance_point;

          balance(r) = balance(node);
          rlink(r) = rlink(node);
          llink(r) = llink(node);
          if (rlink(r)) uplink(rlink(r)) = r;
          if (llink(r)) uplink(llink(r)) = r;
        }
      if (r) {
          AVLMMapNode<K,T>* up = uplink(node);
          uplink(r) = up;
          if (up) {
              if (rlink(up) == node) rlink(up) = r;
              else llink(up) = r;
            }
          if (up == 0) root_ = r;
        }
    }

  // adjust balance won't work if both children are null,
  // so handle this special case here
  if (rebalance_point &&
      llink(rebalance_point) == 0 && rlink(rebalance_point) == 0) {
      balance(rebalance_point) = 0;
      q = rebalance_point;
      rebalance_point = uplink(rebalance_point);
    }
  adjust_balance_remove(rebalance_point, q);

  alloc_.destroy(node);
  alloc_.deallocate(node,1);
}

template <class K, class T, class C, class A>
void
AVLMMap<K,T,C,A>::print(std::ostream &o) const
{
  std::cout << "--------------------------" << std::endl;
  for (AVLMMapNode<K,T>*n=start(); n; next(n)) {
      int d = depth(n) + 1;
      for (int i=0; i<d; i++) o << "     ";
      if (balance(n) == 1) o << " (+) " << key(n) << std::endl;
      else if (balance(n) == -1) o << " (-) " << key(n) << std::endl;
      else if (balance(n) == 0) o << " (.) " << key(n) << std::endl;
      else o << " (" << balance(n) << ") " << key(n) << std::endl;
    }
}

template <class K, class T, class C, class A>
int
AVLMMap<K,T,C,A>::depth(AVLMMapNode<K,T>*node) const
{
  int d = 0;
  while (node) {
      node = uplink(node);
      d++;
    }
  return d;
}

template <class K, class T, class C, class A>
void
AVLMMap<K,T,C,A>::check_node(AVLMMapNode<K,T>*n) const
{
  if (uplink(n) && uplink(n) == n) abort();
  if (llink(n) && llink(n) == n) abort();
  if (rlink(n) && rlink(n) == n) abort();
  if (rlink(n) && rlink(n) == llink(n)) abort();
  if (uplink(n) && uplink(n) == llink(n)) abort();
  if (uplink(n) && uplink(n) == rlink(n)) abort();
  if (uplink(n) && !(llink(uplink(n)) == n
                     || rlink(uplink(n)) == n)) abort();
}

template <class K, class T, class C, class A>
void
AVLMMap<K,T,C,A>::clear(AVLMMapNode<K,T>*n)
{
  if (!n) return;
  clear(llink(n));
  clear(rlink(n));
  alloc_.destroy(n);
  alloc_.deallocate(n,1);
}

template <class K, class T, class C, class A>
int
AVLMMap<K,T,C,A>::height(AVLMMapNode<K,T>* node)
{
  if (!node) return 0;
  int rh = height(rlink(node)) + 1;
  int lh = height(llink(node)) + 1;
  return rh>lh?rh:lh;
}

template <class K, class T, class C, class A>
void
AVLMMap<K,T,C,A>::check()
{
  AVLMMapNode<K,T>* node;
  AVLMMapNode<K,T>* prev=0;
  size_t computed_length = 0;
  for (node = start(); node; next(node)) {
      check_node(node);
      if (prev && compare(key(prev),key(node)) > 0) {
          throw std::runtime_error("nodes out of order");
        }
      prev = node;
      computed_length++;
    }
  for (node = start(); node; next(node)) {
      if (balance(node) != height(rlink(node)) - height(llink(node))) {
          throw std::runtime_error("balance inconsistency");
        }
      if (balance(node) < -1 || balance(node) > 1) {
          throw std::runtime_error("balance out of range");
        }
    }
  if (length_ != computed_length) {
      throw std::runtime_error("length error");
    }
}

template <class K, class T, class C, class A>
void
AVLMMap<K,T,C,A>::next(const AVLMMapNode<K,T>*& node)
{
  const AVLMMapNode<K,T>* r;
  if ((r = rlink(node))) {
      node = r;
      while ((r = llink(node))) node = r;
      return;
    }
  while ((r = uplink(node))) {
      if (node == llink(r)) {
          node = r;
          return;
        }
      node = r;
    }
  node = 0;
}

template <class K, class T, class C, class A>
void
AVLMMap<K,T,C,A>::next(AVLMMapNode<K,T>*& node)
{
  AVLMMapNode<K,T>* r;
  if ((r = rlink(node))) {
      node = r;
      while ((r = llink(node))) node = r;
      return;
    }
  while ((r = uplink(node))) {
      if (node == llink(r)) {
          node = r;
          return;
        }
      node = r;
    }
  node = 0;
}

template <class K, class T, class C, class A>
typename AVLMMap<K,T,C,A>::iterator
AVLMMap<K,T,C,A>::insert_equal(const iterator &hint, const std::pair<K,T> &ipair)
{
  AVLMMapNode<K,T>* n = alloc_.allocate(1);
  new((void*)n)AVLMMapNode<K,T>(ipair);

  llink(n) = 0;
  rlink(n) = 0;
  balance(n) = 0;

  AVLMMapNode<K,T> *x = hint.node;
  if (x == 0) return subtree_insert_equal_(root_,n);
  AVLMMapNode<K,T> *y = uplink(x);
  while (y != 0) {
      int cmp = compare(key(n), key(x));
      if (cmp == 0) return subtree_insert_equal_(x, n);
      int cmp2 = compare(key(n), key(y));
      if (cmp2 == 0) return subtree_insert_equal_(y, n);
      if (cmp == -1 && cmp2 == 1) {
          y = llink(x);
          if (y != 0) return subtree_insert_equal_(y, n);
          else {
              uplink(n) = x;
              llink(x) = n;
              adjust_balance_insert(x,n);
              return iterator(n);
            }
        }
      else if (cmp == 1 && cmp2 == -1) {
          y = rlink(x);
          if (y != 0) return subtree_insert_equal_(y, n);
          else {
              uplink(n) = x;
              rlink(x) = n;
              adjust_balance_insert(x,n);
              return iterator(n);
            }
        }
      x = y;
      cmp = cmp2;
      y = uplink(x);
    }
  return subtree_insert_equal_(x, n);
}

template <class K, class T, class C, class A>
typename AVLMMap<K,T,C,A>::iterator
AVLMMap<K,T,C,A>::insert_new(const iterator &hint, const std::pair<K,T> &ipair)
{
  AVLMMapNode<K,T>* n = alloc_.allocate(1);
  new((void*)n)AVLMMapNode<K,T>(ipair);

  llink(n) = 0;
  rlink(n) = 0;
  balance(n) = 0;

  AVLMMapNode<K,T> *x = hint.node;
  if (x == 0) return subtree_insert_new_(root_,n);
  AVLMMapNode<K,T> *y = uplink(x);
  while (y != 0) {
      bool cmp = key_comp_(key(n), key(x));
      bool cmp2 = key_comp_(key(n), key(y));
      if (cmp && !cmp2) {
          y = llink(x);
          if (y != 0) return subtree_insert_new_(y, n);
          else {
              uplink(n) = x;
              llink(x) = n;
              adjust_balance_insert(x,n);
              return iterator(n);
            }
        }
      else if (!cmp && cmp2) {
          y = rlink(x);
          if (y != 0) return subtree_insert_new_(y, n);
          else {
              uplink(n) = x;
              rlink(x) = n;
              adjust_balance_insert(x,n);
              return iterator(n);
            }
        }
      x = y;
      cmp = cmp2;
      y = uplink(x);
    }
  return subtree_insert_new_(x, n);
}

template <class K, class T, class C, class A>
inline typename AVLMMap<K,T,C,A>::iterator
AVLMMap<K,T,C,A>::insert_equal(const std::pair<K,T> &ipair)
{
  AVLMMapNode<K,T>* n = alloc_.allocate(1);
  new((void*)n)AVLMMapNode<K,T>(ipair);

  return subtree_insert_equal_(root_, n);
}

template <class K, class T, class C, class A>
inline typename AVLMMap<K,T,C,A>::iterator
AVLMMap<K,T,C,A>::insert_unique(const std::pair<K,T> &ipair)
{
  return subtree_insert_unique_(root_, ipair);
}

template <class K, class T, class C, class A>
inline typename AVLMMap<K,T,C,A>::iterator
AVLMMap<K,T,C,A>::insert_new(const std::pair<K,T> &ipair)
{
  AVLMMapNode<K,T>* n = alloc_.allocate(1);
  new((void*)n)AVLMMapNode<K,T>(ipair);

  length_++;

  rlink(n) = 0;
  llink(n) = 0;
  balance(n) = 0;

  return subtree_insert_new_(root_, n);
}

template <class K, class T, class C, class A>
typename AVLMMap<K,T,C,A>::iterator
AVLMMap<K,T,C,A>::subtree_insert_equal_(AVLMMapNode<K,T> *root,
                                        AVLMMapNode<K,T> *n)
{
  length_++;

  rlink(n) = 0;
  llink(n) = 0;
  balance(n) = 0;

  // root is only null if the actual root is null
  if (!root) {
      uplink(n) = 0;
      root_ = n;
      start_ = n;
      return iterator(n);
    }

  // find an insertion point
  AVLMMapNode<K,T>* p = root;
  AVLMMapNode<K,T>* prev_p = 0;
  int cmp;
  while (p) {
      if (p == n) {
          abort();
        }
      prev_p = p;
      cmp = compare(key(n),key(p));
      if (cmp < 0) p = llink(p);
      else {
          p = rlink(p);
        }
    }

  // insert the node
  uplink(n) = prev_p;
  if (prev_p) {
      if (cmp < 0) {
          llink(prev_p) = n;
          if (prev_p == start_) start_ = n;
        }
      else rlink(prev_p) = n;
    }

  // adjust the balance factors
  if (prev_p) {
      adjust_balance_insert(prev_p, n);
    }

  return iterator(n);
}

template <class K, class T, class C, class A>
typename AVLMMap<K,T,C,A>::iterator
AVLMMap<K,T,C,A>::subtree_insert_unique_(AVLMMapNode<K,T> *root,
                                         const std::pair<K,T> &ipair)
{
  // root is only null if the actual root is null
  if (!root) {
      AVLMMapNode<K,T>* n = alloc_.allocate(1);
      new((void*)n)AVLMMapNode<K,T>(ipair);

      length_++;

      rlink(n) = 0;
      llink(n) = 0;
      balance(n) = 0;
      uplink(n) = 0;
      root_ = n;
      start_ = n;
      return iterator(n);
    }

  // find an insertion point
  AVLMMapNode<K,T>* p = root;
  AVLMMapNode<K,T>* prev_p = 0;
  int cmp;
  while (p) {
      prev_p = p;
      cmp = compare(ipair.first,key(p));
      if (cmp < 0) p = llink(p);
      else if (cmp == 0) {
          return iterator(p);
        }
      else {
          p = rlink(p);
        }
    }

  // no match was found; create the node
  AVLMMapNode<K,T>* n = alloc_.allocate(1);
  new((void*)n)AVLMMapNode<K,T>(ipair);

  length_++;

  rlink(n) = 0;
  llink(n) = 0;
  balance(n) = 0;

  // insert the node
  uplink(n) = prev_p;
  if (prev_p) {
      if (cmp < 0) {
          llink(prev_p) = n;
          if (prev_p == start_) start_ = n;
        }
      else rlink(prev_p) = n;
    }

  // adjust the balance factors
  if (prev_p) {
      adjust_balance_insert(prev_p, n);
    }

  return iterator(n);
}

template <class K, class T, class C, class A>
typename AVLMMap<K,T,C,A>::iterator
AVLMMap<K,T,C,A>::subtree_insert_new_(AVLMMapNode<K,T> *root,
                                      AVLMMapNode<K,T> *n)
{
  // root is only null if the actual root is null
  if (!root) {
      uplink(n) = 0;
      root_ = n;
      start_ = n;
      return iterator(n);
    }

  // find an insertion point
  AVLMMapNode<K,T>* p = root;
  AVLMMapNode<K,T>* prev_p = 0;
  bool cmp;
  while (p) {
      prev_p = p;
      cmp = key_comp_(key(n),key(p));
      if (cmp) p = llink(p);
      else {
          p = rlink(p);
        }
    }

  // insert the node
  uplink(n) = prev_p;
  if (prev_p) {
      if (cmp) {
          llink(prev_p) = n;
          if (prev_p == start_) start_ = n;
        }
      else rlink(prev_p) = n;
    }

  // adjust the balance factors
  if (prev_p) {
      adjust_balance_insert(prev_p, n);
    }

  return iterator(n);
}

template <class K, class T, class C, class Alloc>
void
AVLMMap<K,T,C,Alloc>::adjust_balance_insert(AVLMMapNode<K,T>* A, AVLMMapNode<K,T>* child)
{
  if (!A) return;
  int adjustment;
  if (llink(A) == child) adjustment = -1;
  else adjustment = 1;
  int bal = balance(A) + adjustment;
  if (bal == 0) {
      balance(A) = 0;
    }
  else if (bal == -1 || bal == 1) {
      balance(A) = bal;
      adjust_balance_insert(uplink(A), A);
    }
  else if (bal == 2) {
      AVLMMapNode<K,T>* B = rlink(A);
      if (balance(B) == 1) {
          balance(B) = 0;
          balance(A) = 0;
          rlink(A) = llink(B);
          llink(B) = A;
          uplink(B) = uplink(A);
          uplink(A) = B;
          if (rlink(A)) uplink(rlink(A)) = A;
          if (llink(B)) uplink(llink(B)) = B;
          if (uplink(B) == 0) root_ = B;
          else {
              if (rlink(uplink(B)) == A) rlink(uplink(B)) = B;
              else llink(uplink(B)) = B;
            }
        }
      else {
          AVLMMapNode<K,T>* X = llink(B);
          rlink(A) = llink(X);
          llink(B) = rlink(X);
          llink(X) = A;
          rlink(X) = B;
          if (balance(X) == 1) {
              balance(A) = -1;
              balance(B) = 0;
            }
          else if (balance(X) == -1) {
              balance(A) = 0;
              balance(B) = 1;
            }
          else {
              balance(A) = 0;
              balance(B) = 0;
            }
          balance(X) = 0;
          uplink(X) = uplink(A);
          uplink(A) = X;
          uplink(B) = X;
          if (rlink(A)) uplink(rlink(A)) = A;
          if (llink(B)) uplink(llink(B)) = B;
          if (uplink(X) == 0) root_ = X;
          else {
              if (rlink(uplink(X)) == A) rlink(uplink(X)) = X;
              else llink(uplink(X)) = X;
            }
        }
    }
  else if (bal == -2) {
      AVLMMapNode<K,T>* B = llink(A);
      if (balance(B) == -1) {
          balance(B) = 0;
          balance(A) = 0;
          llink(A) = rlink(B);
          rlink(B) = A;
          uplink(B) = uplink(A);
          uplink(A) = B;
          if (llink(A)) uplink(llink(A)) = A;
          if (rlink(B)) uplink(rlink(B)) = B;
          if (uplink(B) == 0) root_ = B;
          else {
              if (llink(uplink(B)) == A) llink(uplink(B)) = B;
              else rlink(uplink(B)) = B;
            }
        }
      else {
          AVLMMapNode<K,T>* X = rlink(B);
          llink(A) = rlink(X);
          rlink(B) = llink(X);
          rlink(X) = A;
          llink(X) = B;
          if (balance(X) == -1) {
              balance(A) = 1;
              balance(B) = 0;
            }
          else if (balance(X) == 1) {
              balance(A) = 0;
              balance(B) = -1;
            }
          else {
              balance(A) = 0;
              balance(B) = 0;
            }
          balance(X) = 0;
          uplink(X) = uplink(A);
          uplink(A) = X;
          uplink(B) = X;
          if (llink(A)) uplink(llink(A)) = A;
          if (rlink(B)) uplink(rlink(B)) = B;
          if (uplink(X) == 0) root_ = X;
          else {
              if (llink(uplink(X)) == A) llink(uplink(X)) = X;
              else rlink(uplink(X)) = X;
            }
        }
    }
}

template <class K, class T, class C, class Alloc>
void
AVLMMap<K,T,C,Alloc>::adjust_balance_remove(AVLMMapNode<K,T>* A, AVLMMapNode<K,T>* child)
{
  if (!A) return;
  int adjustment;
  if (llink(A) == child) adjustment = 1;
  else adjustment = -1;
  int bal = balance(A) + adjustment;
  if (bal == 0) {
      balance(A) = 0;
      adjust_balance_remove(uplink(A), A);
    }
  else if (bal == -1 || bal == 1) {
      balance(A) = bal;
    }
  else if (bal == 2) {
      AVLMMapNode<K,T>* B = rlink(A);
      if (balance(B) == 0) {
          balance(B) = -1;
          balance(A) = 1;
          rlink(A) = llink(B);
          llink(B) = A;
          uplink(B) = uplink(A);
          uplink(A) = B;
          if (rlink(A)) uplink(rlink(A)) = A;
          if (llink(B)) uplink(llink(B)) = B;
          if (uplink(B) == 0) root_ = B;
          else {
              if (rlink(uplink(B)) == A) rlink(uplink(B)) = B;
              else llink(uplink(B)) = B;
            }
        }
      else if (balance(B) == 1) {
          balance(B) = 0;
          balance(A) = 0;
          rlink(A) = llink(B);
          llink(B) = A;
          uplink(B) = uplink(A);
          uplink(A) = B;
          if (rlink(A)) uplink(rlink(A)) = A;
          if (llink(B)) uplink(llink(B)) = B;
          if (uplink(B) == 0) root_ = B;
          else {
              if (rlink(uplink(B)) == A) rlink(uplink(B)) = B;
              else llink(uplink(B)) = B;
            }
          adjust_balance_remove(uplink(B), B);
        }
      else {
          AVLMMapNode<K,T>* X = llink(B);
          rlink(A) = llink(X);
          llink(B) = rlink(X);
          llink(X) = A;
          rlink(X) = B;
          if (balance(X) == 0) {
              balance(A) = 0;
              balance(B) = 0;
            }
          else if (balance(X) == 1) {
              balance(A) = -1;
              balance(B) = 0;
            }
          else {
              balance(A) = 0;
              balance(B) = 1;
            }
          balance(X) = 0;
          uplink(X) = uplink(A);
          uplink(A) = X;
          uplink(B) = X;
          if (rlink(A)) uplink(rlink(A)) = A;
          if (llink(B)) uplink(llink(B)) = B;
          if (uplink(X) == 0) root_ = X;
          else {
              if (rlink(uplink(X)) == A) rlink(uplink(X)) = X;
              else llink(uplink(X)) = X;
            }
          adjust_balance_remove(uplink(X), X);
        }
    }
  else if (bal == -2) {
      AVLMMapNode<K,T>* B = llink(A);
      if (balance(B) == 0) {
          balance(B) = 1;
          balance(A) = -1;
          llink(A) = rlink(B);
          rlink(B) = A;
          uplink(B) = uplink(A);
          uplink(A) = B;
          if (llink(A)) uplink(llink(A)) = A;
          if (rlink(B)) uplink(rlink(B)) = B;
          if (uplink(B) == 0) root_ = B;
          else {
              if (llink(uplink(B)) == A) llink(uplink(B)) = B;
              else rlink(uplink(B)) = B;
            }
        }
      else if (balance(B) == -1) {
          balance(B) = 0;
          balance(A) = 0;
          llink(A) = rlink(B);
          rlink(B) = A;
          uplink(B) = uplink(A);
          uplink(A) = B;
          if (llink(A)) uplink(llink(A)) = A;
          if (rlink(B)) uplink(rlink(B)) = B;
          if (uplink(B) == 0) root_ = B;
          else {
              if (llink(uplink(B)) == A) llink(uplink(B)) = B;
              else rlink(uplink(B)) = B;
            }
          adjust_balance_remove(uplink(B), B);
        }
      else {
          AVLMMapNode<K,T>* X = rlink(B);
          llink(A) = rlink(X);
          rlink(B) = llink(X);
          rlink(X) = A;
          llink(X) = B;
          if (balance(X) == 0) {
              balance(A) = 0;
              balance(B) = 0;
            }
          else if (balance(X) == -1) {
              balance(A) = 1;
              balance(B) = 0;
            }
          else {
              balance(A) = 0;
              balance(B) = -1;
            }
          balance(X) = 0;
          uplink(X) = uplink(A);
          uplink(A) = X;
          uplink(B) = X;
          if (llink(A)) uplink(llink(A)) = A;
          if (rlink(B)) uplink(rlink(B)) = B;
          if (uplink(X) == 0) root_ = X;
          else {
              if (llink(uplink(X)) == A) llink(uplink(X)) = X;
              else rlink(uplink(X)) = X;
            }
          adjust_balance_remove(uplink(X), X);
        }
    }
}

template <class K, class T, class C, class A>
inline
AVLMMap<K,T,C,A>::AVLMMap()
{
  initialize();
}

template <class K, class T, class C, class A>
inline
AVLMMap<K,T,C,A>::AVLMMap(const C &c):
  key_comp_(c)
{
  initialize();
}

template <class K, class T, class C, class A>
inline
AVLMMap<K,T,C,A>::AVLMMap(const AVLMMap<K,T,C,A>&m)
{
  initialize();
  insert(m.begin(), m.end());
}

template <class K, class T, class C, class A>
inline void
AVLMMap<K,T,C,A>::initialize()
{
  root_ = 0;
  start_ = 0;
  length_ = 0;
}

}

#endif

// ///////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
