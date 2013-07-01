//
// eavlmmap.h --- definition for embedded avl multimap class
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
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

#ifndef _util_container_eavlmmap_h
#define _util_container_eavlmmap_h

#ifdef HAVE_CONFIG_H
#  include <mpqc_config.h>
#endif
#include <iostream>
#include <util/misc/exenv.h>
#include <util/container/compare.h>
#include <unistd.h> // for size_t on solaris
#include <stdlib.h>

#ifdef __GNUC__
// gcc typename seems to be broken in some cases
#  define eavl_typename
#else
#  define eavl_typename typename
#endif

namespace sc {

template <class K, class T>
class EAVLMMapNode {
  public:
    K key;
    T* lt;
    T* rt;
    T* up;
    int balance;
  public:
    EAVLMMapNode(const K& k): key(k) {}
};

template <class K, class T>
class EAVLMMap {
  private:
    size_t length_;
    T *root_;
    T *start_;
    EAVLMMapNode<K,T> T::* node_;
  public:
    T*& rlink(T* n) const { return (n->*node_).rt; }
    T* rlink(const T* n) const { return (n->*node_).rt; }
    T*& llink(T* n) const { return (n->*node_).lt; }
    T* llink(const T* n) const { return (n->*node_).lt; }
    T*& uplink(T* n) const { return (n->*node_).up; }
    T* uplink(const T* n) const { return (n->*node_).up; }
    int& balance(T* n) const { return (n->*node_).balance; }
    int balance(const T* n) const { return (n->*node_).balance; }
    K& key(T* n) const { return (n->*node_).key; }
    const K& key(const T* n) const { return (n->*node_).key; }
    int compare(T*n,T*m) const { return sc::compare(key(n), key(m)); }
    int compare(T*n,const K&m) const { return sc::compare(key(n), m); }
  private:
    void adjust_balance_insert(T* A, T* child);
    void adjust_balance_remove(T* A, T* child);
    void clear(T*);
  public:
    class iterator {
      private:
        EAVLMMap<K,T> *map_;
        T *node;
      public:
        iterator(EAVLMMap<K,T> *m, T *n):map_(m),node(n){}
        iterator(const eavl_typename EAVLMMap<K,T>::iterator &i) { map_=i.map_; node=i.node; }
        void operator++() { map_->next(node); }
        void operator++(int) { operator++(); }
        int operator == (const eavl_typename EAVLMMap<K,T>::iterator &i) const
            { return map_ == i.map_ && node == i.node; }
        int operator != (const eavl_typename EAVLMMap<K,T>::iterator &i) const
            { return !operator == (i); }
        void operator = (const eavl_typename EAVLMMap<K,T>::iterator &i)
            { map_ = i.map_; node = i.node; }
        const K &key() const { return map_->key(node); }
        T & operator *() { return *node; }
        T * operator->() { return node; }
    };
  public:
    EAVLMMap();
    EAVLMMap(EAVLMMapNode<K,T> T::* node);
    ~EAVLMMap() { clear(root_); }
    void initialize(EAVLMMapNode<K,T> T::* node);
    void clear_without_delete() { initialize(node_); }
    void clear() { clear(root_); initialize(node_); }
    void insert(T*);
    void remove(T*);
    T* find(const K&) const;

    int height(T* node);
    int height() { return height(root_); }
    void check();
    void check_node(T*) const;

    T* start() const { return start_; }
    void next(const T*&) const;
    void next(T*&) const;

    iterator begin() { return iterator(this,start()); }
    iterator end() { return iterator(this,0); }

    void print() const;
    void detailed_print() const;
    int length() const { return length_; }
    int depth(T*) const;
};

template <class K, class T>
T*
EAVLMMap<K,T>::find(const K& key) const
{
  T* n = root_;

  while (n) {
      int cmp = compare(n, key);
      if (cmp < 0) n = rlink(n);
      else if (cmp > 0) n = llink(n);
      else return n;
    }

  return 0;
}

template <class K, class T>
void
EAVLMMap<K,T>::remove(T* node)
{
  if (!node) return;

  length_--;

  if (node == start_) {
      next(start_);
    }

  T *rebalance_point;
  T *q;

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
      T *r = node;
      next(r);

      if (r == 0 || llink(r) != 0) {
          ExEnv::errn() << "EAVLMMap::remove: inconsistency" << std::endl;
          abort();
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
          T* up = uplink(node);
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
}

template <class K, class T>
void
EAVLMMap<K,T>::print() const
{
  for (T*n=start(); n; next(n)) {
      int d = depth(n) + 1;
      for (int i=0; i<d; i++) ExEnv::out0() << "     ";
      if (balance(n) == 1) ExEnv::out0() << " (+)" << std::endl;
      else if (balance(n) == -1) ExEnv::out0() << " (-)" << std::endl;
      else if (balance(n) == 0) ExEnv::out0() << " (.)" << std::endl;
      else ExEnv::out0() << " (" << balance(n) << ")" << std::endl;
    }
}

template <class K, class T>
void
EAVLMMap<K,T>::detailed_print() const
{
  for (T*n=start(); n; next(n)) {
      int d = depth(n) + 1;
      for (int i=0; i<d; i++) ExEnv::out0() << "     ";
      if (balance(n) == 1) ExEnv::out0() << " (+) "
                                         << key(n) << std::endl;
      else if (balance(n) == -1) ExEnv::out0() << " (-) "
                                               << key(n) << std::endl;
      else if (balance(n) == 0) ExEnv::out0() << " (.) "
                                              << key(n) << std::endl;
      else ExEnv::out0() << " (" << balance(n) << ") "
                         << key(n) << std::endl;
    }
}

template <class K, class T>
int
EAVLMMap<K,T>::depth(T*node) const
{
  int d = 0;
  while (node) {
      node = uplink(node);
      d++;
    }
  return d;
}

template <class K, class T>
void
EAVLMMap<K,T>::check_node(T*n) const
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

template <class K, class T>
void
EAVLMMap<K,T>::clear(T*n)
{
  if (!n) return;
  clear(llink(n));
  clear(rlink(n));
  delete n;
}

template <class K, class T>
int
EAVLMMap<K,T>::height(T* node)
{
  if (!node) return 0;
  int rh = height(rlink(node)) + 1;
  int lh = height(llink(node)) + 1;
  return rh>lh?rh:lh;
}

template <class K, class T>
void
EAVLMMap<K,T>::check()
{
  T* node;
  T* prev=0;
  size_t computed_length = 0;
  for (node = start(); node; next(node)) {
      check_node(node);
      if (prev && compare(prev,node) > 0) {
          ExEnv::errn() << "nodes out of order" << std::endl;
          abort();
        }
      prev = node;
      computed_length++;
    }
  for (node = start(); node; next(node)) {
      if (balance(node) != height(rlink(node)) - height(llink(node))) {
          ExEnv::errn() << "balance inconsistency" << std::endl;
          abort();
        }
      if (balance(node) < -1 || balance(node) > 1) {
          ExEnv::errn() << "balance out of range" << std::endl;
          abort();
        }
    }
  if (length_ != computed_length) {
      ExEnv::errn() << "length error" << std::endl;
      abort();
    }
}

template <class K, class T>
void
EAVLMMap<K,T>::next(const T*& node) const
{
  const T* r = rlink(node);
  if (r) {
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

template <class K, class T>
void
EAVLMMap<K,T>::next(T*& node) const
{
  T* r = rlink(node);
  if (r) {
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

template <class K, class T>
void
EAVLMMap<K,T>::insert(T* n)
{
  if (!n) {
      return;
    }

  length_++;

  rlink(n) = 0;
  llink(n) = 0;
  balance(n) = 0;

  if (!root_) {
      uplink(n) = 0;
      root_ = n;
      start_ = n;
      return;
    }

  // find an insertion point
  T* p = root_;
  T* prev_p = 0;
  int cmp;
  int have_start = 1;
  while (p) {
      if (p == n) {
          abort();
        }
      prev_p = p;
      cmp = compare(n,p);
      if (cmp < 0) p = llink(p);
      else {
          p = rlink(p);
          have_start = 0;
        }
    }

  // insert the node
  uplink(n) = prev_p;
  if (prev_p) {
      if (cmp < 0) llink(prev_p) = n;
      else rlink(prev_p) = n;
    }

  // maybe update the first node in the map
  if (have_start) start_ = n;

  // adjust the balance factors
  if (prev_p) {
      adjust_balance_insert(prev_p, n);
    }
}

template <class K, class T>
void
EAVLMMap<K,T>::adjust_balance_insert(T* A, T* child)
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
      T* B = rlink(A);
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
          T* X = llink(B);
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
      T* B = llink(A);
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
          T* X = rlink(B);
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

template <class K, class T>
void
EAVLMMap<K,T>::adjust_balance_remove(T* A, T* child)
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
      T* B = rlink(A);
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
          T* X = llink(B);
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
      T* B = llink(A);
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
          T* X = rlink(B);
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

template <class K, class T>
inline
EAVLMMap<K,T>::EAVLMMap()
{
  initialize(0);
}

template <class K, class T>
inline
EAVLMMap<K,T>::EAVLMMap(EAVLMMapNode<K,T> T::* node)
{
  initialize(node);
}

template <class K, class T>
inline void
EAVLMMap<K,T>::initialize(EAVLMMapNode<K,T> T::* node)
{
  node_ = node;
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
