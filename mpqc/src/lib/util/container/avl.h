//
// avl.h --- definitions for some list classes
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

#ifndef _util_container_avl_h
#define _util_container_avl_h

#include <iostream.h>
#include <math.h>

template <class T, class K>
class EAVLNode {
  public:
    K key;
    T* lt;
    T* rt;
    T* up;
    int balance;
  public:
    EAVLNode(const K& k): key(k) {}
};

template <class T>
inline int
compare(const T& k1, const T& k2)
{
  return (k1<k2?-1:(k1==k2?0:1));
}

template <class T, class K>
class EAVLList {
  private:
    T *root_;
    T *start_;
    EAVLNode<T,K> T::* node_;
  private:
#if defined(__GNUC__) && defined(__alpha__)
    T*& rlink(T* n) const { T*& r = (n->*node_).rt; return r; }
    T*& llink(T* n) const { T*& r = (n->*node_).lt; return r; }
    T*& uplink(T* n) const { T*& r = (n->*node_).up; return r; }
    int& balance(T* n) const { int& r = (n->*node_).balance; return r; }
    K& key(T* n) const { K& r = (n->*node_).key; return r; }
#else    
    T*& rlink(T* n) const { return (n->*node_).rt; }
    T*& llink(T* n) const { return (n->*node_).lt; }
    T*& uplink(T* n) const { return (n->*node_).up; }
    int& balance(T* n) const { return (n->*node_).balance; }
    K& key(T* n) const { return (n->*node_).key; }
#endif
    int compare(T*n,T*m) const { return ::compare(key(n), key(m)); }
    int compare(T*n,const K&m) const { return ::compare(key(n), m); }

    void adjust_balance_insert(T* A, T* child);
    void adjust_balance_remove(T* A, T* child);
    void clear(T*);
  public:
    EAVLList();
    EAVLList(EAVLNode<T,K> T::* node);
    void initialize(EAVLNode<T,K> T::* node);
    void clear_without_delete() { initialize(node_); }
    void clear() { clear(root_); initialize(node_); }
    void insert(T*);
    void remove(T*);
    T* find(const K&) const;

    int height(T* node);
    int height() { return height(root_); }
    void check();
    void check_node(T*) const;

    int length() const;

    T* start() const { return start_; }
    void next(T*&) const;

    void print();

    int depth(T*);
};

template <class T, class K>
int
EAVLList<T,K>::length() const
{
  int r = 0;
  for (T* i=start(); i; next(i)) r++;
  return r;
}

template <class T, class K>
T*
EAVLList<T,K>::find(const K& key) const
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

template <class T, class K>
void
EAVLList<T,K>::remove(T* node)
{
  if (!node) return;

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
          cerr << "EAVLList::remove: inconsistency" << endl;
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

template <class T, class K>
void
EAVLList<T,K>::print()
{
  for (T*n=start(); n; next(n)) {
      int d = depth(n) + 1;
      for (int i=0; i<d; i++) cout << "     ";
      n->print();
      if (balance(n) == 1) cout << " (+)" << endl;
      else if (balance(n) == -1) cout << " (-)" << endl;
      else if (balance(n) == 0) cout << " (.)" << endl;
      else cout << " (" << balance(n) << ")" << endl;
    }
}

template <class T, class K>
int
EAVLList<T,K>::depth(T*node)
{
  int d = 0;
  while (node) {
      node = uplink(node);
      d++;
    }
  return d;
}

template <class T, class K>
void
EAVLList<T,K>::check_node(T*n) const
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

template <class T, class K>
void
EAVLList<T,K>::clear(T*n)
{
  if (!n) return;
  clear(llink(n));
  clear(rlink(n));
  delete n;
}

template <class T, class K>
int
EAVLList<T,K>::height(T* node)
{
  if (!node) return 0;
  int rh = height(rlink(node)) + 1;
  int lh = height(llink(node)) + 1;
  return rh>lh?rh:lh;
}

template <class T, class K>
void
EAVLList<T,K>::check()
{
  T* node;
  T* prev=0;
  for (node = start(); node; next(node)) {
      check_node(node);
      if (prev && compare(prev,node) > 0) {
          cerr << "nodes out of order" << endl;
          abort();
        }
      prev = node;
    }
  for (node = start(); node; next(node)) {
      if (balance(node) != height(rlink(node)) - height(llink(node))) {
          cerr << "balance inconsistency" << endl;
          abort();
        }
      if (balance(node) < -1 || balance(node) > 1) {
          cerr << "balance out of range" << endl;
          abort();
        }
    }
}

template <class T, class K>
void
EAVLList<T,K>::next(T*& node) const
{
  T* r;
  if (r = rlink(node)) {
      node = r;
      while (r = llink(node)) node = r;
      return;
    }
  while (r = uplink(node)) {
      if (node == llink(r)) {
          node = r;
          return;
        }
      node = r;
    }
  node = 0;
}

template <class T, class K>
void
EAVLList<T,K>::insert(T* n)
{
  if (!n) {
      return;
    }

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

  // maybe update the first node in the list
  if (have_start) start_ = n;

  // adjust the balance factors
  if (prev_p) {
      adjust_balance_insert(prev_p, n);
    }
}

template <class T, class K>
void
EAVLList<T,K>::adjust_balance_insert(T* A, T* child)
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

template <class T, class K>
void
EAVLList<T,K>::adjust_balance_remove(T* A, T* child)
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

template <class T, class K>
inline
EAVLList<T,K>::EAVLList()
{
  initialize(0);
}

template <class T, class K>
inline
EAVLList<T,K>::EAVLList(EAVLNode<T,K> T::* node)
{
  initialize(node);
}

template <class T, class K>
inline void
EAVLList<T,K>::initialize(EAVLNode<T,K> T::* node)
{
  node_ = node;
  root_ = 0;
  start_ = 0;
}


#undef TEST
#define TEST 0
#if TEST

class Data {
  public:
    EAVLNode<Data,int> list1;
    EAVLNode<Data,int> list2;
  public:
    Data(int k1, int k2 = 0): list1(k1), list2(k2) {};
    void print(int indent = 0);
    void change1(int val) { list1.key = val; }
};

void
Data::print(int indent)
{
  for (int i=0; i<indent; i++) cout << " ";
  cout << list1.key;
}

#define TEST1 1
#define TEST2 1
#define TEST3 1
#define TEST4 1
#define TEST5 1
#define TEST6 1
#define TEST7 1
#define TEST8 0
#define TEST9 0

static int Ni = 0;
static int Nr = 0;
static int Nf = 0;

void
testlist(EAVLList<Data, int>& list, Data** data, int n)
{
  for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
          list.insert(data[j]);
          Ni++;
        }
#if 0
      if (i==0) {
          cout << "--------------------------------------------" << endl;
          list.check();
          list.print2();
        }
      cout << "............................. removing ";
      data[i]->print();
      cout << endl;
#endif
      list.remove(data[i]);
      Nr++;
#if 0
      list.print2();
#endif
      list.check();
      list.clear_without_delete();
    }
}

void
rantest(EAVLList<Data, int>&list1, Data** data, int n)
{
  for (int i=0; i<n; i++) {
      Data* d = data[i];
      d->change1(random());
      list1.insert(d);
      Ni++;
    }
  list1.check();
  for (int i=0; i<n; i++) {
      list1.find(i);
      Nf++;
    }
  for (int i=0; i<n; i++) {
      Data* d = data[i];
      list1.remove(d);
      Nr++;
    }
  list1.check();
  list1.clear_without_delete();
}

int
main()
{
  const int maxkey = 9;
  EAVLList<Data,int> list1(&Data::list1);
  Data* data[maxkey][maxkey];
  Data* currentdata[maxkey];
  for (int i=0; i<maxkey; i++) {
      for (int j=0; j<maxkey; j++) {
          data[i][j] = new Data(j);
        }
    }
  int max;

  const int unique = 1;

#if TEST1
  cout << "=================================================" << endl;
  max = 1;
  for (int i=0; i<max; i++) {
      currentdata[0] = data[0][i];

      testlist(list1, currentdata, max);
    }
#endif

#if TEST2
  cout << "=================================================" << endl;
  max = 2;
  for (int i=0; i<max; i++) {
      currentdata[0] = data[0][i];
      for (int j=0; j<max; j++) {
          if (unique && i == j) continue;
          currentdata[1] = data[0][j];
          testlist(list1, currentdata, max);
        }
    }
#endif

#if TEST3
  cout << "=================================================" << endl;
  max = 3;
  for (int i=0; i<max; i++) {
      currentdata[0] = data[0][i];
      for (int j=0; j<max; j++) {
          if (unique && i == j) continue;
          currentdata[1] = data[0][j];
          for (int k=0; k<max; k++) {
              if (unique && k==i || k==j) continue;
              currentdata[2] = data[0][k];

              testlist(list1, currentdata, max);
            }
        }
    }
#endif

#if TEST4
  cout << "=================================================" << endl;
  max = 4;
  for (int i=0; i<max; i++) {
      currentdata[0] = data[0][i];
      for (int j=0; j<max; j++) {
          if (unique && i == j) continue;
          currentdata[1] = data[0][j];
          for (int k=0; k<max; k++) {
              if (unique && k==i || k==j) continue;
              currentdata[2] = data[0][k];
              for (int l=0; l<max; l++) {
                  if (unique && l==i || l==j || l==k) continue;
                  currentdata[3] = data[0][l];

                  testlist(list1, currentdata, max);
                }
            }
        }
    }
#endif

#if TEST5
  cout << "=================================================" << endl;
  max = 5;
  for (int i=0; i<max; i++) {
      currentdata[0] = data[0][i];
      for (int j=0; j<max; j++) {
          if (unique && i == j) continue;
          currentdata[1] = data[0][j];
          for (int k=0; k<max; k++) {
              if (unique && k==i || k==j) continue;
              currentdata[2] = data[0][k];
              for (int l=0; l<max; l++) {
                  if (unique && l==i || l==j || l==k) continue;
                  currentdata[3] = data[0][l];
                  for (int m=0; m<max; m++) {
                  if (unique && m==i || m==j || m==k || m==l) continue;
                  currentdata[4] = data[0][m];

                  testlist(list1, currentdata, max);
                  }
                }
            }
        }
    }
#endif

#if TEST6
  cout << "=================================================" << endl;
  max = 6;
  for (int i=0; i<max; i++) {
      currentdata[0] = data[0][i];
      for (int j=0; j<max; j++) {
          if (unique && i == j) continue;
          currentdata[1] = data[0][j];
          cout << "6: i = " << i << " j = " << j << endl;
          for (int k=0; k<max; k++) {
              if (unique && k==i || k==j) continue;
      currentdata[2] = data[0][k];
              for (int l=0; l<max; l++) {
                  if (unique && l==i || l==j || l==k) continue;
      currentdata[3] = data[0][l];
                  for (int m=0; m<max; m++) {
                  if (unique && m==i || m==j || m==k || m==l) continue;
      currentdata[4] = data[0][m];
                  for (int n=0; n<max; n++) {
                  if (unique && n==i || n==j || n==k || n==l || n==m) continue;
      currentdata[5] = data[0][n];

                  testlist(list1, currentdata, max);

                }
                }
                }
            }
        }
    }
#endif

#if TEST7
  cout << "=================================================" << endl;
  max = 7;
  for (int i=0; i<max; i++) {
      currentdata[0] = data[0][i];
      for (int j=0; j<max; j++) {
          if (unique && i == j) continue;
      currentdata[1] = data[0][j];
          for (int k=0; k<max; k++) {
              if (unique && k==i || k==j) continue;
              currentdata[2] = data[0][k];
              cout << "7: i = " << i << " j = " << j << " k = " << k << endl;
              for (int l=0; l<max; l++) {
                  if (unique && l==i || l==j || l==k) continue;
      currentdata[3] = data[0][l];
                  for (int m=0; m<max; m++) {
                  if (unique && m==i || m==j || m==k || m==l) continue;
      currentdata[4] = data[0][m];
                  for (int n=0; n<max; n++) {
                  if (unique && n==i || n==j || n==k || n==l || n==m) continue;
      currentdata[5] = data[0][n];
                  for (int o=0; o<max; o++) {
                  if (unique && o==i || o==j || o==k || o==l || o==m || o==n) continue;
      currentdata[6] = data[0][o];

                  testlist(list1, currentdata, max);

                }
                }
                }
                }
            }
        }
    }
#endif

#if TEST8
  cout << "=================================================" << endl;
  max = 8;
  for (int i=0; i<max; i++) {
      currentdata[0] = data[0][i];
      for (int j=0; j<max; j++) {
          if (unique && i == j) continue;
      currentdata[1] = data[0][j];
          for (int k=0; k<max; k++) {
              if (unique && k==i || k==j) continue;
      currentdata[2] = data[0][k];
              for (int l=0; l<max; l++) {
                  if (unique && l==i || l==j || l==k) continue;
      currentdata[3] = data[0][l];
                  cout << "7: i = " << i << " j = " << j << " k = " << k
                       << " l = " << l << endl;
                  for (int m=0; m<max; m++) {
                  if (unique && m==i || m==j || m==k || m==l) continue;
      currentdata[4] = data[0][m];
                  for (int n=0; n<max; n++) {
                  if (unique && n==i || n==j || n==k || n==l || n==m) continue;
      currentdata[5] = data[0][n];
                  for (int o=0; o<max; o++) {
                  if (unique && o==i || o==j || o==k || o==l || o==m || o==n) continue;
      currentdata[6] = data[0][o];
                  for (int p=0; p<max; p++) {
                  if (unique && p==i||p==j||p==k||p==l||p==m||p==n||p==o) continue;
      currentdata[7] = data[0][p];

                  testlist(list1, currentdata, max);

                }
                }
                }
                }
                }
            }
        }
    }
#endif

#if TEST9
  cout << "=================================================" << endl;
  max = 9;
  for (int i=0; i<max; i++) {
      currentdata[0] = data[0][i];
      for (int j=0; j<max; j++) {
          if (unique && i == j) continue;
      currentdata[1] = data[0][j];
          for (int k=0; k<max; k++) {
              if (unique && k==i || k==j) continue;
      currentdata[2] = data[0][k];
              for (int l=0; l<max; l++) {
                  if (unique && l==i || l==j || l==k) continue;
      currentdata[3] = data[0][l];
                  for (int m=0; m<max; m++) {
                  if (unique && m==i || m==j || m==k || m==l) continue;
      currentdata[4] = data[0][m];
                  cout << "7: i = " << i << " j = " << j << " k = " << k
                       << " l = " << l << " m = " << m << endl;
                  for (int n=0; n<max; n++) {
                  if (unique && n==i || n==j || n==k || n==l || n==m) continue;
      currentdata[5] = data[0][n];
                  for (int o=0; o<max; o++) {
                  if (unique && o==i || o==j || o==k || o==l || o==m || o==n) continue;
      currentdata[6] = data[0][o];
                  for (int p=0; p<max; p++) {
                  if (unique && p==i||p==j||p==k||p==l||p==m||p==n||p==o) continue;
      currentdata[7] = data[0][p];
                  for (int q=0; q<max; q++) {
                  if (unique && q==i||q==j||q==k||q==l||q==m||q==n||q==o||q==p) continue;
      currentdata[8] = data[0][q];

                  testlist(list1, currentdata, max);

                }
                }
                }
                }
                }
                }
            }
        }
    }
#endif

  cout << "Ni = " << Ni << ", Nr = " << Nr << ", N = " << Ni+Nr << endl;

  const int maxdat2 = 2000;
  Data * data2[maxdat2];
  for (int i=0; i<maxdat2; i++) {
      data2[i] = new Data(i);
    }
  for (int i=0; i<maxdat2; i++) {
      if (i%100 == 0) cout << "-";
    }
  cout << endl;
  for (int i=0; i<maxdat2; i++) {
      if (i%100 == 0) {
          cout << ".";
          fflush(stdout);
        }
      rantest(list1, data2, i);
    }
  cout << endl;

  cout << "Ni = " << Ni << ", Nr = " << Nr << ", Nf = " << Nf
       << ", N = " << Ni+Nr << endl;

  return 0;
}

#endif // TEST

#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
