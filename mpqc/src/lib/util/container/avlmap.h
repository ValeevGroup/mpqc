//
// avlmap.h --- definition for avl map class
//
// Copyright (C) 1998 Limit Point Systems, Inc.
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

#ifndef _util_container_avlmap_h
#define _util_container_avlmap_h

#include <util/container/eavlmmap.h>

namespace sc {
    
template <class K, class T>
class AVLMapNode {
  public:
    T data;
    EAVLMMapNode<K,AVLMapNode<K, T> > node;
  public:
    AVLMapNode(const K& k, const T& d): data(d), node(k) {};
};

template <class K, class T>
class AVLMap {
  public:
    EAVLMMap<K, AVLMapNode<K,T> > map_;
  public:
    class iterator {
      private:
        const EAVLMMap<K, AVLMapNode<K,T> > *map_;
        AVLMapNode<K, T> *node;
      public:
        iterator(): map_(0), node(0) {}
        iterator(const EAVLMMap<K,AVLMapNode<K,T> > *m,
                 AVLMapNode<K,T> *n)
          :map_(m), node(n) {}
        iterator(const eavl_typename AVLMap<K,T>::iterator &i) { map_=i.map_; node=i.node; }
        void operator++() { map_->next(node); }
        void operator++(int) { operator++(); }
        int operator == (const eavl_typename AVLMap<K,T>::iterator &i) const
            { return map_ == i.map_ && node == i.node; }
        int operator != (const eavl_typename AVLMap<K,T>::iterator &i) const
            { return !operator == (i); }
        void operator = (const eavl_typename AVLMap<K,T>::iterator &i)
            { map_ = i.map_; node = i.node; }
        const K &key() const { return node->node.key; }
        T &data() { return node->data; }
    };
  public:
    AVLMap(): map_(&AVLMapNode<K,T>::node) {};
    void clear() { map_.clear(); }
    void insert(const K& key, const T& data);
    void remove(const K& key);
    int contains(const K& k) const { return map_.find(k) != 0; }
    iterator find(const K&) const;
    T &operator[](const K &k);

    int height() { return map_.height(); }
    void check() { map_.check(); }

    int length() const { return map_.length(); }

    iterator begin() const { return iterator(&map_,map_.start()); }
    iterator end() const { return iterator(&map_,0); }

    void print() { map_.print(); }
};

template <class K, class T>
inline void
AVLMap<K,T>::insert(const K& key, const T& data)
{
  AVLMapNode<K,T> *node = map_.find(key);
  if (node) node->data = data;
  else map_.insert(new AVLMapNode<K, T>(key,data));
}

template <class K, class T>
inline void
AVLMap<K,T>::remove(const K& key)
{
  AVLMapNode<K, T> *node = map_.find(key);
  if (node) {
      map_.remove(node);
      delete node;
    }
}

template <class K, class T>
inline typename AVLMap<K,T>::iterator
AVLMap<K,T>::find(const K& k) const
{
  return iterator(&map_,map_.find(k));
}

template <class K, class T>
inline T&
AVLMap<K,T>::operator [](const K& k)
{
  AVLMapNode<K, T> *node = map_.find(k);
  if (node) return node->data;
  insert(k,T());
  node = map_.find(k);
  return node->data;
}

}

#endif

// /////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
