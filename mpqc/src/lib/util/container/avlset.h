//
// avlset.h --- definition for avl set class
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

#ifndef _util_container_avlset_h
#define _util_container_avlset_h

#include <util/container/avlmap.h>

namespace sc {

template <class K>
class AVLSet {
  private:
    AVLMap<K,int> map_;
  public:
    class iterator {
      private:
        const EAVLMMap<K, AVLMapNode<K,int> > *map_;
        const AVLMapNode<K, int> *node;
      public:
        iterator(): map_(0), node(0) {}
        iterator(const EAVLMMap<K,AVLMapNode<K,int> > *m,
                 const AVLMapNode<K,int> *n)
          :map_(m), node(n) {}
        iterator(const eavl_typename AVLSet<K>::iterator &i):map_(i.map_),node(i.node) {}
        void operator++() { map_->next(node); }
        void operator++(int) { operator++(); }
        int operator == (const eavl_typename AVLSet<K>::iterator &i) const
            { return map_ == i.map_ && node == i.node; }
        int operator != (const eavl_typename AVLSet<K>::iterator &i) const
            { return !(map_ == i.map_ && node == i.node); }
        void operator = (const eavl_typename AVLSet<K>::iterator &i)
            { map_ = i.map_; node = i.node; }
        const K &key() const { return node->node.key; }
        const K &operator *() const { return node->node.key; }
        //const K *operator ->() const { return &node->node.key; }
    };
  public:
    AVLSet() {};
    void clear() { map_.clear(); }
    void insert(const K& key) { map_.insert(key,0); }
    void remove(const K& key) { map_.remove(key); }
    int contains(const K& key) const { return map_.contains(key); }
    iterator find(const K& k) const;

    int height() { return map_.height(); }
    void check() { map_.check(); }

    int length() const { return map_.length(); }

    iterator begin() const { return iterator(&map_.map_,map_.map_.start()); }
    iterator end() const { return iterator(&map_.map_,0); }

    void operator -= (const AVLSet<K> &s);
    void operator |= (const AVLSet<K> &s);

    void print() { map_.print(); }
};

template <class K>
void
AVLSet<K>::operator -= (const AVLSet<K> &s)
{
  for (typename AVLSet<K>::iterator i=s.begin(); i!=s.end(); i++) {
      remove(*i);
    }
}

template <class K>
void
AVLSet<K>::operator |= (const AVLSet<K> &s)
{
  for (typename AVLSet<K>::iterator i=s.begin(); i!=s.end(); i++) {
      insert(*i);
    }
}

template <class K>
inline typename AVLSet<K>::iterator
AVLSet<K>::find(const K& k) const
{
  return iterator(&map_.map_,map_.map_.find(k));
}

}

#endif

// ///////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
