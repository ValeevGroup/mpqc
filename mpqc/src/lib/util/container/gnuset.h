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

#ifndef _misc_container_gnuset_h
#define _misc_container_gnuset_h

#ifdef __GNUG__
#pragma interface
#endif

#include <Pix.h>
#include <builtin.h>

template <class T>
class Set
{
  protected:

  public:
    Set() {};
    virtual              ~Set() {};
    
    // current number of items
    virtual int length() = 0;
    int empty() { return length() == 0; }

    // add item; return Pix
    virtual Pix add(T& item) = 0;
    // delete item
    virtual void del(T& item) = 0;
    // is item in set?
    virtual int contains(T& item) {  return seek(item) != 0; }

    // delete all items
    virtual void clear() {
        Pix i = first(); 
        while (i != 0)
          {
            del((*this)(i));
            i = first();
          }
      }


    // Pix of first item or 0
    virtual Pix first() = 0;
    // advance to next or 0
    virtual void next(Pix& i) = 0;
    // access item at i
    virtual T& operator () (Pix i) = 0;

    // is i a valid Pix  ?
    virtual int owns(Pix idx) {
        if (idx == 0) return 0;
        Pix i;
        for (i = first(); i; next(i)) if (i == idx) return 1;
        return 0;
      }

    // Pix of item
    virtual Pix seek(T& item) {
        Pix i;
        for (i = first(); i != 0 && !((*this)(i) == item); next(i));
        return i;
      }
    // add all items in b
    void operator |= (Set<T>& b) {
        if (&b != this) {
            Pix i;
            for (i = b.first(); i; b.next(i)) add(b(i));
          }
      }
    // delete items also in b
    void operator -= (Set<T>& b) {
        if (&b == this)
            clear();
        else {
            Pix i;
            for (i = b.first(); i; b.next(i)) del(b(i));
          }
      }
    // delete items not in b
    void operator &= (Set<T>& b) {
        if (&b != this)
          {
            Pix i = first();
            Pix n = i;
            while (i != 0)
              {
                next(n);
                if (b.seek((*this)(i)) == 0) del((*this)(i));
                i = n;
              }
          }
      }

    int operator == (Set<T>& b) {
        int n = length();
        if (n != b.length()) return 0;
        if (n == 0) return 1;
        Pix i = first();
        Pix j = b.first();
        while (n-- > 0)
          {
            if ((b.seek((*this)(i)) == 0) || (seek(b(j)) == 0)) return 0;
            next(i);
            b.next(j);
          }
        return 1;
      }
    int operator != (Set<T>& b) {
        return !(*this == b);
      }
    int operator <= (Set<T>& b) {
        if (length() > b.length()) return 0;
        if (length() == 0) return 1;
        Pix i;
        for (i = first(); i; next(i))
            if (b.seek((*this)(i)) == 0) return 0;
        return 1;
      }

    void error(const char* msg) {
        (*lib_error_handler)("Set", msg);
      }
    // rep invariant
    virtual int           OK() = 0;
};

#endif
