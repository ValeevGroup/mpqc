//
// settmpl.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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

#ifdef __GNUG__
#pragma interface
#endif

template <class Type>
class Set {
  protected:
    int nelement;
    int array_length;
    enum {array_incr=100 };
    Type *element;
    void range_check(int i) const
    {
      if ((i<0) || (i >= nelement)) {
          cerr << "Set::range_check(" << i << "): nelement=" << nelement
               << endl;
          abort();
        }
    }
    void range_check(Pix i) const
    {
      if (!nelement
          || ((Type*)i<&element[0])
          || ((Type*)i>&element[nelement-1])) {
          cerr << "Set::range_check(0x" << setbase(16) << (void*)element
               << "): start=0x" << setbase(16) << (void*)i
               << " nelement=" << sizeof(Type)
               << " size=" << nelement << endl;
          abort();
        }
    }
    Pix index_to_pix(int i) const
    {
      range_check(i);
      return (Pix)&element[i];
    }
  public:
    int length() const { return nelement; }
    Set():nelement(0),array_length(array_incr),element(new Type[array_incr]) {}
    Set(const Set<Type> & s):nelement(0),array_length(0),element(0)
    {
      this->operator = (s);
    }
    Set<Type>& operator = (const Set<Type> & s)
    {
      nelement = s.nelement;
      if (nelement<array_length) {
          if (element) delete[] element;
          array_length = nelement;
          element = new Type[nelement];
        }
      int i;
      for (i=0; i<nelement; i++) element[i] = s.element[i];
      return *this;
    }
    ~Set() { clear(); }
    void clear()
    {
      if (element) { delete[] element; element = 0; }
      array_length = 0;
      nelement = 0;
    }
    Pix add(const Type & e)
    {
      int i;
      for (i=0; i<nelement; i++) {
          if (e == element[i]) {
              return index_to_pix(i);
            }
        }
      if (nelement == array_length) {
          Type* tmpelement = new Type[array_length + array_incr];
          for (i=0; i<nelement; i++) {
              tmpelement[i] = element[i];
            }
          if (element) delete[] element;
          array_length += array_incr;
          element = tmpelement;
        }
      element[nelement] = e;
      nelement++;
      return index_to_pix(nelement-1);
    }
    void del( Type & e)
    {
      int i;
      for (i=0; i<nelement; i++) {
          if (e == element[i]) {
              nelement--;
              int j;
              for (j=i; j<nelement; j++) {
                  element[j] = element[j+1];
                }
              element[j] = 0;
              break;
            }
        }
    }
    Type& operator()(Pix i)
    {
      range_check(i);
      return *(Type*)i;
    }
    const Type& operator()(const Pix i) const
    {
      range_check(i);
      return *(const Type*)i;
    }
    Set<Type>& operator += (const Set<Type>&s)
    {
      for (int i=0; i<s.nelement; i++) {
          add(s.element[i]);
        }
      return *this;
    }
    Pix seek(const Type &item)
    {
      for (int i=0; i<nelement; i++) {
          if (item == element[i]) return index_to_pix(i);
        }
      return 0;
    }
    int contains(const Type &item)
    {
      return seek(item) != 0;
    }
    int owns(Pix i)
    {
      if (!nelement
          || ((Type*)i<&element[0])
          || ((Type*)i>&element[nelement-1]))
          return 0;
      else return 1;
    }
    Pix first()
    {
      if (nelement) return index_to_pix(0);
      else return 0;
    }
    void next(Pix&i)
    {
      Type* t = (Type*) i;
      if (t < &element[nelement-1]) i = (Pix)&t[1];
      else i = 0;
    }
};

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
