//
// artem.h --- template for the Array class
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

#ifdef __GNUG__
#pragma implementation
#endif

template <class Type>
class Array {
  protected:
    int _length;
    int _managed;
    Type * _array;
  public:
    Array():_length(0),_managed(0),_array(0) {}
    Array(const Array<Type>&a):_length(0),_managed(0),_array(0) {operator=(a);}
    Array(Type* data,int size):_length(size),_managed(0),_array(data){}
    Array(int size):_length(0),_managed(0),_array(0) { set_length(size); }
    ~Array() { clear(); }
    //int length() const { return _length; };
    int size() const { return _length; };
    void clear() { set_length(0); }
    void set_length(int size) {
        if (_array && _managed) delete[] _array;
        _managed = 1;
        if (size) _array = new Type [ size ];
        else _array = 0;
        _length = size;
      }
    void resize(int size) {
        Type*tmp=_array;
        if (size) _array = new Type [ size ];
        else _array = 0;
        int maxi;
        if (size < _length) maxi = size;
        else maxi = _length;
        for (int i=0; i<maxi; i++) _array[i] = tmp[i];
        if (_managed && tmp) delete[] tmp;
        _managed = 1;
        _length = size;
      }
    Array<Type>& operator = (const Array<Type> & s) {
        if (_managed && _array) delete[] _array;
        _managed = 1;
        _length = s._length;
        if (_length) _array = new Type [ _length ];
        else _array = 0;
        for (int i=0; i<_length; i++) {
            _array[i] = s._array[i];
          }
        return (*this);
      }
    Type& operator[] (int i) const {
        if (i<0 || i>=_length) {
            cerr << "Array::operator[](" << i << ") "
                 << "out of range (" << _length << "0" << endl;
            abort();
          };
        return _array[i];
      }
    Type& operator() (int i) const {
        if (i<0 || i>=_length) {
            cerr << "Array::operator()(" << i << ") "
                 << "out of range (" << _length << "0" << endl;
            abort();
          };
        return _array[i];
      }
    void push_back(const Type &d) {
        resize(_length+1);
        _array[_length-1] = d;
    }
    void pop_back() {
        resize(_length-1);
    }
};

// ///////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
