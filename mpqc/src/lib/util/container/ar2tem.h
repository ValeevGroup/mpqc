//
// ar2tem.h --- template for the Array2 class
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

template <class Type>
class Array2 {
  protected:
    int _length0;
    int _length1;
    int _managed;
    Type * _array;
  public:
    Array2():_length0(0),_length1(0),_array(0) {}
    Array2(const Array2<Type> &a): _length0(0),_length1(0),_array(0) {
        operator=(a);
      }
    Array2(Type* data,int size0,int size1):
      _length0(size0),_length1(size1),_managed(0),_array(data) {}
    Array2(int size0,int size1): _length0(0),_length1(0),_array(0) {
        set_lengths(size0,size1);
      }
    ~Array2() { clear(); }
    int length0() const { return _length0; };
    int length1() const { return _length1; };
    void clear() { set_lengths(0,0); }
    void set_lengths(int size0,int size1) {
        if (_managed && _array) delete[] _array;
        _managed = 1;
        if (size0*size1) _array = new Type [ size0*size1 ];
        else _array = 0;
        _length0 = size0;
        _length1 = size1;
      }
    Array2<Type>& operator = (const Array2<Type> & s) {
        if (_managed && _array) delete[] _array;
        _managed = 1;
        _length0 = s._length0;
        _length1 = s._length1;
        if (_length0*_length1) _array = new Type [ _length0*_length1 ];
        else _array = 0;
        for (int i=0; i<_length0*_length1; i++) {
            _array[i] = s._array[i];
          }
        return (*this);
      }
    Type& operator() (int i,int j) {
        if (i<0 || i>=_length0 || j<0 || j>=_length1) {
            cerr << "Array2::operator()(" << i << "," << j << "): "
                 << "out of range (" << _length0 << "," << _length1
                 << ",)" << endl;
            abort();
          };
        return _array[i*_length1+j];
      }
    const Type& operator() (int i,int j) const {
        if (i<0 || i>=_length0 || j<0 || j>=_length1) {
            cerr << "Array2::operator()(" << i << "," << j << "): "
                 << "out of range (" << _length0 << "," << _length1
                 << ",)" << endl;
            abort();
          };
        return _array[i*_length1+j];
      }
};

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
