
#ifndef _util_container_array_h
#define _util_container_array_h

#include <stdio.h>
#include <stdlib.h>

#define ARRAY_dec(Type)							      \
class Array ## Type							      \
{									      \
 private:								      \
  int _length;								      \
  Type * _array;							      \
 public:								      \
  Array ## Type();							      \
  Array ## Type(const Array ## Type &);					      \
  Array ## Type(int size);						      \
  ~Array ## Type();							      \
  void set_length(int size);						      \
  void reset_length(int size);						      \
  Array ## Type& operator = (const Array ## Type & s);			      \
  int length() const;							      \
  void clear();								      \
  Type& operator[] (int i) const;					      \
}

#define TMPARRAY_dec(Type)						      \
class TMPArray ## Type							      \
{									      \
 private:								      \
  int _length;								      \
  Type * _array;							      \
 public:								      \
  TMPArray ## Type(Type* data,int size);				      \
  TMPArray ## Type(const TMPArray ## Type &);				      \
  ~TMPArray ## Type();							      \
  Type& operator[] (int i) const;					      \
}

#define ARRAY2_dec(Type)						      \
class Array2 ## Type							      \
{									      \
 private:								      \
  int _length0;								      \
  int _length1;								      \
  Type * _array;							      \
 public:								      \
  Array2 ## Type();							      \
  Array2 ## Type(const Array2 ## Type &);				      \
  Array2 ## Type(int,int);						      \
  ~Array2 ## Type();							      \
  void set_lengths(int,int);						      \
  Array2 ## Type & operator = (const Array2 ## Type & s);		      \
  int length0() const;							      \
  int length1() const;							      \
  void clear();								      \
  TMPArray ## Type operator[] (int i) const;				      \
}

#define ARRAY_def(Type)							      \
  int Array ## Type::length() const { return _length; };		      \
  void Array ## Type::clear() { set_length(0); }			      \
  Array ## Type::Array ## Type():_array(0),_length(0) {}		      \
  Array ## Type::Array ## Type(const Array ## Type&a):_array(0),_length(0) {  \
    operator=(a);							      \
    }									      \
  Array ## Type::Array ## Type(int size):_array(0),_length(0)		      \
  { set_length(size); }							      \
  Array ## Type::~Array ## Type() { clear(); }				      \
  void Array ## Type::set_length(int size)				      \
  {									      \
    if (_array) delete[] _array;					      \
    if (size) _array = new Type [ size ];				      \
    else _array = 0;							      \
    _length = size;							      \
  }									      \
  void Array ## Type::reset_length(int size)				      \
  {									      \
    Type*tmp=_array;							      \
    if (size) _array = new Type [ size ];				      \
    else _array = 0;							      \
    int maxi;								      \
    if (size < _length) maxi = size;					      \
    else maxi = _length;						      \
    for (int i=0; i<maxi; i++) _array[i] = tmp[i];			      \
    if (tmp) delete[] tmp;						      \
    _length = size;							      \
  }									      \
  Array ## Type& Array ## Type::operator = (const Array ## Type & s)	      \
  {									      \
    if (_array) delete[] _array;					      \
    _length = s._length;						      \
    if (_length) _array = new Type [ _length ];				      \
    else _array = 0;							      \
    for (int i=0; i<_length; i++) {					      \
        _array[i] = s._array[i];					      \
      }									      \
    return (*this);							      \
  }									      \
  Type& Array ## Type::operator[] (int i) const				      \
  {									      \
    if (i<0 || i>=_length) {						      \
        fprintf(stderr,"Array::operator[] out of range: %d (nelement = %d)\n", \
                i,_length);						      \
        abort();							      \
      };								      \
    return _array[i];							      \
  }

#define ARRAY2_def(Type)						      \
  int Array2 ## Type::length0() const { return _length0; };		      \
  int Array2 ## Type::length1() const { return _length1; };		      \
  void Array2 ## Type::clear() { set_lengths(0,0); }			      \
  Array2 ## Type::Array2 ## Type():_array(0),_length0(0),_length1(0) {}	      \
  Array2 ## Type::Array2 ## Type(const Array2 ## Type &a):		      \
    _array(0),_length0(0),_length1(0)					      \
    {									      \
      operator=(a);							      \
    }									      \
  Array2 ## Type::Array2 ## Type(int size0,int size1):			      \
  _array(0),_length0(0),_length1(0)					      \
  { set_lengths(size0,size1); }						      \
  Array2 ## Type::~Array2 ## Type() { clear(); }			      \
  void Array2 ## Type::set_lengths(int size0,int size1)			      \
  {									      \
    if (_array) delete[] _array;					      \
    if (size0*size1) _array = new Type [ size0*size1 ];			      \
    else _array = 0;							      \
    _length0 = size0;							      \
    _length1 = size1;							      \
  }									      \
  Array2 ## Type& Array2 ## Type::operator = (const Array2 ## Type & s)	      \
  {									      \
    if (_array) delete[] _array;					      \
    _length0 = s._length0;						      \
    _length1 = s._length1;						      \
    if (_length0*_length1) _array = new Type [ _length0*_length1 ];	      \
    else _array = 0;							      \
    for (int i=0; i<_length0*_length1; i++) {				      \
        _array[i] = s._array[i];					      \
      }									      \
    return (*this);							      \
  }									      \
  TMPArray ## Type Array2 ## Type::operator[] (int i) const		      \
  {									      \
    if (i<0 || i>=_length0) {						      \
        fprintf(stderr,"Array2::operator[] out of range: %d (nelement = %d)\n", \
                i,_length0);						      \
        abort();							      \
      };								      \
    TMPArray ## Type r(&_array[i*_length1],_length1);			      \
    return r;								      \
  }

#define TMPARRAY_def(Type)						      \
  TMPArray ## Type :: TMPArray ## Type(const TMPArray ## Type &a):	      \
    _length(a._length),_array(a._array) {}				      \
  TMPArray ## Type::TMPArray ## Type(Type* data,int size):		      \
    _length(size),_array(data) {}					      \
  TMPArray ## Type::~TMPArray ## Type() {}				      \
  Type& TMPArray ## Type::operator[] (int i) const			      \
  {									      \
    if (i<0 || i>=_length) {						      \
        fprintf(stderr,"TMPArray::operator[] out of range: %d (nelement = %d)\n", \
                i,_length);						      \
        abort();							      \
      };								      \
    return _array[i];							      \
  }

// declare arrays of the basic types
ARRAY_dec(int);
ARRAY_dec(Arrayint);
TMPARRAY_dec(int);
ARRAY2_dec(int);
ARRAY_dec(double);
ARRAY_dec(Arraydouble);
TMPARRAY_dec(double);
ARRAY2_dec(double);

#endif
