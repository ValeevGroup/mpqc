//
// keyval.h
//
// Copyright (C) 1997 Limit Point Systems, Inc.
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

#ifndef _util_keyval_keyvalval_h
#define _util_keyval_keyvalval_h
#ifdef __GNUG__
#pragma interface
#endif

#include <string>

#include <util/class/class.h>

namespace sc {

class KeyValValue: public RefCount {
  public:
    enum KeyValValueError { OK, WrongType };
  public:
    KeyValValue() {}
    KeyValValue(const KeyValValue&);
    virtual ~KeyValValue();
    // return 1 for success 0, if the datum is of the wrong type
    virtual KeyValValue::KeyValValueError doublevalue(double&) const;
    virtual KeyValValue::KeyValValueError booleanvalue(int&) const;
    virtual KeyValValue::KeyValValueError floatvalue(float&) const;
    virtual KeyValValue::KeyValValueError charvalue(char&) const;
    virtual KeyValValue::KeyValValueError intvalue(int&) const;
    virtual KeyValValue::KeyValValueError sizevalue(size_t&) const;
    virtual KeyValValue::KeyValValueError pcharvalue(const char*&) const;
    virtual KeyValValue::KeyValValueError stringvalue(std::string&) const;
    virtual KeyValValue::KeyValValueError describedclassvalue(Ref<DescribedClass>&) const;
    virtual void print(std::ostream &o=ExEnv::out0()) const;
};
std::ostream& operator<<(std::ostream&,const KeyValValue&);



class KeyValValuedouble: public KeyValValue {
  private:
    double _val;
  public:
    KeyValValuedouble(): _val(0.0) {}
    KeyValValuedouble(double v): _val(v) {}
    KeyValValuedouble(const KeyValValuedouble&);
    ~KeyValValuedouble();
    KeyValValue::KeyValValueError doublevalue(double&) const;
    void print(std::ostream &o=ExEnv::out0()) const;
};

class KeyValValueboolean: public KeyValValue {
  private:
    int _val;
  public:
    KeyValValueboolean(): _val(0) {}
    KeyValValueboolean(int v): _val(v) {}
    KeyValValueboolean(const KeyValValueboolean&);
    ~KeyValValueboolean();
    KeyValValue::KeyValValueError booleanvalue(int&) const;
    void print(std::ostream &o=ExEnv::out0()) const;
};

class KeyValValuefloat: public KeyValValue {
  private:
    float _val;
  public:
    KeyValValuefloat(): _val(0.0) {}
    KeyValValuefloat(float v): _val(v) {}
    KeyValValuefloat(const KeyValValuefloat&);
    ~KeyValValuefloat();
    KeyValValue::KeyValValueError floatvalue(float&) const;
    void print(std::ostream &o=ExEnv::out0()) const;
};

class KeyValValuechar: public KeyValValue {
  private:
    char _val;
  public:
    KeyValValuechar(): _val(0) {}
    KeyValValuechar(char v): _val(v) {}
    KeyValValuechar(const KeyValValuechar&);
    ~KeyValValuechar();
    KeyValValue::KeyValValueError charvalue(char&) const;
    void print(std::ostream &o=ExEnv::out0()) const;
};

class KeyValValueint: public KeyValValue {
  private:
    int _val;
  public:
    KeyValValueint(): _val(0) {}
    KeyValValueint(int v): _val(v) {}
    KeyValValueint(const KeyValValueint&);
    ~KeyValValueint();
    KeyValValue::KeyValValueError intvalue(int&) const;
    void print(std::ostream &o=ExEnv::out0()) const;
};

class KeyValValuesize: public KeyValValue {
  private:
    size_t _val;
  public:
    KeyValValuesize(): _val(0) {}
    KeyValValuesize(int v): _val(v) {}
    KeyValValuesize(const KeyValValuesize&);
    ~KeyValValuesize();
    KeyValValue::KeyValValueError sizevalue(size_t&) const;
    void print(std::ostream &o=ExEnv::out0()) const;
};

class KeyValValuepchar: public KeyValValue {
  private:
    char* _val;
  public:
    KeyValValuepchar(): _val(0) {}
    KeyValValuepchar(const char*);
    KeyValValuepchar(const KeyValValuepchar&);
    ~KeyValValuepchar();
    KeyValValue::KeyValValueError pcharvalue(const char*&) const;
    KeyValValue::KeyValValueError stringvalue(std::string&) const;
    void print(std::ostream &o=ExEnv::out0()) const;
};

class KeyValValuestring: public KeyValValue {
  private:
    std::string _val;
  public:
    KeyValValuestring() {}
    KeyValValuestring(const std::string&);
    KeyValValuestring(const KeyValValuestring&);
    ~KeyValValuestring();
    KeyValValue::KeyValValueError pcharvalue(const char*&) const;
    KeyValValue::KeyValValueError stringvalue(std::string&) const;
    void print(std::ostream &o=ExEnv::out0()) const;
};

class KeyValValueRefDescribedClass: public KeyValValue {
  private:
    Ref<DescribedClass> _val;
  public:
    KeyValValueRefDescribedClass() {}
    KeyValValueRefDescribedClass(const Ref<DescribedClass>& v): _val(v) {}
    KeyValValueRefDescribedClass(const KeyValValueRefDescribedClass&);
    ~KeyValValueRefDescribedClass();
    KeyValValue::KeyValValueError describedclassvalue(Ref<DescribedClass>&) const;
    void print(std::ostream &o=ExEnv::out0()) const;
};

class KeyValValueString: public KeyValValue {
  private:
    const char* _val;
    char *_val_to_delete;
  public:
    // Copy = copy the string data
    // Steal = use the passed pointer and delete it in DTOR
    // Use = use the passed pointer but do not delete it
    enum Storage { Copy, Steal, Use };

    KeyValValueString(const char*,
                      KeyValValueString::Storage s = KeyValValueString::Use);
    KeyValValueString(char*,
                      KeyValValueString::Storage s = KeyValValueString::Use);
    KeyValValueString(const KeyValValueString&);
    ~KeyValValueString();
    KeyValValue::KeyValValueError doublevalue(double&) const;
    KeyValValue::KeyValValueError booleanvalue(int&) const;
    KeyValValue::KeyValValueError floatvalue(float&) const;
    KeyValValue::KeyValValueError charvalue(char&) const;
    KeyValValue::KeyValValueError intvalue(int&) const;
    KeyValValue::KeyValValueError sizevalue(size_t&) const;
    KeyValValue::KeyValValueError pcharvalue(const char*&) const;
    KeyValValue::KeyValValueError stringvalue(std::string&) const;
    void print(std::ostream &o=ExEnv::out0()) const;
};

}

#endif /* _KeyVal_h */

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
