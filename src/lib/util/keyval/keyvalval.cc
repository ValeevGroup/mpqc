//
// keyvalval.cc
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

#include <string.h>
#include <ctype.h>
#include <math.h>

#include <util/keyval/keyvalval.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////////

KeyValValue::KeyValValue(const KeyValValue&)
{
}

KeyValValue::~KeyValValue()
{
}

KeyValValue::KeyValValueError
KeyValValue::doublevalue(double& val) const
{
  KeyValValuedouble def;
  def.doublevalue(val);
  return KeyValValue::WrongType;
}

KeyValValue::KeyValValueError
KeyValValue::booleanvalue(int& val) const
{
  KeyValValueboolean def;
  def.booleanvalue(val);
  return KeyValValue::WrongType;
}

KeyValValue::KeyValValueError
KeyValValue::floatvalue(float& val) const
{
  KeyValValuefloat def;
  def.floatvalue(val);
  return KeyValValue::WrongType;
}

KeyValValue::KeyValValueError
KeyValValue::charvalue(char& val) const
{
  KeyValValuechar def;
  def.charvalue(val);
  return KeyValValue::WrongType;
}

KeyValValue::KeyValValueError
KeyValValue::intvalue(int& val) const
{
  KeyValValueint def;
  def.intvalue(val);
  return KeyValValue::WrongType;
}

KeyValValue::KeyValValueError
KeyValValue::longvalue(long& val) const
{
  KeyValValuelong def;
  def.longvalue(val);
  return KeyValValue::WrongType;
}

KeyValValue::KeyValValueError
KeyValValue::sizevalue(size_t& val) const
{
  KeyValValuesize def;
  def.sizevalue(val);
  return KeyValValue::WrongType;
}

KeyValValue::KeyValValueError
KeyValValue::pcharvalue(const char*& val) const
{
  KeyValValuestring def;
  def.pcharvalue(val);
  return KeyValValue::WrongType;
}

KeyValValue::KeyValValueError
KeyValValue::stringvalue(std::string& val) const
{
  KeyValValuestring def;
  def.stringvalue(val);
  return KeyValValue::WrongType;
}

KeyValValue::KeyValValueError
KeyValValue::describedclassvalue(Ref<DescribedClass>& val) const
{
  KeyValValueRefDescribedClass def;
  def.describedclassvalue(val);
  return KeyValValue::WrongType;
}

void
KeyValValue::print(ostream&o) const
{
  o << "(empty value)";
}

ostream&
sc::operator << (ostream&o, const KeyValValue &val)
{
  val.print(o);
  return o;
}

/////////////////////////////////////////////////////////////////////////

KeyValValuedouble::KeyValValuedouble(const KeyValValuedouble&val):
  _val(val._val)
{
}
KeyValValuedouble::~KeyValValuedouble()
{
}

KeyValValue::KeyValValueError
KeyValValuedouble::doublevalue(double&val) const
{
  val = _val;
  return KeyValValue::OK;
}

void
KeyValValuedouble::print(ostream&o) const
{
  o << _val;
}

/////////////////////////////////////////////////////////////////////////

KeyValValueboolean::KeyValValueboolean(const KeyValValueboolean&val):
  _val(val._val)
{
}
KeyValValueboolean::~KeyValValueboolean()
{
}

KeyValValue::KeyValValueError
KeyValValueboolean::booleanvalue(int&val) const
{
  val = _val;
  return KeyValValue::OK;
}

void
KeyValValueboolean::print(ostream&o) const
{
  o << (_val?"true":"false");
}

/////////////////////////////////////////////////////////////////////////

KeyValValuefloat::KeyValValuefloat(const KeyValValuefloat&val):
  _val(val._val)
{
}
KeyValValuefloat::~KeyValValuefloat()
{
}

KeyValValue::KeyValValueError
KeyValValuefloat::floatvalue(float&val) const
{
  val = _val;
  return KeyValValue::OK;
}

void
KeyValValuefloat::print(ostream&o) const
{
  o << _val;
}

/////////////////////////////////////////////////////////////////////////

KeyValValuechar::KeyValValuechar(const KeyValValuechar&val):
  _val(val._val)
{
}
KeyValValuechar::~KeyValValuechar()
{
}

KeyValValue::KeyValValueError
KeyValValuechar::charvalue(char&val) const
{
  val = _val;
  return KeyValValue::OK;
}

void
KeyValValuechar::print(ostream&o) const
{
  o << _val;
}

/////////////////////////////////////////////////////////////////////////

KeyValValueint::KeyValValueint(const KeyValValueint&val):
  _val(val._val)
{
}
KeyValValueint::~KeyValValueint()
{
}

KeyValValue::KeyValValueError
KeyValValueint::intvalue(int&val) const
{
  val = _val;
  return KeyValValue::OK;
}

void
KeyValValueint::print(ostream&o) const
{
  o << _val;
}

/////////////////////////////////////////////////////////////////////////

KeyValValuelong::KeyValValuelong(const KeyValValuelong&val):
  _val(val._val)
{
}
KeyValValuelong::~KeyValValuelong()
{
}

KeyValValue::KeyValValueError
KeyValValuelong::longvalue(long&val) const
{
  val = _val;
  return KeyValValue::OK;
}

void
KeyValValuelong::print(ostream&o) const
{
  o << _val;
}

/////////////////////////////////////////////////////////////////////////

KeyValValuesize::KeyValValuesize(const KeyValValuesize&val):
  _val(val._val)
{
}
KeyValValuesize::~KeyValValuesize()
{
}

KeyValValue::KeyValValueError
KeyValValuesize::sizevalue(size_t&val) const
{
  val = _val;
  return KeyValValue::OK;
}

void
KeyValValuesize::print(ostream&o) const
{
  o << _val;
}

/////////////////////////////////////////////////////////////////////////

KeyValValuepchar::KeyValValuepchar(const char* val):
  _val(strcpy(new char[strlen(val)+1],val))
{
}
KeyValValuepchar::KeyValValuepchar(const KeyValValuepchar&val):
  _val(strcpy(new char[strlen(val._val)+1],val._val))
{
}
KeyValValuepchar::~KeyValValuepchar()
{
  delete[] _val;
}
KeyValValue::KeyValValueError
KeyValValuepchar::pcharvalue(const char*&val) const
{
  val = _val;
  return KeyValValue::OK;
}
KeyValValue::KeyValValueError
KeyValValuepchar::stringvalue(std::string&val) const
{
  val = _val;
  return KeyValValue::OK;
}

void
KeyValValuepchar::print(ostream&o) const
{
  if (_val == 0) {
      o << "(null)";
    }
  else {
      o << _val;
    }
}

/////////////////////////////////////////////////////////////////////////

KeyValValuestring::KeyValValuestring(const std::string &val):
  _defined(true),
  _val(val)
{
}
KeyValValuestring::KeyValValuestring(const KeyValValuestring&val):
  _defined(true),
  _val(val._val)
{
}
KeyValValuestring::~KeyValValuestring()
{
}
KeyValValue::KeyValValueError
KeyValValuestring::doublevalue(double&val) const
{
  val = atof(_val.c_str());
  return KeyValValue::OK;
}
KeyValValue::KeyValValueError
KeyValValuestring::booleanvalue(int&val) const
{
  char lc_kv[20];
  strncpy(lc_kv,_val.c_str(),20);
  for (int i=0; i<20; i++) {
      if (isupper(lc_kv[i])) lc_kv[i] = tolower(lc_kv[i]);
    }
  if (!strcmp(lc_kv,"yes")) val = 1;
  else if (!strcmp(lc_kv,"true")) val = 1;
  else if (!strcmp(lc_kv,"1")) val = 1;
  else if (!strcmp(lc_kv,"no")) val = 0;
  else if (!strcmp(lc_kv,"false")) val = 0;
  else if (!strcmp(lc_kv,"0")) val = 0;
  else {
      val = 0;
      return KeyValValue::WrongType;
    }

  return KeyValValue::OK;
}
KeyValValue::KeyValValueError
KeyValValuestring::floatvalue(float&val) const
{
  val = (float) atof(_val.c_str());
  return KeyValValue::OK;
}
KeyValValue::KeyValValueError
KeyValValuestring::charvalue(char&val) const
{
  val = _val[0];
  return KeyValValue::OK;
}
KeyValValue::KeyValValueError
KeyValValuestring::intvalue(int&val) const
{
  val = atoi(_val.c_str());
  return KeyValValue::OK;
}
KeyValValue::KeyValValueError
KeyValValuestring::longvalue(long&val) const
{
  val = atol(_val.c_str());
  return KeyValValue::OK;
}
KeyValValue::KeyValValueError
KeyValValuestring::sizevalue(size_t&val) const
{
  int n = ::strlen(_val.c_str());
  int gotdigitspace = 0;
  int gotdigit = 0;
  int gotdecimal = 0;
  int denom = 1;
  double dval = 0;
  for (int i=0; i<n; i++) {
      if (isdigit(_val[i]) && !gotdigitspace) {
          char tmp[2]; tmp[0] = _val[i]; tmp[1] = '\0';
          dval = dval * 10 + atoi(tmp);
          gotdigit = 1;
          if (gotdecimal) denom *= 10;
        }
      else if (_val[i] == '.' && !gotdigitspace && !gotdecimal) {
          gotdecimal = 1;
        }
      else if (_val[i] == ' ') {
          if (gotdigit) gotdigitspace = 1;
        }
      else if (strcmp(&_val[i],"KIB") == 0
               || strcmp(&_val[i],"KiB") == 0
               || strcmp(&_val[i],"kIB") == 0
               || strcmp(&_val[i],"kiB") == 0) {
          dval *= 1024;
          i+=2;
        }
      else if (strcmp(&_val[i],"MIB") == 0
               || strcmp(&_val[i],"MiB") == 0) {
          dval *= 1048576;
          i+=2;
        }
      else if (strcmp(&_val[i],"GIB") == 0
               ||strcmp(&_val[i],"GiB") == 0) {
          dval *= 1073741824;
          i+=2;
        }
      else if (strcmp(&_val[i],"KB") == 0
               || strcmp(&_val[i],"kB") == 0) {
          dval *= 1000;
          i++;
        }
      else if (strcmp(&_val[i],"MB") == 0) {
          dval *= 1000000;
          i++;
        }
      else if (strcmp(&_val[i],"GB") == 0) {
          dval *= 1000000000;
          i++;
        }
      else {
          val = 0;
          return KeyValValue::WrongType;
        }
    }
  dval /= denom;
  val = size_t(dval);
  return KeyValValue::OK;
}
KeyValValue::KeyValValueError
KeyValValuestring::pcharvalue(const char*&val) const
{
  if (!_defined) val = 0;
  else val = _val.c_str();
  return KeyValValue::OK;
}
KeyValValue::KeyValValueError
KeyValValuestring::stringvalue(std::string&val) const
{
  val = _val;
  return KeyValValue::OK;
}

void
KeyValValuestring::print(ostream&o) const
{
  o << _val;
}

/////////////////////////////////////////////////////////////////////////

KeyValValueRefDescribedClass::
  KeyValValueRefDescribedClass(const KeyValValueRefDescribedClass& val):
  _val(val._val)
{
}
KeyValValueRefDescribedClass::
  ~KeyValValueRefDescribedClass()
{
}
KeyValValue::KeyValValueError
KeyValValueRefDescribedClass::describedclassvalue(Ref<DescribedClass>&val) const
{
  val = _val;
  return KeyValValue::OK;
}

void
KeyValValueRefDescribedClass::print(ostream&o) const
{
  if (_val.nonnull()) {
      o << "<" << _val->class_name()
        << " " << _val->identifier()
        << ">";
    }
  else {
      o << "<empty object>";
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
