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

#include <iostream.h>
#include <util/keyval/keyval.h>

REF_def(KeyValValue)

KeyValValue::KeyValValue()
{
}

KeyValValue::KeyValValue(KeyValValue&)
{
}

KeyValValue::~KeyValValue()
{
}

KeyVal::KeyValError
KeyValValue::doublevalue(double& val)
{
  val = KeyVal::Defaultdouble();
  return KeyVal::WrongType;
}

KeyVal::KeyValError
KeyValValue::booleanvalue(int& val)
{
  val = KeyVal::Defaultboolean();
  return KeyVal::WrongType;
}

KeyVal::KeyValError
KeyValValue::floatvalue(float& val)
{
  val = KeyVal::Defaultfloat();
  return KeyVal::WrongType;
}

KeyVal::KeyValError
KeyValValue::charvalue(char& val)
{
  val = KeyVal::Defaultchar();
  return KeyVal::WrongType;
}

KeyVal::KeyValError
KeyValValue::intvalue(int& val)
{
  val = KeyVal::Defaultint();
  return KeyVal::WrongType;
}

KeyVal::KeyValError
KeyValValue::pcharvalue(const char*& val)
{
  val = KeyVal::Defaultpchar();
  return KeyVal::WrongType;
}

KeyVal::KeyValError
KeyValValue::describedclassvalue(RefDescribedClass& val)
{
  val = KeyVal::DefaultRefDescribedClass();
  return KeyVal::WrongType;
}

KeyValValuedouble::KeyValValuedouble(double val):
  _val(val)
{
}
KeyValValuedouble::KeyValValuedouble(const KeyValValuedouble&val):
  _val(val._val)
{
}
KeyValValuedouble::~KeyValValuedouble()
{
}

KeyVal::KeyValError
KeyValValuedouble::doublevalue(double&val)
{
  val = _val;
  return KeyVal::OK;
}

KeyValValueboolean::KeyValValueboolean(int val):
  _val(val)
{
}
KeyValValueboolean::KeyValValueboolean(const KeyValValueboolean&val):
  _val(val._val)
{
}
KeyValValueboolean::~KeyValValueboolean()
{
}

KeyVal::KeyValError
KeyValValueboolean::booleanvalue(int&val)
{
  val = _val;
  return KeyVal::OK;
}

KeyValValuefloat::KeyValValuefloat(float val):
  _val(val)
{
}
KeyValValuefloat::KeyValValuefloat(const KeyValValuefloat&val):
  _val(val._val)
{
}
KeyValValuefloat::~KeyValValuefloat()
{
}

KeyVal::KeyValError
KeyValValuefloat::floatvalue(float&val)
{
  val = _val;
  return KeyVal::OK;
}

KeyValValuechar::KeyValValuechar(char val):
  _val(val)
{
}
KeyValValuechar::KeyValValuechar(const KeyValValuechar&val):
  _val(val._val)
{
}
KeyValValuechar::~KeyValValuechar()
{
}

KeyVal::KeyValError
KeyValValuechar::charvalue(char&val)
{
  val = _val;
  return KeyVal::OK;
}

KeyValValueint::KeyValValueint(int val):
  _val(val)
{
}
KeyValValueint::KeyValValueint(const KeyValValueint&val):
  _val(val._val)
{
}
KeyValValueint::~KeyValValueint()
{
}

KeyVal::KeyValError
KeyValValueint::intvalue(int&val)
{
  val = _val;
  return KeyVal::OK;
}

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
KeyVal::KeyValError
KeyValValuepchar::pcharvalue(const char*&val)
{
  val = _val;
  return KeyVal::OK;
}

KeyValValueRefDescribedClass::
  KeyValValueRefDescribedClass(const RefDescribedClass& val):
  _val(val)
{
}
KeyValValueRefDescribedClass::
  KeyValValueRefDescribedClass(const KeyValValueRefDescribedClass& val):
  _val(val._val)
{
}
KeyValValueRefDescribedClass::
  ~KeyValValueRefDescribedClass()
{
}
KeyVal::KeyValError
KeyValValueRefDescribedClass::describedclassvalue(RefDescribedClass&val)
{
  val = _val;
  return KeyVal::OK;
}

KeyValValueString::KeyValValueString(const char* val):
  _val(val)
{
}
KeyValValueString::KeyValValueString(const KeyValValueString&val):
  _val(val._val)
{
}
KeyValValueString::~KeyValValueString()
{
}
KeyVal::KeyValError
KeyValValueString::doublevalue(double&val)
{
  val = atof(_val);
  return KeyVal::OK;
}
KeyVal::KeyValError
KeyValValueString::booleanvalue(int&val)
{
  char lc_kv[20];
  strncpy(lc_kv,_val,20);
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
      return KeyVal::WrongType;
    }

  return KeyVal::OK;
}
KeyVal::KeyValError
KeyValValueString::floatvalue(float&val)
{
  val = (float) atof(_val);
  return KeyVal::OK;
}
KeyVal::KeyValError
KeyValValueString::charvalue(char&val)
{
  val = _val[0];
  return KeyVal::OK;
}
KeyVal::KeyValError
KeyValValueString::intvalue(int&val)
{
  val = atoi(_val);
  return KeyVal::OK;
}
KeyVal::KeyValError
KeyValValueString::pcharvalue(const char*&val)
{
  val = _val;
  return KeyVal::OK;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
