//
// keyval.cc
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

extern "C" {
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
}

#include <iostream>

#include <util/class/proxy.h>
#include <util/misc/formio.h>
#include <util/keyval/keyval.h>

using namespace std;
using namespace sc;

////////////////////////////////////////////////////////////////////////
// KeyVal

KeyVal::KeyVal() :
  errcod(OK),
  verbose_(0)
{
}

KeyVal::~KeyVal()
{
}

const char* KeyVal::errormsg(KeyValError err)
  {
  const char* msg1 = "No problem.";
  const char* msg2 = "The keyword was not found.";
  const char* msg3 = "The requested operation failed.";
  const char* msg4 = "The datum is not of the appropiate type.";
  const char* msg5 = "The keyword has no value.";
  const char* invalid = "The KeyValError is invalid.";
  if      (err == OK             ) return msg1;
  else if (err == UnknownKeyword ) return msg2;
  else if (err == OperationFailed) return msg3;
  else if (err == WrongType      ) return msg4;
  else if (err == HasNoValue     ) return msg5;
  else return invalid;
  }
int KeyVal::key_count(const char* key)
  {
  int i=0;
  while(exists(key,i)) i++;
  if (i!=0) seterror(OK);
  return i;
  }

double
KeyVal::key_doublevalue(const char* key, const KeyValValue& def)
{
  Ref<KeyValValue> val(key_value(key,def));
  double result;
  if (val.nonnull()) {
      seterror(val->doublevalue(result));
    }
  else {
      KeyValValue::KeyValValueError err = def.doublevalue(result);
      if (error() == OK) seterror(err);
    }
  return result;
}
int
KeyVal::key_booleanvalue(const char* key, const KeyValValue& def)
{
  Ref<KeyValValue> val(key_value(key,def));
  int result;
  if (val.nonnull()) {
      seterror(val->booleanvalue(result));
    }
  else {
      KeyValValue::KeyValValueError err = def.booleanvalue(result);
      if (error() == OK) seterror(err);
    }
  return result;
}
int
KeyVal::key_intvalue(const char* key, const KeyValValue& def)
{
  Ref<KeyValValue> val(key_value(key,def));
  int result;
  if (val.nonnull()) {
      seterror(val->intvalue(result));
    }
  else {
      KeyValValue::KeyValValueError err = def.intvalue(result);
      if (error() == OK) seterror(err);
    }
  return result;
}
long
KeyVal::key_longvalue(const char* key, const KeyValValue& def)
{
  Ref<KeyValValue> val(key_value(key,def));
  long result;
  if (val.nonnull()) {
      seterror(val->longvalue(result));
    }
  else {
      KeyValValue::KeyValValueError err = def.longvalue(result);
      if (error() == OK) seterror(err);
    }
  return result;
}
size_t
KeyVal::key_sizevalue(const char* key, const KeyValValue& def)
{
  Ref<KeyValValue> val(key_value(key,def));
  size_t result;
  if (val.nonnull()) {
      seterror(val->sizevalue(result));
    }
  else {
      KeyValValue::KeyValValueError err = def.sizevalue(result);
      if (error() == OK) seterror(err);
    }
  return result;
}
float
KeyVal::key_floatvalue(const char* key, const KeyValValue& def)
{
  Ref<KeyValValue> val(key_value(key,def));
  float result;
  if (val.nonnull()) {
      seterror(val->floatvalue(result));
    }
  else {
      KeyValValue::KeyValValueError err = def.floatvalue(result);
      if (error() == OK) seterror(err);
    }
  return result;
}
char
KeyVal::key_charvalue(const char* key, const KeyValValue& def)
{
  Ref<KeyValValue> val(key_value(key,def));
  char result;
  if (val.nonnull()) {
      seterror(val->charvalue(result));
    }
  else {
      KeyValValue::KeyValValueError err = def.charvalue(result);
      if (error() == OK) seterror(err);
    }
  return result;
}
char*
KeyVal::key_pcharvalue(const char* key, const KeyValValue& def)
{
  Ref<KeyValValue> val(key_value(key,def));
  const char* result;
  if (val.nonnull()) {
      seterror(val->pcharvalue(result));
    }
  else {
      KeyValValue::KeyValValueError err = def.pcharvalue(result);
      if (error() == OK) seterror(err);
    }
  if (result) return strcpy(new char[strlen(result)+1], result);
  else return 0;
}
std::string
KeyVal::key_stringvalue(const char* key, const KeyValValue& def)
{
  Ref<KeyValValue> val(key_value(key,def));
  std::string result;
  if (val.nonnull()) {
      seterror(val->stringvalue(result));
    }
  else {
      KeyValValue::KeyValValueError err = def.stringvalue(result);
      if (error() == OK) seterror(err);
    }
  return result;
}
Ref<DescribedClass>
KeyVal::key_describedclassvalue(const char* key, const KeyValValue& def)
{
  Ref<KeyValValue> val(key_value(key,def));
  Ref<DescribedClass> result;
  if (val.nonnull()) {
      seterror(val->describedclassvalue(result));
      val = 0; // fix for gcc 2.7.0 bug
    }
  else {
      KeyValValue::KeyValValueError err = def.describedclassvalue(result);
      if (error() == OK) seterror(err);
    }
  return result;
}

int
KeyVal::exists(const char*key)
{
  return key_exists(key);
}
int
KeyVal::count(const char*key)
{
  return key_count(key);
}
Ref<KeyValValue>
KeyVal::value(const char*key,const KeyValValue &def)
{
  return key_value(key,def);
}
int
KeyVal::booleanvalue(const char*key,const KeyValValue& def)
{
  return key_booleanvalue(key,def);
}
double
KeyVal::doublevalue(const char*key,const KeyValValue& def)
{
  return key_doublevalue(key,def);
}
float
KeyVal::floatvalue(const char*key,const KeyValValue& def)
{
  return key_floatvalue(key,def);
}
char
KeyVal::charvalue(const char*key,const KeyValValue& def)
{
  return key_charvalue(key,def);
}
int
KeyVal::intvalue(const char*key,const KeyValValue& def)
{
  return key_intvalue(key,def);
}
long
KeyVal::longvalue(const char*key,const KeyValValue& def)
{
  return key_longvalue(key,def);
}
size_t
KeyVal::sizevalue(const char*key,const KeyValValue& def)
{
  return key_sizevalue(key,def);
}
char*
KeyVal::pcharvalue(const char*key,const KeyValValue& def)
{
  return key_pcharvalue(key,def);
}
std::string
KeyVal::stringvalue(const char*key,const KeyValValue& def)
{
  return key_stringvalue(key,def);
}
Ref<DescribedClass>
KeyVal::describedclassvalue(const char*key,const KeyValValue& def)
{
  return key_describedclassvalue(key,def);
}

const char*
KeyVal::classname(const char * key)
{
  return 0;
}

Ref<DescribedClass>
KeyVal::describedclass(const char* classname)
{
  const ClassDesc* cd = ClassDesc::name_to_class_desc(classname);
  if (!cd) {
      ClassDesc::load_class(classname);
      cd = ClassDesc::name_to_class_desc(classname);

      if (cd == 0) {
        std::ostringstream oss;
        oss << "KeyVal::describedclass is asked to construct an object of unknown type \"" << classname << "\"";
        throw InputError(oss.str().c_str(),
                         __FILE__, __LINE__, "", 0);
      }
    }
  // the original error status must be saved
  KeyValError original_error = error();
  Ref<DescribedClass> newdc(cd->create(this));
  if (newdc.nonnull()) {
    DescribedClassProxy *proxy
    = dynamic_cast<DescribedClassProxy*>(newdc.pointer());
    if (proxy) {
      newdc = proxy->object();
    }
  }
  seterror(original_error);

  return newdc;
}

static void getnewkey(char*newkey,const char*key,int n1)
  {
  if (key) sprintf(newkey,"%s:%d",key,n1);
  else sprintf(newkey,"%d",n1);
  }

static void getnewkey(char*newkey,const char*key,int n1,int n2)
  {
  if (key) sprintf(newkey,"%s:%d:%d",key,n1,n2);
  else  sprintf(newkey,"%d:%d",n1,n2);
  }

//static void getnewkey(char*newkey,const char*key,int n1,int n2,int n3)
//  {
//  if (key) sprintf(newkey,"%s:%d:%d:%d",key,n1,n2,n3);
//  else sprintf(newkey,"%d:%d:%d",n1,n2,n3);
//  }

//static void getnewkey(char*newkey,const char*key,int n1,int n2,int n3,int n4)
//  {
//  if (key) sprintf(newkey,"%s:%d:%d:%d:%d",key,n1,n2,n3,n4);
//  else  sprintf(newkey,"%d:%d:%d:%d",n1,n2,n3,n4);
//  }

// For vectors:
int KeyVal::exists(const char* key,int n1)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1);
  return key_exists(newkey);
  }
int KeyVal::count(const char* key,int n1)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1);
  return key_count(newkey);
  }
double KeyVal::doublevalue(const char* key,int n1,const KeyValValue& def)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1);
  return key_doublevalue(newkey,def);
  }
float KeyVal::floatvalue(const char* key,int n1,const KeyValValue& def)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1);
  return key_floatvalue(newkey,def);
  }
char KeyVal::charvalue(const char* key,int n1,const KeyValValue& def)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1);
  return key_charvalue(newkey,def);
  }
int KeyVal::intvalue(const char* key,int n1,const KeyValValue& def)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1);
  return key_intvalue(newkey,def);
  }
long KeyVal::longvalue(const char* key,int n1,const KeyValValue& def)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1);
  return key_longvalue(newkey,def);
  }
size_t KeyVal::sizevalue(const char* key,int n1,const KeyValValue& def)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1);
  return key_sizevalue(newkey,def);
  }
int KeyVal::booleanvalue(const char* key,int n1,const KeyValValue& def)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1);
  return key_booleanvalue(newkey,def);
  }
char* KeyVal::pcharvalue(const char* key,int n1,const KeyValValue& def)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1);
  return key_pcharvalue(newkey,def);
  }
std::string KeyVal::stringvalue(const char* key,int n1,const KeyValValue& def)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1);
  return key_stringvalue(newkey,def);
  }
Ref<DescribedClass> KeyVal::describedclassvalue(const char* key,int n1,
                                              const KeyValValue& def)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1);
  return key_describedclassvalue(newkey,def);
  }

// For arrays:
int KeyVal::exists(const char* key,int n1,int n2)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1,n2);
  return key_exists(newkey);
  }
int KeyVal::count(const char* key,int n1,int n2)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1,n2);
  return key_count(newkey);
  }
double KeyVal::doublevalue(const char* key,int n1,int n2,
                           const KeyValValue& def)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1,n2);
  return key_doublevalue(newkey,def);
  }
float KeyVal::floatvalue(const char* key,int n1,int n2,
                         const KeyValValue& def)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1,n2);
  return key_floatvalue(newkey,def);
  }
char KeyVal::charvalue(const char* key,int n1,int n2,
                       const KeyValValue& def)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1,n2);
  return key_charvalue(newkey,def);
  }
int KeyVal::intvalue(const char* key,int n1,int n2,
                     const KeyValValue& def)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1,n2);
  return key_intvalue(newkey,def);
  }
long KeyVal::longvalue(const char* key,int n1,int n2,
                       const KeyValValue& def)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1,n2);
  return key_longvalue(newkey,def);
  }
size_t KeyVal::sizevalue(const char* key,int n1,int n2,
                         const KeyValValue& def)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1,n2);
  return key_sizevalue(newkey,def);
  }
int KeyVal::booleanvalue(const char* key,int n1,int n2,
                         const KeyValValue& def)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1,n2);
  return key_booleanvalue(newkey,def);
  }
char* KeyVal::pcharvalue(const char* key,int n1,int n2,
                         const KeyValValue& def)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1,n2);
  return key_pcharvalue(newkey,def);
  }
std::string KeyVal::stringvalue(const char* key,int n1,int n2,
                                const KeyValValue& def)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1,n2);
  return key_stringvalue(newkey,def);
  }
Ref<DescribedClass> KeyVal::describedclassvalue(const char* key,int n1,int n2,
                                              const KeyValValue& def)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1,n2);
  return key_describedclassvalue(newkey,def);
  }

// new and improved for the intel, we can once again use va_arg(), so the
// 12 arg limit is gone.  Unfortunately we can't use new here, so the vals
// array is hardwired to 80.  That should suffice in the foreseeable future
//
#define getnewvakey(newkey,key,narg) \
  strcpy(newkey,key); \
  if(narg!=0) { \
    int vals[80]; \
    if(narg > 80) { \
        ExEnv::errn() << "keyval.cc: getnewvakey: too many varargs...sorry" << endl; \
        exit(1); \
      } \
    va_start(args,narg); \
    int i; \
    for(i=0; i < narg; i++) \
      vals[i] = va_arg(args,int); \
    va_end(args); \
    for(i=0; i < narg; i++)  \
      sprintf((newkey+strlen(newkey)),":%d",vals[i]); \
    }

// For all else:
int KeyVal::Va_exists(const char* key,int narg,...)
  {
  va_list args;
  char newkey[MaxKeywordLength];
  getnewvakey(newkey,key,narg);
  return key_exists(newkey);
  }
int KeyVal::Va_count(const char* key,int narg,...)
  {
  va_list args;
  char newkey[MaxKeywordLength];
  getnewvakey(newkey,key,narg);
  return key_count(newkey);
  }
double KeyVal::Va_doublevalue(const char* key,int narg,...)
  {
  va_list args;
  char newkey[MaxKeywordLength];
  getnewvakey(newkey,key,narg);
  return key_doublevalue(newkey,KeyValValuedouble());
  }
float KeyVal::Va_floatvalue(const char* key,int narg,...)
  {
  va_list args;
  char newkey[MaxKeywordLength];
  getnewvakey(newkey,key,narg);
  return key_floatvalue(newkey,KeyValValuefloat());
  }
char KeyVal::Va_charvalue(const char* key,int narg,...)
  {
  va_list args;
  char newkey[MaxKeywordLength];
  getnewvakey(newkey,key,narg);
  return key_charvalue(newkey,KeyValValuechar());
  }
int KeyVal::Va_intvalue(const char* key,int narg,...)
  {
  va_list args;
  char newkey[MaxKeywordLength];
  getnewvakey(newkey,key,narg);
  return key_intvalue(newkey,KeyValValueint());
  }
long KeyVal::Va_longvalue(const char* key,int narg,...)
  {
  va_list args;
  char newkey[MaxKeywordLength];
  getnewvakey(newkey,key,narg);
  return key_longvalue(newkey,KeyValValuelong());
  }
size_t KeyVal::Va_sizevalue(const char* key,int narg,...)
  {
  va_list args;
  char newkey[MaxKeywordLength];
  getnewvakey(newkey,key,narg);
  return key_sizevalue(newkey,KeyValValuesize());
  }
char* KeyVal::Va_pcharvalue(const char* key,int narg,...)
  {
  va_list args;
  char newkey[MaxKeywordLength];
  getnewvakey(newkey,key,narg);
  return key_pcharvalue(newkey,KeyValValuestring());
  }
std::string KeyVal::Va_stringvalue(const char* key,int narg,...)
  {
  va_list args;
  char newkey[MaxKeywordLength];
  getnewvakey(newkey,key,narg);
  return key_stringvalue(newkey,KeyValValuestring());
  }
Ref<DescribedClass> KeyVal::Va_describedclassvalue(const char* key,int narg,...)
  {
  va_list args;
  char newkey[MaxKeywordLength];
  getnewvakey(newkey,key,narg);
  return key_describedclassvalue(newkey,KeyValValueRefDescribedClass());
  }

void KeyVal::errortrace(ostream&fp)
{
  fp << indent << "KeyVal: error: \"" << errormsg() << "\"" << endl;
}

void KeyVal::dump(ostream&fp)
{
  fp << indent << "KeyVal: error: \"" << errormsg() << "\"" << endl;
}

void KeyVal::print_unseen(ostream&fp)
{
  fp << indent << "(this keyval does not record unread variables)" << endl;
}

int KeyVal::have_unseen()
{
  return -1;
}

void
KeyVal::seterror(KeyValValue::KeyValValueError e)
{
  if (e == KeyValValue::OK) {
      seterror(KeyVal::OK);
    }
  else if (e == KeyValValue::WrongType) {
      seterror(KeyVal::WrongType);
    }
  else {
      // shouldn't get here
      seterror(KeyVal::OperationFailed);
    }
}

// here are some inline candidates that are here for now because
// they were making executables big
void
KeyVal::seterror(KeyValError err)
{
  errcod = err;
}

int
KeyVal::exists(int i)
{
  return exists((const char*)0,i);
}

int
KeyVal::count(int i)
{
  return count((const char*)0,i);
}

int
KeyVal::booleanvalue(int i,const KeyValValue& def)
{
  return booleanvalue((const char*)0,i,def);
}

double
KeyVal::doublevalue(int i,const KeyValValue& def)
{
  return doublevalue((const char*)0,i,def);
}

float
KeyVal::floatvalue(int i,const KeyValValue& def)
{
  return floatvalue((const char*)0,i,def);
}

char
KeyVal::charvalue(int i,const KeyValValue& def)
{
  return charvalue((const char*)0,i,def);
}

int
KeyVal::intvalue(int i,const KeyValValue& def)
{
  return intvalue((const char*)0,i,def);
}

long
KeyVal::longvalue(int i,const KeyValValue& def)
{
  return longvalue((const char*)0,i,def);
}

size_t
KeyVal::sizevalue(int i,const KeyValValue& def)
{
  return sizevalue((const char*)0,i,def);
}

char*
KeyVal::pcharvalue(int i,const KeyValValue& def)
{
  return pcharvalue((const char*)0,i,def);
}

std::string
KeyVal::stringvalue(int i,const KeyValValue& def)
{
  return stringvalue((const char*)0,i,def);
}

Ref<DescribedClass>
KeyVal::describedclassvalue(int i,const KeyValValue& def)
{
  return describedclassvalue((const char*)0,i,def);
}

int
KeyVal::exists(int i,int j)
{
  return exists((const char*)0,i,j);
}

int
KeyVal::count(int i,int j)
{
  return count((const char*)0,i,j);
}

int
KeyVal::booleanvalue(int i,int j,const KeyValValue& def)
{
  return booleanvalue((const char*)0,i,j,def);
}

double
KeyVal::doublevalue(int i,int j,const KeyValValue& def)
{
  return doublevalue((const char*)0,i,j,def);
}

float
KeyVal::floatvalue(int i,int j,const KeyValValue& def)
{
  return floatvalue((const char*)0,i,j,def);
}

char
KeyVal::charvalue(int i,int j,const KeyValValue& def)
{
  return charvalue((const char*)0,i,j,def);
}

int
KeyVal::intvalue(int i,int j,const KeyValValue& def)
{
  return intvalue((const char*)0,i,j,def);
}

long
KeyVal::longvalue(int i,int j,const KeyValValue& def)
{
  return longvalue((const char*)0,i,j,def);
}

size_t
KeyVal::sizevalue(int i,int j,const KeyValValue& def)
{
  return sizevalue((const char*)0,i,j,def);
}

char*
KeyVal::pcharvalue(int i,int j,const KeyValValue& def)
{
  return pcharvalue((const char*)0,i,j,def);
}

std::string
KeyVal::stringvalue(int i,int j,const KeyValValue& def)
{
  return stringvalue((const char*)0,i,j,def);
}

Ref<DescribedClass>
KeyVal::describedclassvalue(int i,int j,const KeyValValue& def)
{
  return describedclassvalue((const char*)0,i,j,def);
}

KeyVal::KeyValError
KeyVal::error()
{
  return errcod;
}

const char*
KeyVal::errormsg()
{
  return errormsg(errcod);
}
