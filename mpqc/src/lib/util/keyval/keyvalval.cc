
extern "C" {
#include <string.h>
#include <ctype.h>
#include <math.h>
};

#include <iostream.h>
#include <util/keyval/keyval.h>

REF_def(KeyValValue)

KeyValValue::~KeyValValue()
{
}

KeyVal::KeyValError
KeyValValue::value(double& val)
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
KeyValValue::value(float& val)
{
  val = KeyVal::Defaultfloat();
  return KeyVal::WrongType;
}

KeyVal::KeyValError
KeyValValue::value(char& val)
{
  val = KeyVal::Defaultchar();
  return KeyVal::WrongType;
}

KeyVal::KeyValError
KeyValValue::value(int& val)
{
  val = KeyVal::Defaultint();
  return KeyVal::WrongType;
}

KeyVal::KeyValError
KeyValValue::value(char*& val)
{
  val = KeyVal::Defaultpchar();
  return KeyVal::WrongType;
}

KeyVal::KeyValError
KeyValValue::value(RefDescribedClass& val)
{
  val = KeyVal::DefaultRefDescribedClass();
  return KeyVal::WrongType;
}

KeyValValuedouble::KeyValValuedouble(double val):
  _val(val)
{
}

KeyVal::KeyValError
KeyValValuedouble::value(double&val)
{
  val = _val;
  return KeyVal::OK;
}

KeyValValueboolean::KeyValValueboolean(int val):
  _val(val)
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

KeyVal::KeyValError
KeyValValuefloat::value(float&val)
{
  val = _val;
  return KeyVal::OK;
}

KeyValValuechar::KeyValValuechar(char val):
  _val(val)
{
}

KeyVal::KeyValError
KeyValValuechar::value(char&val)
{
  val = _val;
  return KeyVal::OK;
}

KeyValValueint::KeyValValueint(int val):
  _val(val)
{
}

KeyVal::KeyValError
KeyValValueint::value(int&val)
{
  val = _val;
  return KeyVal::OK;
}

KeyValValuepchar::KeyValValuepchar(const char* val):
  _val(strcpy(new char[strlen(val)+1],val))
{
}
KeyValValuepchar::~KeyValValuepchar()
{
  delete[] _val;
}
KeyVal::KeyValError
KeyValValuepchar::value(char*&val)
{
  val = _val;
  return KeyVal::OK;
}

KeyValValueRefDescribedClass::
  KeyValValueRefDescribedClass(RefDescribedClass& val):
  _val(val)
{
}
KeyVal::KeyValError
KeyValValueRefDescribedClass::value(RefDescribedClass&val)
{
  val = _val;
  return KeyVal::OK;
}

KeyValValueString::KeyValValueString(const char* val):
  _val(val)
{
}
KeyVal::KeyValError
KeyValValueString::value(double&val)
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
KeyValValueString::value(float&val)
{
  val = (float) atof(_val);
  return KeyVal::OK;
}
KeyVal::KeyValError
KeyValValueString::value(char&val)
{
  val = _val[0];
  return KeyVal::OK;
}
KeyVal::KeyValError
KeyValValueString::value(int&val)
{
  val = atoi(_val);
  return KeyVal::OK;
}
KeyVal::KeyValError
KeyValValueString::value(char*&val)
{
  val = strcpy(new char[strlen(_val)+1],_val);
  return KeyVal::OK;
}
