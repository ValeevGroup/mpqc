
#include <util/keyval/keyval.h>
#include "keyvalImplMap.h"

AssignedKeyVal::AssignedKeyVal()
{
  _map = new MAPCTOR;
}

AssignedKeyVal::~AssignedKeyVal()
{
  delete _map;
}

int
AssignedKeyVal::key_exists(const char * key)
{
  KeyValKeyword k(key); 
  int result = _map->contains(k);
  if (!result) {
      seterror(UnknownKeyword);
    }
  else {
      seterror(OK);
    }
  return result;
}

RefKeyValValue
AssignedKeyVal::key_value(const char * key)
{
  KeyValKeyword k(key); 
  if (exists(key)) {
      seterror(OK);
      return _map->operator[](k);
    }
  else {
      seterror(UnknownKeyword);
      return 0;
    }
}

void
AssignedKeyVal::assign(const char*key,const RefKeyValValue& val)
{
  KeyValKeyword k(key);
  _map->operator[](k) = val;
}
void
AssignedKeyVal::assign(const char*key,double val)
{
  assign(key,new KeyValValuedouble(val));
}
void
AssignedKeyVal::assignboolean(const char*key,int val)
{
  assign(key,new KeyValValueboolean(val));
}
void
AssignedKeyVal::assign(const char*key,float val)
{
  assign(key,new KeyValValuefloat(val));
}

void
AssignedKeyVal::assign(const char*key,char val)
{
  assign(key,new KeyValValuechar(val));
}
void
AssignedKeyVal::assign(const char*key,int val)
{
  assign(key,new KeyValValueint(val));
}
void
AssignedKeyVal::assign(const char*key,const char* val)
{
  assign(key,new KeyValValuepchar(val));
}
void
AssignedKeyVal::assign(const char*key,const RefDescribedClass&val)
{
  assign(key,new KeyValValueRefDescribedClass(val));
}

void
AssignedKeyVal::clear()
{
  _map->clear();
}
