
extern "C" {
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>

}
#include <iostream.h>
#include "keyval.h"

#if defined(I860) && !defined(PARAGON)
#include <util/unix/cct_cprot.h>
#endif

////////////////////////////////////////////////////////////////////////
// KeyVal

KeyVal::KeyVal() :
  errcod(OK)
{
}

KeyVal::~KeyVal()
{
}

char* KeyVal::errormsg(KeyValError err)
  {
  char* msg1 = "No problem.";
  char* msg2 = "The keyword was not found.";
  char* msg3 = "The requested operation failed.";
  char* msg4 = "The datum is not of the appropiate type.";
  char* msg5 = "The keyword has no value.";
  char* invalid = "The KeyValError is invalid.";
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
KeyVal::key_doublevalue(const char* key)
{
  RefKeyValValue val(value(key));
  if (val.nonnull()) {
      double result;
      seterror(val->value(result));
      return result;
    }
  return Defaultdouble();
}
int
KeyVal::key_booleanvalue(const char* key)
{
  RefKeyValValue val(value(key));
  if (val.nonnull()) {
      int result;
      seterror(val->booleanvalue(result));
      return result;
    }
  return Defaultboolean();
}
int
KeyVal::key_intvalue(const char* key)
{
  RefKeyValValue val(value(key));
  if (val.nonnull()) {
      int result;
      seterror(val->value(result));
      return result;
    }
  return Defaultint();
}
float
KeyVal::key_floatvalue(const char* key)
{
  RefKeyValValue val(value(key));
  if (val.nonnull()) {
      float result;
      seterror(val->value(result));
      return result;
    }
  return Defaultfloat();
}
char
KeyVal::key_charvalue(const char* key)
{
  RefKeyValValue val(value(key));
  if (val.nonnull()) {
      char result;
      seterror(val->value(result));
      return result;
    }
  return Defaultchar();
}
char*
KeyVal::key_pcharvalue(const char* key)
{
  RefKeyValValue val(value(key));
  if (val.nonnull()) {
      char* result;
      seterror(val->value(result));
      return result;
    }
  return Defaultpchar();
}
RefDescribedClass
KeyVal::key_describedclassvalue(const char* key)
{
  RefKeyValValue val(value(key));
  if (val.nonnull()) {
      RefDescribedClass result;
      seterror(val->value(result));
      return result;
    }
  return DefaultRefDescribedClass();
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
RefKeyValValue
KeyVal::value(const char*key)
{
  return key_value(key);
}
int
KeyVal::booleanvalue(const char*key)
{
  return key_booleanvalue(key);
}
double
KeyVal::doublevalue(const char*key)
{
  return key_doublevalue(key);
}
float
KeyVal::floatvalue(const char*key)
{
  return key_floatvalue(key);
}
char
KeyVal::charvalue(const char*key)
{
  return key_charvalue(key);
}
int
KeyVal::intvalue(const char*key)
{
  return key_intvalue(key);
}
char*
KeyVal::pcharvalue(const char*key)
{
  return key_pcharvalue(key);
}
RefDescribedClass
KeyVal::describedclassvalue(const char*key)
{
  return key_describedclassvalue(key);
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

static void getnewkey(char*newkey,const char*key,int n1,int n2,int n3)
  {
  if (key) sprintf(newkey,"%s:%d:%d:%d",key,n1,n2,n3);
  else sprintf(newkey,"%d:%d:%d",n1,n2,n3);
  }

static void getnewkey(char*newkey,const char*key,int n1,int n2,int n3,int n4)
  {
  if (key) sprintf(newkey,"%s:%d:%d:%d:%d",key,n1,n2,n3,n4);
  else  sprintf(newkey,"%d:%d:%d:%d",n1,n2,n3,n4);
  }

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
double KeyVal::doublevalue(const char* key,int n1)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1);
  return key_doublevalue(newkey);
  }
float KeyVal::floatvalue(const char* key,int n1)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1);
  return key_floatvalue(newkey);
  }
char KeyVal::charvalue(const char* key,int n1)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1);
  return key_charvalue(newkey);
  }
int KeyVal::intvalue(const char* key,int n1)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1);
  return key_intvalue(newkey);
  }
int KeyVal::booleanvalue(const char* key,int n1)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1);
  return key_booleanvalue(newkey);
  }
char* KeyVal::pcharvalue(const char* key,int n1)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1);
  return key_pcharvalue(newkey);
  }
RefDescribedClass KeyVal::describedclassvalue(const char* key,int n1)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1);
  return key_describedclassvalue(newkey);
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
double KeyVal::doublevalue(const char* key,int n1,int n2)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1,n2);
  return key_doublevalue(newkey);
  }
float KeyVal::floatvalue(const char* key,int n1,int n2)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1,n2);
  return key_floatvalue(newkey);
  }
char KeyVal::charvalue(const char* key,int n1,int n2)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1,n2);
  return key_charvalue(newkey);
  }
int KeyVal::intvalue(const char* key,int n1,int n2)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1,n2);
  return key_intvalue(newkey);
  }
int KeyVal::booleanvalue(const char* key,int n1,int n2)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1,n2);
  return key_booleanvalue(newkey);
  }
char* KeyVal::pcharvalue(const char* key,int n1,int n2)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1,n2);
  return key_pcharvalue(newkey);
  }
RefDescribedClass KeyVal::describedclassvalue(const char* key,int n1,int n2)
  {
  char newkey[MaxKeywordLength];
  getnewkey(newkey,key,n1,n2);
  return key_describedclassvalue(newkey);
  }

#if !defined(I860)  || defined(PARAGON)
#define getnewvakey(newkey,key,narg) \
  strcpy(newkey,key); \
  if(narg!=0) { \
    va_start(args,narg); \
    for(int i=0; i < narg; i++)  \
      sprintf((newkey+strlen(newkey)),":%d",va_arg(args,int)); \
    va_end(args); \
    }
#else
// new and improved for the intel, we can once again use va_arg(), so the
// 12 arg limit is gone.  Unfortunately we can't use new here, so the vals
// array is hardwired to 80.  That should suffice in the foreseeable future
//
#define getnewvakey(newkey,key,narg) \
  strcpy(newkey,key); \
  if(narg!=0) { \
    int vals[80]; \
    if(narg > 80) err_quit("getnewvakey: too many varargs for intel...sorry"); \
    va_start(args,narg); \
    for(int i=0; i < narg; i++) \
      vals[i] = va_arg(args,int); \
    va_end(args); \
    for(i=0; i < narg; i++)  \
      sprintf((newkey+strlen(newkey)),":%d",vals[i]); \
    }
#endif

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
  return key_doublevalue(newkey);
  }
float KeyVal::Va_floatvalue(const char* key,int narg,...)
  {
  va_list args;
  char newkey[MaxKeywordLength];
  getnewvakey(newkey,key,narg);
  return key_floatvalue(newkey);
  }
char KeyVal::Va_charvalue(const char* key,int narg,...)
  {
  va_list args;
  char newkey[MaxKeywordLength];
  getnewvakey(newkey,key,narg);
  return key_charvalue(newkey);
  }
int KeyVal::Va_intvalue(const char* key,int narg,...)
  {
  va_list args;
  char newkey[MaxKeywordLength];
  getnewvakey(newkey,key,narg);
  return key_intvalue(newkey);
  }
char* KeyVal::Va_pcharvalue(const char* key,int narg,...)
  {
  va_list args;
  char newkey[MaxKeywordLength];
  getnewvakey(newkey,key,narg);
  return key_pcharvalue(newkey);
  }
RefDescribedClass KeyVal::Va_describedclassvalue(const char* key,int narg,...)
  {
  va_list args;
  char newkey[MaxKeywordLength];
  getnewvakey(newkey,key,narg);
  return key_describedclassvalue(newkey);
  }

void KeyVal::offset(FILE*fp,int n)
{
  for (int i=0; i<n; i++) {
      fprintf(fp," ");
    }
}

void KeyVal::errortrace(FILE*fp,int n)
{
  offset(fp,n); fprintf(fp,"KeyVal: error: \"%s\"\n",errormsg());
}

void KeyVal::dump(FILE*fp,int n)
{
  offset(fp,n); fprintf(fp,"KeyVal: error: \"%s\"\n",errormsg());
}
