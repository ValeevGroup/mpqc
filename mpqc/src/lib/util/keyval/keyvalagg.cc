
extern "C" {
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
}
#include "keyval.h"

/////////////////////////////////////////////////////////////////////
// AggregateKeyVal

AggregateKeyVal::AggregateKeyVal(KeyVal&kv0)
{
  kv[0] = &kv0;
  kv[1] = 0;
  kv[2] = 0;
  kv[3] = 0;
}

AggregateKeyVal::AggregateKeyVal(KeyVal&kv0,KeyVal&kv1)
{
  kv[0] = &kv0;
  kv[1] = &kv1;
  kv[2] = 0;
  kv[3] = 0;
}

AggregateKeyVal::AggregateKeyVal(KeyVal&kv0,KeyVal&kv1,KeyVal&kv2)
{
  kv[0] = &kv0;
  kv[1] = &kv1;
  kv[2] = &kv2;
  kv[3] = 0;
}

AggregateKeyVal::AggregateKeyVal(KeyVal&kv0,KeyVal&kv1,KeyVal&kv2,KeyVal&kv3)
{
  kv[0] = &kv0;
  kv[1] = &kv1;
  kv[2] = &kv2;
  kv[3] = &kv3;
}

AggregateKeyVal::~AggregateKeyVal()
{
}

KeyVal*
AggregateKeyVal::getkeyval(const char* keyword)
{
  for (int i=0; i<MaxKeyVal && kv[i]; i++) {
      kv[i]->exists(keyword);
      seterror(kv[i]->error());
      if (error() != KeyVal::UnknownKeyword) return kv[i];
    }
  return 0;
}

RefKeyValValue
AggregateKeyVal::key_value(const char*arg)
{
  KeyVal* kval = getkeyval(arg);
  if (kval) return kval->value(arg);
  else return 0;
}

int
AggregateKeyVal::key_exists(const char* key)
{
  KeyVal* kval = getkeyval(key);
  if (kval) return kval->exists(key);
  else return 0;
}

void
AggregateKeyVal::errortrace(FILE*fp,int n)
{
  offset(fp,n); fprintf(fp,"AggregateKeyVal: error: \"%s\"\n",errormsg());
  for (int i = 0; i<4; i++) {
      if (kv[i]) {
          offset(fp,n); fprintf(fp,"  KeyVal #%d:\n",i);
          kv[i]->errortrace(fp,n+OffsetDelta);
        }
    }
}

void
AggregateKeyVal::dump(FILE*fp,int n)
{
  offset(fp,n); fprintf(fp,"AggregateKeyVal: error: \"%s\"\n",errormsg());
  for (int i = 0; i<4; i++) {
      if (kv[i]) {
          offset(fp,n); fprintf(fp,"  KeyVal #%d:\n",i);
          kv[i]->dump(fp,n+OffsetDelta);
        }
    }
}
