
extern "C" {
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
}
#include <util/keyval/keyval.h>

/////////////////////////////////////////////////////////////////////
// AggregateKeyVal

AggregateKeyVal::AggregateKeyVal(const RefKeyVal&kv0)
{
  kv[0] = kv0;
  kv[1] = 0;
  kv[2] = 0;
  kv[3] = 0;
}

AggregateKeyVal::AggregateKeyVal(const RefKeyVal&kv0,const RefKeyVal&kv1)
{
  kv[0] = kv0;
  kv[1] = kv1;
  kv[2] = 0;
  kv[3] = 0;
}

AggregateKeyVal::AggregateKeyVal(const RefKeyVal&kv0,const RefKeyVal&kv1,
                                 const RefKeyVal&kv2)
{
  kv[0] = kv0;
  kv[1] = kv1;
  kv[2] = kv2;
  kv[3] = 0;
}

AggregateKeyVal::AggregateKeyVal(const RefKeyVal&kv0,const RefKeyVal&kv1,
                                 const RefKeyVal&kv2,const RefKeyVal&kv3)
{
  kv[0] = kv0;
  kv[1] = kv1;
  kv[2] = kv2;
  kv[3] = kv3;
}

AggregateKeyVal::~AggregateKeyVal()
{
}

RefKeyVal
AggregateKeyVal::getkeyval(const char* keyword)
{
  for (int i=0; i<MaxKeyVal && kv[i].nonnull(); i++) {
      kv[i]->exists(keyword);
      seterror(kv[i]->error());
      if (error() != KeyVal::UnknownKeyword) return kv[i];
    }
  return 0;
}

RefKeyValValue
AggregateKeyVal::key_value(const char*arg)
{
  KeyVal* kval = getkeyval(arg).pointer();
  if (kval) return kval->value(arg);
  else return 0;
}

int
AggregateKeyVal::key_exists(const char* key)
{
  KeyVal* kval = getkeyval(key).pointer();
  if (kval) return kval->exists(key);
  else return 0;
}

void
AggregateKeyVal::errortrace(ostream&fp,int n)
{
  offset(fp,n);
  fp << "AggregateKeyVal: error: \"" << errormsg() << "\"" << endl;
  for (int i = 0; i<4; i++) {
      if (kv[i].nonnull()) {
          offset(fp,n); fp << "  KeyVal #" << i << ":" << endl;
          kv[i]->errortrace(fp,n+OffsetDelta);
        }
    }
}

void
AggregateKeyVal::dump(ostream&fp,int n)
{
  offset(fp,n);
  fp << "AggregateKeyVal: error: \"" << errormsg() << "\"" << endl;
  for (int i = 0; i<4; i++) {
      if (kv[i].nonnull()) {
          offset(fp,n); fp << "  KeyVal #" << i << ":" << endl;
          kv[i]->dump(fp,n+OffsetDelta);
        }
    }
}
