
extern "C" {
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
}
#include "keyval.h"

//////////////////////////////////////////////////////////////////////
// utility functions

static void getnewkey(char*newkey,const char*key,int n1)
  {
  sprintf(newkey,"%s:%d",key,n1);
  }

static void getnewkey(char*newkey,const char*key,int n1,int n2)
  {
  sprintf(newkey,"%s:%d:%d",key,n1,n2);
  }

static void getnewkey(char*newkey,const char*key,int n1,int n2,int n3)
  {
  sprintf(newkey,"%s:%d:%d:%d",key,n1,n2,n3);
  }

static void getnewkey(char*newkey,const char*key,int n1,int n2,int n3,int n4)
  {
  sprintf(newkey,"%s:%d:%d:%d:%d",key,n1,n2,n3,n4);
  }

///////////////////////////////////////////////////////////////////////
// PrefixKeyVal

PrefixKeyVal::PrefixKeyVal(const char *prefix,KeyVal&kv,int n1,int n2,int n3,int n4):
keyval(&kv)
{
  setup(prefix,4,n1,n2,n3,n4);
}

PrefixKeyVal::PrefixKeyVal(const char *prefix,KeyVal&kv,int n1,int n2,int n3):
keyval(&kv)
{
  setup(prefix,3,n1,n2,n3,0);
}

PrefixKeyVal::PrefixKeyVal(const char *prefix,KeyVal&kv,int n1,int n2):
keyval(&kv)
{
  setup(prefix,2,n1,n2,0,0);
}

PrefixKeyVal::PrefixKeyVal(const char *prefix,KeyVal&kv,int n1):
keyval(&kv)
{
  setup(prefix,1,n1,0,0,0);
}

PrefixKeyVal::PrefixKeyVal(const char *prefix,KeyVal&kv):
keyval(&kv)
{
  setup(prefix,0,0,0,0,0);
}

void PrefixKeyVal::setup(const char*prefix,int n_dim,int n1,int n2,int n3,int n4)
{
  char* prefix_ = strdup(prefix);
  char* token;
  int i;
  for (i=0,token = strtok(prefix_," "); token; i++,token=strtok(0," "));
  nprefix = i;
  free(prefix_);
  if (!nprefix) {
      prefices = 0;
    }
  else {
      prefix_ = strdup(prefix);
      prefices = new char*[nprefix];
      for (i=0,token = strtok(prefix_," "); token; i++,token=strtok(0," ")) {
        char newtoken[MaxKeywordLength];
        if (n_dim == 0) strcpy(newtoken,token);
        else if (n_dim == 1) getnewkey(newtoken,token,n1);
        else if (n_dim == 2) getnewkey(newtoken,token,n1,n2);
        else if (n_dim == 3) getnewkey(newtoken,token,n1,n2,n3);
        else if (n_dim == 4) getnewkey(newtoken,token,n1,n2,n3,n4);
        prefices[i] = new char[strlen(newtoken)+1];
        strcpy(prefices[i],newtoken);
        }
      free(prefix_);
    }
  return;
}

PrefixKeyVal::~PrefixKeyVal()
{
  for (int i=0; i<nprefix; i++) {
      free(prefices[i]);
    }
  delete[] prefices;
}

void PrefixKeyVal::errortrace(FILE*fp,int n)
{
  offset(fp,n); fprintf(fp,"PrefixKeyVal: error: \"%s\"\n",errormsg());
  offset(fp,n); fprintf(fp,"  prefixes:\n");
  for (int i=0; i<nprefix; i++) {
      offset(fp,n); fprintf(fp,"    \"%s\"\n",prefices[i]);
    }
  offset(fp,n); fprintf(fp,"  keyval:\n");
  keyval->errortrace(fp,n + OffsetDelta);
}

void PrefixKeyVal::dump(FILE*fp,int n)
{
  offset(fp,n); fprintf(fp,"PrefixKeyVal: error: \"%s\"\n",errormsg());
  offset(fp,n); fprintf(fp,"  prefixes:\n");
  for (int i=0; i<nprefix; i++) {
      offset(fp,n); fprintf(fp,"    \"%s\"\n",prefices[i]);
    }
  offset(fp,n); fprintf(fp,"  keyval:\n");
  keyval->dump(fp,n + OffsetDelta);
}

int PrefixKeyVal::getnewprefixkey(const char*key,char*newkey)
{
  int i;
  int result;

  if (key[0] == ':') {
    strcpy(newkey,key);
    result = keyval->exists(key);
    seterror(keyval->error());
    }
  else {
      for (i=0; i<nprefix; i++) {
          sprintf(newkey,"%s:%s",prefices[i],key);
          result = keyval->exists(newkey);
          seterror(keyval->error());
          if (error() != KeyVal::UnknownKeyword) break;
        }
    }
  return result;
}

RefKeyValValue
PrefixKeyVal::key_value(const char * arg)
{
  char newkey[MaxKeywordLength];
  getnewprefixkey(arg,newkey);
  RefKeyValValue result(keyval->value(newkey));
  seterror(keyval->error());
  return result;
}

int PrefixKeyVal::key_exists(const char* key)
{
  char newkey[MaxKeywordLength];
  return getnewprefixkey(key,newkey);
}
