//
// keyvalpre.cc
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
#include <stdlib.h>
#include <stdio.h>
}
#include <util/misc/formio.h>
#include <util/keyval/keyval.h>

using namespace std;
using namespace sc;

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

PrefixKeyVal::PrefixKeyVal(const Ref<KeyVal>&kv,const char *prefix_a,
                           int n1,int n2,int n3,int n4):
keyval(kv)
{
  setup(prefix_a,4,n1,n2,n3,n4);
}

PrefixKeyVal::PrefixKeyVal(const Ref<KeyVal>&kv,const char *prefix_a,
                           int n1,int n2,int n3):
keyval(kv)
{
  setup(prefix_a,3,n1,n2,n3,0);
}

PrefixKeyVal::PrefixKeyVal(const Ref<KeyVal>&kv,const char *prefix_a,
                           int n1,int n2):
keyval(kv)
{
  setup(prefix_a,2,n1,n2,0,0);
}

PrefixKeyVal::PrefixKeyVal(const Ref<KeyVal>&kv,const char *prefix_a,int n1):
keyval(kv)
{
  setup(prefix_a,1,n1,0,0,0);
}

PrefixKeyVal::PrefixKeyVal(const Ref<KeyVal>&kv,const char *prefix_a):
keyval(kv)
{
  setup(prefix_a,0,0,0,0,0);
}

PrefixKeyVal::PrefixKeyVal(const Ref<KeyVal>&kv,int n1):
keyval(kv)
{
  setup(0,1,n1,0,0,0);
}

PrefixKeyVal::PrefixKeyVal(const Ref<KeyVal>&kv,int n1,int n2):
keyval(kv)
{
  setup(0,2,n1,n2,0,0);
}

PrefixKeyVal::PrefixKeyVal(const Ref<KeyVal>&kv,int n1,int n2,int n3):
keyval(kv)
{
  setup(0,3,n1,n2,n3,0);
}

PrefixKeyVal::PrefixKeyVal(const Ref<KeyVal>&kv,int n1,int n2,int n3,int n4):
keyval(kv)
{
  setup(0,4,n1,n2,n3,n4);
}

void PrefixKeyVal::setup(const char*pref,int n_dim,int n1,int n2,int n3,int n4)
{
  if (!pref) {
      prefix = 0;
    }
  else {
      char newtoken[MaxKeywordLength];
      if (n_dim == 0) strcpy(newtoken,pref);
      else if (n_dim == 1) getnewkey(newtoken,pref,n1);
      else if (n_dim == 2) getnewkey(newtoken,pref,n1,n2);
      else if (n_dim == 3) getnewkey(newtoken,pref,n1,n2,n3);
      else if (n_dim == 4) getnewkey(newtoken,pref,n1,n2,n3,n4);
      prefix = new char[strlen(newtoken)+1];
      strcpy(prefix,newtoken);
    }
  return;
}

PrefixKeyVal::~PrefixKeyVal()
{
  if (prefix) {
      delete[] prefix;
      prefix=0;
    }
}

void PrefixKeyVal::errortrace(ostream&fp)
{
  fp << indent << "PrefixKeyVal: error: \"" << errormsg() << "\"" << endl;
  fp << indent << "  prefix:" << endl;
  fp << indent << "    \"" << prefix << "\"" << endl;
  fp << indent << "  keyval:" << endl;
  fp << incindent;
  keyval->errortrace(fp);
  fp << decindent;
}

void PrefixKeyVal::dump(ostream&fp)
{
  fp << indent << "PrefixKeyVal: error: \"" << errormsg() << "\"" << endl;
  fp << indent << "  prefixes:" << endl;
  fp << indent << "    \"" << prefix << "\"" << endl;
  fp << indent << "  keyval:" << endl;
  fp << incindent;
  keyval->dump(fp);
  fp << decindent;
}

int PrefixKeyVal::getnewprefixkey(const char*key,char*newkey)
{
  int result=0;

  if (key[0] == ':') {
      strcpy(newkey,key);
      result = keyval->exists(key);
      seterror(keyval->error());
    }
  else {
      sprintf(newkey,"%s:%s",prefix,key);
      result = keyval->exists(newkey);
      seterror(keyval->error());
    }
  return result;
}

Ref<KeyValValue>
PrefixKeyVal::key_value(const char * arg, const KeyValValue &def)
{
  char newkey[MaxKeywordLength];
  getnewprefixkey(arg,newkey);
  Ref<KeyValValue> result(keyval->key_value(newkey,def));
  seterror(keyval->error());
  return result;
}

int PrefixKeyVal::key_exists(const char* key)
{
  char newkey[MaxKeywordLength];
  return getnewprefixkey(key,newkey);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
