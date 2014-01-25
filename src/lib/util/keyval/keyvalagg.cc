//
// keyvalagg.cc
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
}
#include <util/misc/formio.h>
#include <util/keyval/keyval.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////
// AggregateKeyVal

AggregateKeyVal::AggregateKeyVal(const Ref<KeyVal>&kv0)
{
  kv[0] = kv0;
  kv[1] = 0;
  kv[2] = 0;
  kv[3] = 0;
}

AggregateKeyVal::AggregateKeyVal(const Ref<KeyVal>&kv0,const Ref<KeyVal>&kv1)
{
  kv[0] = kv0;
  kv[1] = kv1;
  kv[2] = 0;
  kv[3] = 0;
}

AggregateKeyVal::AggregateKeyVal(const Ref<KeyVal>&kv0,const Ref<KeyVal>&kv1,
                                 const Ref<KeyVal>&kv2)
{
  kv[0] = kv0;
  kv[1] = kv1;
  kv[2] = kv2;
  kv[3] = 0;
}

AggregateKeyVal::AggregateKeyVal(const Ref<KeyVal>&kv0,const Ref<KeyVal>&kv1,
                                 const Ref<KeyVal>&kv2,const Ref<KeyVal>&kv3)
{
  kv[0] = kv0;
  kv[1] = kv1;
  kv[2] = kv2;
  kv[3] = kv3;
}

AggregateKeyVal::~AggregateKeyVal()
{
}

Ref<KeyVal>
AggregateKeyVal::getkeyval(const char* keyword)
{
  Ref<KeyVal> lastkeyval;
  for (int i=0; i<MaxKeyVal && kv[i].nonnull(); i++) {
      kv[i]->exists(keyword);
      seterror(kv[i]->error());
      if (error() != KeyVal::UnknownKeyword) return kv[i];
      lastkeyval = kv[i];
    }
  // The last keyval in the list is used to lookup the value
  // if the keyword is not found.  This only affects printing
  // in verbose keyvals.
  return lastkeyval;
}

Ref<KeyValValue>
AggregateKeyVal::key_value(const char*arg, const KeyValValue &def)
{
  Ref<KeyVal> kval = getkeyval(arg);
  if (kval.nonnull()) return kval->key_value(arg,def);
  else return 0;
}

int
AggregateKeyVal::key_exists(const char* key)
{
  Ref<KeyVal> kval = getkeyval(key);
  if (kval.nonnull()) return kval->exists(key);
  else return 0;
}

const char*
AggregateKeyVal::classname(const char * key)
{
  Ref<KeyVal> kval = getkeyval(key);
  if (kval.nonnull()) return kval->classname(key);
  else return 0;
}

void
AggregateKeyVal::errortrace(ostream&fp)
{
  fp << indent << "AggregateKeyVal: error: \"" << errormsg() << "\"" << endl;
  for (int i = 0; i<4; i++) {
      if (kv[i].nonnull()) {
          fp << indent << "  KeyVal #" << i << ":" << endl;
          fp << incindent;
          kv[i]->errortrace(fp);
          fp << decindent;
        }
    }
}

void
AggregateKeyVal::dump(ostream&fp)
{
  fp << indent << "AggregateKeyVal: error: \"" << errormsg() << "\"" << endl;
  for (int i = 0; i<4; i++) {
      if (kv[i].nonnull()) {
          fp << indent << "  KeyVal #" << i << ":" << endl;
          fp << incindent;
          kv[i]->dump(fp);
          fp << decindent;
        }
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
