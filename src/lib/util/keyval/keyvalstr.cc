//
// keyvalstr.cc
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

#include <iostream>

#include <util/misc/formio.h>
#include <util/keyval/keyval.h>
#include <util/class/proxy.h>

using namespace std;
using namespace sc;

////////////////////////////////////////////////////////////////////////
// StringKeyVal



StringKeyVal::StringKeyVal()
{

}

StringKeyVal::~StringKeyVal()
{
}

const char*
StringKeyVal::classname(const char * key)
{
  return 0;
}

const char*
StringKeyVal::truekeyword(const char * key)
{
  return key;
}

// This does not cause objects to be constructed.
int
StringKeyVal::key_exists(const char* key)
{
  stringrep(key);
  if (error() == OK || error() == HasNoValue) {
      return 1;
    }
  return 0;
}

Ref<KeyValValue>
StringKeyVal::key_value(const char* key, const KeyValValue &def)
{
  Ref<KeyValValue> result;

  if (!key) key = "TOP";

  // convert the key to the true key so variable assignments in the
  // input will effectively be done by reference
  // check to see if the datum is a described class
  const char* tkw = truekeyword(key);
  //ExEnv::outn() << "truekeyword = "<< tkw << '\n';
  if (!tkw) {
      result = 0;
    }
  else {
      // if a classname exists then read in the datum as an object
      const char* classn = classname(tkw);
      //ExEnv::outn() << "classname = " << classn << '\n';
      if (classn) {
          std::string truekey(tkw);
          if (_map.find(truekey) != _map.end()) {
              result = _map[truekey];
            }
          else {
              // create a new instance of this datum
              Ref<KeyVal> pkv = new PrefixKeyVal(this, tkw);
              const ClassDesc* cd = ClassDesc::name_to_class_desc(classn);
              if (!cd) {
                  ClassDesc::load_class(classn);
                  cd = ClassDesc::name_to_class_desc(classn);

                  if (cd == 0) {
                    std::ostringstream oss;
                    oss << "StringKeyVal is asked to construct an object of unknown type \"" << classn << "\"";
                    throw InputError(oss.str().c_str(),
                                     __FILE__, __LINE__, tkw, 0);
                  }
                }
              // the original error status must be saved
              KeyValError original_error = error();
              Ref<DescribedClass> newdc(cd->create(pkv));
              if (newdc == 0) {
                  ExEnv::errn() << "StringKeyVal::value: create failed for:" << endl
                       << " keyword = \"" << tkw << "\" class = \"" << classn
                       << "\"" << endl
                       << " either the KeyVal create operator doesn't" << endl
                       << " exist or memory was exhausted" << endl;
                }
              DescribedClassProxy *proxy
                  = dynamic_cast<DescribedClassProxy*>(newdc.pointer());
              if (proxy) {
                  newdc = proxy->object();
                }
              seterror(original_error);
              KeyValValueRefDescribedClass* keyvalvalue
                  = new KeyValValueRefDescribedClass(newdc);
              _map[truekey] = keyvalvalue;
              result = keyvalvalue;
            }
        }
      else {
          std::string str = stringrep(tkw);
          if (error() != OK) {
              result = 0;
            }
          else {
              result = new KeyValValuestring(str);
            }
        }
    }

  if (verbose_) {
      ExEnv::out0() << indent << key << " = ";
      if (result == 0) {
          ExEnv::out0() << def << " (default)";
        }
      else ExEnv::out0() << *result.pointer();
      ExEnv::out0() << endl;
    }

  return result;
}

void
StringKeyVal::errortrace(ostream&o)
{
  o << indent << "StringKeyVal: error: \"" << errormsg() << "\"" << endl;
}

void
StringKeyVal::dump(ostream&o)
{
  o << indent << "StringKeyVal: error: \"" << errormsg() << "\"" << endl;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
