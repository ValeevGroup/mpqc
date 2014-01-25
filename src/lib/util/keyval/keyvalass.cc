//
// keyvalass.cc
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

#include <mpqc_config.h>
#include <util/keyval/keyval.h>

using namespace sc;

AssignedKeyVal::AssignedKeyVal()
{
}

AssignedKeyVal::~AssignedKeyVal()
{
}

int
AssignedKeyVal::key_exists(const char * key)
{
  std::string k(key); 
  int result = (_map.find(k) != _map.end());
  if (!result) {
      // An exact match does not exist in the map. However, it is
      // possible that key is part of a longer segment that does
      // exist. Test all keys that are not less than k, contain k as a substring, followed by at least one ':'
      std::map<std::string,Ref<KeyValValue> >::iterator trial = _map.lower_bound(k);
      while (trial != _map.end() && trial->first.find(k,0) != std::string::npos) {
          std::string tk(trial->first);
          if (tk.size() > k.size()
              && tk.find(':',k.size()) != std::string::npos
              && tk.substr(0,k.size()) == k) {
              seterror(HasNoValue);
              return 1;
          }
          ++trial; // try next key ...
      }
      seterror(UnknownKeyword);
    }
  else {
      seterror(OK);
    }
  return result;
}

Ref<KeyValValue>
AssignedKeyVal::key_value(const char * key, const KeyValValue &def)
{
  std::string k(key); 
  if (exists(key)) {
      seterror(OK);
      return _map[k];
    }
  else {
      seterror(UnknownKeyword);
      return 0;
    }
}

void
AssignedKeyVal::assign(const char*key,const Ref<KeyValValue>& val)
{
  std::string k(key);
  _map[k] = val;
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
AssignedKeyVal::assign(const char*key,long val)
{
  assign(key,new KeyValValuelong(val));
}
void
AssignedKeyVal::assign(const char*key,const char* val)
{
  assign(key,new KeyValValuestring(val));
}
void
AssignedKeyVal::assign(const char*key,const std::string& val)
{
  assign(key,new KeyValValuestring(val));
}
void
AssignedKeyVal::assign(const char*key,const Ref<DescribedClass>&val)
{
  assign(key,new KeyValValueRefDescribedClass(val));
}

const char*
AssignedKeyVal::classname(const char * key)
{
  Ref<KeyValValueRefDescribedClass> kv_dc;
  kv_dc << this->key_value(key, KeyValValueRefDescribedClass());
  const char* result = 0;
  if (kv_dc.nonnull()) {
    Ref<DescribedClass> dc;
    kv_dc->describedclassvalue(dc);
    result = dc->class_name();
  }
  return result;
}

void
AssignedKeyVal::clear()
{
  _map.clear();
}

void
AssignedKeyVal::print(std::ostream & os) const
{
  // sort keys
  typedef std::set<std::string> keys_t;
  typedef keys_t::const_iterator keys_iter;
  keys_t keys;
  typedef _map_t::const_iterator mapiter;
  for(mapiter i = _map.begin(); i != _map.end(); ++i)
    keys.insert(i->first);

  // print
  os << indent << "AssignedKeyVal = (" << incindent << std::endl;
  for(keys_iter ki = keys.begin(); ki != keys.end(); ++ki) {
    os << indent << *ki << " = ";
    _map.find(*ki)->second->print(os);
    os << std::endl;
  }
  os << decindent << indent << ")" << std::endl;
}

#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
template class EAVLMMapNode<std::string, AVLMapNode<std::string,Ref<KeyValValue> > >;
template class EAVLMMap<std::string, AVLMapNode<std::string,Ref<KeyValValue> > >;
template class AVLMapNode<std::string,Ref<KeyValValue> >;
template class AVLMap<std::string,Ref<KeyValValue> >;
#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
