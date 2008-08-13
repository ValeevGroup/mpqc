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

#ifdef HAVE_CONFIG_H
#include <scconfig.h>
#endif
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
      // exist. Test for this.
      std::map<std::string,Ref<KeyValValue> >::iterator lb = _map.lower_bound(k);
      if (lb != _map.end()) {
          std::string lbk(lb->first);
          if (lbk.size() > k.size()
              && lbk[k.size()] == ':'
              && lbk.substr(0,k.size()) == k) {
              seterror(HasNoValue);
              return 1;
          }
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
AssignedKeyVal::assign(const char*key,const char* val)
{
  assign(key,new KeyValValuestring(val));
}
void
AssignedKeyVal::assign(const char*key,const Ref<DescribedClass>&val)
{
  assign(key,new KeyValValueRefDescribedClass(val));
}

void
AssignedKeyVal::clear()
{
  _map.clear();
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
