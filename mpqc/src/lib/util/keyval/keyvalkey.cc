//
// keyvalkey.cc
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

#include <util/keyval/keyval.h>

KeyValKeyword::KeyValKeyword() :
  keyword_(0)
{
}

KeyValKeyword::KeyValKeyword(const char* name)
{
  // get rid of the leading ':' character
  if (name[0] == ':') {
      keyword_ = ::strcpy(new char[strlen(name)],&name[1]);
    }
  else {
      keyword_ = ::strcpy(new char[strlen(name)+1],name);
    }
}

KeyValKeyword::KeyValKeyword(const KeyValKeyword& key):
  keyword_(::strcpy(new char[strlen(key.keyword_)+1],key.keyword_))
{
}

KeyValKeyword::~KeyValKeyword()
{
  if (keyword_) delete[] keyword_;
}

KeyValKeyword& KeyValKeyword::operator=(const KeyValKeyword& key)
{
  if (keyword_ && keyword_ != key.keyword_) {
      delete[] keyword_;
      keyword_ = ::strcpy(new char[strlen(key.keyword_)+1],key.keyword_);
    }
  return *this;
}

int KeyValKeyword::operator==(const KeyValKeyword& ck) const
{
  return cmp(ck) == 0;
}

int KeyValKeyword::operator<(const KeyValKeyword& ck) const
{
  return cmp(ck)<0;
}

int
KeyValKeyword::hash() const
{
  int r=0;
  size_t i;

  // Even numbered bytes make up the lower part of the hash index
  for (i=0; i < ::strlen(keyword_); i+=2) {
      r ^= keyword_[i];
    }

  // Odd numbered bytes make up the upper part of the hash index
  for (i=1; i < ::strlen(keyword_); i+=2) {
      r ^= keyword_[i]<<8;
    }

  return r;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
