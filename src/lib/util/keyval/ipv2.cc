//
// ipv2.cc
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

#include <stdlib.h>
#include <string.h>

#include <util/keyval/ipv2.h>
#include <util/keyval/ipv2_parse.h>

using namespace std;
using namespace sc;

IPV2::IPV2():
filename_(0),
table_keywords(0),
current_table_keyword(0),
table_sub_tree(0),
table_array_depth(0),
karray_indices(0),
sub_tree(0),
cwkstack(0),
ip_initialized(0),
ip_in(0),
ip_out(0),
ip_tree(0),
ip_cwk(0),
ip_keyword(0)
{
  lastkeyword[0] = '\0';
  lexer = new IPV2FlexLexer;
}

IPV2::~IPV2()
{
  if (ip_tree) ip_free_keyword_tree(ip_tree);
  ip_tree = 0;
  delete lexer;
  delete[] filename_;
  free_keyword_tree_list(ip_cwk);
}

void IPV2::read(istream&in,ostream&out,const char *filename)
{
  delete[] filename_;
  if (filename) {
      filename_ = strcpy(new char[strlen(filename)+1], filename);
    }
  else filename_ = 0;
  if (ip_initialized) {
    ip_append(in,out);
    }
  else {
    ip_initialize(in,out);
    cwk_root();
    }
}

void IPV2::yerror(const char* s)
{
  error(s);
}

/* Show position. */
void
IPV2::showpos()
{
  ExEnv::errn() << "error occurred at line number "
       << lexer->lineno() << " (roughly)" << endl;
  if (filename_) {
      ExEnv::errn() << "in file \"" << filename_ << "\"";
    }
  ExEnv::errn() << endl;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
