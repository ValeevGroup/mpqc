//
// ipv2_alloc.cc
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
#include <util/keyval/ipv2.h>

using namespace sc;

ip_keyword_tree_t *
IPV2::ip_alloc_keyword_tree()
{
  ip_keyword_tree_t *result;

  result = (ip_keyword_tree_t *) malloc(sizeof(ip_keyword_tree_t));
  if (!result) {
    ExEnv::errn() << "ip_alloc_keyword_tree: malloc failed";
    error(0);
    }

  result->up = 0;
  result->down = 0;
  result->across = 0;
  result->keyword = 0;
  result->classname = 0;
  result->truename = 0;
  result->value = 0;
  result->variable = 0;
  result->seen = 0;

  return result;
  }

void
IPV2::ip_free_keyword_tree(ip_keyword_tree_t* tree)
{
  ip_keyword_tree_t *I,*start,*nextI;

  if (!tree) return;

  /* Convert the circular list into a standard linked list (to
   * avoid saber-c error messages) */
  start = tree->across;
  tree->across = 0;
  for (I=start; I!=0; I=nextI) {
    ip_free_keyword_tree(I->down);
    if (I->keyword) free(I->keyword);
    if (I->classname) free(I->classname);
    if (I->truename) free(I->truename);
    if (I->value) free(I->value);
    if (I->variable) free(I->variable);
    nextI = I->across;
    free(I);
    }
  }
