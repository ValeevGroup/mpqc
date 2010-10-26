//
// ipv2_print.cc
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

#include <iostream>
#include <util/misc/formio.h>
#include <util/keyval/ipv2.h>

using namespace std;
using namespace sc;

void
IPV2::print_keyword(ostream&fp,ip_keyword_tree_t*st)
{
  if (st) {
    if (st->up) print_keyword(fp,st->up);
    fp << st->keyword << ":";
    }
  }

/* This prints out a keyword tree, tree.  If tree is NULL then ip_tree
 * is printed out. */
void
IPV2::print_tree(ostream&fp,ip_keyword_tree_t*tree)
{
  if (!tree) tree = ip_tree;
  if (!tree) return;

  print_tree_(fp,tree);
  }


/* This prints out a keyword tree, tree.  If tree is NULL then ip_tree
 * is printed out. */
void
IPV2::print_tree_(ostream&fp,ip_keyword_tree_t*tree)
{
  ip_keyword_tree_t *I;

  I=tree;
  do {
    //if (I->value && I->down) {
    //  warn("print_tree: tree has both value and subtrees - can't print");
    //  warn("keyword is %s, value is %s, subtree key is %s\n",
    //       I->keyword,I->value,I->down->keyword);
    //  }

    if (!I->keyword) {
      warn("print_tree: tree has no keyword - impossible");
      }

    fp << indent;
    if (ip_special_characters(I->keyword)) {
      fp << "\"" << I->keyword << "\"" << endl;
      }
    else {
      fp << I->keyword << endl;
      }

    if (I->classname) {
      fp << "<" << I->keyword << ">" << endl;
      }

    if (!(I->value || I->down || I->variable)) {
      fp << ": (" << endl;
      }

    if (I->variable) {
      fp << " = $" << I->variable << endl;
      }
    if (I->truename) {
      fp << "\"" << I->truename << "\"";
      }

    if (I->value) {
      if (I->down) fp << " (= " << I->value << ")";
      else fp << " = " << I->value << endl;
      }
    if (I->down) {
      fp << ": (" << endl;
      fp << incindent;
      print_tree_(fp,I->down);
      fp << decindent;
      fp << indent << ")" << endl;
      }

    } while ((I = I->across) != tree);

  }

/* This prints out a keyword tree, tree.  If tree is NULL then ip_tree
 * is printed out. */
void
IPV2::print_unseen(ostream&fp,ip_keyword_tree_t*I)
{
  if (!I) I = ip_tree;
  if (!I) return;
  ip_keyword_tree_t *start = I;
  do {
      if (!I->seen) {
          fp << indent;
          print_keyword(fp,I->up);
          fp << I->keyword << endl;
        }
      else if (I->down) {
          print_unseen(fp,I->down);
        }
    } while ((I = I->across) != start);
}

/* This prints out a keyword tree, tree.  If tree is NULL then ip_tree
 * is printed out. */
int
IPV2::have_unseen(ip_keyword_tree_t*I)
{
  if (!I) I = ip_tree;
  if (!I) return 0;
  ip_keyword_tree_t *start = I;
  do {
      if (!I->seen) {
          return 1;
        }
      else if (I->down) {
          if (have_unseen(I->down)) return 1;
        }
    } while ((I = I->across) != start);
  return 0;
}

int
IPV2::ip_special_characters(char*keyword)
{
  char *ch=keyword;

  if (!keyword) return 0;
  while (*ch) {
    if (!(  (*ch >= 'a' && *ch <= 'z')
          ||(*ch >= 'A' && *ch <= 'Z')
          ||(*ch >= '0' && *ch <= '9')
          ||(*ch == '<')
          ||(*ch == '>')
          ||(*ch == '+')
          ||(*ch == '-')
          ||(*ch == '.')
          ||(*ch == '_'))) return 1;

    ch++;
    }

  return 0;
  }
