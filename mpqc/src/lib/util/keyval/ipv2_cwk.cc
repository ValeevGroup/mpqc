//
// ipv2_cwk.cc
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

/* These routines manipulate the current working keyword.  This
 * is an ordered list of keyword_tree's.  When a relative
 * keyword is searched for we start looking under the first keyword
 * tree in the current working keyword list and if it is not found
 * continue looking under successive members of the list. */

#include <stdlib.h>
#include <string.h>
#include <util/keyval/ipv2.h>

using namespace std;
using namespace sc;

/* This sets up the current working keyword path to the declaration
 * list. */
void
IPV2::cwk_root()
{
  free_keyword_tree_list(ip_cwk);
  ip_cwk = splice_keyword_tree_list(ip_tree,NULL);
  }

/* This sets up the current working keyword path to NULL
 * list. */
void
IPV2::cwk_clear()
{
  free_keyword_tree_list(ip_cwk);
  ip_cwk = NULL;
  }

/* This adds a keyword tree to the keyword path. */
void
IPV2::ip_cwk_add_kt(ip_keyword_tree_t *kt)
{
  ip_cwk = splice_keyword_tree_list(kt,ip_cwk);
  }

/* This adds a keyword to the keyword path. */
/* NOTE: the last path to be searched must be added first. */
void
IPV2::cwk_add(const char* keyword)
{
  ip_keyword_tree_t *kt;
  ip_keyword_tree_list_t *I,*old_cwk;

  old_cwk = ip_cwk;

  /* Initialize the new cwk list. */
  ip_cwk = NULL;

  /* See if the keyword we were given is NULL.
   * If so, just copy the cwk. */
  if (!keyword) {
    ip_cwk = old_cwk;
    }
  /* See if we have been given an absolute path. */
  else if (keyword[0] == ':') {
    /* Copy the old keyword tree list to the new ip_cwk global. */
    for (I=old_cwk; I!=NULL; I=I->p) {
      ip_cwk = splice_keyword_tree_list(I->kt,ip_cwk);
      }
    /* Add the keyword into the keyword list. */
    kt = ip_descend_tree(ip_tree,&(keyword[1]));
    if (kt) ip_cwk = splice_keyword_tree_list(kt->down,ip_cwk);
    free_keyword_tree_list(old_cwk);
    }
  else {
    /* For an relative path append the keyword to each of the keyword
     * paths in the current working keyword list. */
    ip_keyword_tree_t *kt;
    for (I=old_cwk; I!=NULL; I=I->p) {
      kt = ip_descend_tree(I->kt,keyword);
      if (kt) {
        kt = ip_descend_tree(I->kt,keyword);
        ip_cwk = splice_keyword_tree_list(kt->down,ip_cwk);
        }
      }
    free_keyword_tree_list(old_cwk);
    }

  if (ip_keyword) {
    *ip_out << "IP_KEYWORDS from IPV2::cwk_add (" << keyword << "): {"
            << endl;
    for (I=ip_cwk; I!=NULL; I=I->p) {
        *ip_out << "  ";
        print_keyword(*ip_out,I->kt);
        *ip_out << endl;
      }
    *ip_out << "  }" << endl;
    }

  }

/* This pushes the old cwk list without modifying the current cwk list. */
void
IPV2::cwk_push()
{
  ip_keyword_tree_list_t *I;

  /* Allocate a stack slot to hold the old cwk. */
  if (!cwkstack) {
    cwkstack = (ip_cwk_stack_t *) malloc(sizeof(ip_cwk_stack_t));
    cwkstack->p = NULL;
    }
  else {
    ip_cwk_stack_t *tmp = cwkstack;
    cwkstack = (ip_cwk_stack_t *) malloc(sizeof(ip_cwk_stack_t));
    cwkstack->p = tmp;
    }

  /* Push the previous cwk list onto the stack. */
  cwkstack->ktl = ip_cwk;

  /* Copy the old keyword tree list to the ip_cwk global. */
  ip_cwk = NULL;
  for (I=cwkstack->ktl; I!=NULL; I=I->p) {
    ip_cwk = splice_keyword_tree_list(I->kt,ip_cwk);
    }
  }

/* This moves up the keyword tree for each member of the cwk list.
 * If a cwk is already at the top of the tree, then that cwk list entry
 * will be deleted. */
void
IPV2::cwk_pop()
{
  ip_cwk_stack_t *tmp;
  if (!cwkstack) {
    error("IPV2::cwk_pop: tried to pop above the top");
    }
  free_keyword_tree_list(ip_cwk);
  ip_cwk = cwkstack->ktl;
  tmp = cwkstack;
  cwkstack = tmp->p;
  free(tmp);
  }

/* Descend the keyword tree using the cwk and obtain a new keyword tree. */
ip_keyword_tree_t *
IPV2::ip_cwk_descend_tree(const char* keyword)
{
  ip_keyword_tree_list_t *I;
  ip_keyword_tree_t *kt=NULL;

  /* If the keyword is NULL, then the first value in the cwk list is returned.*/
  if (keyword[0] == '\0') {
    if (ip_cwk) kt = ip_cwk->kt;
    else kt = NULL;
    }
  /* Is the keyword an absolute path? */
  else if (keyword[0] != ':') {
    /* See if we can descend to this keyword in any of the cwk's */
    for (I=ip_cwk; I!=NULL; I=I->p) {
      if ((kt = ip_descend_tree(I->kt,keyword)) != NULL) break;
      }
    }
  else {
    kt = ip_descend_tree(ip_tree,&(keyword[1]));
    }

  return kt;
  }

////////////////////////////////////////////////////////////////////////
// IPV2StrTok provides strtok functionality, but allows multiple strings
// to be processed at a time.

class IPV2StrTok {
  private:
    char* str;
    const char* delim;
    int ndelim;
  public:
    IPV2StrTok(char* s, const char*d): str(s), delim(d), ndelim(strlen(d)) {}
    char* tok();
    int is_white(char c);
};

int
IPV2StrTok::is_white(char c)
{
  for (int i=0; i<ndelim; i++) {
      if (c == delim[i]) {
          return 1;
        }
    }
  return 0;
}

char*
IPV2StrTok::tok()
{
  // move str past the white space
  while (*str && is_white(*str)) str++;
  char *ret = str;

  // put 0 at the end of the string and advance str
  while (*str && !is_white(*str)) str++;
  if (*str) {
      *str = '\0';
      str++;
    }

  if (*ret) return ret;
  else return NULL;
}

////////////////////////////////////////////////////////////////////////

/* Descend the given keyword tree using the info in the passed string.
 * The new keyword tree or NULL, if it is not found, will be returned. */
ip_keyword_tree_t *
IPV2::ip_descend_tree(ip_keyword_tree_t* kt,const char* keyword)
{
  ip_keyword_tree_t *I,*r;
  char ch[KEYWORD_LENGTH];
  char *token;
  int found;

  if (!keyword) return kt;

  if (strlen(keyword)+1 > KEYWORD_LENGTH) {
      error("ip_descend_tree: maximum KEYWORD_LENGTH has been exceeded");
      }

  if (keyword[0] == ':') {
      kt = ip_tree;
      strcpy(ch,&keyword[1]);
    }
  else {
      strcpy(ch,keyword);
    }

  r = kt;
  //IPV2StrTok tok(ch, ": \t");
  IPV2StrTok tok(ch, ":");
  token = tok.tok();
  while ((r != NULL) && (token != NULL)) {
    /* Transverse the circular list. */
    found = 0;
    I = r;
    do {
      if (!strcmp(token,"..")) {
        r = I->up;
        token = tok.tok();
        if (token == NULL) return I;
        found = 1;
        break;
        }
      else if (! I->keyword) {
          return NULL;
        }
      else if (!strcmp(token,I->keyword)) {
        I->seen = 1;
        if (I->variable) I = ip_descend_tree(I,I->variable);
        token = tok.tok();
        if (token == NULL) return I;
        r = I->down;
        if (!r) {
            return NULL;
          }
        found = 1;
        break;
        }
      } while ((I = I->across) != r);
    if (!found) {
        return NULL;
      }
    }

  if (r && ip_keyword) {
    *ip_out << "IP_KEYWORD from ip_descend_tree: ";
    print_keyword(*ip_out,r);
    *ip_out << endl;
    }

  return r;
  }

/* Return the value of the given keyword. */
char *
IPV2::ip_key_value(const char* keyword)
{
  ip_keyword_tree_t *kt;

  kt = ip_cwk_descend_tree(keyword);

  if (kt && ip_keyword) {
    *ip_out << "IP_KEYWORD from ip_key_value: ";
    print_keyword(*ip_out,kt);
    *ip_out << endl;
    }

  if (kt) return kt->value;
  else return NULL;
  }

/* Free memory for a keyword tree list. */
void
IPV2::free_keyword_tree_list(ip_keyword_tree_list_t *ktl)
{
  if (!ktl) return;
  free_keyword_tree_list(ktl->p);
  free(ktl);
  }

/* Splice a new keyword tree into a keyword tree list. */
ip_keyword_tree_list_t*
IPV2::splice_keyword_tree_list(ip_keyword_tree_t*kt,ip_keyword_tree_list_t*p)
{
  ip_keyword_tree_list_t *r;

  if (kt==NULL) return p;
  r = (ip_keyword_tree_list_t *) malloc(sizeof(ip_keyword_tree_list_t));
  r->kt = kt;
  r->p = p;
  return r;
  }
