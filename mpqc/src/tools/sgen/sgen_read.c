
/* This contains the routines that YACC needs when it reads in the input
 * file. */

/* $Log$
 * Revision 1.3  1996/03/23 02:39:06  cljanss
 * Everything can now be configured with autoconf.
 *
 * Revision 1.2  1994/10/18 23:04:03  etseidl
 * fix many warnings, use memset rather than bzero
 *
 * Revision 1.1.1.1  1993/12/29  12:53:58  etseidl
 * SC source tree 0.1
 *
 * Revision 1.2  1992/03/30  23:10:03  seidl
 * merge in sandia changes
 *
 * Revision 1.4  91/09/28  16:40:40  cljanss
 * new naming convention is used to keep names <= 14 characters.
 * 
 * Revision 1.3  91/07/19  18:21:37  cljanss
 * The default_modules global is now properly initialized.
 * 
 * Revision 1.2  1991/07/19  14:42:33  cljanss
 * Changed the way that generation modules are selected.
 *
 * Revision 1.1  1991/06/15  21:13:57  janssen
 * Initial revision
 * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "types.h"
#include "global.h"

#include "sgen_read.gbl"
#include "sgen_read.lcl"

#include "error.gbl"

GLOBAL_FUNCTION void
init_declarations()
{
  free_name_list(default_modules);
  default_modules = NULL;

  free_declaration(dl);
  dl = NULL;
  }

LOCAL_FUNCTION void
free_declaration(dl)
declaration_t *dl;
{
  if (!dl) return;
  free_declaration(dl->p);
  free_member_list(dl->members);
  free(dl->name);
  free(dl);
  }

LOCAL_FUNCTION void
free_member_list(members)
member_list_t *members;
{
  if (!members) return;
  free_member_list(members->p);
  free_member(members->member);
  free(members);
  }

LOCAL_FUNCTION void
free_name_list(names)
name_list_t *names;
{
  if (!names) return;
  free_name_list(names->p);
  free(names->name);
  free(names);
  }

LOCAL_FUNCTION void
free_member(member)
member_t *member;
{
  if (!member) return;
  free(member->type);
  free(member->name);
  free_index_list(member->indices);
  }

LOCAL_FUNCTION void
free_index_list(indices)
index_list_t *indices;
{
  if (!indices) return;
  free_index_list(indices->p);
  if (indices->index.type == TYPE_VARIABLE) free(indices->index.v.v);
  free(indices);
  }

LOCAL_FUNCTION declaration_t *
alloc_declaration()
{
  declaration_t *r = (declaration_t *) malloc(sizeof(declaration_t));
  if (!r) error("malloc failed");
  return r;
  }

LOCAL_FUNCTION member_t *
alloc_member()
{
  member_t *r = (member_t *) malloc(sizeof(member_t));
  if (!r) error("malloc failed");
  return r;
  }

LOCAL_FUNCTION name_list_t *
alloc_name_list()
{
  name_list_t *r = (name_list_t *) malloc(sizeof(name_list_t));
  if (!r) error("malloc failed");
  return r;
  }

LOCAL_FUNCTION member_list_t *
alloc_member_list()
{
  member_list_t *r = (member_list_t *) malloc(sizeof(member_list_t));
  if (!r) error("malloc failed");
  return r;
  }

LOCAL_FUNCTION index_list_t *
alloc_index_list()
{
  index_list_t *r = (index_list_t *) malloc(sizeof(index_list_t));
  if (!r) error("malloc failed");
  return r;
  }

/* Add a declaration to the end of the global declaration list. */
GLOBAL_FUNCTION void
declaration(name,modules,members)
char *name;
name_list_t *modules;
member_list_t *members;
{
  declaration_t *I;

  if (dl) {
    for (I=dl; I->p!=NULL; I=I->p);
    I->p = alloc_declaration();
    I = I->p;
    }
  else {
    dl = alloc_declaration();
    I = dl;
    }

  I->p = NULL;
  I->name = name;
  I->members = members;
  I->modules = modules;
  }

/* This sets up the global name list that contains the default
 * modules. */
GLOBAL_FUNCTION void
set_default_modules(nl)
name_list_t *nl;
{
  if (default_modules) free_name_list(default_modules);
  default_modules = nl;
  }

/* This adds a module to the name list (this is the same
 * as add_name. */
GLOBAL_FUNCTION name_list_t *
add_module(names,name)
name_list_t *names;
char *name;
{
  return add_name(names,name);
  }

/* This adds a module to be excluded the name list (this is the same
 * as add_name. */
GLOBAL_FUNCTION name_list_t *
exclude_module(names,name)
name_list_t *names;
char *name;
{
  char *xname;

  xname = (char *) malloc(strlen(name)+2);
  strcpy(xname,"!");
  strcat(xname,name);
  free(name);
  return add_name(names,xname);
  }

/* This puts a name into the end of the name list. */
GLOBAL_FUNCTION name_list_t *
add_name(names,name)
name_list_t *names;
char *name;
{
  name_list_t *I;

  if (names) {
    for (I=names; I->p!=NULL; I=I->p);
    I->p = alloc_name_list();
    I = I->p;
    }
  else {
    names = alloc_name_list();
    I = names;
    }

  I->name = name;
  I->p = NULL;

  return names;
  }

/* This puts members into the end of the member list. */
GLOBAL_FUNCTION member_list_t *
add_member(members,member)
member_list_t *members;
member_t *member;
{
  member_list_t *I;

  if (members) {
    for (I=members; I->p!=NULL; I=I->p);
    I->p = alloc_member_list();
    I = I->p;
    }
  else {
    members = alloc_member_list();
    I = members;
    }

  I->member = member;
  I->p = NULL;

  return members;
  }

/* This puts members into the end of the member list. */
GLOBAL_FUNCTION member_list_t *
merge_members(m1,m2)
member_list_t *m1;
member_list_t *m2;
{
  member_list_t *I;

  if (!m1) return m2;
  if (!m2) return m1;

  for (I=m1; I->p!=NULL; I=I->p);
  I->p = m2;

  return m1;
  }

GLOBAL_FUNCTION member_list_t *
setup_union(uselname,uname,members)
char *uselname;
char *uname;
member_list_t *members;
{
  member_list_t *I;

  for (I=members; I!=NULL; I=I->p) {
    I->member->uselname = uselname;
    I->member->uname = uname;
    }

  return members;
  }

GLOBAL_FUNCTION member_t *
create_member(qualifier,type,pointer,name,indices,uselname,uselval,uname)
int qualifier;
char *type;
int pointer;
char *name;
index_list_t *indices;
char *uselname;
char *uselval;
char *uname;
{
  member_t *r = alloc_member();
  r->qualifier = qualifier;
  r->type = type;
  r->pointer = pointer;
  r->name = name;
  r->indices = indices;
  r->uselname = uselname;
  r->uselval = uselval;
  r->uname = uname;
  return r;
  }

GLOBAL_FUNCTION index_list_t *
add_index(indices,index)
index_list_t *indices;
char *index;
{
  index_list_t *I;

  if (indices) {
    for (I=indices; I->p!=NULL; I=I->p);
    I->p = alloc_index_list();
    I = I->p;
    }
  else {
    indices = alloc_index_list();
    I = indices;
    }

  /* Determine whether or not this index is a constant or a variable. */
  if (isdigit(index[0])) {
    I->index.type = TYPE_CONSTANT;
    I->index.v.c = atoi(index);
    if (!I->index.v.c) {
      fprintf(stderr,"%s: bad index: %s\n",progname,index);
      error(NULL);
      }
    }
  else {
    I->index.type = TYPE_VARIABLE;
    I->index.v.v = index;
    }

  I->p = NULL;
  return indices;
  }

