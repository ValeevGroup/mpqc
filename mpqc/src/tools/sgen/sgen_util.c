
/* Some convenient functions are provided here. */

/* $Log$
 * Revision 1.1  1993/12/29 12:53:58  etseidl
 * Initial revision
 *
 * Revision 1.2  1992/06/17  23:07:47  jannsen
 * modified to generate clean code
 *
 * Revision 1.1.1.1  1992/03/17  18:11:08  seidl
 * Struct GENerator 2.0
 *
 * Revision 1.1  1992/03/17  18:11:07  seidl
 * Initial revision
 *
 * Revision 1.2  1991/12/23  13:46:48  seidl
 * add functions member_name2 and index_dimension2 for use with binary
 * sgen functions
 *
 * Revision 1.1  1991/12/20  16:18:23  seidl
 * Initial revision
 *
 * Revision 1.4  91/09/28  16:40:41  cljanss
 * new naming convention is used to keep names <= 14 characters.
 * 
 * Revision 1.3  91/07/19  14:42:33  cljanss
 * Changed the way that generation modules are selected.
 * 
 * Revision 1.2  1991/06/16  15:11:15  janssen
 * added support for boolean
 *
 * Revision 1.1  1991/06/15  21:13:57  janssen
 * Initial revision
 * */
static char *rcsid = "$Id$";

#include <stdio.h>
#include <tmpl.h>
#include "types.h"
#include "global.h"

#include "sgen_util.gbl"
#include "sgen_util.lcl"

/* Convert an index and a structure name to a dimension name. */
GLOBAL_FUNCTION char *
index_dimension(structname,index)
char *structname;
index_t *index;
{
  static char dim[STRING_LENGTH];

  if (index->type == TYPE_CONSTANT) sprintf(dim,"%d",index->v.c);
  else sprintf(dim,"_%s->%s",structname,index->v.v);

  return dim;
  }

/* Given a member find the name for the member. */
GLOBAL_FUNCTION char *
member_name(structname,member)
char *structname;
member_t *member;
{
  static char name[STRING_LENGTH];

  if (member->uname) {
    sprintf(name,"_%s->%s.%s",structname,member->uname,member->name);
    }
  else {
    sprintf(name,"_%s->%s",structname,member->name);
    }

  return name;
  }

/* Given a member find the name for the type. */
GLOBAL_FUNCTION char *
type_name(member)
member_t *member;
{
  static char name[STRING_LENGTH];

  if (basic_type(member->type)) {
    if (member->qualifier == Q_UNSIGNED)
      sprintf(name,"unsigned %s",member->type);
    else if (member->qualifier == Q_SIGNED)
      sprintf(name,"signed %s",member->type);
    else
      sprintf(name,"%s",member->type);
    if (!strcmp(member->type,"string"))
      sprintf(name,"char *");
    else if (!strcmp(member->type,"boolean"))
      sprintf(name,"int");
    }
  else {
    sprintf(name,"%s_t",member->type);
    }

  return name;
  }

/* This routines 1 of type is one of the basic types. */
GLOBAL_FUNCTION int
basic_type(type)
char *type;
{
  if (!strcmp(type,"int")) return 1;
  if (!strcmp(type,"char")) return 1;
  if (!strcmp(type,"double")) return 1;
  if (!strcmp(type,"float")) return 1;
  if (!strcmp(type,"long")) return 1;
  if (!strcmp(type,"string")) return 1;
  if (!strcmp(type,"boolean")) return 1;

  return 0;
  }

/* This finds out how many indices are needed to loop thru the elements
 * of an array and writes out a declaration of these indices to the
 * global output file.  The structname is needed incase an error occurs. */
GLOBAL_FUNCTION VOID
declare_indices(members,structname)
member_list_t *members;
char *structname;
{
  int max = 0;
  int i,n;
  char *maxname;
  member_list_t *I;
  index_list_t *J;

  /* Go thru all of the members and find the maximum number of indices on
   * a given member. */
  for (I=members; I!=NULL; I=I->p) {
    for (n=0,J=I->member->indices; J!=NULL; n++,J=J->p);
    if (n>max) { max = n; maxname = I->member->name; }
    }

  if (max > 'z' - 'i') {
    fprintf(stderr,"Problem generating a loop thru the elements of the\n");
    fprintf(stderr,"member %s of the structure %s:\n", maxname, structname);
    fprintf(stderr,"too many indices\n");
    error("declare_indices: aborting");
    }

  if (!max) return;

  fprintf(output,"  int ");

  for (i=0; i<max; i++) {
    if (i>0) fprintf(output,",");
    fprintf(output,"%c",'i'+i);
    }

  fprintf(output,";\n");
  }

/* This is like declare_indices, but one less index is declared. */
GLOBAL_FUNCTION VOID
declare_indices_less1(members,structname)
member_list_t *members;
char *structname;
{
  int max = 0;
  int i,n;
  char *maxname;
  member_list_t *I;
  index_list_t *J;

  /* Go thru all of the members and find the maximum number of indices on
   * a given member. */
  for (I=members; I!=NULL; I=I->p) {
    for (n=0,J=I->member->indices; J!=NULL; n++,J=J->p);
    if (n>max) { max = n; maxname = I->member->name; }
    }

  if (max) max--;

  if (max > 'z' - 'i') {
    fprintf(stderr,"Problem generating a loop thru the elements of the\n");
    fprintf(stderr,"member %s of the structure %s:\n", maxname, structname);
    fprintf(stderr,"too many indices\n");
    error("declare_indices: aborting");
    }

  if (!max) return;

  fprintf(output,"  int ");

  for (i=0; i<max; i++) {
    if (i>0) fprintf(output,",");
    fprintf(output,"%c",'i'+i);
    }

  fprintf(output,";\n");
  }

/* This is like declare_indices_less1, but is only for basic types. */
GLOBAL_FUNCTION VOID
declare_indices_for_nonpointer(members,structname)
member_list_t *members;
char *structname;
{
  int max = 0;
  int i,n;
  char *maxname;
  member_list_t *I;
  index_list_t *J;

  /* Go thru all of the members and find the maximum number of indices on
   * a given member. */
  for (I=members; I!=NULL; I=I->p) {
    for (n=0,J=I->member->indices; J!=NULL; n++,J=J->p);
    if (I->member->pointer) continue;
    if (n>max) { max = n; maxname = I->member->name; }
    }

  if (max > 'z' - 'i') {
    fprintf(stderr,"Problem generating a loop thru the elements of the\n");
    fprintf(stderr,"member %s of the structure %s:\n", maxname, structname);
    fprintf(stderr,"too many indices\n");
    error("declare_indices: aborting");
    }

  if (!max) return;

  fprintf(output,"  int ");

  for (i=0; i<max; i++) {
    if (i>0) fprintf(output,",");
    fprintf(output,"%c",'i'+i);
    }

  fprintf(output,";\n");
  }

/* This is like declare_indices_less1, but is only for basic types. */
GLOBAL_FUNCTION VOID
declare_indices_less1basic(members,structname)
member_list_t *members;
char *structname;
{
  int max = 0;
  int i,n;
  char *maxname;
  member_list_t *I;
  index_list_t *J;

  /* Go thru all of the members and find the maximum number of indices on
   * a given member. */
  for (I=members; I!=NULL; I=I->p) {
    for (n=0,J=I->member->indices; J!=NULL; n++,J=J->p);
    if (basic_type(I->member->type) && n>0 && I->member->pointer==0) n--;
    if (n>max) { max = n; maxname = I->member->name; }
    }

  if (max > 'z' - 'i') {
    fprintf(stderr,"Problem generating a loop thru the elements of the\n");
    fprintf(stderr,"member %s of the structure %s:\n", maxname, structname);
    fprintf(stderr,"too many indices\n");
    error("declare_indices: aborting");
    }

  if (!max) return;

  fprintf(output,"  int ");

  for (i=0; i<max; i++) {
    if (i>0) fprintf(output,",");
    fprintf(output,"%c",'i'+i);
    }

  fprintf(output,";\n");
  }

/* This sees if any of the declarations in a declaration list need
 * a given module. */
GLOBAL_FUNCTION int
is_entirely_excluded(module)
char *module;
{
  declaration_t *I;

  for (I=dl; I!=NULL; I=I->p) {
    if (!is_excluded(module,I)) return 0;
    }
  return 1;
  }


/* This is given the name of a routine generation module and a declaration
 * struct.  If the name is to be excluded 1 is returned, otherwise 0
 * is returned. */
GLOBAL_FUNCTION int
is_excluded(module,dec)
char *module;
declaration_t *dec;
{
  name_list_t *I;
  char notmodule[STRING_LENGTH];

  strcpy(notmodule,"!");
  strcat(notmodule,module);

  /* Check through the option list for this declaration. */
  for (I=dec->modules; I!=NULL; I=I->p) {
    if (!strcmp(I->name,module)) return 0;
    if (!strcmp(I->name,notmodule)) return 1;
    }

  /* Check the list of default_modules. */
  for (I=default_modules; I!=NULL; I=I->p) {
    if (!strcmp(I->name,module)) return 0;
    if (!strcmp(I->name,notmodule)) return 1;
    }
  /* The default default is to exclude the module. */
  return 1;
  }

/* This returns true if the member is an array index for another
 * member of the structure. */
GLOBAL_FUNCTION int
is_array_index(member,members)
member_t *member;
member_list_t *members;
{
  member_list_t *I;
  index_list_t *J;

  for (I=members; I!=NULL; I=I->p) {
    for (J=I->member->indices; J!=NULL; J=J->p) {
      if (  (J->index.type == TYPE_VARIABLE)
          &&(!strcmp(member->name,J->index.v.v))) return 1;
      }
    }
  return 0;
  }

/* This returns true if the member is an array index for another
 * member of the structure. */
GLOBAL_FUNCTION int
is_union_selector(member,members)
member_t *member;
member_list_t *members;
{
  member_list_t *I;
  index_list_t *J;

  for (I=members; I!=NULL; I=I->p) {
    if (I->member->uselname && (!strcmp(member->name,I->member->uselname)))
      return 1;
    }
  return 0;
  }

/* the following functions added by ets for binary routines */
/* Convert an index and a structure name to a dimension name with _2
 * appended. For use with binary routines iseq and asgn */
GLOBAL_FUNCTION char *
index_dimension2(structname,ind,index)
char *structname;
int ind;
index_t *index;
{
  static char dim[STRING_LENGTH];

  if (index->type == TYPE_CONSTANT) sprintf(dim,"%d",index->v.c);
  else sprintf(dim,"_%s_%d->%s",structname,ind,index->v.v);

  return dim;
  }

/* Given a member find the name for the member. */
GLOBAL_FUNCTION char *
member_name2(structname,ind,member)
char *structname;
int ind;
member_t *member;
{
  static char name[STRING_LENGTH];

  if (member->uname) {
    sprintf(name,"_%s_%d->%s.%s",structname,ind,member->uname,member->name);
    }
  else {
    sprintf(name,"_%s_%d->%s",structname,ind,member->name);
    }

  return name;
  }
