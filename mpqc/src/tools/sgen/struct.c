
/* $Log$
 * Revision 1.3  1996/03/23 02:39:09  cljanss
 * Everything can now be configured with autoconf.
 *
 * Revision 1.2  1994/10/18 23:04:06  etseidl
 * fix many warnings, use memset rather than bzero
 *
 * Revision 1.1.1.1  1993/12/29  12:53:59  etseidl
 * SC source tree 0.1
 *
 * Revision 1.3  1993/04/28  17:30:50  jannsen
 * added struct_close() and moved close call from struct.c to sgen.c
 *
 * Revision 1.2  1992/03/30  23:10:15  seidl
 * merge in sandia changes
 *
 * Revision 1.3  91/09/28  16:40:42  cljanss
 * new naming convention is used to keep names <= 14 characters.
 * 
 * Revision 1.2  91/06/16  15:11:15  janssen
 * added support for boolean
 * 
 * Revision 1.1  1991/06/15  21:13:57  janssen
 * Initial revision
 * */

#include <stdio.h>
#include <string.h>
#include "types.h"
#include "global.h"

#include "struct.gbl"
#include "struct.lcl"

#include "sgen.gbl"
#include "error.gbl"
#include "sgen_util.gbl"

GLOBAL_FUNCTION void
struct_open()
{
  char outfile[FILENAME_MAX];

  /* Convert filename to the name of the output file. */
  strcpy(outfile,"");
  strcat(outfile,BaseName);
  strcat(outfile,".h");

  /* Open the output file. */
  output = fopen(outfile,"w");
  if (!output) {
    fprintf(stderr,"%s: couldn't open the output file: %s\n",progname,outfile);
    error(NULL);
    }
  includeoutput = output;
  }

GLOBAL_FUNCTION void
struct_close()
{
  fclose(includeoutput);
  includeoutput = 0;
}

GLOBAL_FUNCTION void
struct_gen()
{
  declaration_t *I;
  member_list_t *J;

  /* Go thru the list of declarations and generate the structure definitions
   * for the data. */
  for (I=dl; I!=NULL; I=I->p) {
    fprintf(output,"\n/* Generated structure for %s: */\n",I->name);
    fprintf(output,"struct struct_%s {\n",I->name);
    /* Reset the marks. */
    for (J=I->members; J!=NULL; J=J->p) J->member->mark = 0;
    /* Loop the the members of the structure. */
    for (J=I->members; J!=NULL; J=J->p) {
      /* Check to see if we've got a union. */
      if (J->member->uname == NULL) member_gen(J->member,I->members,I->name);
      else if (J->member->mark == 0) union_gen(J->member,I->members,I->name);
      }
    fprintf(output,"  };\n");
    /* Set up a typedef for this struct. */
    fprintf(output,"typedef struct struct_%s %s_t;\n\n",I->name,I->name);
    }

  }

LOCAL_FUNCTION void
union_gen(member,members,structname)
member_t *member;
member_list_t *members;
char *structname;
{
  member_list_t *I;

  /* This union member's name must be unique for the reading and printing
   * routines. */
  for (I=members; I->member!=member; I=I->p) {
    if (!strcmp(I->member->name,member->name)) {
      fprintf(stderr, 
              "the structure \"%s\" has two members with the name \"%s\"\n",
              structname,member->name);
      error("union_gen: aborting");
      }
    }

  fprintf(output,"  union {\n");
  for (I=members; I!=NULL; I=I->p) {
    if (!I->member->uname) continue;
    /* Make sure that the selector name has been properly declared. */
    find_variable_index(I->member->uselname,
                        member,members,structname,member->name);
    if (!strcmp(member->uname,I->member->uname)) {
      fprintf(output,"  ");
      member_gen(I->member,members,structname);
      I->member->mark = 1;
      }
    }
  fprintf(output,"    } %s;\n",member->uname);
  }

/* This needs not only the member of the structure to be written out,
 * but it also needs the entire member list for the current structure
 * so that it can check to make sure that all variable indices have
 * been declared as integers within the structure. */
LOCAL_FUNCTION void
member_gen(member,members,structname)
member_t *member;
member_list_t *members;
char *structname;
{
  int i;
  int n_pointer;
  int n_index = 0;
  int n_constant_index = 0;
  int n_variable_index = 0;
  index_list_t *I;
  char outline[STRING_LENGTH];
  char tmp[STRING_LENGTH];

  if (!member) {
    warn("no member to member_gen");
    return;
    }

  /* Count the number of indices. */
  for (I=member->indices; I!=NULL; I=I->p) {
    n_index++;
    if (I->index.type == TYPE_CONSTANT) n_constant_index++;
    else if (I->index.type == TYPE_VARIABLE) {
      find_variable_index(I->index.v.v,member,members,structname,member->name);
      n_variable_index++;
      }
    else {
      error("member_gen: internal error: bad index.type");
      }
    }

  if ((n_constant_index != 0)&&(n_variable_index != 0)) {
    fprintf(stderr,"Problem for member %s of struct %s:\n",
            structname,member->name);
    fprintf(stderr,"constant and variable indices cannot\n");
    fprintf(stderr,"be mixed in the same array.\n");
    error("member_gen: aborting");
    }

  n_pointer = member->pointer + n_variable_index + n_constant_index;

  /* Build up the output string. */
  outline[0] = '\0';
  /* Is there a qualifier to print out? */
  if (member->qualifier == Q_STRUCT) { 
    strcat(outline,"struct ");
    }
  else if (member->qualifier == Q_UNSIGNED) { 
    strcat(outline,"unsigned ");
    }
  else if (member->qualifier == Q_SIGNED) { 
    strcat(outline,"signed ");
    }
  else if (member->qualifier != Q_NONE) {
    error("member_gen: internal error: bad qualifier");
    }
  /* Is this member one of the basic types?  If not, then,
   * if no qualifier has been given, insert the struct qualifier
   * and modify the name. */
  if ((!basic_type(member->type))&&(member->qualifier == Q_NONE)) {
    strcat(outline,"struct struct_");
    strcat(outline,member->type);
    }
  else if (!strcmp(member->type,"string")) {
    if (member->qualifier != Q_NONE) {
      fprintf(stderr,"cannot qualify a string\n");
      exit(1);
      }
    strcat(outline,"char *");
    }
  else if (!strcmp(member->type,"boolean")) {
    if (member->qualifier != Q_NONE) {
      fprintf(stderr,"cannot qualify a boolean\n");
      exit(1);
      }
    strcat(outline,"int");
    }
  else {
    strcat(outline,member->type);
    }
  strcat(outline," ");
  for (i=0; i<n_pointer; i++) strcat(outline,"*");
  strcat(outline,member->name);
  for (I=member->indices; I!=NULL; I=I->p) {
      if (I->index.type == TYPE_CONSTANT)
        sprintf(outline,"%s /*[%d]*/",strcpy(tmp,outline),I->index.v.c);
      else
        sprintf(outline,"%s /*[%s]*/",strcpy(tmp,outline),I->index.v.v);
    }
  if (!member->uname) fprintf(output,"  %s;\n",outline);
  else fprintf(output,"  %s; /* %s == %s */\n",
               outline,member->uselname,member->uselval);

  }

/* This checks to make sure that all variable indices and union
 * selectors have been declared as integers within the structure. */
LOCAL_FUNCTION void
find_variable_index(indexname,member,members,structname,membername)
char *indexname;
member_t *member;
member_list_t *members;
char *structname;
char *membername;
{
  member_list_t *I;

  for (I=members; I->member!=member; I=I->p) {
    if (!strcmp(I->member->name,indexname)) {
      if (  (I->member->pointer != 0)
          ||(I->member->indices != NULL)
          ||(I->member->uname != NULL)
          ||(strcmp(I->member->type,"int") != 0)
         ) {
        fprintf(stderr,"Problem for the variable %s in member %s of struct %s:\n",
                indexname, membername, structname);
        fprintf(stderr,"The variable was not properly declared.\n");
        fprintf(stderr,"Restrictions apply to varibles which are used\n");
        fprintf(stderr,"as dimensions or as union selectors:\n");
        fprintf(stderr,"They must be of type \"int\" and cannot be part of");
        fprintf(stderr,"union.\n");
        error("find_variable_index: aborting");
        }
      /* Everying is OK, so we can return control. */
      return;
      }
    }

  fprintf(stderr,"Problem for the index %s in member %s of struct %s:\n",
          indexname, membername, structname);
  fprintf(stderr,"the variable index was not properly declared:\n");
  fprintf(stderr,"it must be declared within the structure and before\n");
  fprintf(stderr,"the array uses it is declared\n");
  error("find_variable_index: aborting");
  }
