
/* This is copied from ets's bread and slighted modified to do more
 * general read type operations. The bread routine and others
 * now call this routine with the appropiate arguments to generate
 * code. */

/* $Log$
 * Revision 1.5  1994/10/31 17:55:34  etseidl
 * include string.h in generated c file
 *
 * Revision 1.4  1994/10/21  20:39:49  cljanss
 * Work arounds for IRIX 6.0 IDO/C++ bugs.
 *
 * Revision 1.3  1994/10/18  23:03:49  etseidl
 * fix many warnings, use memset rather than bzero
 *
 * Revision 1.2  1994/10/14  18:28:28  etseidl
 * replace bzero with memset
 *
 * Revision 1.1.1.1  1993/12/29  12:53:58  etseidl
 * SC source tree 0.1
 *
 * Revision 1.5  1993/04/29  00:38:22  jannsen
 * fixed potential overwrite and cleaned up generated code
 *
 * Revision 1.4  1992/07/20  18:37:33  seidl
 * add support for string arrays
 *
 * Revision 1.3  1992/06/17  23:07:23  jannsen
 * modified to generate clean code
 *
 * Revision 1.2  1992/03/30  23:05:58  seidl
 * merge in sandia changes
 *
 * Revision 1.9  1991/11/18  17:22:36  cljanss
 * fixed some null pointer initialization problems
 *
 * Revision 1.8  91/09/28  20:39:21  cljanss
 * fixed the way pointers are nulled
 * 
 * Revision 1.7  91/09/28  18:13:35  cljanss
 * Output files now have shorter names.  Input file rootnames
 * of 8 or fewer characters will have output filenames of 14
 * or fewer characters.
 * 
 * Revision 1.6  91/09/28  16:40:28  cljanss
 * new naming convention is used to keep names <= 14 characters.
 * 
 * Revision 1.5  91/09/28  16:10:50  cljanss
 * fixed a bug
 * 
 * Revision 1.4  1991/08/08  22:29:00  cljanss
 * now memory is always allocated, even if a pointer is nonnull.  This is
 * done to make it easier to parallize, since it seems that memory in
 * nodes (in the Linda impl. of PICL) is dirty when the process starts.
 *
 * Revision 1.3  1991/08/08  22:21:06  cljanss
 * basicname and read or write name can now be different
 *
 * Revision 1.2  1991/07/19  19:06:24  cljanss
 * Fixed handling of NULL tpointerargs.
 *
 * Revision 1.1  1991/07/19  17:49:39  cljanss
 * Initial revision
 * */

/* a blatant copy of clj's print.c hacked to do binary writes */

/* The old log from bread.c:
 * Revision 1.3  1991/06/20  16:33:45  seidl
 * reads in all pointers. if a pointer is non-null then what it points
 * to is read in.  allocates storage for the buffer you're reading
 * into if it is null
 *
 * Revision 1.2  1991/06/17  18:33:04  seidl
 * fix how things like double_matrix are handled
 *
 * Revision 1.1  1991/06/17  17:24:15  seidl
 * Initial revision
 * */

#include <stdio.h>
#include <string.h>
#include <tmpl.h>
#include "types.h"
#include "global.h"

#include "gen_read.gbl"
#include "gen_read.lcl"

#include "sgen_util.gbl"
#include "error.gbl"

GLOBAL_FUNCTION VOID
general_read_gen(suffix,readname,basicname,
                 protoargsfmt,funcargsfmt,funcdecsfmt,elemargsfmt,tpointerargs,
                 strpointerargs,useoffsets,extra_includes)
char *suffix;
char *readname;
char *basicname;
char *protoargsfmt;
char *funcargsfmt;
char *funcdecsfmt;
char *elemargsfmt;
char *tpointerargs;
char *strpointerargs;
int useoffsets;
char *extra_includes;
{
  declaration_t *I;
  member_list_t *J;
  char outfile[FILENAME_MAX];
  FILE *include;
  char protoargs[STRING_LENGTH];
  char funcargs[STRING_LENGTH];
  char funcdecs[STRING_LENGTH];
  char *functiontype;

  if (useoffsets) functiontype = "int";
  else            functiontype = "void";

  /* Convert filename to the name of the output file. */
  sprintf(outfile,"%s%s.c",BaseName,suffix);

  /* Open the output file. */
  output = fopen(outfile,"w");
  if (!output) {
    fprintf(stderr,"%s: couldn't open the output file: %s\n",progname,outfile);
    error(NULL);
    }

  /* Convert filename to the name of an include file. */
  sprintf(outfile,"%s%s.h",BaseName,suffix);

  /* Open the output file. */
  include = fopen(outfile,"w");
  if (!include) {
    fprintf(stderr,"%s: couldn't open the include file: %s\n",progname,outfile);
    error(NULL);
    }

  /* Include the following files. */
  fprintf(output,"\n");
  fprintf(output,"#include <stdio.h>\n");
  fprintf(output,"#include <stdlib.h>\n");
  fprintf(output,"#include <string.h>\n");
  fprintf(output,"#include <util/sgen/sgen.h>\n");
  if (extra_includes) fprintf(output,"%s",extra_includes);
  fprintf(output,"#include \"%s.h\"\n",BaseName);
  fprintf(output,"#include \"%s%s.h\"\n",BaseName,suffix);

  /* Go thru the list of declarations and generate the general_read functions
   * for the data. */
  for (I=dl; I!=NULL; I=I->p) {
    /* Set up the argument lists and declarations. */
    sprintf(protoargs,protoargsfmt,I->name,I->name);
    sprintf(funcargs,funcargsfmt,I->name,I->name);
    sprintf(funcdecs,funcdecsfmt,I->name,I->name);

    /* Write to the C include file. */
    fprintf(include,"#if defined(NO_PROTO)\n");
    fprintf(include,"  %s %s_%s();\n",functiontype,readname,I->name);
    fprintf(include,"#else\n");
    fprintf(include,"  %s %s_%s(%s);\n",
            functiontype,readname,I->name,protoargs);
    fprintf(include,"#endif\n");

    if (is_excluded(readname,I)) continue;
    /* Write to the C source file. */
    fprintf(output,"\n/* Generated %s function for %s: */\n",readname,I->name);
    fprintf(output,"%s\n",functiontype);
    fprintf(output,"%s_%s(%s)\n",readname,I->name,funcargs);
    fprintf(output,"%s",funcdecs);
    fprintf(output,"{\n");
    fprintf(output,"  typedef int boolean;\n");
    fprintf(output,"  typedef char * string;\n");
    fprintf(output,"  int * %s_test_pointer();\n",readname);
    declare_indices_less1basic(I->members,I->name);
    if (useoffsets) fprintf(output,"  int _offset_init= *_offset;\n");
#if 0
    fprintf(output,
     "  char *error_msg=\"%s_%s: error: read through NULL pointer\";\n",
        readname,I->name);
#endif

  /* Loop over the members of the structure. */
    for (J=I->members; J!=NULL; J=J->p) {
      general_read_member(readname,basicname,
                          elemargsfmt,tpointerargs,strpointerargs,
                          J->member,I->name);
      }
    if (useoffsets) fprintf(output,"  return(*_offset-_offset_init);\n");
    else            fprintf(output,"  return;\n");
    fprintf(output,"  }\n");
    }

  fclose(output);

  }

LOCAL_FUNCTION VOID
general_read_member(readname,basicname, elemargsfmt,tpointerargs,
    strpointerargs,member,structname)
char *readname;
char *basicname;
char *elemargsfmt;
char *tpointerargs;
char *strpointerargs;
member_t *member;
char *structname;
{
  char indices[STRING_LENGTH];
  char spaces[STRING_LENGTH];
  char range[STRING_LENGTH];
  char stars[STRING_LENGTH];
  index_list_t *I;
  int i;
  int array_type;
  char tpointerargswithcomma[STRING_LENGTH];

  if (tpointerargs) {
    strcpy(tpointerargswithcomma,tpointerargs);
    strcat(tpointerargswithcomma,",");
    }
  else {
    tpointerargswithcomma[0] = '\0';
    }

  spaces[0] = '\0';

  /* If this is a union then we conditional execute the following code. */
  if (member->uname) {
    fprintf(output,"  if (_%s->%s==%s) {\n",
            structname,member->uselname,member->uselval);
    strcat(spaces,"  ");
    }

  indices[0] = '\0';
  stars[0]='\0';
  array_type=0;
  for (I=member->indices,i=0;  I!=NULL ; I=I->p,i++) {
    array_type++;
    strcat(stars,"*");
    }
  for(i=0; i<member->pointer ; i++) strcat(stars,"*");

  for (I=member->indices,i=0;  I!=NULL ; I=I->p,i++) {
    fprintf(output,"%s  if(%s!=0) {\n",spaces
                                ,index_dimension(structname,&I->index));
    strcat(spaces,"  ");
    fprintf(output,"%s  if(%s_test_pointer(%s",
            spaces,basicname,tpointerargswithcomma);
    if(basic_type(member->type)) {
      fprintf(output,"sizeof(%s %s))!=NULL) {\n",member->type,stars);
      }
    else {
      fprintf(output,"sizeof(%s_t %s))!=NULL) {\n",member->type,stars);
      }
    strcat(spaces,"  ");
    if(basic_type(member->type)) {
      fprintf(output,"%s  %s%s = (%s %s) malloc(sizeof(%s",spaces,
        member_name(structname,member),
        indices,member->type,stars,member->type);
      stars[strlen(stars)-1]='\0';
      fprintf(output," %s)*%s);\n",stars,
                                 index_dimension(structname,&I->index));
      fprintf(output,"%s  memset(%s%s,0,sizeof(%s",spaces,
        member_name(structname,member),
        indices,member->type);
      fprintf(output," %s)*%s);\n",stars,
                                 index_dimension(structname,&I->index));
      }
    else {
      fprintf(output,"%s  %s%s = (%s_t %s) malloc(sizeof(%s_t",spaces,
        member_name(structname,member),
        indices,member->type,stars,member->type);
      stars[strlen(stars)-1]='\0';
      fprintf(output," %s)*%s);\n",stars,
                                 index_dimension(structname,&I->index));
      fprintf(output,"%s  memset(%s%s,0,sizeof(%s_t",spaces,
        member_name(structname,member),
        indices,member->type);
      fprintf(output," %s)*%s);\n",stars,
                                 index_dimension(structname,&I->index));
      }
    fprintf(output,"%s  sgen_chkmalloc(%s%s);\n",
            spaces,member_name(structname,member),indices);

    if(I->p != NULL || !basic_type(member->type) || member->pointer ||
         !strcmp(member->type,"string")) {
      strcat(indices,"[");
      indices[strlen(indices)+1] = '\0';
      indices[strlen(indices)] = 'i' + i;
      strcat(indices,"]");
      strcat(spaces,"  ");
      fprintf(output,"%sfor (%c=0; %c<%s; %c++)",spaces,
            'i'+i,'i'+i,index_dimension(structname,&I->index),'i'+i);
      fprintf(output,"  {\n"); /* } */
      }
    strcpy(range,index_dimension(structname,&I->index));
    }

 /* if this is pointer, only read in what's at the pointer if non-NULL */
  if(member->pointer) {
    for(i=0; i < member->pointer ; i++) {
      fprintf(output,"%s  if(%s_test_pointer(%s",
              spaces,basicname,tpointerargswithcomma);
      if(basic_type(member->type)) {
        fprintf(output,"sizeof(%s %s))!=NULL) {\n",member->type,stars);
        }
      else {
        fprintf(output,"sizeof(%s_t %s))!=NULL) {\n",member->type,stars);
        }
      strcat(spaces,"  ");
      fprintf(output,"%s  if (",spaces);
      fprintf(output,"(%s%s == NULL)",
           member_name(structname,member),indices);
      fprintf(output,") {\n");
      if(basic_type(member->type)) {
        fprintf(output,"%s    %s%s = (%s %s) malloc(sizeof(%s",spaces,
          member_name(structname,member),
          indices,member->type,stars,member->type);
        }
      else {
        fprintf(output,"%s    %s%s = (%s_t %s) malloc(sizeof(%s_t",spaces,
          member_name(structname,member),
          indices,member->type,stars,member->type);
        }
      stars[strlen(stars)-1]='\0';
      fprintf(output," %s));\n",stars);
      fprintf(output,"%s    sgen_chkmalloc(%s%s);\n",
              spaces,member_name(structname,member),indices);
      fprintf(output,"%s    }\n",spaces);
      strcat(indices,"[0]");
      }
    }

  general_read_elementary(readname,basicname,elemargsfmt,strpointerargs,
                          member,structname,indices,spaces,array_type,range);

  if(member->pointer) {
    for(i=0; i < member->pointer ; i++) {
      fprintf(output,"%s  }\n",spaces);
      spaces[strlen(spaces)-2] = '\0';
      }
    }

  if (member->indices && member->pointer) {
    char *lastbrack;

    /* Cleave off the last two indices. */
    if (basic_type(member->type)) {
      lastbrack = strrchr(indices,'[');
      if (lastbrack) *lastbrack = '\0';
      lastbrack = strrchr(indices,'[');
      if (lastbrack) *lastbrack = '\0';
      }
    }

  for (I=member->indices;  I!=NULL ; I=I->p) {
    char *lastbrack;
    if(I->p != NULL || !basic_type(member->type) || member->pointer ||
         !strcmp(member->type,"string")) {
      fprintf(output,"%s  }\n",spaces); /* Close for loop (for ptrs). */
      spaces[strlen(spaces)-2] = '\0';
      }
    fprintf(output,"%s  }\n",spaces); /* Close test_pointer (for ptrs). */
    spaces[strlen(spaces)-2] = '\0';
    /* Cleave off the last index. */
    if (!basic_type(member->type)) {
      lastbrack = strrchr(indices,'[');
      if (lastbrack) *lastbrack = '\0';
      }
    fprintf(output,"%s  else %s%s = NULL;\n",spaces,
            member_name(structname,member),indices);
    fprintf(output,"%s  }\n",spaces); /* Close dimension test (for ptrs)*/
    spaces[strlen(spaces)-2] = '\0';
    if (member->pointer || (I->p == NULL)) {
      fprintf(output,"%s  else %s%s = NULL;\n",spaces,
            member_name(structname,member),indices);
      }
     else {
     fprintf(output,"%s  /* Skipped: else %s%s = NULL;*/\n",spaces,
           member_name(structname,member),indices);
     }
    /* Cleave off the last index. */
    if (basic_type(member->type)) {
      lastbrack = strrchr(indices,'[');
      if (lastbrack) *lastbrack = '\0';
      }
    }

  if (member->uname) {
    fprintf(output,"%s  }\n",spaces);
    }

  }

LOCAL_FUNCTION VOID
general_read_elementary(readname,basicname,elemargsfmt,strpointerargs,
                        member,structname,indices,spaces,array_type,range)
char *readname;
char *basicname;
char *elemargsfmt;
char *strpointerargs;
member_t *member;
char *structname;
char *indices;
char *spaces;
int array_type;
char *range;
{
  char elemargs[STRING_LENGTH];
  char elemargsfmt2[STRING_LENGTH];
  char selemargsfmt2[STRING_LENGTH];

  sprintf(elemargsfmt2,elemargsfmt,"%s",member_name(structname,member),indices);
  sprintf(selemargsfmt2,strpointerargs,"%s",member_name(structname,member),
                                                                       indices);
  if(member->pointer) {
    if(basic_type(member->type)) {
      sprintf(elemargs,elemargsfmt2,"&");
      fprintf(output,"%s  %s_%s(%s,",spaces,basicname,member->type,elemargs);
      fprintf(output,"sizeof(%s));\n",member->type);
      }
    else {
      sprintf(elemargs,elemargsfmt2,"&");
      fprintf(output,"%s  %s_%s(%s);\n",spaces,readname,member->type,elemargs);
      }
    }
  else if(array_type) {
    if(!strcmp(member->type,"string")) {
      sprintf(elemargs,selemargsfmt2,"");
      fprintf(output,"%s  %s_%s(%s,",spaces,basicname,member->type,elemargs);
      fprintf(output,"sizeof(%s)*1);\n",member->type);
      }
    else if(basic_type(member->type)) {
      sprintf(elemargs,elemargsfmt2,"");
      fprintf(output,"%s  %s_%s(%s,",spaces,basicname,member->type,elemargs);
      fprintf(output,"sizeof(%s)*%s);\n",member->type,range);
      }
    else {
      sprintf(elemargs,elemargsfmt2,"&");
      fprintf(output,"%s  %s_%s(%s);\n",spaces,readname,member->type,elemargs);
      }
    }
  else if(!basic_type(member->type)) {
    sprintf(elemargs,elemargsfmt2,"&");
    fprintf(output,"%s  %s_%s(%s);\n",spaces,readname,member->type,elemargs);
    }
  else {
    sprintf(elemargs,elemargsfmt2,"&");
    fprintf(output,"%s  %s_%s(%s,",spaces,basicname,member->type,elemargs);
    fprintf(output,"sizeof(%s)*%s);\n",member->type,"1");
    }
  }
