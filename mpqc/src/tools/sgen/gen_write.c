
/* This is copied from ets's bwrite and slighted modified to do more
 * general read type operations. The bread routine and others
 * now call this routine with the appropiate arguments to generate
 * code. */

/* $Log$
 * Revision 1.2  1994/10/18 23:03:51  etseidl
 * fix many warnings, use memset rather than bzero
 *
 * Revision 1.1.1.1  1993/12/29  12:53:57  etseidl
 * SC source tree 0.1
 *
 * Revision 1.3  1992/07/20  18:37:35  seidl
 * add support for string arrays
 *
 * Revision 1.2  1992/06/17  23:07:25  jannsen
 * modified to generate clean code
 *
 * Revision 1.1.1.1  1992/03/17  18:10:28  seidl
 * Struct GENerator 2.0
 *
 * Revision 1.1  1992/03/17  18:10:27  seidl
 * Initial revision
 *
 * Revision 1.1  1991/12/20  16:18:23  seidl
 * Initial revision
 *
 * Revision 1.5  91/09/28  18:13:36  cljanss
 * Output files now have shorter names.  Input file rootnames
 * of 8 or fewer characters will have output filenames of 14
 * or fewer characters.
 * 
 * Revision 1.4  91/09/28  16:40:29  cljanss
 * new naming convention is used to keep names <= 14 characters.
 * 
 * Revision 1.3  91/09/28  16:10:33  cljanss
 * fixed a bug
 * 
 * Revision 1.2  1991/08/08  22:21:06  cljanss
 * basicname and read or write name can now be different
 *
 * Revision 1.1  1991/07/19  17:49:39  cljanss
 * Initial revision
 * */

/* a blatant copy of clj's print.c hacked to do binary writes */

/* The old log from bwrite.c:
 * Revision 1.5  1991/06/22  06:08:58  seidl
 * write out pointers, no longer exit if a pointer is null, just skip
 * writing that element
 *
 * Revision 1.4  1991/06/17  18:33:24  seidl
 * fix how things like double_matrix are handled
 *
 * Revision 1.3  1991/06/17  17:24:22  seidl
 * add ability to recursively write structs of structs
 *
 * Revision 1.2  1991/06/16  20:31:59  seidl
 * first working version
 *
 * Revision 1.1  1991/06/16  16:07:13  seidl
 * Initial revision
 * */

#include <stdio.h>
#include <string.h>
#include <tmpl.h>
#include "types.h"
#include "global.h"

#include "gen_write.gbl"
#include "gen_write.lcl"

#include "sgen_util.gbl"
#include "error.gbl"

GLOBAL_FUNCTION VOID
general_write_gen(suffix,writename,basicname,
                 protoargsfmt,funcargsfmt,funcdecsfmt,elemargsfmt,
                 selemargsfmt,useoffsets,extra_includes)
char *suffix;
char *writename;
char *basicname;
char *protoargsfmt;
char *funcargsfmt;
char *funcdecsfmt;
char *elemargsfmt;
char *selemargsfmt;
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
  fprintf(output,"#include <util/sgen/sgen.h>\n");
  fprintf(output,"#include \"%s.h\"\n",BaseName);
  fprintf(output,"#include \"%s%s.h\"\n",BaseName,suffix);

  /* Go thru the list of declarations and generate the general_write functions
   * for the data. */
  for (I=dl; I!=NULL; I=I->p) {
    /* Set up the argument lists and declarations. */
    sprintf(protoargs,protoargsfmt,I->name,I->name);
    sprintf(funcargs,funcargsfmt,I->name,I->name);
    sprintf(funcdecs,funcdecsfmt,I->name,I->name);

    /* Write to the C include file. */
    fprintf(include,"#if defined(NO_PROTO)\n");
    fprintf(include,"  %s %s_%s();\n",functiontype,writename,I->name);
    fprintf(include,"#else\n");
    fprintf(include,"  %s %s_%s(%s);\n",
            functiontype,writename,I->name,protoargs);
    fprintf(include,"#endif\n");

    if (is_excluded(writename,I)) continue;
    /* Write to the C source file. */
    fprintf(output,"\n/* Generated %s function for %s: */\n",writename,I->name);
    fprintf(output,"%s\n",functiontype);
    fprintf(output,"%s_%s(%s)\n",writename,I->name,funcargs);
    fprintf(output,"%s",funcdecs);
    fprintf(output,"{\n");
    fprintf(output,"  typedef int boolean;\n");
    fprintf(output,"  typedef char * string;\n");
    declare_indices_less1basic(I->members,I->name);
    if (useoffsets) fprintf(output,"  int _offset_init= *_offset;\n");
#if 0
    fprintf(output,
     "  char *error_msg=\"%s_%s: error: write through NULL pointer\";\n",
        writename,I->name);
#endif

  /* Loop over the members of the structure. */
    for (J=I->members; J!=NULL; J=J->p) {
      general_write_member(writename,basicname,elemargsfmt,selemargsfmt,
                           J->member,I->name);
      }
    if (useoffsets) fprintf(output,"  return(*_offset-_offset_init);\n");
    else            fprintf(output,"  return;\n");
    fprintf(output,"  }\n");
    }

  fclose(output);

  }

LOCAL_FUNCTION VOID
general_write_member(writename,basicname,elemargsfmt,selemargsfmt,
                   member,structname)
char *writename;
char *basicname;
char *elemargsfmt;
char *selemargsfmt;
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
  char elemargs[STRING_LENGTH];

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
    sprintf(elemargs,elemargsfmt,"&",member_name(structname,member),indices);
    fprintf(output,"%s  %s_pointer(%s,",spaces,
         basicname,elemargs);
    if(basic_type(member->type)) {
      fprintf(output,"sizeof(%s %s));\n",member->type,stars);
      }
    else {
      fprintf(output,"sizeof(%s_t %s));\n",member->type,stars);
      }
    stars[strlen(stars)-1]='\0';
    fprintf(output,"%s  if(%s%s!=NULL) {\n",spaces,
                         member_name(structname,member),indices);
    strcat(spaces,"  ");

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

 /* if this is pointer, only print out what's at the pointer if non-NULL */
  if(member->pointer) {
    for(i=0; i < member->pointer ; i++) {
      sprintf(elemargs,elemargsfmt,"&",member_name(structname,member),indices);
      fprintf(output,"%s  %s_pointer(%s,\n",
         spaces,basicname,elemargs);
      if(basic_type(member->type)) {
        fprintf(output,"%s    sizeof(%s %s));\n",spaces,
           member->type,stars);
        }
      else {
        fprintf(output,"%s    sizeof(%s_t %s));\n",spaces,
           member->type,stars);
        }
      strcat(indices,"[0]");
      fprintf(output,"%s  if (",spaces);
      fprintf(output,"(&(%s%s) != NULL)",
           member_name(structname,member),indices);
      fprintf(output,") {\n");
      strcat(spaces,"  ");
      stars[strlen(stars)-1]='\0';
      }
    }

  general_write_elementary(writename,basicname,elemargsfmt,selemargsfmt,
                           member,structname,indices,spaces,array_type,range);

  if(member->pointer) {
    for(i=0; i < member->pointer ; i++) {
      fprintf(output,"%s  }\n",spaces);
      spaces[strlen(spaces)-2] = '\0';
      }
    }

  for (I=member->indices;  I!=NULL ; I=I->p) {
    fprintf(output,"%s  }\n",spaces);
    spaces[strlen(spaces)-2] = '\0';
    fprintf(output,"%s  }\n",spaces);
    spaces[strlen(spaces)-2] = '\0';
    if(I->p != NULL || !basic_type(member->type) || member->pointer ||
        !strcmp(member->type,"string")) {
      fprintf(output,"%s  }\n",spaces);
      spaces[strlen(spaces)-2] = '\0';
      }
    }

  if (member->uname) {
    fprintf(output,"%s  }\n",spaces);
    }

  }

LOCAL_FUNCTION VOID
general_write_elementary(writename,basicname,elemargsfmt,selemargsfmt,
                         member,structname,indices,spaces,array_type,range)
char *writename;
char *basicname;
char *elemargsfmt;
char *selemargsfmt;
member_t *member;
char *structname;
char *indices;
char *spaces;
int array_type;
char *range;
{
  char elemargsfmt2[STRING_LENGTH];
  char selemargsfmt2[STRING_LENGTH];
  char elemargs[STRING_LENGTH];

  sprintf(elemargsfmt2,elemargsfmt,"%s",member_name(structname,member),indices);
  sprintf(selemargsfmt2,selemargsfmt,
        "%s",member_name(structname,member),indices);
  if(member->pointer) {
    if(basic_type(member->type)) {
      sprintf(elemargs,elemargsfmt2,"&");
      fprintf(output,"%s  %s_%s(%s,",spaces,basicname,member->type,elemargs);
      fprintf(output,"sizeof(%s));\n",member->type);
      }
    else {
      sprintf(elemargs,elemargsfmt2,"&");
      fprintf(output,"%s  %s_%s(%s);\n",spaces,writename,member->type,elemargs);
      }
    }
  else if(array_type) {
    if(!strcmp(member->type,"string")) {
      sprintf(elemargs,selemargsfmt2,"");
      fprintf(output,"%s  %s_%s(%s,", spaces,basicname,member->type,elemargs);
      fprintf(output,"sizeof(%s)*1);\n",member->type);
      }
    else if(basic_type(member->type)) {
      sprintf(elemargs,elemargsfmt2,"");
      fprintf(output,"%s  %s_%s(%s,", spaces,basicname,member->type,elemargs);
      fprintf(output,"sizeof(%s)*%s);\n",member->type,range);
      }
    else {
      sprintf(elemargs,elemargsfmt2,"&");
      fprintf(output,"%s  %s_%s(%s);\n",spaces,writename,member->type,elemargs);
      }
    }
  else if(!basic_type(member->type)) {
    sprintf(elemargs,elemargsfmt2,"&");
    fprintf(output,"%s  %s_%s(%s);\n",spaces,writename,member->type,elemargs);
    }
  else {
    sprintf(elemargs,elemargsfmt2,"&");
    fprintf(output,"%s  %s_%s(%s,",spaces,basicname,member->type,elemargs);
    fprintf(output,"sizeof(%s)*%s);\n",member->type,"1");
    }
  }
