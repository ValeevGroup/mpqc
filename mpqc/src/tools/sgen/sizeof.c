/* yet another manifestation of the print routines.
 * routine sizeof_x will return (obviously) the size of x
 * in bytes, including arrays, etc. */

/* $Log$
 * Revision 1.1  1993/12/29 12:53:58  etseidl
 * Initial revision
 *
 * Revision 1.2  1992/06/17  23:07:49  jannsen
 * modified to generate clean code
 *
 * Revision 1.1.1.1  1992/03/17  18:11:10  seidl
 * Struct GENerator 2.0
 *
 * Revision 1.1  1992/03/17  18:11:09  seidl
 * Initial revision
 *
 * Revision 1.1  1991/12/23  22:42:57  seidl
 * Initial revision
 * */
 
static char rcsid[] = "$Id$";

#include <stdio.h>
#include <tmpl.h>
#include <string.h>
#include "types.h"
#include "global.h"

#include "sizeof.gbl"
#include "sizeof.lcl"

GLOBAL_FUNCTION VOID
sizeof_gen()
{
  declaration_t *I;
  member_list_t *J;
  char outfile[FILENAME_MAX];
  FILE *include;

  /* Convert filename to the name of the output file. */
  strcpy(outfile,BaseName);
  strcat(outfile,"size");
  strcat(outfile,".c");

  /* Open the output file. */
  output = fopen(outfile,"w");
  if (!output) {
    fprintf(stderr,"%s: couldn't open the output file: %s\n",progname,outfile);
    error(NULL);
    }

  /* Convert filename to the name of an include file. */
  strcpy(outfile,BaseName);
  strcat(outfile,"size");
  strcat(outfile,".h");

  /* Open the output file. */
  include = fopen(outfile,"w");
  if (!include) {
    fprintf(stderr,"%s: couldn't open the include file: %s\n",progname,outfile);
    error(NULL);
    }

  /* Include the following files. */
  fprintf(output,"\n");
  fprintf(output,"#include <stdio.h>\n",BaseName);
  fprintf(output,"#include <string.h>\n",BaseName);
  fprintf(output,"#include <util/sgen/sgen.h>\n");
  fprintf(output,"#include \"%s.h\"\n",BaseName);
  fprintf(output,"#include \"%ssize.h\"\n",BaseName);

  /* Go thru the list of declarations and generate the sizeof functions
   * for the data. */
  for (I=dl; I!=NULL; I=I->p) {
    /* Write to the C include file. */
    fprintf(include,"#if defined(NO_PROTO)\n");
    fprintf(include,"  int sizeof_%s();\n",I->name);
    fprintf(include,"#else\n");
    fprintf(include,"  int sizeof_%s(%s_t *%s);\n",
            I->name,I->name,I->name,I->name,I->name);
    fprintf(include,"#endif\n");

    if (is_excluded("sizeof",I)) continue;
    /* Write to the C source file. */
    fprintf(output,"\n/* Generated sizeof function for %s: */\n",I->name);
    fprintf(output,"int\n");
    fprintf(output,"sizeof_%s(_%s)\n",I->name,I->name);
    fprintf(output,"%s_t *_%s;\n",I->name,I->name);
    fprintf(output,"{\n");
    declare_indices(I->members,I->name);
    fprintf(output,"  int sizeof_%s=0;\n",I->name);
    fprintf(output,"  typedef int boolean;\n");
    fprintf(output,"  typedef char * string;\n");
    /* Loop over the members of the structure. */
    for (J=I->members; J!=NULL; J=J->p) {
      sizeof_member(J->member,I->name);
      }
    fprintf(output,"\n  return sizeof_%s;\n",I->name);
    fprintf(output,"  }\n");
    }

  fclose(output);

  }

LOCAL_FUNCTION VOID
sizeof_member(member,structname)
member_t *member;
char *structname;
{
  char indices[STRING_LENGTH];
  char spaces[STRING_LENGTH];
  char range[STRING_LENGTH];
  char stars[STRING_LENGTH];
  index_list_t *I;
  int i,len,array_type;

  spaces[0] = '\0';

  /* If this is a union then we conditional execute the following code. */
  if (member->uname) {
    fprintf(output,"  if (_%s->%s==%s) {\n",
            structname,member->uselname,member->uselval);
    strcat(spaces,"  ");
    }

  indices[0] = '\0';
  stars[0] = '\0';
  array_type=0;
  for (I=member->indices,i=0;  I!=NULL; I=I->p,i++) {
    array_type++;
    strcat(stars,"*");
    }
  for(i=0; i<member->pointer ; i++) strcat(stars,"*");

  for (I=member->indices,i=0;  I!=NULL; I=I->p,i++) {
    fprintf(output,"%s  if(%s!=0) {\n",spaces /* } */
                                ,index_dimension(structname,&I->index));
    strcat(spaces,"  ");
    if(basic_type(member->type)) {
      fprintf(output,"%s  sizeof_%s += sizeof(%s %s);\n",spaces,structname,
                  member->type,stars);
      }
    else {
      fprintf(output,"%s  sizeof_%s += sizeof(%s_t %s);\n",spaces,structname,
                  member->type,stars);
      }
    if(len=strlen(stars)) stars[len-1]='\0';

    fprintf(output,"%s  if(%s%s!=NULL) {\n",spaces, /* } */
                  member_name(structname,member),indices);
    strcat(spaces,"  ");
    if(I->p != NULL || !basic_type(member->type) || member->pointer) {
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

  /* Print out a test so that NULL data is not printed out. */
  if (member->pointer) fprintf(output,"%s  if (",spaces);
  for (i=0; i<member->pointer; i++) {
     fprintf(output,"(%s%s != NULL)",
             member_name(structname,member),indices);
     }
  if (member->pointer) {
    fprintf(output,") {\n"); /* } */
    strcat(spaces,"  ");
    }

  sizeof_elementary(member,structname,indices,spaces,array_type,range);

  if (member->pointer) {
    for(i=0; i < member->pointer ; i++) {
    /* { */
      fprintf(output,"%s  }\n",spaces);
      spaces[strlen(spaces)-2] = '\0';
      }
    }

  for (I=member->indices;  I!=NULL ; I=I->p) {
    char *lastbrack;
    if(I->p != NULL || !basic_type(member->type) || member->pointer) {
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
    fprintf(output,"%s  }\n",spaces); /* Close dimension test (for ptrs)*/
    spaces[strlen(spaces)-2] = '\0';

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
sizeof_elementary(member,structname,indices,spaces,array_type,range)
member_t *member;
char *structname;
char *indices;
char *spaces;
int array_type;
char *range;
{

  if (!strcmp(member->type,"string")) {
    fprintf(output,"%s  sizeof_%s += strlen(%s%s);\n",spaces,structname,
            member_name(structname,member),indices);
    }
  else if (basic_type(member->type)) {
    if(!array_type) {
      fprintf(output,"%s  sizeof_%s += sizeof(%s);\n",spaces,structname,
            member->type);
      }
    else if(member->pointer) {
      fprintf(output,"%s  sizeof_%s += sizeof(%s *)*%s;\n",spaces,structname,
            member->type,range);
      }
    else {
      fprintf(output,"%s  sizeof_%s += sizeof(%s)*%s;\n",spaces,structname,
            member->type,range);
      }
    }
  else {
    fprintf(output,"%s  sizeof_%s += sizeof_%s(&(%s%s));\n",
            spaces,structname,member->type,
            member_name(structname,member),indices);
    }
 }
