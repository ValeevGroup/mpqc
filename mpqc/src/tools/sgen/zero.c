/* yet another modification of the print routines.
 * routine zero_x will zero out arrays in x */

/* $Log$
 * Revision 1.1  1993/12/29 12:53:58  etseidl
 * Initial revision
 *
 * Revision 1.3  1992/07/20  18:37:46  seidl
 * add support for string arrays
 *
 * Revision 1.2  1992/06/17  23:07:52  jannsen
 * modified to generate clean code
 *
 * Revision 1.1.1.1  1992/03/17  18:11:20  seidl
 * Struct GENerator 2.0
 *
 * Revision 1.1  1992/03/17  18:11:19  seidl
 * Initial revision
 *
 * Revision 1.2  1992/01/07  18:01:43  seidl
 * include sgen_util.gbl
 *
 * Revision 1.1  1991/12/23  20:50:35  seidl
 * Initial revision
 * */
 
static char rcsid[] = "$Id$";

#include <stdio.h>
#include <tmpl.h>
#include <string.h>
#include "types.h"
#include "global.h"

#include "sgen_util.gbl"

#include "zero.gbl"
#include "zero.lcl"

GLOBAL_FUNCTION VOID
zero_gen()
{
  declaration_t *I;
  member_list_t *J;
  char outfile[FILENAME_MAX];
  FILE *include;

  /* Convert filename to the name of the output file. */
  strcpy(outfile,BaseName);
  strcat(outfile,"zero");
  strcat(outfile,".c");

  /* Open the output file. */
  output = fopen(outfile,"w");
  if (!output) {
    fprintf(stderr,"%s: couldn't open the output file: %s\n",progname,outfile);
    error(NULL);
    }

  /* Convert filename to the name of an include file. */
  strcpy(outfile,BaseName);
  strcat(outfile,"zero");
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
  fprintf(output,"#include \"%szero.h\"\n",BaseName);

  /* Go thru the list of declarations and generate the zero functions
   * for the data. */
  for (I=dl; I!=NULL; I=I->p) {
    /* Write to the C include file. */
    fprintf(include,"#if defined(NO_PROTO)\n");
    fprintf(include,"  void zero_%s();\n",I->name);
    fprintf(include,"#else\n");
    fprintf(include,"  void zero_%s(%s_t *%s);\n",
            I->name,I->name,I->name,I->name,I->name);
    fprintf(include,"#endif\n");

    if (is_excluded("zero",I)) continue;
    /* Write to the C source file. */
    fprintf(output,"\n/* Generated zero function for %s: */\n",I->name);
    fprintf(output,"void\n");
    fprintf(output,"zero_%s(_%s)\n",I->name,I->name);
    fprintf(output,"%s_t *_%s;\n",I->name,I->name);
    fprintf(output,"{\n");
    fprintf(output,"  typedef int boolean;\n");
    fprintf(output,"  typedef char * string;\n");
    declare_indices_less1basic(I->members,I->name);
    /* Loop over the members of the structure. */
    for (J=I->members; J!=NULL; J=J->p) {
      zero_member(J->member,I->name);
      }
    fprintf(output,"  }\n");
    }

  fclose(output);

  }

LOCAL_FUNCTION VOID
zero_member(member,structname)
member_t *member;
char *structname;
{
  char indices[STRING_LENGTH];
  char spaces[STRING_LENGTH];
  char range[STRING_LENGTH];
  index_list_t *I;
  int i,array_type;

  spaces[0] = '\0';

  /* If this is a union then we conditional execute the following code. */
  if (member->uname) {
    fprintf(output,"  if (_%s->%s==%s) {\n",
            structname,member->uselname,member->uselval);
    strcat(spaces,"  ");
    }

  indices[0] = '\0';
  array_type=0;
  for (I=member->indices,i=0;  I!=NULL; I=I->p,i++) {
    array_type++;
    }

  for (I=member->indices,i=0;  I!=NULL; I=I->p,i++) {
    fprintf(output,"%s  if(%s!=0) {\n",spaces /* } */
                                ,index_dimension(structname,&I->index));
    strcat(spaces,"  ");
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

  zero_elementary(member,structname,indices,spaces,array_type,range);

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
zero_elementary(member,structname,indices,spaces,array_type,range)
member_t *member;
char *structname;
char *indices;
char *spaces;
int array_type;
char *range;
{

  if (basic_type(member->type) && !member->pointer) {
    if(array_type) {
      fprintf(output,"%s  bzero(%s%s,sizeof(%s)*%s);\n",
        spaces,member_name(structname,member),indices,member->type,range);
      }
    }
  else if(!member->pointer) {
    fprintf(output,"%s  zero_%s(&(%s%s));\n",
            spaces,member->type,member_name(structname,member),indices); 
    }
  else if(member->pointer) {
    fprintf(output,"%s  %s%s=NULL;\n",spaces,
                     member_name(structname,member),indices);
    }
 }
