/* yet another modification of the print routines.
 * routine assign_x will return set x1 equal to x2 */

/* $Log$
 * Revision 1.1  1993/12/29 12:53:56  etseidl
 * Initial revision
 *
 * Revision 1.6  1992/07/20  18:37:30  seidl
 * add support for string arrays
 *
 * Revision 1.5  1992/06/17  23:07:10  jannsen
 * modified to generate clean code
 *
 * Revision 1.4  1992/04/01  17:28:09  seidl
 * fix initialization of ind
 *
 * Revision 1.3  1992/03/30  23:01:37  seidl
 * merge in sandia changes
 *
 * Revision 1.2  1992/01/09  15:47:18  cljanss
 * now sgen_util.gbl is included
 *
 * Revision 1.1  1992/01/09  15:42:36  cljanss
 * Initial revision
 *
 * Revision 1.2  1991/12/23  20:48:12  seidl
 * if the member is a pointer, allocate memory for it, but do not
 * assign anything to it
 *
 * Revision 1.1  1991/12/23  13:49:52  seidl
 * Initial revision
 * */
 
static char rcsid[] = "$Id$";

#include <stdio.h>
#include <tmpl.h>
#include <string.h>
#include "types.h"
#include "global.h"

#include "assign.gbl"
#include "assign.lcl"

#include "sgen_util.gbl"

GLOBAL_FUNCTION VOID
assign_gen()
{
  declaration_t *I;
  member_list_t *J;
  char outfile[FILENAME_MAX];
  FILE *include;

  /* Convert filename to the name of the output file. */
  strcpy(outfile,BaseName);
  strcat(outfile,"asgn");
  strcat(outfile,".c");

  /* Open the output file. */
  output = fopen(outfile,"w");
  if (!output) {
    fprintf(stderr,"%s: couldn't open the output file: %s\n",progname,outfile);
    error(NULL);
    }

  /* Convert filename to the name of an include file. */
  strcpy(outfile,BaseName);
  strcat(outfile,"asgn");
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
  fprintf(output,"#include <stdlib.h>\n",BaseName);
  fprintf(output,"#include <string.h>\n",BaseName);
  fprintf(output,"#include <util/sgen/sgen.h>\n",BaseName);
  fprintf(output,"#include \"%s.h\"\n",BaseName);
  fprintf(output,"#include \"%sasgn.h\"\n",BaseName);

  /* Go thru the list of declarations and generate the assign functions
   * for the data. */
  for (I=dl; I!=NULL; I=I->p) {
    /* Write to the C include file. */
    fprintf(include,"#if defined(NO_PROTO)\n");
    fprintf(include,"  int assign_%s();\n",I->name);
    fprintf(include,"#else\n");
    fprintf(include,"  int assign_%s(%s_t *%s_1,%s_t *%s_2);\n",
            I->name,I->name,I->name,I->name,I->name);
    fprintf(include,"#endif\n");

    if (is_excluded("assign",I)) continue;
    /* Write to the C source file. */
    fprintf(output,"\n/* Generated assign function for %s: */\n",I->name);
    fprintf(output,"int\n");
    fprintf(output,"assign_%s(_%s_1,_%s_2)\n",I->name,I->name,I->name);
    fprintf(output,"%s_t *_%s_1;\n",I->name,I->name);
    fprintf(output,"%s_t *_%s_2;\n",I->name,I->name);
    fprintf(output,"{\n");
    fprintf(output,"  typedef int boolean;\n");
    fprintf(output,"  typedef char * string;\n");
    declare_indices(I->members,I->name);
    /* Loop over the members of the structure. */
    for (J=I->members; J!=NULL; J=J->p) {
      assign_member(J->member,I->name);
      }
    fprintf(output,"\n  return 0;\n");
    fprintf(output,"  }\n");
    }

  fclose(output);

  }

LOCAL_FUNCTION VOID
assign_member(member,structname)
member_t *member;
char *structname;
{
  char indices[STRING_LENGTH];
  char spaces[STRING_LENGTH];
  char range[STRING_LENGTH];
  char stars[STRING_LENGTH];
  index_list_t *I;
  int i,array_type;

  spaces[0] = '\0';

  /* If this is a union then we conditional execute the following code. */
  if (member->uname) {
    fprintf(output,"  if (_%s_2->%s==%s) {\n",
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
                                ,index_dimension2(structname,2,&I->index));
    strcat(spaces,"  ");
    fprintf(output,"%s  if(%s%s!=NULL) {\n",spaces, /* } */
                  member_name2(structname,2,member),indices);
    strcat(spaces,"  ");
    if(basic_type(member->type)) {
      fprintf(output,"%s  %s%s = (%s %s) malloc(sizeof(%s",spaces,
        member_name2(structname,1,member),
        indices,member->type,stars,member->type);
      stars[strlen(stars)-1]='\0';
      fprintf(output," %s)*%s);\n",stars,
                                 index_dimension2(structname,2,&I->index));
      fprintf(output,"%s  if(%s%s==NULL) return(-1);\n",spaces,
        member_name2(structname,1,member),indices);
      }
    else {
      fprintf(output,"%s  %s%s = (%s_t %s) malloc(sizeof(%s_t",spaces,
        member_name2(structname,1,member),
        indices,member->type,stars,member->type);
      stars[strlen(stars)-1]='\0';
      fprintf(output," %s)*%s);\n",stars,
                                 index_dimension2(structname,2,&I->index));
      fprintf(output,"%s  if(%s%s==NULL) return(-1);\n",spaces,
        member_name2(structname,1,member),indices);
      }
    if(I->p != NULL || !basic_type(member->type) || member->pointer) {
      strcat(indices,"[");
      indices[strlen(indices)+1] = '\0';
      indices[strlen(indices)] = 'i' + i;
      strcat(indices,"]");
      strcat(spaces,"  ");
      fprintf(output,"%sfor (%c=0; %c<%s; %c++)",spaces,
            'i'+i,'i'+i,index_dimension2(structname,2,&I->index),'i'+i);
      fprintf(output,"  {\n"); /* } */
      }
    strcpy(range,index_dimension2(structname,1,&I->index));
    }

  /* Print out a test so that NULL data is not printed out. */
  if (member->pointer) fprintf(output,"%s  if (",spaces);
  for (i=0; i<member->pointer; i++) {
     fprintf(output,"(%s%s != NULL)",
             member_name2(structname,2,member),indices);
     }
  if (member->pointer) {
    fprintf(output,") {\n"); /* } */
    strcat(spaces,"  ");
    }

  if(array_type && !member->pointer && basic_type(member->type)) {
    strcat(indices,"[");
    indices[strlen(indices)+1] = '\0';
    indices[strlen(indices)] = 'i' + array_type - 1;
    strcat(indices,"]");
    }

  assign_elementary(member,structname,indices,spaces,array_type,range);

  if(array_type) indices[strlen(indices)-3] = '\0';

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
    fprintf(output,"%s  else %s%s = NULL;\n",spaces,
            member_name2(structname,1,member),indices);
    fprintf(output,"%s  }\n",spaces); /* Close dimension test (for ptrs)*/
    spaces[strlen(spaces)-2] = '\0';
    if (member->pointer || (I->p == NULL)) {
      fprintf(output,"%s  else %s%s = NULL; /* DT. */\n",spaces,
            member_name2(structname,1,member),indices);
      }
     else {
     fprintf(output,"%s  /* Skipped: else %s%s = NULL;*/ /* DT. */\n",spaces,
           member_name2(structname,1,member),indices);
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
assign_elementary(member,structname,indices,spaces,array_type,range)
member_t *member;
char *structname;
char *indices;
char *spaces;
int array_type;
char *range;
{
  char n1[80];
  char n2[80];
  char ind[2];

  ind[0]='i'+array_type-1;
  ind[1]='\0';

  strcpy(n1,member_name2(structname,1,member));
  strcpy(n2,member_name2(structname,2,member));

  if(member->pointer) {
    }
  else if (!strcmp(member->type,"string")) {
    if(!array_type) {
      fprintf(output,"%s  if(%s%s!=NULL) {\n",spaces,n2,indices);
      fprintf(output,"%s    %s%s = (char *) malloc(strlen(%s%s)+1);\n",spaces,
                                    n1,indices,n2,indices);
      fprintf(output,"%s    if(%s%s==NULL) return(-1);\n",spaces,n1,indices);
      fprintf(output,"%s    strcpy(%s%s,%s%s);\n",spaces,
                                    n1,indices,n2,indices);
      fprintf(output,"%s    }\n",spaces);
      }
    else {
      fprintf(output,"%s  for(%s=0;%s<%s;%s++) {\n",spaces,ind,ind,range,ind);
      fprintf(output,"%s    if(%s%s!=NULL) {\n",spaces,n2,indices);
      fprintf(output,"%s      %s%s = (char *) malloc(strlen(%s%s)+1);\n",
                       spaces, n1,indices,n2,indices);
      fprintf(output,"%s      if(%s%s==NULL) return(-1);\n",spaces,n1,indices);
      fprintf(output,"%s      strcpy(%s%s,%s%s);\n",spaces,
                                  n1,indices,n2,indices);
      fprintf(output,"%s      }\n",spaces);
      fprintf(output,"%s    }\n",spaces);
      }
    }
  else if (basic_type(member->type)) {
    if(!array_type) {
      fprintf(output,"%s  %s%s=%s%s;\n",spaces,n1,indices,n2,indices);
      }
    else {
      fprintf(output,"%s  for(%s=0;%s<%s;%s++)\n",spaces,ind,ind,range,ind);
      fprintf(output,"%s    %s%s=%s%s;\n",spaces,n1,indices, n2,indices);
      }
    }
  else {
    fprintf(output,"%s  assign_%s(&(%s%s),&(%s%s));\n",
            spaces,member->type,n1,indices,
            member_name2(structname,2,member),indices);
    }
 }
