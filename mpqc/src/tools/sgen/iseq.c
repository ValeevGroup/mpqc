/* yet another modification of the print routines, this time to compare
 * two data types.  routine iseq_x will return 1 if the two x's are
 * the same, 0 otherwise */

/* $Log$
 * Revision 1.2  1994/10/18 23:03:55  etseidl
 * fix many warnings, use memset rather than bzero
 *
 * Revision 1.1.1.1  1993/12/29  12:53:58  etseidl
 * SC source tree 0.1
 *
 * Revision 1.4  1992/06/17  23:07:31  jannsen
 * modified to generate clean code
 *
 * Revision 1.3  1992/03/30  23:08:23  seidl
 * merge in sandia changes
 *
 * Revision 1.1  1992/01/09  15:42:36  cljanss
 * Initial revision
 *
 * Revision 1.3  1991/12/23  20:49:13  seidl
 * do not perform the comparison if the member is a pointer
 *
 * Revision 1.2  1991/12/23  13:44:25  seidl
 * add a test to make reasonably sure you're not going to try reading through
 * a null pointer
 *
 * Revision 1.1  1991/12/21  21:39:12  seidl
 * Initial revision
 * */
 
#include <stdio.h>
#include <string.h>
#include <tmpl.h>
#include "types.h"
#include "global.h"

#include "iseq.gbl"
#include "iseq.lcl"

#include "error.gbl"
#include "sgen_util.gbl"

GLOBAL_FUNCTION VOID
iseq_gen()
{
  declaration_t *I;
  member_list_t *J;
  char outfile[FILENAME_MAX];
  FILE *include;

  /* Convert filename to the name of the output file. */
  strcpy(outfile,BaseName);
  strcat(outfile,"iseq");
  strcat(outfile,".c");

  /* Open the output file. */
  output = fopen(outfile,"w");
  if (!output) {
    fprintf(stderr,"%s: couldn't open the output file: %s\n",progname,outfile);
    error(NULL);
    }

  /* Convert filename to the name of an include file. */
  strcpy(outfile,BaseName);
  strcat(outfile,"iseq");
  strcat(outfile,".h");

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
  fprintf(output,"#include \"%siseq.h\"\n",BaseName);

  /* Go thru the list of declarations and generate the iseq functions
   * for the data. */
  for (I=dl; I!=NULL; I=I->p) {
    /* Write to the C include file. */
    fprintf(include,"#if defined(NO_PROTO)\n");
    fprintf(include,"  int iseq_%s();\n",I->name);
    fprintf(include,"#else\n");
    fprintf(include,"  int iseq_%s(%s_t *%s_1,%s_t *%s_2);\n",
            I->name,I->name,I->name,I->name,I->name);
    fprintf(include,"#endif\n");

    if (is_excluded("iseq",I)) continue;
    /* Write to the C source file. */
    fprintf(output,"\n/* Generated iseq function for %s: */\n",I->name);
    fprintf(output,"int\n");
    fprintf(output,"iseq_%s(_%s_1,_%s_2)\n",I->name,I->name,I->name);
    fprintf(output,"%s_t *_%s_1;\n",I->name,I->name);
    fprintf(output,"%s_t *_%s_2;\n",I->name,I->name);
    fprintf(output,"{\n");
    declare_indices_for_nonpointer(I->members,I->name);
    /* Loop over the members of the structure. */
    for (J=I->members; J!=NULL; J=J->p) {
      if(!J->member->pointer) comp_member(J->member,I->name);
      }
    fprintf(output,"  return 1;\n");
    fprintf(output,"  }\n");
    }

  fclose(output);

  }

LOCAL_FUNCTION VOID
comp_member(member,structname)
member_t *member;
char *structname;
{
  char indices[STRING_LENGTH];
  char spaces[STRING_LENGTH];
  index_list_t *I;
  int i;

  spaces[0] = '\0';

  /* If this is a union then we conditional execute the following code. */
  if (member->uname) {
    fprintf(output,"  if (_%s_2->%s==%s) {\n",
            structname,member->uselname,member->uselval);
    strcat(spaces,"  ");
    }

  indices[0] = '\0';
  for (I=member->indices,i=0;  I!=NULL; I=I->p,i++) {
    strcat(spaces,"  ");
    fprintf(output,"%sif (%s && %s%s==NULL &&",spaces,
             index_dimension2(structname,1,&I->index),
             member_name2(structname,1,member),indices);
    fprintf(output,"%s%s!=NULL) return 0;\n",
             member_name2(structname,2,member),indices);
    strcat(indices,"[");
    indices[strlen(indices)+1] = '\0';
    indices[strlen(indices)] = 'i' + i;
    strcat(indices,"]");
    fprintf(output,"%sfor (%c=0; %c<%s; %c++)",spaces,
            'i'+i,'i'+i,index_dimension2(structname,2,&I->index),'i'+i);
    fprintf(output,"  {\n"); /* } */
    }

  /* Print out a test so that NULL data is not printed out. */
  if (member->pointer) fprintf(output,"%s  if (",spaces);
  for (i=0; i<member->pointer; i++) {
     if (i!=0) fprintf(output,"\n%s      &&",spaces);
     else fprintf(output,"  ");
     fprintf(output,"(%s%s != NULL)",
             member_name2(structname,2,member),indices);
     strcat(indices,"[0]");
     }
  if (member->pointer) {
    fprintf(output,") {\n"); /* } */
    strcat(spaces,"  ");
    }

  iseq_elementary(member,structname,indices,spaces);

  if (member->pointer) {
    /* { */
    fprintf(output,"%s  }\n",spaces);
    spaces[strlen(spaces)-2] = '\0';
    }

  for (I=member->indices,i=0;  I!=NULL; I=I->p,i++) {
    fprintf(output,"%s  }\n",spaces);
    spaces[strlen(spaces)-2] = '\0';
    }


  if (member->uname) {
    fprintf(output,"%s  }\n",spaces);
    }

  }

LOCAL_FUNCTION VOID
iseq_elementary(member,structname,indices,spaces)
member_t *member;
char *structname;
char *indices;
char *spaces;
{
  char n1[80];
  strcpy(n1,member_name2(structname,1,member));

  if (basic_type(member->type)) {
    fprintf(output,"%s  if(!iseq_%s(%s%s,%s%s)) return 0;\n",spaces,
            member->type,n1,indices,
            member_name2(structname,2,member),indices);
    }
  else {
    fprintf(output,"%s  if(!iseq_%s(&(%s%s),&(%s%s))) return 0;\n",
            spaces,member->type,n1,indices,
            member_name2(structname,2,member),indices);
    }
 }
