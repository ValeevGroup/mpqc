
/* $Log$
 * Revision 1.2  1994/10/18 23:03:58  etseidl
 * fix many warnings, use memset rather than bzero
 *
 * Revision 1.1.1.1  1993/12/29  12:53:58  etseidl
 * SC source tree 0.1
 *
 * Revision 1.3  1992/06/17  23:07:33  jannsen
 * modified to generate clean code
 *
 * Revision 1.2  1992/03/30  23:08:46  seidl
 * merge in sandia changes
 *
 * Revision 1.4  91/09/28  18:13:38  cljanss
 * Output files now have shorter names.  Input file rootnames
 * of 8 or fewer characters will have output filenames of 14
 * or fewer characters.
 * 
 * Revision 1.3  91/09/28  16:40:34  cljanss
 * new naming convention is used to keep names <= 14 characters.
 * 
 * Revision 1.2  91/06/15  21:46:20  seidl
 * search: cets061591
 * added a few lines to pretty-print long arrays
 * 
 * Revision 1.1  1991/06/15  21:13:57  janssen
 * Initial revision
 * */

#include <stdio.h>
#include <string.h>
#include <tmpl.h>
#include "types.h"
#include "global.h"

#include "print.gbl"
#include "print.lcl"

#include "error.gbl"
#include "sgen_util.gbl"

GLOBAL_FUNCTION VOID
print_gen()
{
  declaration_t *I;
  member_list_t *J;
  char outfile[FILENAME_MAX];
  FILE *include;

  /* Convert filename to the name of the output file. */
  strcpy(outfile,BaseName);
  strcat(outfile,"prnt");
  strcat(outfile,".c");

  /* Open the output file. */
  output = fopen(outfile,"w");
  if (!output) {
    fprintf(stderr,"%s: couldn't open the output file: %s\n",progname,outfile);
    error(NULL);
    }

  /* Convert filename to the name of an include file. */
  strcpy(outfile,BaseName);
  strcat(outfile,"prnt");
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
  fprintf(output,"#include \"%sprnt.h\"\n",BaseName);
  fprintf(output,"extern int sgen_print_nindent;\n");
  fprintf(output,"#define SPI sgen_print_indent(fp)\n");

  /* Go thru the list of declarations and generate the print functions
   * for the data. */
  for (I=dl; I!=NULL; I=I->p) {
    /* Write to the C include file. */
    fprintf(include,"#if defined(NO_PROTO)\n");
    fprintf(include,"  void print_%s();\n",I->name);
    fprintf(include,"#else\n");
    fprintf(include,"  void print_%s(FILE *fp,%s_t *%s);\n",
            I->name,I->name,I->name);
    fprintf(include,"#endif\n");

    if (is_excluded("print",I)) continue;
    /* Write to the C source file. */
    fprintf(output,"\n/* Generated print function for %s: */\n",I->name);
    fprintf(output,"void\n");
    fprintf(output,"print_%s(fp,_%s)\n",I->name,I->name);
    fprintf(output,"FILE *fp;\n");
    fprintf(output,"%s_t *_%s;\n",I->name,I->name);
    fprintf(output,"{\n");
    declare_indices(I->members,I->name);
    fprintf(output,"  int orig_indent = sgen_print_nindent;\n");
    /* Loop over the members of the structure. */
    for (J=I->members; J!=NULL; J=J->p) {
      if (J!=I->members) fprintf(output,"  SPI;\n");
      print_member(J->member,I->name);
      }
    fprintf(output,"  }\n");
    }

  fclose(output);

  }

LOCAL_FUNCTION VOID
print_member(member,structname)
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
    fprintf(output,"  if (_%s->%s==%s) {\n",
            structname,member->uselname,member->uselval);
    strcat(spaces,"  ");
    }

  if (basic_type(member->type)) {
    fprintf(output,"%s  fprintf(fp,\"%s = \");\n",spaces,member->name);
    adjust_indent(strlen(member->name)+3,NULL);
    }
  else {
    fprintf(output,"%s  fprintf(fp,\"%s:\");\n",spaces,member->name);
    adjust_indent(strlen(member->name)+1,NULL);
    }

  indices[0] = '\0';
  for (I=member->indices,i=0;  I!=NULL; I=I->p,i++) {
    strcat(indices,"[");
    indices[strlen(indices)+1] = '\0';
    indices[strlen(indices)] = 'i' + i;
    strcat(indices,"]");
    fprintf(output,"%s  fprintf(fp,\"[\");\n",spaces);
    adjust_indent(1,spaces);
    strcat(spaces,"  ");
    fprintf(output,"%sfor (%c=0; %c<%s; %c++)",spaces,
            'i'+i,'i'+i,index_dimension(structname,&I->index),'i'+i);
    fprintf(output,"  {\n");
    if ((!basic_type(member->type)) || I->p != NULL) {
      fprintf(output,"%s  if (%c!=0) SPI;\n",spaces,'i'+i);
      }

    }

  /* Print out a test so that NULL data is not printed out. */
  if (member->pointer) fprintf(output,"%s  if (",spaces);
  for (i=0; i<member->pointer; i++) {
     if (i!=0) fprintf(output,"\n%s      &&",spaces);
     else fprintf(output,"  ");
     fprintf(output,"(%s%s != NULL)",
             member_name(structname,member),indices);
     strcat(indices,"[0]");
     }
  if (member->pointer) {
    fprintf(output,") {\n");
    strcat(spaces,"  ");
    }

  print_elementary(member,structname,indices,spaces);

  if (member->pointer) {
    fprintf(output,"%s  }\n",spaces);
    spaces[strlen(spaces)-2] = '\0';
    }

  for (I=member->indices,i=0;  I!=NULL; I=I->p,i++) {

/* added by e seidl Sat Jun 15 17:45:07 EDT 1991 (cets061591) */
    int j=strlen(indices);

    if(!i) {
      char ets_ind = indices[j-2];
      fprintf(output,
       "%s  if(!((%c+1)%c8) && %c) {\n%s    fprintf(fp,\"\\n\");\n%s    SPI;\
        \n%s    }\n",
        spaces,ets_ind,'%',ets_ind,spaces,spaces,spaces);
      }
/* end cets061591 */

    fprintf(output,"%s  }\n",spaces);
    spaces[strlen(spaces)-2] = '\0';
    if ((!basic_type(member->type)) || I != member->indices) {
      fprintf(output,"%s  SPI;\n",spaces);
      }
    fprintf(output,"%s  fprintf(fp,\"]\\n\");\n",spaces);
    adjust_indent(-1,spaces);
    }


  if (basic_type(member->type)) {
    if (member->indices == NULL) {
      fprintf(output,"%s  fprintf(fp,\"\\n\");\n",spaces);
      }
    }

  if (member->uname) {
    fprintf(output,"  }\n");
    }

  reset_indent();
  }

LOCAL_FUNCTION VOID
print_elementary(member,structname,indices,spaces)
member_t *member;
char *structname;
char *indices;
char *spaces;
{
#define NO_BASIC_TYPES_DIRECT
#ifdef BASIC_TYPES_DIRECT
  if (!strcmp(member->type,"double")) {
    if (member->qualifier != Q_NONE) error("bad qualifier for double");
    fprintf(output,"%s  fprintf(fp,\" %%f\",%s%s);\n",spaces,
            member_name(structname,member),indices);
    }
  else if (!strcmp(member->type,"float")) {
    if (member->qualifier != Q_NONE) error("bad qualifier for float");
    fprintf(output,"%s  fprintf(fp,\" %%f\",(double)%s%s);\n",spaces,
            member_name(structname,member),indices);
    }
  else if (!strcmp(member->type,"int")) {
    if (member->qualifier == Q_SIGNED || member->qualifier == Q_NONE) {
      fprintf(output,"%s  fprintf(fp,\" %%d\",%s%s);\n",spaces,
              member_name(structname,member),indices);
      }
    else if (member->qualifier == Q_UNSIGNED) {
      fprintf(output,"%s  fprintf(fp,\" %%u\",%s%s);\n",spaces,
              member_name(structname,member),indices);
      }
    else error("bad qualifier for int");
    }
  else if (!strcmp(member->type,"long")) {
    if (member->qualifier == Q_SIGNED || member->qualifier == Q_NONE) {
      fprintf(output,"%s  fprintf(fp,\" %%ld\",%s%s);\n",spaces,
              member_name(structname,member),indices);
      }
    else if (member->qualifier == Q_UNSIGNED) {
      fprintf(output,"%s  fprintf(fp,\" %%lu\",%s%s);\n",spaces,
              member_name(structname,member),indices);
      }
    else error("bad qualifier for long");
    }
  else if (!strcmp(member->type,"char")) {
    fprintf(output,"%s  fprintf(fp,\"%%c\",%s%s);\n",spaces,
            member_name(structname,member),indices);
    }
  else {
    fprintf(output,"%s  fprintf(fp,\"(\");\n",spaces);
    adjust_indent(1,spaces);
    fprintf(output,"%s  print_%s(fp,&(%s%s));\n",
            spaces,member->type,member_name(structname,member),indices);
    fprintf(output,"%s  SPI;\n",spaces);
    fprintf(output,"%s  fprintf(fp,\")\\n\");\n",spaces);
    adjust_indent(-1,spaces);
    }
#else  /* BASIC_TYPES_DIRECT */
  if (basic_type(member->type)) {
    fprintf(output,"%s  print_%s(fp,&(%s%s));\n",
            spaces,member->type,member_name(structname,member),indices);
    }
  else {
    fprintf(output,"%s  fprintf(fp,\"(\");\n",spaces);
    adjust_indent(1,spaces);
    fprintf(output,"%s  print_%s(fp,&(%s%s));\n",
            spaces,member->type,member_name(structname,member),indices);
    fprintf(output,"%s  SPI;\n",spaces);
    fprintf(output,"%s  fprintf(fp,\")\\n\");\n",spaces);
    adjust_indent(-1,spaces);
    }
#endif /* BASIC_TYPES_DIRECT */
 }

LOCAL_FUNCTION VOID
adjust_indent(n,spaces)
int n;
char *spaces;
{
  if (spaces) fprintf(output,"%s  sgen_print_nindent += %d;\n",spaces,n);
  else fprintf(output,"  sgen_print_nindent += %d;\n",n);
  }

LOCAL_FUNCTION VOID
reset_indent()
{
  fprintf(output,"  sgen_print_nindent = orig_indent;\n");
  }
