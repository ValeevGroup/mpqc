
/* $Log$
 * Revision 1.2  1994/10/18 23:03:48  etseidl
 * fix many warnings, use memset rather than bzero
 *
 * Revision 1.1.1.1  1993/12/29  12:53:57  etseidl
 * SC source tree 0.1
 *
 * Revision 1.6  1993/04/29  00:37:08  jannsen
 * fixed memory leak in generated code and went back to declare_indices_less1
 *
 * Revision 1.5  1992/09/18  10:48:11  seidl
 * fix freeing of string arrays
 *
 * Revision 1.4  1992/07/20  18:37:32  seidl
 * add support for string arrays
 *
 * Revision 1.3  1992/06/17  23:07:19  jannsen
 * modified to generate clean code
 *
 * Revision 1.2  1992/03/30  23:02:58  seidl
 * merge in sandia changes
 *
 * Revision 1.3  91/09/28  18:13:34  cljanss
 * Output files now have shorter names.  Input file rootnames
 * of 8 or fewer characters will have output filenames of 14
 * or fewer characters.
 * 
 * Revision 1.2  91/09/28  16:40:27  cljanss
 * new naming convention is used to keep names <= 14 characters.
 * 
 * Revision 1.1  91/06/15  21:13:57  janssen
 * Initial revision
 *  */

#include <stdio.h>
#include <string.h>
#include <tmpl.h>
#include "types.h"
#include "global.h"

#include "free.gbl"
#include "free.lcl"

#include "sgen_util.gbl"
#include "error.gbl"

GLOBAL_FUNCTION VOID
free_gen()
{
  declaration_t *I;
  char outfile[FILENAME_MAX];
  FILE *include;

  /* Convert filename to the name of the output file. */
  strcpy(outfile,BaseName);
  strcat(outfile,"free");
  strcat(outfile,".c");

  /* Open the output file. */
  output = fopen(outfile,"w");
  if (!output) {
    fprintf(stderr,"%s: couldn't open the output file: %s\n",progname,outfile);
    error(NULL);
    }

  /* Convert filename to the name of an include file. */
  strcpy(outfile,BaseName);
  strcat(outfile,"free");
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
  fprintf(output,"#include \"%sfree.h\"\n",BaseName);
  fprintf(output,"#include \"%sinit.h\"\n",BaseName);
  fprintf(output,"#define NAMES_LENGTH 256\n");

  /* Go thru the list of declarations and generate the ip functions
   * for the data. */
  for (I=dl; I!=NULL; I=I->p) {
    /* Write to the C include file. */
    fprintf(include,"#if defined(NO_PROTO)\n");
    fprintf(include,"  void free_%s();\n",I->name);
    fprintf(include,"#else\n");
    fprintf(include,"  void free_%s(%s_t *_%s);\n",
            I->name,I->name,I->name);
    fprintf(include,"#endif\n");

    if (is_excluded("free",I)) continue;

    free_data(I);
    }

  fclose(output);
  }

/* This generates a function that takes a pointer to a
 * datum and frees storage. */
LOCAL_FUNCTION VOID
free_data(dec)
declaration_t *dec;
{
  member_list_t *J;

  fprintf(output,"\n/* Generated free function for %s: */\n",dec->name);
  fprintf(output,"void\n");

  fprintf(output,"free_%s(_%s)\n",dec->name,dec->name);
  fprintf(output,"%s_t *_%s;\n",dec->name,dec->name);

  fprintf(output,"{\n");
  declare_indices_less1(dec->members,dec->name);

  /* Loop over the members of the structure. */
  for (J=dec->members; J!=NULL; J=J->p) {
    free_member(J->member,dec->name);
    }

  /* If we get here then return that everything is OK. */
  fprintf(output,"  }\n");
  }

LOCAL_FUNCTION VOID
free_member(member,structname)
member_t *member;
char *structname;
{
  char indices[STRING_LENGTH];
  char spaces[STRING_LENGTH];
  char pointer[STRING_LENGTH];
  char pointerfree[STRING_LENGTH];
  index_list_t *I;
  int i;
  int n_indices;

  spaces[0] = '\0';

  /* If this is a union then we conditionally execute the following code. */
  if (member->uname) {
    fprintf(output,"  if (_%s->%s==%s) {\n",
            structname,member->uselname,member->uselval);
    strcat(spaces,"  ");
    }

  /* Compute the number of indices on the member. */
  for ((I=member->indices),(n_indices=0);  I!=NULL; (I=I->p),(n_indices++));

  /* We only deal with the member if all of its indices have nonzero
   * dimension. */
  if (member->indices) {
    fprintf(output,"%s  if (",spaces);
    for (I=member->indices; I!=NULL; I=I->p) {
      fprintf(output,"%s",index_dimension(structname,&I->index));
      if (I->p) fprintf(output," && ");
      }
    fprintf(output,") {\n");
    strcat(spaces,"  ");
    }
  else if(!strcmp(member->type,"string")) {
    fprintf(output,"%s  if(%s!=NULL) free(%s);\n",spaces,
             member_name(structname,member),member_name(structname,member));
    }

  /* free storage for contained types */
  if (!member->indices && !basic_type(member->type)) {
      fprintf(output,"%s  free_%s(&%s);\n",
              spaces,
              member->type,
              member_name(structname,member));
    }

  /* Generate code to free storage for the data. */
  indices[0] = '\0';
  strcpy(pointerfree," ");
  for (i=0; i<n_indices+member->pointer-1; i++) strcat(pointerfree,"*");
  strcpy(pointer,pointerfree);
  strcat(pointer,"*");
  for (I=member->indices,i=0;  I!=NULL; I=I->p,i++) {
    strcat(spaces,"  ");
    if (strlen(pointer)==1) pointer[0] = '\0';
    if (strlen(pointerfree)==1) pointerfree[0] = '\0';
    if (strlen(pointer)) pointer[strlen(pointer)-1] = '\0';
    if (strlen(pointerfree)) pointerfree[strlen(pointerfree)-1] = '\0';
    strcat(indices,"[");
    indices[strlen(indices)+1] = '\0';
    indices[strlen(indices)] = 'i' + i;
    strcat(indices,"]");
    if(strcmp(member->type,"string")) {
      if (I->p!=NULL) fprintf(output,"%sfor (%c=0; %c<%s; %c++)",spaces,
		  'i'+i,'i'+i,
		  index_dimension(structname,&I->index),
		  'i'+i);
      if (I->p!=NULL) fprintf(output,"  {\n");
      }
    else {
      fprintf(output,"%sfor (%c=0; %c<%s; %c++)",spaces, 'i'+i,'i'+i,
		  index_dimension(structname,&I->index), 'i'+i);
      fprintf(output,"  {\n");
      }
    }

  if(!strcmp(member->type,"string")) {
    for (I=member->indices; I!=NULL; I=I->p) {
      fprintf(output,"%s  if(%s%s!=NULL) free(%s%s);\n", spaces,
            member_name(structname,member), indices,
            member_name(structname,member), indices);
      fprintf(output,"%s  }\n",spaces);
      indices[strlen(indices)-3] = '\0';
      spaces[strlen(spaces)-2] = '\0';
      }
    }
  else {
    for (I=member->indices; I!=NULL; I=I->p) {
      if (I!=member->indices) fprintf(output,"%s  }\n",spaces);
      indices[strlen(indices)-3] = '\0';
      if (!basic_type(member->type)) {
        fprintf(output,"%sint ii;\n",spaces);
        fprintf(output,"%sfor (ii=0; ii<%s; ii++) {\n",
                spaces,index_dimension(structname,&I->index));
        fprintf(output,"%s  free_%s(&(%s%s[ii]));\n",
                spaces,
                member->type,
                member_name(structname,member),
                indices);
        fprintf(output,"%s  }\n",spaces);
        }
      fprintf(output,"%sfree(%s%s);\n",
            spaces,
            member_name(structname,member),
            indices);
      spaces[strlen(spaces)-2] = '\0';
      }
    }

  /* Finish the test for nonzero dimensions. */
  if (member->indices) {
    spaces[strlen(spaces)-2] = '\0';
    if(!strcmp(member->type,"string"))
      fprintf(output,"%s    free(%s);\n",
        spaces,member_name(structname,member));
    fprintf(output,"%s    }\n",spaces);
    }

  if (member->uname) {
    fprintf(output,"    }\n");
    }

  }

