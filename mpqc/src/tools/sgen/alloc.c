
/* $Log$
 * Revision 1.2  1994/10/18 23:03:39  etseidl
 * fix many warnings, use memset rather than bzero
 *
 * Revision 1.1.1.1  1993/12/29  12:53:56  etseidl
 * SC source tree 0.1
 *
 * Revision 1.3  1992/10/09  19:37:21  seidl
 * alloc memory for strings
 *
 * Revision 1.2  1992/06/17  23:07:07  jannsen
 * modified to generate clean code
 *
 * Revision 1.1.1.1  1992/03/17  18:09:59  seidl
 * Struct GENerator 2.0
 *
 * Revision 1.1  1992/03/17  18:09:56  seidl
 * Initial revision
 *
 * Revision 1.4  1992/01/15  18:00:08  seidl
 * no longer include strings.h as this is not defined by ANSI
 *
 * Revision 1.3  1991/12/23  13:49:35  seidl
 * fix minor bug for unions
 *
 * Revision 1.2  1991/12/20  16:44:18  seidl
 * have generated routines include string.h.  this is where ANSI
 * says strtok should be
 *
 * Revision 1.1  1991/12/20  16:18:23  seidl
 * Initial revision
 *
 * Revision 1.3  91/09/28  18:11:55  cljanss
 * Output files now have shorter names.  Input file rootnames
 * of 8 or fewer characters will have output filenames of 14
 * or fewer characters.
 * 
 * Revision 1.2  91/09/28  16:40:21  cljanss
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

#include "alloc.gbl"
#include "alloc.lcl"

#include "error.gbl"
#include "sgen_util.gbl"

GLOBAL_FUNCTION VOID
alloc_gen()
{
  declaration_t *I;
  char outfile[FILENAME_MAX];
  FILE *include;

  /* Convert filename to the name of the output file. */
  strcpy(outfile,BaseName);
  strcat(outfile,"allc");
  strcat(outfile,".c");

  /* Open the output file. */
  output = fopen(outfile,"w");
  if (!output) {
    fprintf(stderr,"%s: couldn't open the output file: %s\n",progname,outfile);
    error(NULL);
    }

  /* Convert filename to the name of an include file. */
  strcpy(outfile,BaseName);
  strcat(outfile,"allc");
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
  fprintf(output,"#include <stdlib.h>\n");
  fprintf(output,"#include <stdarg.h>\n");
  fprintf(output,"#include <string.h>\n");
  fprintf(output,"#include <util/sgen/sgen.h>\n");
  fprintf(output,"#include \"%s.h\"\n",BaseName);
  fprintf(output,"#include \"%sallc.h\"\n",BaseName);
  fprintf(output,"#include \"%sinit.h\"\n",BaseName);
  fprintf(output,"#define NAMES_LENGTH 256\n");

  /* Go thru the list of declarations and generate the ip functions
   * for the data. */
  for (I=dl; I!=NULL; I=I->p) {
    /* Write to the C include file. */
    fprintf(include,"#if defined(NO_PROTO)\n");
    fprintf(include,"  int allocbn_%s();\n",I->name);
    fprintf(include,"  int alloc_%s();\n",I->name);
    fprintf(include,"#else\n");
    fprintf(include,"  int allocbn_%s(%s_t *_%s, char *names, ...);\n",
            I->name,I->name,I->name);
    fprintf(include,"  int alloc_%s(%s_t *_%s);\n",
            I->name,I->name,I->name);
    fprintf(include,"#endif\n");

    if (is_excluded("alloc",I)) continue;
    /* Generate the vararg version of the alloc function. */
    alloc_by_name(I);

    alloc_initialized(I);
    }

  fclose(output);
  }

/* This generates a function to allocate storage using the names of
 * the dimension members. */
LOCAL_FUNCTION VOID
alloc_by_name(dec)
declaration_t *dec;
{
  char *elsestring;
  member_list_t *I;

#define F fprintf
#define O output
  /* The RS/6000 xlc compiler requires the following distinction for
   * var arg function declarations. */
  F(O,"int\n");
  F(O,"#if defined(NO_PROTO)\n");
    F(O,"allocbn_%s(_%s, names)\n",dec->name,dec->name);
    F(O,"%s_t *_%s;\n",dec->name,dec->name);
    F(O,"char *names;\n");
  F(O,"#else\n");
    F(O,"allocbn_%s(",dec->name);
    F(O,"%s_t *_%s,",dec->name,dec->name);
    F(O," char *names,");
    F(O," ...)\n");
  F(O,"#endif\n");

  F(O,"{\n");

  F(O,"  va_list args;\n");
  F(O,"  char tmp[NAMES_LENGTH],*tok;\n");

  F(O,"  if (strlen(names)+1 >= NAMES_LENGTH) {\n");
  F(O,"    fprintf(stderr,\"ERROR: allocbn_%s: names(=\\\"%%s\\\") too long\",names);\n",dec->name);
  F(O,"    return -1;\n");
  F(O,"    }\n");

  /* Initilize the passed data. */
  F(O,"  init_%s(_%s);\n",dec->name,dec->name);

  /* Find each token in the names list. */
  F(O,"  strcpy(tmp,names);\n");
  F(O,"  va_start(args, names);\n");
  F(O,"  tok = strtok(tmp,\" \\t\\n\");\n");
  F(O,"  while (tok) {\n");

  /* Generate code to go thru each member and compare the name to member.
   * If the names match we get the corresponding va_arg.
   * NOTE: If types don't match there will be problems. */
  elsestring = "";
  for (I=dec->members; I!=NULL; I=I->p) {
    if ( I->member->pointer || I->member->indices ) continue;
    F(O,"    %sif (!strcmp(tok,\"%s\")) {\n",elsestring,I->member->name);
    if(I->member->uname) 
      F(O,"      _%s->%s.%s = va_arg(args,%s);\n",
        dec->name,I->member->uname,I->member->name,type_name(I->member));
    else if(!strcmp(I->member->type,"string")) {
      F(O,"      char *str=va_arg(args,char*);\n");
      F(O,"      _%s->%s = malloc(strlen(str)+1);\n",dec->name,I->member->name);
      F(O,"      if(_%s->%s==NULL) return -1;\n",dec->name,I->member->name);
      F(O,"      strcpy(_%s->%s,str);\n",dec->name,I->member->name);
      }
    else
      F(O,"      _%s->%s = va_arg(args,%s);\n",
        dec->name,I->member->name,type_name(I->member));
    F(O,"      }\n");
    elsestring = "else ";
    }

  if (dec->members) {
    F(O,"    else {\n");
    F(O,"      fprintf(stderr,\"allocbn_%s: cannot init \\\"%%s\\\"\\n\",tok);\n", dec->name);
    F(O,"      return -1;\n");
    F(O,"      }\n");
    }

  F(O,"    tok = strtok(NULL,\" \\t\\n\");\n");
  F(O,"    }\n");
  F(O,"  va_end(args);\n");
  F(O,"  return alloc_%s(_%s);\n",dec->name,dec->name);
  F(O,"  }\n");

  }

/* This generates a function that takes a pointer to a preinitialized
 * datum and allocates storage consistent with the dimensions given
 * in the preinitialized datum. */
LOCAL_FUNCTION VOID
alloc_initialized(dec)
declaration_t *dec;
{
  member_list_t *J;

  fprintf(output,"\n/* Generated alloc function for %s: */\n",dec->name);
  fprintf(output,"int\n");

  fprintf(output,"alloc_%s(_%s)\n",dec->name,dec->name);
  fprintf(output,"%s_t *_%s;\n",dec->name,dec->name);

  fprintf(output,"{\n");
  declare_indices_less1(dec->members,dec->name);

  /* Loop over the members of the structure. */
  for (J=dec->members; J!=NULL; J=J->p) {
    alloc_member(J->member,dec->name);
    }

  /* If we get here then return that everything is OK. */
  fprintf(output,"  return 0;\n");
  fprintf(output,"  }\n");
  }

LOCAL_FUNCTION VOID
alloc_member(member,structname)
member_t *member;
char *structname;
{
  char indices[STRING_LENGTH];
  char spaces[STRING_LENGTH];
  char pointer[STRING_LENGTH];
  char pointermalloc[STRING_LENGTH];
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

  /* Generate code to allocate storage for the data. */
  indices[0] = '\0';
  strcpy(pointermalloc," ");
  for (i=0; i<n_indices+member->pointer-1; i++) strcat(pointermalloc,"*");
  strcpy(pointer,pointermalloc);
  strcat(pointer,"*");
  for (I=member->indices,i=0;  I!=NULL; I=I->p,i++) {
    strcat(spaces,"  ");
    if (strlen(pointer)==1) pointer[0] = '\0';
    if (strlen(pointermalloc)==1) pointermalloc[0] = '\0';
    fprintf(output,"%s%s%s = (%s%s)malloc(sizeof(%s%s)*%s);\n",
            spaces,
            member_name(structname,member),
            indices,
            type_name(member),
            pointer,
            type_name(member),
            pointermalloc,
            index_dimension(structname,&I->index));
    fprintf(output,"%sif (!%s%s) return -1;\n",
            spaces,
            member_name(structname,member),
            indices);
    if (strlen(pointer)) pointer[strlen(pointer)-1] = '\0';
    if (strlen(pointermalloc)) pointermalloc[strlen(pointermalloc)-1] = '\0';
    strcat(indices,"[");
    indices[strlen(indices)+1] = '\0';
    indices[strlen(indices)] = 'i' + i;
    strcat(indices,"]");
    if (I->p!=NULL) fprintf(output,"%sfor (%c=0; %c<%s; %c++)",spaces,
                        'i'+i,'i'+i,
                        index_dimension(structname,&I->index),
                        'i'+i);
    if (I->p!=NULL) fprintf(output,"  {\n");
    }

  for (I=member->indices;  I!=NULL; I=I->p) {
    if (I!=member->indices) fprintf(output,"%s  }\n",spaces);
    spaces[strlen(spaces)-2] = '\0';
    }

  /* Finish the test for nonzero dimensions. */
  if (member->indices) {
    spaces[strlen(spaces)-2] = '\0';
    fprintf(output,"%s    }\n",spaces);
    }

  if (member->uname) {
    fprintf(output,"    }\n");
    }

  }

