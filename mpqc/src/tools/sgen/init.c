
/* $Log$
 * Revision 1.3  1996/03/23 02:38:59  cljanss
 * Everything can now be configured with autoconf.
 *
 * Revision 1.2  1994/10/18 23:03:52  etseidl
 * fix many warnings, use memset rather than bzero
 *
 * Revision 1.1.1.1  1993/12/29  12:53:58  etseidl
 * SC source tree 0.1
 *
 * Revision 1.3  1992/06/17  23:07:27  jannsen
 * modified to generate clean code
 *
 * Revision 1.2  1992/03/30  23:08:12  seidl
 * merge in sandia changes
 *
 * Revision 1.3  91/09/28  18:13:37  cljanss
 * Output files now have shorter names.  Input file rootnames
 * of 8 or fewer characters will have output filenames of 14
 * or fewer characters.
 * 
 * Revision 1.2  91/09/28  16:40:30  cljanss
 * new naming convention is used to keep names <= 14 characters.
 * 
 * Revision 1.1  91/06/15  21:13:57  janssen
 * Initial revision
 *  */

#include <stdio.h>
#include <string.h>
#include "types.h"
#include "global.h"

#include "init.gbl"
#include "init.lcl"

#include "sgen_util.gbl"
#include "error.gbl"

GLOBAL_FUNCTION void
init_gen()
{
  declaration_t *I;
  char outfile[FILENAME_MAX];
  FILE *include;

  /* Convert filename to the name of the output file. */
  strcpy(outfile,BaseName);
  strcat(outfile,"init");
  strcat(outfile,".c");

  /* Open the output file. */
  output = fopen(outfile,"w");
  if (!output) {
    fprintf(stderr,"%s: couldn't open the output file: %s\n",progname,outfile);
    error(NULL);
    }

  /* Convert filename to the name of an include file. */
  strcpy(outfile,BaseName);
  strcat(outfile,"init");
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
  fprintf(output,"#include \"%sinit.h\"\n",BaseName);

  /* Go thru the list of declarations and generate the ip functions
   * for the data. */
  for (I=dl; I!=NULL; I=I->p) {
    /* Write to the C include file. */
    fprintf(include,"#if defined(NO_PROTO)\n");
    fprintf(include,"  void init_%s();\n",I->name);
    fprintf(include,"#else\n");
    fprintf(include,"  void init_%s(%s_t *%s);\n",
            I->name,I->name,I->name);
    fprintf(include,"#endif\n");

    if (is_excluded("init",I)) continue;
    /* Write to the C source file. */
    init_declaration(I);
    }

  fclose(output);
  }

LOCAL_FUNCTION void
init_declaration(dec)
declaration_t *dec;
{
  member_list_t *J;

  fprintf(output,"\n/* Generated init function for %s: */\n",dec->name);
  fprintf(output,"void\n");
  fprintf(output,"init_%s(%s)\n",dec->name,dec->name);
  fprintf(output,"%s_t *%s;\n",dec->name,dec->name);
  fprintf(output,"{\n");

  /* Loop over the members of the structure. */
  for (J=dec->members; J!=NULL; J=J->p) {
    init_member(J->member,dec->name);
    }

  fprintf(output,"  }\n");
  }

LOCAL_FUNCTION void
init_member(member,structname)
member_t *member;
char *structname;
{
  index_list_t *I;
  int n_indices;

  /* Cannot initialize unions. */
  if (member->uname) return;

  /* Compute the number of indices on the member. */
  for ((I=member->indices),(n_indices=0);  I!=NULL; (I=I->p),(n_indices++));

  /* If this is a pointer make it NULL. */
  if (member->pointer || n_indices || member->indices) {
    fprintf(output,"  %s->%s = NULL;\n",structname,member->name);
    }
  /* If this is a basic data type, make it 0. */
  else if (basic_type(member->type)) {
    fprintf(output,"  %s->%s = 0;\n",structname,member->name);
    }
  /* Otherwise call the appropiate initialization function. */
  else {
    fprintf(output,"  init_%s(&(%s->%s));\n",
            member->type,structname,member->name);
    }

  }

