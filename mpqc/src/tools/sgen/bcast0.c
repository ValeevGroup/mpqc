
/* $Log$
 * Revision 1.3  1995/03/18 00:09:56  cljanss
 * Using util/group to provide picl support.  Deleted the comm directory.
 *
 * Revision 1.2  1994/10/18  23:03:43  etseidl
 * fix many warnings, use memset rather than bzero
 *
 * Revision 1.1.1.1  1993/12/29  12:53:57  etseidl
 * SC source tree 0.1
 *
 * Revision 1.3  1992/06/17  23:07:12  jannsen
 * modified to generate clean code
 *
 * Revision 1.2  1992/03/30  23:01:59  seidl
 * merge in sandia changes
 *
 * Revision 1.3  91/09/28  18:13:30  cljanss
 * Output files now have shorter names.  Input file rootnames
 * of 8 or fewer characters will have output filenames of 14
 * or fewer characters.
 * 
 * Revision 1.2  91/09/28  16:40:22  cljanss
 * new naming convention is used to keep names <= 14 characters.
 * 
 * Revision 1.1  91/07/19  19:08:06  cljanss
 * Initial revision
 *  */

#include <stdio.h>
#include <string.h>
#include <tmpl.h>
#include "types.h"
#include "global.h"

#include "bcast0.gbl"
#include "bcast0.lcl"

#include "sgen_util.gbl"
#include "error.gbl"

GLOBAL_FUNCTION VOID
bcast0_gen()
{
  declaration_t *I;
  char outfile[FILENAME_MAX];
  FILE *include;

  /* Convert filename to the name of the output file. */
  strcpy(outfile,BaseName);
  strcat(outfile,"bc0");
  strcat(outfile,".c");

  /* Open the output file. */
  output = fopen(outfile,"w");
  if (!output) {
    fprintf(stderr,"%s: couldn't open the output file: %s\n",progname,outfile);
    error(NULL);
    }

  /* Convert filename to the name of an include file. */
  strcpy(outfile,BaseName);
  strcat(outfile,"bc0");
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
  fprintf(output,"#include <util/group/picl.h>\n");
  fprintf(output,"#include \"%s.h\"\n",BaseName);
  fprintf(output,"#include \"%sbc0.h\"\n",BaseName);
  fprintf(output,"#include \"%ssbc0.h\"\n",BaseName);
  fprintf(output,"#include \"%srbc0.h\"\n",BaseName);

  /* Go thru the list of declarations and generate the ip functions
   * for the data. */
  for (I=dl; I!=NULL; I=I->p) {
    /* Write to the C include file. */
    fprintf(include,"#if defined(NO_PROTO)\n");
    fprintf(include,"  void bcast0_%s();\n",I->name);
    fprintf(include,"#else\n");
    fprintf(include,"  void bcast0_%s(%s_t *%s,int type,int root);\n",
            I->name,I->name,I->name);
    fprintf(include,"#endif\n");

    if (is_excluded("bcast0",I)) continue;
    /* Write to the C source file. */
    bcast0_declaration(I);
    }

  fclose(output);
  }

LOCAL_FUNCTION VOID
bcast0_declaration(dec)
declaration_t *dec;
{

  fprintf(output,"\n/* Generated bcast0 function for %s: */\n",dec->name);
  fprintf(output,"void\n");
  fprintf(output,"bcast0_%s(%s,type,root)\n",dec->name,dec->name);
  fprintf(output,"%s_t *%s;\n",dec->name,dec->name);
  fprintf(output,"int type;\n");
  fprintf(output,"int root;\n");
  fprintf(output,"{\n");
  fprintf(output,"  int numproc,me,host;\n");

  fprintf(output,"  who0(&numproc,&me,&host);\n");

  fprintf(output,"  if (root==me) sbcast0_%s(%s,type,root);\n",
          dec->name,dec->name);
  fprintf(output,"  else          rbcast0_%s(%s,type,root);\n",
          dec->name,dec->name);

  fprintf(output,"  }\n");
  }

