
/* $Log$
 * Revision 1.2  1996/03/23 02:38:56  cljanss
 * Everything can now be configured with autoconf.
 *
 * Revision 1.1.1.1  1993/12/29 12:53:57  etseidl
 * SC source tree 0.1
 *
 * Revision 1.4  1992/07/20  18:37:42  seidl
 * add support for string arrays
 *
 * Revision 1.3  1992/06/17  23:07:21  jannsen
 * modified to generate clean code
 *
 * Revision 1.2  1992/03/30  23:05:45  seidl
 * merge in sandia changes
 *
 * Revision 1.1  91/11/18  17:24:22  cljanss
 * Initial revision
 * 
 * */

static char *rcsid = "$Id$";

#include <stdio.h>

#include "fwrite.gbl"
#include "fwrite.lcl"

#include "gen_write.gbl"

GLOBAL_FUNCTION void
fwrite_gen()
{
  char *protoargsfmt = "FILE *_fp,%s_t *_%s,int *_offset";
  char *funcargsfmt  = "_fp,_%s,_offset";
  char *funcdecsfmt  = "FILE *_fp;\n%s_t *_%s;\nint *_offset;\n";
  char *elemargsfmt  = "_fp,%s(%s%s),_offset";
  char *selemargsfmt  = "_fp,%s&(%s%s),_offset";
  int useoffsets = 1;

  general_write_gen("fwr","fwrite","fwrite",
                 protoargsfmt,funcargsfmt,funcdecsfmt,elemargsfmt,
                 selemargsfmt,useoffsets,NULL);
 }

