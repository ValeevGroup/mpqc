
/* $Log$
 * Revision 1.2  1996/03/23 02:38:55  cljanss
 * Everything can now be configured with autoconf.
 *
 * Revision 1.1.1.1  1993/12/29 12:53:57  etseidl
 * SC source tree 0.1
 *
 * Revision 1.4  1992/07/20  18:37:39  seidl
 * add support for string arrays
 *
 * Revision 1.3  1992/06/17  23:07:17  jannsen
 * modified to generate clean code
 *
 * Revision 1.2  1992/03/30  23:02:44  seidl
 * merge in sandia changes
 *
 * Revision 1.1  91/11/18  17:24:20  cljanss
 * Initial revision
 * 
 * */

static char *rcsid = "$Id$";

#include <stdio.h>

#include "fread.gbl"
#include "fread.lcl"

#include "gen_read.gbl"

GLOBAL_FUNCTION void
fread_gen()
{
  char *protoargsfmt = "FILE *_fp,%s_t *_%s, int *_offset";
  char *funcargsfmt  = "_fp,_%s,_offset";
  char *funcdecsfmt  = "FILE *_fp;\n%s_t *_%s;\nint *_offset;\n";
  char *elemargsfmt  = "_fp,%s(%s%s),_offset";
  char *tpointerargs = "_fp,_offset";
  char *strpointerargs  = "_fp,%s&(%s%s),_offset";
  int useoffsets = 1;

  general_read_gen("frd","fread","fread",
                 protoargsfmt,funcargsfmt,funcdecsfmt,elemargsfmt,tpointerargs,
                 strpointerargs,useoffsets,NULL);
 }

