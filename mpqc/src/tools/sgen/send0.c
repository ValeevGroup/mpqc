
/* $Log$
 * Revision 1.2  1995/03/18 00:09:59  cljanss
 * Using util/group to provide picl support.  Deleted the comm directory.
 *
 * Revision 1.1.1.1  1993/12/29  12:53:58  etseidl
 * SC source tree 0.1
 *
 * Revision 1.3  1992/07/20  18:37:44  seidl
 * add support for string arrays
 *
 * Revision 1.2  1992/06/17  23:07:45  jannsen
 * modified to generate clean code
 *
 * Revision 1.1.1.1  1992/03/17  18:11:00  seidl
 * Struct GENerator 2.0
 *
 * Revision 1.1  1992/03/17  18:11:00  seidl
 * Initial revision
 *
 * Revision 1.1  1992/01/09  12:15:34  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#include <stdio.h>
#include <tmpl.h>

#include "send0.gbl"
#include "send0.lcl"

#include "gen_write.gbl"

GLOBAL_FUNCTION VOID
send0_gen()
{
  char *protoargsfmt = "%s_t *_%s,int _type,int _dest";
  char *funcargsfmt  = "_%s,_type,_dest";
  char *funcdecsfmt  = "%s_t *_%s;\nint _type;\nint _dest;\n";
  char *elemargsfmt  = "%s(%s%s),_type,_dest";
  char *selemargsfmt  = "%s&(%s%s),_type,_dest";
  int useoffsets = 0;

  general_write_gen("snd0","send0","send0",
                 protoargsfmt,funcargsfmt,funcdecsfmt,elemargsfmt,
                 selemargsfmt,useoffsets,"#include <util/group/picl.h>\n");
 }
