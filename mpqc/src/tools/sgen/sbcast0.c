
/* $Log$
 * Revision 1.1  1993/12/29 12:53:57  etseidl
 * Initial revision
 *
 * Revision 1.4  1992/07/20  18:37:45  seidl
 * add support for string arrays
 *
 * Revision 1.3  1992/06/17  23:07:43  jannsen
 * modified to generate clean code
 *
 * Revision 1.2  1992/03/30  23:09:15  seidl
 * merge in sandia changes
 *
 * Revision 1.4  91/09/28  18:13:40  cljanss
 * Output files now have shorter names.  Input file rootnames
 * of 8 or fewer characters will have output filenames of 14
 * or fewer characters.
 * 
 * Revision 1.3  91/09/28  16:40:36  cljanss
 * new naming convention is used to keep names <= 14 characters.
 * 
 * Revision 1.2  91/08/08  22:22:22  cljanss
 * updated for new argument lists in the general read and write routines
 * 
 * Revision 1.1  1991/07/19  19:08:06  cljanss
 * Initial revision
 * */

static char *rcsid = "$Id:";

#include <stdio.h>
#include <tmpl.h>

#include "sbcast0.gbl"
#include "sbcast0.lcl"

#include "gen_read.gbl"

GLOBAL_FUNCTION VOID
sbcast0_gen()
{
  char *protoargsfmt = "%s_t *_%s,int _type,int _root";
  char *funcargsfmt  = "_%s,_type,_root";
  char *funcdecsfmt  = "%s_t *_%s;\nint _type;\nint _root;\n";
  char *elemargsfmt  = "%s(%s%s),_type,_root";
  char *selemargsfmt  = "%s&(%s%s),_type,_root";
  int useoffsets = 0;

  general_write_gen("sbc0","sbcast0","sbcast0",
                 protoargsfmt,funcargsfmt,funcdecsfmt,elemargsfmt,
                 selemargsfmt,useoffsets,"#include <comm/picl/picl.h>\n");
 }
