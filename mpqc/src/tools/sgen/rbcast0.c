
/* $Log$
 * Revision 1.3  1996/03/23 02:39:03  cljanss
 * Everything can now be configured with autoconf.
 *
 * Revision 1.2  1995/03/18 00:09:57  cljanss
 * Using util/group to provide picl support.  Deleted the comm directory.
 *
 * Revision 1.1.1.1  1993/12/29  12:53:58  etseidl
 * SC source tree 0.1
 *
 * Revision 1.4  1992/07/20  18:37:38  seidl
 * add support for string arrays
 *
 * Revision 1.3  1992/06/17  23:07:35  jannsen
 * modified to generate clean code
 *
 * Revision 1.2  1992/03/30  23:08:59  seidl
 * merge in sandia changes
 *
 * Revision 1.4  91/09/28  18:13:39  cljanss
 * Output files now have shorter names.  Input file rootnames
 * of 8 or fewer characters will have output filenames of 14
 * or fewer characters.
 * 
 * Revision 1.3  91/09/28  16:40:35  cljanss
 * new naming convention is used to keep names <= 14 characters.
 * 
 * Revision 1.2  91/08/08  22:22:22  cljanss
 * updated for new argument lists in the general read and write routines
 * 
 * Revision 1.1  1991/07/19  19:06:56  cljanss
 * Initial revision
 * */

static char *rcsid = "$Id:";

#include <stdio.h>

#include "rbcast0.gbl"
#include "rbcast0.lcl"

#include "gen_read.gbl"

GLOBAL_FUNCTION void
rbcast0_gen()
{
  char *protoargsfmt = "%s_t *_%s,int _type,int _root";
  char *funcargsfmt  = "_%s,_type,_root";
  char *funcdecsfmt  = "%s_t *_%s;\nint _type;\nint _root;\n";
  char *elemargsfmt  = "%s(%s%s),_type,_root";
  char *tpointerargs = "_type,_root";
  char *strpointerargs  = "%s&(%s%s),_type,_root";
  int useoffsets = 0;

  general_read_gen("rbc0","rbcast0","rbcast0",
                 protoargsfmt,funcargsfmt,funcdecsfmt,elemargsfmt,tpointerargs,
                 strpointerargs,useoffsets,"#include <util/group/picl.h>\n");
 }
