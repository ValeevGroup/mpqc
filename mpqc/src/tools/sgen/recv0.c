
/* $Log$
 * Revision 1.1  1993/12/29 12:53:57  etseidl
 * Initial revision
 *
 * Revision 1.3  1992/07/20  18:37:36  seidl
 * add support for string arrays
 *
 * Revision 1.2  1992/06/17  23:07:40  jannsen
 * modified to generate clean code
 *
 * Revision 1.1.1.1  1992/03/17  18:10:53  seidl
 * Struct GENerator 2.0
 *
 * Revision 1.1  1992/03/17  18:10:52  seidl
 * Initial revision
 *
 * Revision 1.1  1992/01/09  12:15:30  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#include <stdio.h>
#include <tmpl.h>

#include "recv0.gbl"
#include "recv0.lcl"

#include "gen_read.gbl"

GLOBAL_FUNCTION VOID
recv0_gen()
{
  char *protoargsfmt = "%s_t *_%s,int _type,int _from";
  char *funcargsfmt  = "_%s,_type,_from";
  char *funcdecsfmt  = "%s_t *_%s;\nint _type;\nint _from;\n";
  char *elemargsfmt  = "%s(%s%s),_type,_from";
  char *tpointerargs = "_type,_from";
  char *strpointerargs = "%s&(%s%s),_type,_from";
  int useoffsets = 0;

  general_read_gen("rcv0","recv0","recv0",
                 protoargsfmt,funcargsfmt,funcdecsfmt,elemargsfmt,tpointerargs,
                 strpointerargs,useoffsets,"#include <comm/picl/picl.h>\n");
 }
