
/* $Log$
 * Revision 1.1  1993/12/29 12:53:56  etseidl
 * Initial revision
 *
 * Revision 1.4  1992/07/20  18:37:40  seidl
 * add support for string arrays
 *
 * Revision 1.3  1992/06/17  23:07:14  jannsen
 * modified to generate clean code
 *
 * Revision 1.2  1992/03/30  23:02:16  seidl
 * merge in sandia changes
 *
 * Revision 1.7  91/09/28  18:13:32  cljanss
 * Output files now have shorter names.  Input file rootnames
 * of 8 or fewer characters will have output filenames of 14
 * or fewer characters.
 * 
 * Revision 1.6  91/09/28  16:40:23  cljanss
 * new naming convention is used to keep names <= 14 characters.
 * 
 * Revision 1.5  91/08/08  22:22:22  cljanss
 * updated for new argument lists in the general read and write routines
 * 
 * Revision 1.4  1991/07/19  17:51:50  cljanss
 * The initial general read and write routines have been put in.
 * The contents of previous versions of this file are now in
 * general_{read,write}.c
 *
 * Revision 1.3  1991/06/20  16:33:45  seidl
 * reads in all pointers. if a pointer is non-null then what it points
 * to is read in.  allocates storage for the buffer you're reading
 * into if it is null
 *
 * Revision 1.2  1991/06/17  18:33:04  seidl
 * fix how things like double_matrix are handled
 *
 * Revision 1.1  1991/06/17  17:24:15  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#include <stdio.h>
#include <tmpl.h>

#include "bread.gbl"
#include "bread.lcl"

#include "gen_read.gbl"

GLOBAL_FUNCTION VOID
bread_gen()
{
  char *protoargsfmt = "int _unit,%s_t *_%s,int *_offset";
  char *funcargsfmt  = "_unit,_%s,_offset";
  char *funcdecsfmt  = "int _unit;\n%s_t *_%s;\nint *_offset;\n";
  char *elemargsfmt  = "_unit,%s(%s%s),_offset";
  char *tpointerargs = "_unit,_offset";
  char *strpointerargs  = "_unit,%s&(%s%s),_offset";
  int useoffsets = 1;

  general_read_gen("brd","bread","bread",
                 protoargsfmt,funcargsfmt,funcdecsfmt,elemargsfmt,tpointerargs,
                 strpointerargs,useoffsets,NULL);
 }

