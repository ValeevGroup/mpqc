
/* $Log$
 * Revision 1.1  1993/12/29 12:53:56  etseidl
 * Initial revision
 *
 * Revision 1.4  1992/07/20  18:37:41  seidl
 * add support for string arrays
 *
 * Revision 1.3  1992/06/17  23:07:16  jannsen
 * modified to generate clean code
 *
 * Revision 1.2  1992/03/30  23:02:29  seidl
 * merge in sandia changes
 *
 * Revision 1.9  91/09/28  18:13:33  cljanss
 * Output files now have shorter names.  Input file rootnames
 * of 8 or fewer characters will have output filenames of 14
 * or fewer characters.
 * 
 * Revision 1.8  91/09/28  16:40:25  cljanss
 * new naming convention is used to keep names <= 14 characters.
 * 
 * Revision 1.7  91/08/08  22:22:22  cljanss
 * updated for new argument lists in the general read and write routines
 * 
 * Revision 1.6  1991/07/19  17:51:50  cljanss
 * The initial general read and write routines have been put in.
 * The contents of previous versions of this file are now in
 * general_{read,write}.c
 *
 * Revision 1.5  1991/06/22  06:08:58  seidl
 * write out pointers, no longer exit if a pointer is null, just skip
 * writing that element
 *
 * Revision 1.4  1991/06/17  18:33:24  seidl
 * fix how things like double_matrix are handled
 *
 * Revision 1.3  1991/06/17  17:24:22  seidl
 * add ability to recursively write structs of structs
 *
 * Revision 1.2  1991/06/16  20:31:59  seidl
 * first working version
 *
 * Revision 1.1  1991/06/16  16:07:13  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#include <stdio.h>
#include <tmpl.h>

#include "bwrite.gbl"
#include "bwrite.lcl"

#include "gen_write.gbl"

GLOBAL_FUNCTION VOID
bwrite_gen()
{
  char *protoargsfmt = "int _unit,%s_t *_%s,int *_offset";
  char *funcargsfmt  = "_unit,_%s,_offset";
  char *funcdecsfmt  = "int _unit;\n%s_t *_%s;\nint *_offset;\n";
  char *elemargsfmt  = "_unit,%s(%s%s),_offset";
  char *selemargsfmt  = "_unit,%s&(%s%s),_offset";
  int useoffsets = 1;

  general_write_gen("bwr","bwrite","bwrite",
                 protoargsfmt,funcargsfmt,funcdecsfmt,elemargsfmt,
                 selemargsfmt,useoffsets,NULL);
 }

