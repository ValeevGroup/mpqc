
/* These utility routines assist in reading in the basis functions. */

/* $Log$
 * Revision 1.5  1994/08/26 22:45:09  etseidl
 * fix a bunch of warnings, get rid of rcs id's, get rid of bread/bwrite and
 * fread/fwrite modules
 *
 * Revision 1.4  1994/08/24  16:05:25  etseidl
 * get rid of ip functions, they have been replaced by keyval equivalents
 *
 * Revision 1.3  1994/05/27  23:44:12  cljanss
 * Changed some char* to const char*.  Included ipv2 interface from
 * keyval/ipv2c.h.  Fixed a bug with basis->n not being initialized.
 *
 * Revision 1.2  1993/12/30  13:32:43  etseidl
 * mostly rcs id stuff
 *
 * Revision 1.4  1993/04/28  00:30:13  jannsen
 * added int_read_basis global function
 *
 * Revision 1.3  1992/06/17  22:04:23  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.2  1992/03/31  01:21:20  jannsen
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.1  1991/12/14  00:19:36  cljanss
 * Initial revision
 * */
/*
 * Revision 1.6  91/10/31  14:32:24  cljanss
 * The basis set name is kept in the basis structure.
 * 
 * Revision 1.5  91/09/28  19:26:42  cljanss
 * Switch to new file naming scheme
 * 
 * Revision 1.4  91/09/10  19:32:34  cljanss
 * If the basis name is "nothing", an empty basis set
 * will be returned.
 * 
 * Revision 1.3  1991/08/09  16:41:23  cljanss
 * fixed basis set reader, it returned incorrect return codes
 *
 * Revision 1.2  1991/07/16  17:55:27  cljanss
 * slight change to make a character initialization compile on the SGI
 *
 * Revision 1.1  1991/06/16  16:40:07  janssen
 * Initial revision
 * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <util/sgen/sgen.h>

#include "atoms.h"

/* This will print out a shell type.  The integer am is converted to
 * a character and if puream is nonnull, then the number of functions
 * in the shell is appended to the name.
 */
void
print_shell_type(FILE *fp, shell_type_t *shell_type)
{
  char *amnames = "spdfghijklmnoqrtuvwxyz";

  fprintf(fp," am = %c",amnames[shell_type->am]);
  if (shell_type->puream) fprintf(fp,"%d",2*shell_type->am + 1);

  sgen_print_suppress_indent();
}
