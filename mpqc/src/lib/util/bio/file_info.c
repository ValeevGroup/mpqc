
/* $Log$
 * Revision 1.2  1994/06/08 01:17:11  cljanss
 * Many changes.  These include: newmat7 and nihmatrix -> scmat
 * and mpqcic -> MPSCF and updated optimize stuff.
 *
 * Revision 1.1.1.1  1993/12/29  12:53:41  etseidl
 * SC source tree 0.1
 *
 * Revision 1.2  1992/03/30  22:34:23  jannsen
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.3  91/12/02  17:19:28  cljanss
 * clean up use of VOID and VOID_PTR for old compilers
 * 
 * Revision 1.2  91/09/28  21:05:32  cljanss
 * switch to new naming convention
 * 
 * Revision 1.1  91/06/16  20:28:41  seidl
 * Initial revision
 *  */

static char *rcsid = "$Id$";

#include <stdio.h>
#include <tmpl.h>
#include <util/keyval/ipv2c.h>
#include "param.h"

#include "file_info.gbl"
#include "file_info.lcl"

GLOBAL_FUNCTION int
get_file_info(name,token,format,val)
char *name;
char *token;
char *format;
VOID_PTR val;
{
  int errcod;
  char ip_token[MAX_STRING];

  sprintf(ip_token,":files:%s:%s",name,token);
  errcod = ip_data(ip_token,format,val,0);
  if(errcod == IPE_OK) return(0);

  sprintf(ip_token,":files:default:%s",token);
  errcod = ip_data(ip_token,format,val,0);
  if(errcod == IPE_OK) return(0);

  sprintf(ip_token,":default:files:%s:%s",name,token);
  errcod = ip_data(ip_token,format,val,0);
  if(errcod == IPE_OK) return(0);

  sprintf(ip_token,":default:files:default:%s",token);
  errcod = ip_data(ip_token,format,val,0);
  if(errcod == IPE_OK) return(0);

  return(-1);
  }
