
/* $Log$
 * Revision 1.4  1995/03/18 00:11:27  cljanss
 * Using util/group to provide picl support.  Deleted the comm directory.
 *
 * Revision 1.3  1995/03/17  01:51:38  cljanss
 * Removed -I. and -I$(SRCDIR) from the default include path in
 * GlobalMakefile to avoid name conflicts with system include files.
 * Modified files under src.lib to include all files relative to src.lib.
 * Makefiles under src.bin need to add the -I. and -I$(SRCDIR) back onto
 * INCLUDE and CXXINCLUDE or make other arrangements.
 *
 * Revision 1.2  1994/08/25  22:48:11  etseidl
 * remove rcsids and fix some warnings
 *
 * Revision 1.1.1.1  1993/12/29  12:53:40  etseidl
 * SC source tree 0.1
 *
 * Revision 1.2  1992/06/17  22:16:44  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.1.1.1  1992/03/17  17:09:34  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:09:33  seidl
 * Initial revision
 *
 * Revision 1.2  1992/01/09  12:20:18  seidl
 * add function sgen_reset_bcast0
 *
 * Revision 1.1  1991/12/20  16:10:32  seidl
 * Initial revision
 *
 * Revision 1.1  91/11/18  18:16:23  cljanss
 * Initial revision
 *  */

#include <stdio.h>
#include <util/group/picl.h>
#include <util/sgen/sgen.h>

#define ALLOC_GLOBALS
#include <util/sgen/bcast0.h>
#undef ALLOC_GLOBALS

/* This is used to set up the print options for the sgen interface
 * to the bcast0 routines.
 *  send = print level for sending
 *  receive = print level for receiving
 *  host = host for which printing is done (-1 = all)
 */
void
sgen_bcast0_set_print(send,receive,host)
int send;
int receive;
int host;
{
  sgen_bcast0_print_send = send;
  sgen_bcast0_print_receive = receive;
  sgen_bcast0_print_host = host;
  }

/* This prints out information when debugging is on. */
void
sgen_bcast_print(sorr,type,root,size)
int sorr;
int type;
int root;
int size;
{
  int numproc,me,host;
  char *endline;
  static char *val = ", value = ";
  static char *new = "\n";

  who0(&numproc,&me,&host);

  if ((me == sgen_bcast0_print_host) || (sgen_bcast0_print_host == -1)) {
    if ((sgen_bcast0_print_send & 2) && (sorr == 's')) endline = val;
    else if ((sgen_bcast0_print_receive & 2) && (sorr == 'r')) endline = val;
    else endline = new;
    printf("%cbcast0: me = %d, type = %d, root = %d, size = %d%s",
           sorr,me,type,root,size,endline);
    }
  }


void
bcast0_sync()
{
  int i;
  int numproc,me,host;
  int junk=0;

  who0(&numproc,&me,&host);

  /* printf("synchronizing ...\n"); */

  if (me == 0) {
    for (i=1; i<numproc; i++) {
      /* printf("sending sync_start to %d\n",i); */
      send0(&junk,sizeof(junk),SYNC_START,i);
      }
    }
  else {
    /* printf("receiving sync_start at %d\n",me); */
    recv0(&junk,sizeof(junk),SYNC_START);
    }

  if (me == 0) {
    for (i=1; i<numproc; i++) {
      /* printf("receiving sync_end\n"); */
      recv0(&junk,sizeof(junk),SYNC_END);
      }
    }
  else {
    /* printf("sending sync_end to %d\n",0); */
    send0(&junk,sizeof(junk),SYNC_END,0);
    }
  }

/* this routine will reset the sgen_bcast0_type counter
 * this is used to ensure that the host and node programs 
 * are in sync */
void
sgen_reset_bcast0()
{
  sgen_bcast0_type = MIN_TYPE;
  }

