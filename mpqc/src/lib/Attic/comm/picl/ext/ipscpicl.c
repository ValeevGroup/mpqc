/* This library includes IPSC to PICL conversions for a small
   subset of the IPSC library routines */
/*      All questions and criticisms should be directed to
        Michael Colvin  Center for Computational Engineering
        Sandia National Laboratory
        Livermore CA 94551                         */
/* $Log$
 * Revision 1.1  1993/12/29 12:53:47  etseidl
 * Initial revision
 *
 * Revision 1.6  1992/06/24  14:40:26  seidl
 * use intel functions rather than picl
 *
 * Revision 1.5  1992/06/23  19:59:47  seidl
 * add gdcomb and gdcomb_, fast global sum for the intel
 *
 * Revision 1.4  1992/06/17  22:25:28  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.3  1992/04/22  16:07:55  seidl
 * don't include mynode0() on sun
 *
 * Revision 1.2  1992/04/06  12:59:04  seidl
 * merge in sandia changes
 *
 * Revision 1.1.1.1  1992/03/17  17:12:16  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:12:15  seidl
 * Initial revision
 *
 * Revision 1.2  1992/02/20  12:07:29  seidl
 * *** empty log message ***
 *
 * Revision 1.1  1991/12/20  16:15:25  seidl
 * Initial revision
 *
 * Revision 1.2  1991/10/29  22:31:50  cljanss
 * Some of the routines were not needed on the NCUBE.
 *
 * Revision 1.1  91/08/21  00:17:47  cljanss
 * Initial revision
 * 
 * Revision 1.5  1991/07/17  20:44:40  cljanss
 * Checked in Mike Colvin's latest version after I made the copy.
 *
 * Revision 1.4  1991/03/11  22:37:05  colvin
 * Added routines without '0' appended to their names.
 *
 * Revision 1.3  1991/03/06  00:14:27  colvin
 * Fixed some nasty bugs in cubedim.
 *
 * Revision 1.2  91/03/05  09:42:23  colvin
 * Added infocount() [which may or may not be an actual IPSC routine,
 * but Bob Whiteside used it in his MPSCF program.
 * 
 * Revision 1.1  1991/03/04  23:47:29  colvin
 * Initial revision
 * */
static char rcsid[]="$Id$";

#include <comm/picl/picl.h>
#include "piclext.h"

#if defined(PARAGON)
#include <nx.h>
#elif defined(I860)
#include <cube.h>
#endif

#if !defined(I860)
int cubedim_() {
  int n,p,me, host,i=0;
  who0(&p,&me,&host);
  n=p;
  while (n > 0) {n >>= 1; i++;}
  return (i-1);
}
#ifndef NCUBE
int cubedim() {
  int n,p,me, host,i=0;
  who0(&p,&me,&host);
  n=p;
  while (n > 0) {n >>= 1; i++;}
  return (i-1);
}
int mynode() {
  int p, me, host;
  who0(&p, &me, &host);
  return (me);
}
int mynode_() {
  int p, me, host;
  who0(&p, &me, &host);
  return (me);
}
#endif /*NCUBE*/
int numnodes_() {
  int p, me, host;
  who0(&p, &me, &host);
  return (p);
}
int numnodes() {
  int p, me, host;
  who0(&p, &me, &host);
  return (p);
}

int infocount() {
  int bytes,type,source;
  recvinfo0(&bytes,&type,&source);
  return(bytes);
}
#endif /* !I860 */

/* Here are many of the above routines with 0 appended to their names 
   Bob Whiteside wrote his matrix lib using these names; I don't know why */
#if !defined(I860)
int mynode0() {
  int p, me, host;
  who0(&p, &me, &host);
  return (me);
}
int mynode0_() {
  int p, me, host;
  who0(&p, &me, &host);
  return (me);
}

int cubedim0_() {
  int n,p,me, host,i=0;
  who0(&p,&me,&host);
  n=p;
  while (n > 0) {n >>= 1; i++;}
  return (i-1);
}
int cubedim0() {
  int n,p,me, host,i=0;
  who0(&p,&me,&host);
  n=p;
  while (n > 0) {n >>= 1; i++;}
  return (i-1);
}
int numnodes0_() {
  int p, me, host;
  who0(&p, &me, &host);
  return (p);
}
int numnodes0() {
  int p, me, host;
  who0(&p, &me, &host);
  return (p);
}
#else /* I860 */
int mynode0()  { return(mynode()); }
int mynode0_() { return(mynode()); }
int cubedim0_() { return(nodedim()); }
int cubedim0()  { return(nodedim()); }
int numnodes0_() { return (numnodes()); }
int numnodes0()  { return (numnodes()); }
#endif /* !I860 */

