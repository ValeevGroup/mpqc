
/* This file contains a global double sum routine for use only on the iPSC */

/* $Log$
 * Revision 1.2  1994/08/26 23:54:17  etseidl
 * remove rcs ids and fix a warning
 *
 * Revision 1.1.1.1  1993/12/29  12:53:47  etseidl
 * SC source tree 0.1
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

#include <comm/picl/picl.h>
#include "piclext.h"

#include <cube.h>

/* the following is a fast global sum for long vectors.
 * It was written by Stan Erwin, an Intel employee at the NIH.
 *
 * FORCE_TYPE is defined in <cube.h>
 *
 * modifications of van de Geijn's gdcomp to
 * 1. vectorize the add (via a daxpy call)
 * 2. change the three trip protocol to a 2 trip protocal using
 *    force types and control messages
 */

void
gdcomb(n, x, y, dim,idim)
int n;
double *x;
double *y;
int dim;
int idim;
{
  int i,ibit, l1, l2 , me, tempdim,msgid,msgid2,dummy,ione;
  double done;
  ione = 1;
  done = 1.0;
  me = mynode();

  if (dim == idim) return;

  l1 = n/2;
  l2 = n-l1;
  ibit = 1<<dim;
  if ((me&ibit) == 0) {
    msgid = _irecv(FORCE_TYPE + me^ibit,&y[l1], l2*sizeof(double));
    _csend(me,dummy, 0, me^ibit,0);
    _crecv(me^ibit,dummy, 0);
    msgid2 = _isend(FORCE_TYPE + me,x, l1*sizeof(double), me^ibit,0);
    _msgwait(msgid);

    daxpy_(&l2, &done, &(y[l1]), &ione, &(x[l1]), &ione);

    tempdim = dim + 1;
    _msgwait(msgid2);

    gdcomb(l2, &x[l1], &y[l1],tempdim,idim);

    msgid = _irecv(FORCE_TYPE + me^ibit,x, l1*sizeof(double));
    _csend(me,dummy, 0, me^ibit,0);
    _crecv(me^ibit,dummy, 0);
    _csend(FORCE_TYPE + me,&x[l1], l2*sizeof(double), me^ibit,0);
    _msgwait(msgid);
    }
  else {
    msgid = _irecv(FORCE_TYPE + me^ibit,y,l1*sizeof(double));
    _csend(me,dummy, 0, me^ibit,0);
    _crecv(me^ibit,dummy, 0);
    msgid2 =_isend(FORCE_TYPE + me,&x[l1], l2*sizeof(double), me^ibit,0);
    _msgwait(msgid);

    daxpy_(&l1, &done, y, &ione, x, &ione);

    tempdim = dim + 1;
    _msgwait(msgid2);

    gdcomb(l1, x, y,tempdim,idim);

    msgid = _irecv(FORCE_TYPE + me^ibit,&x[l1], l2*sizeof(double));
    _csend(me,dummy, 0, me^ibit,0);
    _crecv(me^ibit,dummy, 0);
    _csend(FORCE_TYPE + me,x, l1*sizeof(double), me^ibit,0);
    _msgwait(msgid);
    }
  }

