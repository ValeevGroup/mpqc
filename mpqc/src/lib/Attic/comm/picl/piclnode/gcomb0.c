/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *  PICL source code                                               *
 *                                                                 *
 *  We welcome questions, comments, and bug reports, and request   *
 *  that you share any modifications with us.                      *
 *                                                                 *
 *  Patrick Worley                                                 *
 *  Oak Ridge National Laboratory                                  *
 *  worley@msr.epm.ornl.gov                                        *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "portglobals.h"

gcomb0(buf, items, datatype, msgtype, root, comb)
char *buf ;
int items, datatype, msgtype, root;
void (*comb)() ;

/*
 *  Componentwise combination of a vector over all processors,
 *  using given topology.
 *
 *  buf        array of data to be combined.
 *  items      number of items in array buf.
 *  datatype   code number for type of data:
 *  datatype = 0 : char
 *             1 : short   (int*2 in Fortran)
 *             2 : int
 *             3 : long    (int*4 in Fortran)
 *             4 : float   (real*4 in Fortran)
 *             5 : double  (real*8 in Fortran)
 *  msgtype    user-defined id to distinguish messages.
 *  root       processor in which final result will reside.
 *  comb       name of user-supplied function to be applied in
 *             combining data.  The operation defined by comb
 *             must be associative and commutative in order for
 *             result to be well defined independent of topology
 *             and order.  Examples are max, min, +, *, &, |, ^.
 *
 *  Caution: buf is overwritten.
 */
{
  int p, me, host, top, ord, dir, root0, bytes, datasize0() ;

  traceblockbegin(GCOMB, msgtype, root);

  /**** calculate the gcomb ****/
  who0(&p, &me, &host) ;
  getarc0(&p, &top, &ord, &dir) ;

  root0 = (root == host) ? 0 : root ;

  if (p > 1) switch (top) {
    case 1: gcomb_cube(buf, items, datatype, msgtype, root0, me, p, comb) ;
            break ;
    case 2: gcomb_full(buf, items, datatype, msgtype, root0, me, p, comb) ;
            break ;
    case 3: gcomb_ring(buf, items, datatype, msgtype, root0, me, p, dir, comb) ;
            break ;
    case 4: gcomb_bidi(buf, items, datatype, msgtype, root0, me, p, comb) ;
            break ;
    };

  if ((root == host) && (me == 0)) {
    bytes = items*datasize0(datatype) ;
    send0(buf, bytes, msgtype, host) ;
    };
  /*******************************/

  traceblockend(GCOMB, msgtype, root);

}

gcomb_cube ( buf, items, datatype, msgtype, root, myname, p, comb )
char *buf ;
int items, datatype, msgtype, root, myname, p;
void (*comb)() ;
/*
 *  Componentwise combination of a vector over all processors,
 *  using spanning tree.  Caution: buf is overwritten.
 */
{
	char *malloc0() ;
	int bytes, i, j, me, newroot, datasize0() ;
	char *t ;

	bytes = items*datasize0(datatype) ;
	t = (char *)malloc0(bytes) ;

	newroot = root < envHLF ? root : root - envHLF ;
	if (myname < envHLF) {
                if (envTWN < p && envTWN != root) {
			recv0(t, bytes, msgtype) ;
			(*comb)(buf, t, items, datatype) ;
		}
                me = myname^newroot ;
                for (i = envHLF/2 ; i > me ; i /= 2) {
			recv0(t, bytes, msgtype) ;
			(*comb)(buf, t, items, datatype) ;
		}
		if (myname != newroot)send0(buf,bytes,msgtype,(me-i)^newroot);
		else if (myname != root) send0(buf, bytes, msgtype, root) ; 
	}
	else {
		if (myname == root) {
			recv0(t, bytes, msgtype) ;
			(*comb)(buf, t, items, datatype) ;
		}
		else send0(buf, bytes, msgtype, envTWN) ;
	}
	free (t) ;
}

gcomb_full ( buf, items, datatype, msgtype, root, myname, p, comb )
char *buf ;
int items, datatype, msgtype, root, myname, p;
void (*comb)() ;
/*
 *  Componentwise combination of a vector over all processors,
 *  using full interconnections.  Caution: buf is overwritten.
 */
{
	char *malloc0() ;
	int i, k, bytes, datasize0() ;
	char *t ;

	bytes = items*datasize0(datatype) ;
	t = (char *)malloc0(bytes) ;

	if (root != myname) send0(buf, bytes, msgtype, root) ;
	else for ( i = 1 ; i < p ; i++ ) {
		recv0(t, bytes, msgtype) ;
		(*comb)(buf, t, items, datatype) ;
	}
	free (t) ;
}

gcomb_ring ( buf, items, datatype, msgtype, root, myname, p, dir, comb )
char *buf ;
int items, datatype, msgtype, root, myname, p, dir;
void (*comb)() ;
/*
 *  Componentwise combination of a vector over all processors,
 *  using unidirectional ring.  Caution: buf is overwritten.
 */
{
	char *malloc0() ;
	int neighbor0() ;
	int forward, back, i, bytes, datasize0() ;
	char *t ;

	bytes = items*datasize0(datatype) ;
	t = (char *)malloc0(bytes) ;

	forward = neighbor0(1, myname) ;
	back = neighbor0(-1, myname) ;

	if (dir == 1) {
		if (root != back) {
			recv0(t, bytes, msgtype) ;
			(*comb)(buf, t, items, datatype) ;
		}
		if (myname != root) send0(buf, bytes, msgtype, forward) ;
	}
	else {
		if (root != forward) {
			recv0(t, bytes, msgtype) ;
			(*comb)(buf, t, items, datatype) ;
		}
		if (myname != root) send0(buf, bytes, msgtype, back) ;
	}
	free (t) ;
}

gcomb_bidi ( buf, items, datatype, msgtype, root, myname, p, comb )
char *buf ;
int items, datatype, msgtype, root, myname, p;
void (*comb)() ;
/*
 *  Componentwise combination of a vector over all processors,
 *  using a bidirectional ring.  Caution: buf is overwritten.
 */
{
	char *malloc0() ;
	int neighbor0(), antipode0() ;
	int forward, back, antipode, node, bytes, i, datasize0() ;
	char *t ;

	bytes = items*datasize0(datatype) ;
	t = (char *)malloc0(bytes) ;

	forward = neighbor0(1, myname) ;
	back = neighbor0(-1, myname) ;
	antipode = antipode0(root) ;

	if (myname == root) {
		recv0(t, bytes, msgtype) ;
		(*comb)(buf, t, items, datatype) ;
		if (p > 2) {
			recv0(t, bytes, msgtype) ;
			(*comb)(buf, t, items, datatype) ;
		}
	}
	else if (myname == antipode) {
		send0(buf, bytes, msgtype, forward) ;
	}
	else if (forward == antipode && p > 2) {
		send0(buf, bytes, msgtype, back) ;
	}
	else {
		recv0(t, bytes, msgtype) ;
		recvinfo0(&bytes, &msgtype, &node) ;
		(*comb)(buf, t, items, datatype) ;
		if (node == back) {
			send0(buf, bytes, msgtype, forward) ;
		}
		else {
			send0(buf, bytes, msgtype, back) ;
		}
	}
	free (t) ;
}

