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

bcast0(buf, bytes, type, root)
char *buf ;
int bytes, type, root ;

/*
 *  Broadcast message buf, of given length and type, to
 *  all processors, using given interconnection topology.
 *  This version meant for synchronous applications in
 *  which all processors communicate at about the same time.
 */
{
  int p, me, host, top, ord, dir, node ;

  traceblockbegin(BCAST0, type, root);

  /**** calculate the bcast ****/
  who0(&p, &me, &host) ;
  getarc0(&p, &top, &ord, &dir) ;

  if (root == host) {
    if (me == 0) {
      recv0(buf, bytes, type) ;
      recvinfo0(&bytes, &type, &node) ;
      };
    root = 0 ;
    }

  if (p > 1) switch (top) {
    case 1: bcast_cube(buf, bytes, type, root, me, p) ;
            break ;
    case 2: bcast_full(buf, bytes, type, root, me, p, dir) ;
            break ;
    case 3: bcast_ring(buf, bytes, type, root, me, p, dir) ;
            break ;
    case 4: bcast_bidi(buf, bytes, type, root, me, p, dir) ;
            break ;
    };
  /*******************************/

  traceblockend(BCAST0, type, root);

}

bcast_cube ( buf, bytes, type, root, myname, p )
	char *buf ;
	int bytes, type, root, myname, p ;
/*
 *  Broadcast message buf of length bytes to all processors
 *  using a minimum spanning tree with given root.
 */
{
	int i, me, newroot, node ;

	newroot = root < envHLF ? root : root - envHLF ;
	if (myname != root) {
		recv0(buf, bytes, type) ;
		recvinfo0(&bytes, &type, &node) ;
	}
	else if (myname != newroot) send0(buf, bytes, type, newroot) ;
	if (myname < envHLF) {
		me = myname^newroot ;
		for ( i = 1 ; i < envHLF ; i *= 2 )
			if (i > me){
                                   send0(buf, bytes, type, (me+i)^newroot) ;
};
		if (envTWN < p && envTWN != root)
			send0(buf, bytes, type, envTWN) ;
	}
}

bcast_full ( buf, bytes, type, root, myname, p, dir )
	char *buf ;
	int bytes, type, root, myname, p, dir ;
/*
 *  Broadcast message buf of length bytes to all processors
 *  using full interconnections (direct sends).
 */
{
	int neighbor0() ;
	int node, i ;

	if (myname != root) {
		recv0(buf, bytes, type) ;
		recvinfo0(&bytes, &type, &node) ;
	}
	else for ( i = 1 ; i < p ; i++ )
		send0(buf, bytes, type, neighbor0(i*dir, myname)) ;
}

bcast_ring ( buf, bytes, type, root, myname, p, dir )
	char *buf ;
	int bytes, type, root, myname, p, dir ;
/*
 *  Broadcast message buf of length bytes to all processors
 *  using a unidirectional ring.
 */
{
	int neighbor0() ;
	int forward, back, node ;

	forward = neighbor0(1, myname) ;
	back = neighbor0(-1, myname) ;
	if (root != myname) {
		recv0(buf, bytes, type) ;
		recvinfo0(&bytes, &type, &node) ;
	}
	if (dir == 1) {
		if (forward != root) send0(buf, bytes, type, forward) ;
	}
	else {
		if (back != root) send0(buf, bytes, type, back) ;
	}
}

bcast_bidi ( buf, bytes, type, root, myname, p, dir )
	char *buf ;
	int bytes, type, root, myname, p, dir ;
/*
 *  Broadcast message buf of length bytes to all processors
 *  using a bidirectional ring.
 */
{
	int neighbor0(), antipode0() ;
	int forward, back, antipode, node ;

	forward = neighbor0(1, myname) ;
	back = neighbor0(-1, myname) ;

	if (root == myname) {
		if (dir == 1) {
			send0(buf, bytes, type, forward) ;
			if (p > 2) send0(buf, bytes, type, back) ;
		}
		else {
			send0(buf, bytes, type, back) ;
			if (p > 2) send0(buf, bytes, type, forward) ;
		}
	}
	else {
		recv0(buf, bytes, type) ;
		recvinfo0(&bytes, &type, &node) ;
		antipode = antipode0(root) ;
		if (myname != antipode) {
			if (node == forward) {
				send0(buf, bytes, type, back);
			}
			else if (forward != antipode) {
				send0(buf, bytes, type, forward);
			}
		}
	}
}

