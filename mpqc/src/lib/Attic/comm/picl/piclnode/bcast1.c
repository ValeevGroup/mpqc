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

bcast1(buf, bytes, type, root)
char *buf;
int bytes, type, root;

/*
 *  Broadcast message buf, of given length and type, to
 *  all processors, using given interconnection topology.
 *  This version is meant for asynchronous applications in
 *  which communication and computation are overlapped by
 *  pipelining.
 */
{
  int p, me, host, top, ord, dir, node ;

  traceblockbegin(BCAST1, type, root);

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
    case 1: relay_cube(buf, bytes, type, root, me, p) ;
            break ;
    case 2: relay_full(buf, bytes, type, root, me, p, dir) ;
            break ;
    case 3: relay_ring(buf, bytes, type, root, me, p, dir) ;
            break ;
    case 4: relay_bidi(buf, bytes, type, root, me, p, dir) ;
            break ;
    };
  /*******************************/

  traceblockend(BCAST1, type, root);

}

relay_cube ( buf, bytes, type, root, myname, p )
	char *buf ;
	int bytes, type, root, myname, p ;
/*
 *  Relay message buf of length bytes to all sons
 *  in a minimum spanning tree with given root.
 */
{
	int i, me, newroot, node ;

	newroot = root < envHLF ? root : root - envHLF ;
	if (myname == root && myname != newroot)
		 send0(buf, bytes, type, newroot) ;
	if (myname < envHLF) {
		me = myname^newroot ;
		for ( i = 1 ; i < envHLF ; i *= 2 )
			if (i > me) send0(buf, bytes, type, (me+i)^newroot) ;
		if (envTWN < p && envTWN != root)
			send0(buf, bytes, type, envTWN) ;
	}
}

relay_full ( buf, bytes, type, root, myname, p, dir )
	char *buf ;
	int bytes, type, root, myname, p, dir ;
/*
 *  If root, send message buf of length bytes to all processors,
 *  using full interconnections (direct sends); otherwise do nothing.
 */
{
	int neighbor0() ;
	int i ;

	if (myname == root) for ( i = 1 ; i < p ; i++ )
		send0(buf, bytes, type, neighbor0(i*dir, myname)) ;
}

relay_ring ( buf, bytes, type, root, myname, p, dir )
	char *buf ;
	int bytes, type, root, myname, p, dir ;
/*
 *  Relay message buf of length bytes to next processor
 *  in a unidirectional ring.
 */
{
	int neighbor0() ;
	int forward, back, node ;

	forward = neighbor0(1, myname) ;
	back = neighbor0(-1, myname) ;

	if (dir == 1) {
		if (forward != root) send0(buf, bytes, type, forward) ;
	}
	else {
		if (back != root) send0(buf, bytes, type, back) ;
	}
}

relay_bidi ( buf, bytes, type, root, myname, p, dir )
	char *buf ;
	int bytes, type, root, myname, p, dir ;
/*
 *  Relay message buf of length bytes to neighboring processors
 *  in a bidirectional ring.
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
	else
	{
		recvinfo0(&bytes, &type, &node) ;
		antipode = antipode0(root) ;
		if (myname != antipode) {
			if (node == forward)
				send0(buf, bytes, type, back);
			else if (forward != antipode)
				send0(buf, bytes, type, forward);
		}
	}
}

