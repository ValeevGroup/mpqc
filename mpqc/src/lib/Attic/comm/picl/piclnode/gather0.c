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

gather0(vec, n, type, root)
float *vec ;
int n, type, root ;

/*
 *  Gather components of a vector into root from all processors
 *  using given topology.  Caution: vec is overwritten.
 */
{
  int p, me, host, top, ord, dir, root0 ;

  traceblockbegin(GATHER, type, root);

  /**** calculate the gather ****/
  who0(&p, &me, &host) ;
  getarc0(&p, &top, &ord, &dir) ;

  root0 = (root == host) ? 0 : root ;

  if (p > 1) switch (top) {
    case 1: gather_cube(vec, n, type, root0, me, p) ;
            break ;
    case 2: gather_full(vec, n, type, root0, me, p) ;
            break ;
    case 3: gather_ring(vec, n, type, root0, me, p) ;
            break ;
    case 4: gather_bidi(vec, n, type, root0, me, p) ;
            break ;
    };

  if ((root == host) && (me == 0)) send0(vec, n*sizeof(float), type, host) ;
  /*******************************/

  traceblockend(GATHER, type, root);

}

gather_cube ( vec, n, type, root, myname, p )
	float *vec ;
	int n, type, root, myname, p ;
/*
 *  Gather components of a vector into root from all processors,
 *  using spanning tree.  Caution: vec is overwritten.
 */
{
	char *malloc0() ;
	int i, me, bytes, cnt, node ;
	float *t ;

	bytes = n*sizeof(float) ;
	t = (float *)malloc0(bytes) ;
	me = myname^root ;
	p /= 2 ;
	if (me < p) {
		recv0(t, bytes, type) ;
		for ( i = 0 ; i < n ; i++ )
			vec[i] = vec[i] != 0.0 ? vec[i] : t[i] ;
		if (p != 1) gather_cube(vec, n, type, root, myname, p ) ;
	}
	else send0(vec, bytes, type, (me-p)^root) ;
	free (t) ;
}

gather_full ( vec, n, type, root, myname, p )
	float *vec ;
	int n, type, root, myname, p ;
/*
 *  Gather components of a vector into root
 *  from all processors, using full interconnections
 *  Caution: vec is overwritten.
 */
{
	char *malloc0() ;
	int i, j, bytes, cnt, node ;
	float *t ;

	bytes = n*sizeof(float) ;

	if (root != myname) send0(vec, bytes, type, root) ;
	else {
		t = (float *)malloc0(bytes) ;
		for ( i = 1 ; i < p ; i++ ) {
			recv0(t, bytes, type) ;
			for ( j = 0 ; j < n ; j++ )
				vec[j] = vec[j] != 0.0 ? vec[j] : t[j] ;
		}
		free (t) ;
	}
}

gather_ring ( vec, n, type, root, myname, p )
	float *vec ;
	int n, type, root, myname, p ;
/*
 *  Gather components of a vector into root
 *  from all processors, using a unidirectional ring
 *  Caution: vec is overwritten.
 */
{
	char *malloc0() ;
	int neighbor0() ;
	int i, bytes, cnt, node ;
	float *t ;

	bytes = n*sizeof(float) ;
	t = (float *)malloc0(bytes) ;
	if (root != neighbor0(-1, myname)) {
		recv0(t, bytes, type) ;
		for ( i = 0 ; i < n ; i++ )
			vec[i] = vec[i] != 0.0 ? vec[i] : t[i] ;
	}
	if (myname != root) send0(vec, bytes, type, neighbor0(1, myname)) ;
	free (t) ;
}

gather_bidi ( vec, n, type, root, myname, p )
	float *vec ;
	int n, type, root, myname, p ;
/*
 *  Gather components of a vector into root
 *  from all processors, using a bidirectional ring
 *  Caution: vec is overwritten.
 */
{
	char *malloc0() ;
	int neighbor0(), antipode0() ;
	int i, bytes, cnt, node, forward, back, antipode ;
	float *t ;

	forward = neighbor0(1, myname) ;
	back = neighbor0(-1, myname) ;
	antipode = antipode0(root) ;

	bytes = n*sizeof(float) ;
	t = (float *)malloc0(bytes) ;

	if (myname == root) {
		recv0(t, bytes, type) ;
		for ( i = 0 ; i < n ; i++ )
			vec[i] = vec[i] != 0.0 ? vec[i] : t[i] ;
		if (p > 2) {
			recv0(t, bytes, type) ;
			for ( i = 0 ; i < n ; i++ )
				vec[i] = vec[i] != 0.0 ? vec[i] : t[i] ;
		}
	}
	else if (myname == antipode) {
		send0(vec, bytes, type, forward) ;
	}
	else if (forward == antipode && p > 2) {
		send0(vec, bytes, type, back) ;
	}
	else {
		recv0(t, bytes, type) ;
		for ( i = 0 ; i < n ; i++ )
			vec[i] = vec[i] != 0.0 ? vec[i] : t[i] ;
		if (node == back) {
			send0(vec, bytes, type, forward) ;
		}
		else {
			send0(vec, bytes, type, back) ;
		}
	}
	free (t) ;
}

