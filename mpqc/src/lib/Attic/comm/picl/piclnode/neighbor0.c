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


int neighbor0 ( i, node )
	int i, node ;

/*  Compute the neighbor i links away from given node in a ring topology.
 *
 *  A positive value for i gives i-th neighbor in forward direction,
 *  a negative value gives i-th neighbor in backward direction.  Values
 *  of 1 and -1 give immediate neighbors, larger values, 1 < i < p, give
 *  more distant "neighbors".
 *
 */
{
	int gray0(), ginv0() ;
	int p, top, ord, dir, retval ;

	getarc0(&p, &top, &ord, &dir) ;
	retval = (ord == 1) ? (gray0((ginv0(node)+p+i)%p)) : ((node+p+i)%p) ;
	return( retval ) ;
}

