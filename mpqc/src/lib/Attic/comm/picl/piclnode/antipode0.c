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

int antipode0 ( root )
	int root ;
/*
 *  Compute antipodal node from root in ring topology.
 *  Used in bidirectional ring communication.
 */
{
	int gray0(), ginv0() ;
	int p, top, ord, dir, retval ;

	getarc0(&p, &top, &ord, &dir) ;
	retval = (ord == 1) ? (gray0((ginv0(root)+(p+1)/2)%p))
		: ((root+(p+1)/2)%p) ;
	return( retval ) ;
}

