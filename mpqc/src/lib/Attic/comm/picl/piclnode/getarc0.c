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

#include "globals.h"
#include "portglobals.h"

/*  Get architectural parameters to be used.
 *
 *  nprocs = actual number of processors to be used (may be smaller than
 *           the number allocated and/or loaded with load0).
 *
 *  top = 1 : hypercube
 *        2 : full connectivity (i.e., direct sends between all nodes)
 *        3 : unidirectional ring
 *        4 : bidirectional ring
 *
 *  ord = 0 : natural ordering of node numbers
 *        1 : Gray code ordering of nodes
 *
 *  dir =-1 : backward (direction used for unidirectional ring topology and
 *        1 : forward   for order in which "full" broadcasts are made)
 *
 *  Note: Use of hypercube topology and/or Gray code ordering requires
 *        that the number of processors be a power of two.
 */

getarc0 ( nprocs, top, ord, dir )
int *nprocs, *top, *ord, *dir ;
{
	*nprocs = envNPA ;
	*top = envTOP ;
	*ord = envORD ;
	*dir = envDIR ;
}

