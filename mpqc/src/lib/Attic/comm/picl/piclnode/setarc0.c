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
#include "envglobals.h"
#include "portglobaldefs.h"

/*  Set architectural parameters to be used.
 *
 *  nprocs = actual number of processors to be used (may be smaller than
 *           the number allocated and/or loaded with load0).
 *
 *  top = 1 : hypercube
 *        2 : full connectivity (i.e., direct sends between all nodes)
 *        3 : unidirectional ring
 *        4 : bidirectional ring
 *
 *  dir =-1 : backward (direction used for unidirectional ring topology and
 *        1 : forward   for order in which "full" broadcasts are made)
 *
 *  ord = 0 : natural ordering of node numbers
 *        1 : Gray code ordering of nodes
 *
 *  Note: Use of Gray code ordering requires that the number of processors
 *        be a power of two.
 */

setarc0(nprocs, top, ord, dir)
int *nprocs, *top, *ord, *dir ;
{
  int info[4], i, host, me, p ;

  if (envHOSTED == 1){

    who0(&p, &me, &host) ;
    recv0(info, 4*sizeof(int), PORT_INIT) ;
    if (envNPH <= 0) envNPH = envNPA;
    if ((info[0] > 0) && (info[0] <= envNPH)){

      /* unnormalize the clock if the number of processors */
      /* is changing */
      if (info[0] != envNPA){
        envCLKSYNC = 0;
        envCLKSTART = 0.0;
        };

      /* define virtual machine */
      *nprocs = envNPA = info[0] ;
      *top = envTOP = info[1] ;
      *ord = envORD = info[2] ;
      *dir = envDIR = info[3] ;
      for (i = 1 ; i < envNPA ; i *= 2) ;
      envHLF = i == envNPA ? i : i/2 ;
      envTWN = me < envHLF ? me + envHLF : me - envHLF ;

      }
      else{

      /* unnormalize the clock if the number of processors */
      /* is changing */
      if (envNPH != envNPA){
        envCLKSYNC = 0;
        envCLKSTART = 0.0;
        };

      /* exiting setarc0 "logic", so reset defaults */
      *nprocs = 0 ;
      envNPA = envNPH;
      envNPH = -1;
      *top = envTOP = 2 ;
      *ord = envORD = 0 ;
      *dir = envDIR = 1 ;
      for (i = 1 ; i < envNPA ; i *= 2) ;
      envHLF = i == envNPA ? i : i/2 ;
      envTWN = me < envHLF ? me + envHLF : me - envHLF ;

      };

    /* rebroadcast info to other nodes */
    if (me == 0) send0(info, 4*sizeof(int), PORT_INIT, -1) ;

    }
    else getarc0(nprocs, top, ord, dir);

}

