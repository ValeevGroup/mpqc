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

/* node version of clock - returns time in seconds */

double clock0()
{
  double t, _dclock();

  if (envCLKSYNC == 0) return(_dclock());
    else{
    t = _dclock(); 
     /* t + envCLKOFFSET is the "node 0" time if clocks don't drift */
     /* t - envCLKREF is the time since envCLKOFFSET was */
     /* calculated, and therefore is the time over which the clock will */
     /* have drifted. */
    return(t+envCLKOFFSET+envCLKDRIFT*(t-envCLKREF)-envCLKSTART);
    };

}

