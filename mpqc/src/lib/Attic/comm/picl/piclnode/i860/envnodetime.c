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

/* translates "raw" node timestamps to double precision value in seconds */

double envNODETIME(time1, time2)
unsigned long time1, time2;
{
  double t;

  t = NSFACTOR1*time1 + NSFACTOR2*time2;
  /* t + envNOF is the "node 0" time if clocks don't drift */
  /* t - envNRF is the time since envNOF was */
  /*  calculated, and therefore is the time over which the clock will */
  /*  have drifted. */
  return(t + envCLKOFFSET + envCLKDRIFT*(t - envCLKREF) - envCLKSTART);

}
