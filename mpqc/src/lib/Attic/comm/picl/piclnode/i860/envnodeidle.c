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

/* translates "raw" node idle times to double precision value in seconds */

double envNODEIDLE(time1, time2)
unsigned long time1, time2;
{
  double t;

  t = NSFACTOR1*time1 + NSFACTOR2*time2;
  return(t + envCLKDRIFT*t);

}
