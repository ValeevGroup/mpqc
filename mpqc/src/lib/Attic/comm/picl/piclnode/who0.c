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

/* node version of who */

who0(numproc, me, host)
int *numproc, *me, *host;
{

  if (envOPN == 1){

    *numproc = envNPA;
    *me = envME;
    *host = HOST0;

    }
    else envERRMSG("who0: communication channel not open - exiting");

}


