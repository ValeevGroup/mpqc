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
#include "envprobe.h"

/* routine to test whether a message of the specified type is waiting */
/* to be read */
/* if message not found, 0 is returned, else 1 is returned */

int probe0(type)
int type;
{

  /* check whether open statement has been executed */
  if (envOPN == 1){

    /* update trace info */
    envPRB++;
 
    /* test for message */
    return( ENVPROBE(type, envCHECKING) );

    }
    else envERRMSG("probe0: communication channel not open - exiting");

}
