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
#include "envrecvinfo.h"

/* node version of received message info */

recvinfo0(bytes, type, node)
int *bytes, *type, *node;
{

  if (envOPN == 1){
    ENVRECVINFO(bytes, type, node);
    }
    else envERRMSG("recvinfo0: communication channel not open - exiting");

}
