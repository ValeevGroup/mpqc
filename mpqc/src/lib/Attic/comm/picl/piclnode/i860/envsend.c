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
#include <stdio.h>

/* node version of send a message. */

envSEND(buf, bytes, type, node)
int bytes, type, node;
char *buf;
{
  int i;
  long numnodes();

  /* check whether destination is legal */
  if ((node < -1) || (node >= envNPA)){
    if (node != HOST0){
      sprintf(envMBUF, 
"send0:  node %d is not a legal destination on iPSC/860 - exiting",node);
      envERRMSG(envMBUF);
      };
    };

  /* check whether type is legal */
  if ((type < 0) || (type >= TYPELIMIT1)){
    if ((type <= TYPELIMIT2) || (type >= TYPELIMIT3)){
      sprintf(envMBUF, 
"send0:  message type %d is prohibited on iPSC/860 - exiting",type);
      envERRMSG(envMBUF);
      };
    };

  /* check whether the message is going to the host and is short enough */
  if (node == HOST0){
    if (bytes >= HMAXLENGTH){
      sprintf(envMBUF, 
"send0:  %d byte message too large to send to iPSC/860 host - exiting", bytes);
      envERRMSG(envMBUF);
      };
    };

  if (node != -1){

    /* do the send (with the true destination) */
    if (node == HOST0) csend(type, buf, bytes, envHST, DPID);
      else csend(type, buf, bytes, node, DPID);

    }
    else{

    if (numnodes() != envNPA){

      /* simulate a broadcast */
      for (i=0; i<envNPA; i++)
        if (i != envME) csend(type, buf, bytes, i, DPID);

      }
      else csend(type, buf, bytes, node, DPID);

    };

}
