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

/* node version of "begin sending a message". */

envSENDBEG(buf, bytes, type, node, checking)
int bytes, type, node, checking;
char *buf;
{
  int i;
  long isend(), _isend(), numnodes();

  /* try to add a cell to the linked list */
  if (envISFREE != NULL){

    /* add cell to end */
    if (envISHEAD == NULL) envISHEAD = envISFREE;
      else envISTAIL->next = envISFREE;
    envISTAIL = envISFREE;
    envISFREE = envISFREE->next;
    envISTAIL->msgtype = type;
    envISTAIL->next = NULL;

    /* check the legality of the arguments */
    if (checking == 1){

      /* check whether destination is legal */
      if ((node < -1) || (node >= envNPA)){
        if (node != HOST0){
          sprintf(envMBUF, 
"sendbegin0:  node %d is not a legal destination on iPSC/860 - exiting",node);
          envERRMSG(envMBUF);
          };
        };

      /* check whether type is legal */
      if ((type < 0) || (type >= TYPELIMIT1)){
        if ((type <= TYPELIMIT2) || (type >= TYPELIMIT3)){
          sprintf(envMBUF, 
"sendbegin0:  message type %d is prohibited on iPSC/860 - exiting",type);
          envERRMSG(envMBUF);
          };
        };

      /* check whether the message is going to the host and is short enough */
      if (node == HOST0){
        if (bytes >= HMAXLENGTH){
          sprintf(envMBUF, 
"sendbegin0:  %d byte message too large to send to iPSC/860 host - exiting", 
                bytes);
          envERRMSG(envMBUF);
          };
        };

      if (node != -1){

        /* do the send (with the true destination) */
        if (node == HOST0) 
          envISTAIL->msgid = isend(type, buf, bytes, envHST, DPID);
          else 
          envISTAIL->msgid = isend(type, buf, bytes, node, DPID);

        }
        else{

        if (numnodes() != envNPA){

          /* simulate a broadcast */
          for (i=0; i<envNPA; i++){
            if (i != envME) csend(type, buf, bytes, i, DPID);
            envISTAIL->msgid = -1;
            };

          }
          else 
          envISTAIL->msgid = isend(type, buf, bytes, node, DPID);

        };

      }
      else{

      if (node != -1){

        /* do the send (with the true destination) */
        if (node == HOST0) 
          envISTAIL->msgid = _isend(type, buf, bytes, envHST, DPID);
          else 
          envISTAIL->msgid = _isend(type, buf, bytes, node, DPID);

        }
        else{

        if (numnodes() != envNPA){

          /* simulate a broadcast */
          for (i=0; i<envNPA; i++){
            if (i != envME) _csend(type, buf, bytes, i, DPID);
            envISTAIL->msgid = -1;
            };

          }
          else 
          envISTAIL->msgid = _isend(type, buf, bytes, node, DPID);

        };

      };

    }
    else{
    sprintf(envMBUF, 
"sendbegin0:  %d outstanding sendbegins is too many for iPSC/860 - exiting", 
            ISNDLIMIT);
    envERRMSG(envMBUF);
    };

}
