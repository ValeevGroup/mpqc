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

/* node version of "begin receiving a message" */

envRECVBEG(buf, bytes, type, checking)
int bytes, type, checking;
char *buf;
{
  long irecv();

  /* try to add a cell to the linked list */
  if (envIRFREE != NULL){

    /* add cell to end */
    if (envIRHEAD == NULL) envIRHEAD = envIRFREE;
      else envIRTAIL->next = envIRFREE;
    envIRTAIL = envIRFREE;
    envIRFREE = envIRFREE->next;
    envIRTAIL->msgtype = type;
    envIRTAIL->next = NULL;

    /* check for errors */
    if (checking == 1){

      /* check whether type is legal */
      if (type >= TYPELIMIT1){
        if ((type <= TYPELIMIT2) || (type >= TYPELIMIT3)){
          sprintf(envMBUF, 
"recvbegin0:  message type %d is prohibited on iPSC/860 - exiting",type);
          envERRMSG(envMBUF);
          };
        };

      /* initiate the receive */
      envIRTAIL->msgid = irecv(type, buf, bytes);

      }
      else envIRTAIL->msgid = irecv(type, buf, bytes);

    }
    else{
    sprintf(envMBUF, 
"recvbegin0:  %d outstanding recvbegins is too many for iPSC/860 - exiting", 
            IRCVLIMIT);
    envERRMSG(envMBUF);
    };

}

