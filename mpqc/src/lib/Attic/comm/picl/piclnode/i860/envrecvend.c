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

/* node version of "finish receiving a message" */

long envRECVEND(type, checking)
int type, checking;
{
  long found, ftype, msgdone();
  struct cell *i, *j;

  if (envIRHEAD != NULL){

    /* look for first cell with desired type */
    i = envIRHEAD;
    j = i;
    while ((j != envIRTAIL) && (j->msgtype != type)){
      i = j;
      j = j->next;
      };

    if (j->msgtype == type){
      /* if found, wait for message to arrive */
      found = msgdone(j->msgid);
      if (found == 0) msgwait(j->msgid);

      /* and update linked list */
      if (i != j){
        /* desired cell not at the front of the list */
        i->next = j->next;
        j->next = envIRFREE;
        envIRFREE = j;
        if (envIRTAIL == j) envIRTAIL = i;
        }
        else{
        /* removing header cell */
        envIRHEAD = j->next;
        if (envIRTAIL == j) envIRTAIL = NULL;
        j->next = envIRFREE;
        envIRFREE = j;
        };

      }
      else{
      sprintf(envMBUF, 
"recvend0: no outstanding recvbegins of type %d - exiting", type);
      envERRMSG(envMBUF);
      };

    /* check for errors */
    if (checking == 1){

      /* check for a promiscuous receive mistakenly finding a prohibited value */
      if (type < 0){

        /* check type of message received */
        ftype = infotype();

        if (ftype >= TYPELIMIT1){
          if ((ftype <= TYPELIMIT2) || (ftype >= TYPELIMIT3)){

            sprintf(envMBUF,
"recvend0:  receive type %d matched a prohibited type - exiting",type);
            envERRMSG(envMBUF);

            };
    
          };

        };

      };

    }
    else envERRMSG("recvend0: no outstanding recvbegins - exiting");

  if (found == 0) return(1);
    else return(0);

}

