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

/* node version of "finish sending a message". */

envSENDEND(type)
int type;
{
  struct cell *i, *j;

  if (envISHEAD != NULL){

    /* look for first cell with desired type */
    i = envISHEAD;
    j = i;
    while ((j != envISTAIL) && (j->msgtype != type)){
      i = j;
      j = j->next;
      };

    if (j->msgtype == type){
      /* if found, wait for completion of the send */
      if (j->msgid != -1) msgwait(j->msgid);

      /* and update linked list */
      if (i != j){
        /* desired cell not at the front of the list */
        i->next = j->next;
        j->next = envISFREE;
        envISFREE = j;
        if (envISTAIL == j) envISTAIL = i;
        }
        else{
        /* removing header cell */
        envISHEAD = j->next;
        if (envISTAIL == j) envISTAIL = NULL;
        j->next = envISFREE;
        envISFREE = j;
        };

      }
      else{
      sprintf(envMBUF, 
"sendend0: no outstanding sendbegins of type %d - exiting", type);
      envERRMSG(envMBUF);
      };

    }
    else envERRMSG("sendend0: no outstanding sendbegins - exiting");
}
