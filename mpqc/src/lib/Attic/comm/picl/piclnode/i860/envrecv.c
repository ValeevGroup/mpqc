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

/* node version of receive a message. */

envRECV(buf, bytes, type)
int bytes, type;
char *buf;
{
  long ftype, infotype();

  /* check whether type is legal */
  if (type >= TYPELIMIT1){
    if ((type <= TYPELIMIT2) || (type >= TYPELIMIT3)){
      sprintf(envMBUF, 
"recv0:  message type %d is prohibited on iPSC/860 - exiting",type);
      envERRMSG(envMBUF);
      };
    };

  /* check for a promiscuous receive mistakenly finding a prohibited value */
  if (type < 0){

    /* wait until a message is available */
    cprobe(type);

    ftype = infotype();

    if (ftype >= TYPELIMIT1){
      if ((ftype <= TYPELIMIT2) || (ftype >= TYPELIMIT3)){

        sprintf(envMBUF,
"recv0:  receive type %d matched a prohibited type - exiting",type);
        envERRMSG(envMBUF);

        };
    
      };

    };

  /* read in the message */
  crecv(type, buf, bytes);

}

