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
#include "envclock.h"
#include "envrecvbeg.h"
#include "envupdatetime.h"

/* node version of "begin receiving a message". */

recvbegin0(buf, bytes, type)
int bytes, type;
char *buf;
{
  long oldctim1, oldctim2, start1, start2, stop1, stop2;

  /* check whether open statement has been executed */
  if (envOPN == 1){

    /* collect trace info? */
    if (envTOPN == 1){

      /* save previous idle time values */
      oldctim1 = envCTIM1;
      oldctim2 = envCTIM2;

      /* record start time */
      ENVCLOCK(&start1, &start2);

      /* initiate the receive */
      ENVRECVBEG(buf, bytes, type, envCHECKING);

      /* record stop time */
      ENVCLOCK(&stop1, &stop2);

      /* update running total */
      ENVUPDATETIME(&envCTIM1, &envCTIM2, start1, start2, stop1, stop2);

      /* write out trace info */
      if (envTCP > 1){
    
        /* write out info preceding sendend */
        envTDATA[envNXT+KEY] = COMPSTATS;
        envTDATA[envNXT+CLOCK1] = start1;
        envTDATA[envNXT+CLOCK2] = start2;
        envTDATA[envNXT+IDLE1] = oldctim1;
        envTDATA[envNXT+IDLE2] = oldctim2;
        envNXT = envNXT + COMPSTATSIZE;
       
        /* write out info following sendend */
        envTDATA[envNXT+KEY] = COMPSTATS;
        envTDATA[envNXT+CLOCK1] = stop1;
        envTDATA[envNXT+CLOCK2] = stop2;
        envTDATA[envNXT+IDLE1] = envCTIM1;
        envTDATA[envNXT+IDLE2] = envCTIM2;
        envNXT = envNXT + COMPSTATSIZE;
      
        if (envNXT > envMLTH) envFLUSH(stop1, stop2);
     
        };
    
      }
      else{
      ENVRECVBEG(buf, bytes, type, envCHECKING);
      };

    }
    else envERRMSG("recvbegin0: communication channel not open - exiting");

}

