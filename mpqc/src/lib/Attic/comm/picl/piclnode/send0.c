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
#include "envsend.h"
#include "envclock.h"
#include "envupdatetime.h"

/* node version of send a message. The time the message is sent */
/* is computed and the trace line is written. */

send0(buf, bytes, type, node)
int bytes, type, node;
char *buf;
{
  int i;
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

      /* do the send */
      ENVSEND(buf, bytes, type, node, envCHECKING);

      /* record stop time */
      ENVCLOCK(&stop1, &stop2);

      /* update running totals */
      ENVUPDATETIME(&envCTIM1, &envCTIM2, start1, start2, stop1, stop2);

      /* write out trace info */
      if (envTRA > 1){
  
        /* write out trace info */
        envTDATA[envNXT+KEY] = SEND;
        envTDATA[envNXT+CLOCK1] = start1;
        envTDATA[envNXT+CLOCK2] = start2;
        envTDATA[envNXT+TYPE] = type;
        envTDATA[envNXT+TO] = node;
        envTDATA[envNXT+LENGTH] = bytes;
        envNXT = envNXT + SENDSIZE;
    
        if (envNXT > envMLTH) envFLUSH(start1, start2);
      
        };
    
      if (envTCP > 1){
    
        /* write out info preceding send */
        envTDATA[envNXT+KEY] = COMPSTATS;
        envTDATA[envNXT+CLOCK1] = start1;
        envTDATA[envNXT+CLOCK2] = start2;
        envTDATA[envNXT+IDLE1] = oldctim1;
        envTDATA[envNXT+IDLE2] = oldctim2;
        envNXT = envNXT + COMPSTATSIZE;
       
        /* write out info following send */
        envTDATA[envNXT+KEY] = COMPSTATS;
        envTDATA[envNXT+CLOCK1] = stop1;
        envTDATA[envNXT+CLOCK2] = stop2;
        envTDATA[envNXT+IDLE1] = envCTIM1;
        envTDATA[envNXT+IDLE2] = envCTIM2;
        envNXT = envNXT + COMPSTATSIZE;
      
        if (envNXT > envMLTH) envFLUSH(stop1, stop2);
     
        };
    
      if (envTCM > 1){
      
        /* write out info preceding send */
        envTDATA[envNXT+KEY] = COMMSTATS;
        envTDATA[envNXT+CLOCK1] = start1;
        envTDATA[envNXT+CLOCK2] = start2;
        envTDATA[envNXT+MRECV] = envRCV;
        envTDATA[envNXT+MSENT] = envSNT;
        envTDATA[envNXT+MRVOL] = envRVOL;
        envTDATA[envNXT+MSVOL] = envSVOL;
        envTDATA[envNXT+MPROB] = envPRB;
        envNXT = envNXT + COMMSTATSIZE;
    
        /* update send info */
        envSNT++;
        envSVOL = envSVOL + bytes;

        /* write out info following send */
        envTDATA[envNXT+KEY] = COMMSTATS;
        envTDATA[envNXT+CLOCK1] = stop1;
        envTDATA[envNXT+CLOCK2] = stop2;
        envTDATA[envNXT+MRECV] = envRCV;
        envTDATA[envNXT+MSENT] = envSNT;
        envTDATA[envNXT+MRVOL] = envRVOL;
        envTDATA[envNXT+MSVOL] = envSVOL;
        envTDATA[envNXT+MPROB] = envPRB;
        envNXT = envNXT + COMMSTATSIZE;
    
        if (envNXT > envMLTH) envFLUSH(stop1, stop2);
   
        }
        else{
        /* update send info */
        envSNT++;
        envSVOL = envSVOL + bytes;
        };
     
      }
      else{
      ENVSEND(buf, bytes, type, node, envCHECKING);
      };

    }
    else envERRMSG("send0: communication channel not open - exiting");

}
