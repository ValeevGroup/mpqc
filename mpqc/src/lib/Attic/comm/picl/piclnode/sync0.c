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
#include "envsync.h"
#include "envupdatetime.h"

/* sync0 - barrier synchronization using dimensional exchange */
/* written in machine dependent routines for efficiency */

sync0()
{
  long cnt, oldctim1, oldctim2, start1, start2, stop1, stop2;

  /* check whether open statement has been exeuted */
  if (envOPN == 1){

    /* collect trace info */
    if (envTOPN == 1){

      /* save previous idle time values */
      oldctim1 = envCTIM1;
      oldctim2 = envCTIM2;

      /* record start time */
      ENVCLOCK(&start1, &start2);

      /* execute barrier */
      cnt = ENVSYNC();

      /* record stop time */
      ENVCLOCK(&stop1, &stop2);

      /* update running total */
      ENVUPDATETIME(&envCTIM1, &envCTIM2, start1, start2, stop1, stop2);

      /* write out trace info */
      if (envTRA > 1){
  
        envTDATA[envNXT+KEY] = SYNC;
        envTDATA[envNXT+CLOCK1] = start1;
        envTDATA[envNXT+CLOCK2] = start2;
        envNXT = envNXT + SYNCSIZE;
  
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

        /* update sync info */  
        envSNT = envSNT + cnt;
        envRCV = envRCV + cnt;
        envSVOL = envSVOL + cnt;
        envRVOL = envRVOL + cnt;

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
        /* update sync info */  
        envSNT = envSNT + cnt;
        envRCV = envRCV + cnt;
        envSVOL = envSVOL + cnt;
        envRVOL = envRVOL + cnt;
        };

      }
      else{
      ENVSYNC();
      };

    } 
    else envERRMSG("sync0: communication channel not open - exiting");

}

