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
#include "envupdatetime.h"
#include "envrecv.h"
#include "envsend.h"
#include "envsync.h"

/* node clock synchronization routine */
/* (calls envCLOCKSYNC to calculate clock offset and drift) */
/* (then synchronizes the processors) */
/* */
/* if first time called, then define the synchronization time to be the clock */
/* start time and normalize all trace record timestamps and clock0 accordingly*/
/* Otherwise, don't change the current clock start time, and use envSYNC */
/* to approximately resynchronize */
/* Note: setarc0 resets synchronization time and normalization parameters */
/* to zero. */
clocksync0()
{
  int times, first;
  long cnt, oldctim1, oldctim2, start1, start2, stop1, stop2;
  double t1, t2, clock0();

  /* check whether open statement has been exeuted */
  if (envOPN == 1){

    /* collect trace info */
    if (envTOPN == 1){

      /* save previous idle time values */
      oldctim1 = envCTIM1;
      oldctim2 = envCTIM2;

      /* record start time */
      ENVCLOCK(&start1, &start2);

      /* calculate clock offset and drift */
      /* after making sure that the network is quiescent */
      first = ( envCLKSYNC == 0 ? 1 : 0 );
      cnt = ENVSYNC();
      envCLOCKSYNC();

      /* resynchronize processors */

      /* execute barrier */
      ENVSYNC();

      if (first == 1){

        if (envNPA > 1){

          /* check how long a broadcast and a sync take */
          /* (second sync is to make sure that the estimate is conservative) */
          if (envME == 0){

            t1 = clock0();
            ENVSEND(&t1, sizeof(double), CLOCK_PREP, -1, 0) ;  
            ENVSYNC();
            ENVSYNC();
            t2 = clock0();
            envWAIT = 2.0*(t2-t1);

            }
            else{

            /* wait for message from node 0 */      
            ENVRECV(&t1, sizeof(double), CLOCK_PREP, 0);  
            /* and wait until all processors have received it */
            ENVSYNC();
            ENVSYNC();

            };

          /* wait until an agreed upon time */
          if (envME == 0){
            /* pick a start time sufficiently far in the future */
            t2 = clock0() + envWAIT;
            /* send the "start" time to the other nodes */
            ENVSEND(&t2, sizeof(double), CLOCK_SYNC, -1, 0) ;  
            }
            else{
            /* wait for start time from node 0 */      
            ENVRECV(&t2, sizeof(double), CLOCK_SYNC, 0);  
            };

          if (clock0() < t2){
            while (clock0() < t2) ;
            }
            else message0("clocksync0: synchronization algorithm failed\n");

          }
          else t2 = clock0();

        /* record clock start time */
        ENVCLOCK(&envCLKSTART1,&envCLKSTART2);
        envCLKSTART = t2;

        };

      /* record stop time */
      ENVCLOCK(&stop1, &stop2);

      /* update running time counter */
      ENVUPDATETIME(&envCTIM1, &envCTIM2, start1, start2, stop1, stop2);

      /* write out trace info */
      if (envFUL != 1){
  
        envTDATA[envNXT+KEY] = CLOCKSYNC;
        envTDATA[envNXT+CLOCK1] = start1;
        envTDATA[envNXT+CLOCK2] = start2;
        envTDATA[envNXT+CLKSYNC1] = envCLKSTART1;
        envTDATA[envNXT+CLKSYNC2] = envCLKSTART2;
        envTDATA[envNXT+OFFSET1] = envCLKOFFSET;
        envTDATA[envNXT+OFFSET2] = 
         1000000*(envCLKOFFSET - envTDATA[envNXT+OFFSET1]);
        envTDATA[envNXT+DRIFT1] = envCLKDRIFT;
        envTDATA[envNXT+DRIFT2] = 
         1000000*(envCLKDRIFT - envTDATA[envNXT+DRIFT1]);
        envTDATA[envNXT+CLKREF1] = envCLKREF;
        envTDATA[envNXT+CLKREF2] = 
         1000000*(envCLKREF - envTDATA[envNXT+CLKREF1]);
        envNXT = envNXT + CLKSYNCSIZE;

        if (envNXT > envMLTH) envFLUSH(start1, start2);
    
        };
  
      if (envTCP > 1){
  
        /* write out info preceding clocksync */
        envTDATA[envNXT+KEY] = COMPSTATS;
        envTDATA[envNXT+CLOCK1] = start1;
        envTDATA[envNXT+CLOCK2] = start2;
        envTDATA[envNXT+IDLE1] = oldctim1;
        envTDATA[envNXT+IDLE2] = oldctim2;
        envNXT = envNXT + COMPSTATSIZE;
     
        /* write out info following clocksync */
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

        /* update clocksync info */  
        if (envNPA > 1){
          times = 0;
          if ((first == 1) && (CALCOFFSET == 1)) times = 1;
          if (CALCDRIFT == 1) times = times + 1;
          }
          else{
          times = 0;
          cnt = 0;
          };

        if (envME == 0){
          envSNT = envSNT + 2*(1+first)*cnt + times*envNPA*(1 + REPS);
          envRCV = envRCV + 2*(1+first)*cnt + times*envNPA*REPS;
          envSVOL = envSVOL + 2*(1+first)*cnt + 
                              times*envNPA*(1 + REPS)*sizeof(double);
          envRVOL = envRVOL + 2*(1+first)*cnt + 
                              times*envNPA*REPS*sizeof(double);
          }
          else{
          envSNT = envSNT + 2*(1+first)*cnt + times*envNPA*REPS;
          envRCV = envRCV + 2*(1+first)*cnt + times*envNPA*(1 + REPS);
          envSVOL = envSVOL + 2*(1+first)*cnt + 
                              times*envNPA*REPS*sizeof(double);
          envRVOL = envRVOL + 2*(1+first)*cnt + 
                              times*envNPA*(1 + REPS)*sizeof(double);
          };

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

        /* update clocksync info */  
        if (envNPA > 1){
          times = 0;
          if ((first == 1) && (CALCOFFSET == 1)) times = 1;
          if (CALCDRIFT == 1) times = times + 1;
          }
          else{
          times = 0;
          cnt = 0;
          };

        if (envME == 0){
          envSNT = envSNT + 2*(1+first)*cnt + times*envNPA*(1 + REPS);
          envRCV = envRCV + 2*(1+first)*cnt + times*envNPA*REPS;
          envSVOL = envSVOL + 2*(1+first)*cnt + 
                              times*envNPA*(1 + REPS)*sizeof(double);
          envRVOL = envRVOL + 2*(1+first)*cnt + 
                              times*envNPA*REPS*sizeof(double);
          }
          else{
          envSNT = envSNT + 2*(1+first)*cnt + times*envNPA*REPS;
          envRCV = envRCV + 2*(1+first)*cnt + times*envNPA*(1 + REPS);
          envSVOL = envSVOL + 2*(1+first)*cnt + 
                              times*envNPA*REPS*sizeof(double);
          envRVOL = envRVOL + 2*(1+first)*cnt + 
                              times*envNPA*(1 + REPS)*sizeof(double);
          };

        };

      }
      else{

      /* calculate clock offset and drift */
      /* after making sure that the network is quiescent */
      first = ( envCLKSYNC == 0 ? 1 : 0 );
      ENVSYNC();
      envCLOCKSYNC();

      /* resynchronize processors */

      /* execute barrier */
      ENVSYNC();

      if (first == 1){

        if (envNPA > 1){

          /* check how long a broadcast and a sync take */
          /* (second sync is to make sure that the estimate is conservative) */

          if (envME == 0){

            t1 = clock0();
            ENVSEND(&t1, sizeof(double), CLOCK_PREP, -1, 0) ;  
            ENVSYNC();
            ENVSYNC();
            t2 = clock0();
            envWAIT = 2.0*(t2-t1);

            }
            else{

            /* wait for message from node 0 */      
            ENVRECV(&t1, sizeof(double), CLOCK_PREP, 0);  
            /* and wait until all processors have received it */
            ENVSYNC();
            ENVSYNC();

            };

          /* wait until an agreed upon time */
          if (envME == 0){
            /* pick a trace start time sufficiently far in the future */
            t2 = clock0() + envWAIT;
            /* send the "start" time to the other nodes */
            ENVSEND(&t2, sizeof(double), CLOCK_SYNC, -1, 0) ;  
            }
            else{
            /* wait for start time from node 0 */      
            ENVRECV(&t2, sizeof(double), CLOCK_SYNC, 0);  
            };

          if (clock0() < t2){
            while (clock0() < t2) ;
            }
            else message0("clocksync0: synchronization algorithm failed\n");

          }
          else t2 = clock0();

        /* record clock start time (if first time) */
        ENVCLOCK(&envCLKSTART1,&envCLKSTART2);
        envCLKSTART = t2;

        };

      };

    } 
    else envERRMSG("clocksync0: communication channel not open - exiting");

}




