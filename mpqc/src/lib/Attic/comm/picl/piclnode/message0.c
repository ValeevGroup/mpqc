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
#include "envprintf.h"
#include "envsend.h"
#include "envupdatetime.h"

/* send a message of type MSG0 back to stadard out. If there is a host, */
/* then send the message via the host, else send it directly. */
/* The message should be a string of length less than 80 characters. */

message0(message)
char *message;
{
  int i, j, lth, lth2, strlen();
  long oldctim1, oldctim2, start1, start2, stop1, stop2;
  double time, clock0();
  char mesbuf[120];

  /* check whether a communication channel has been opened */
  if (envOPN == 1){

    /* check whether tracing enabled */
    if (envTOPN == 1){

      /* save previous idle time values */
      oldctim1 = envCTIM1;
      oldctim2 = envCTIM2;

      /* record start time */
      ENVCLOCK(&start1, &start2);

      /* prepare the message */         
      time = clock0();
      sprintf(mesbuf, "message  node %d clock %f ", envME, (float) time);
      lth = strlen(mesbuf);
      lth2 = lth + strlen(message);
      if (lth2 > 118) lth2 = 118;
      j=0;
      for (i=lth; i<lth2; i++, j++) mesbuf[i] = message[j];
      mesbuf[i] = '\n'; 
      mesbuf[i+1] = '\0';
      lth = strlen(mesbuf) + 1;

      /* send the message */
      if (envHOSTED == 1){ /*if hosted, send to host to print on standard out */
        ENVSEND(mesbuf, lth, MSG0, HOST0, 0);
        }
        else{ /* else, send to standard out directly */
        ENVPRINTF("%s",mesbuf);
        };

      /* record stop time */
      ENVCLOCK(&stop1, &stop2);

      /* update running totals */
      ENVUPDATETIME(&envCTIM1, &envCTIM2, start1, start2, stop1, stop2);

      /* write out trace records */
      if (envTRA > 1){

        envTDATA[envNXT+KEY] = MESSAGE;
        envTDATA[envNXT+CLOCKSTART1] = start1;
        envTDATA[envNXT+CLOCKSTART2] = start2;
        envNXT = envNXT + MESSAGESIZE;

        if (envNXT > envMLTH) envFLUSH(stop1, stop2);

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

        /* update message info */    
        envSNT++;
        envSVOL = envSVOL + lth;

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
        /* update message info */    
        envSNT++;
        envSVOL = envSVOL + lth;
        };
    
      }
      else{

      /* prepare the message */         
      time = clock0();
      sprintf(mesbuf, "message  node %d clock %f ", envME, (float) time);
      lth = strlen(mesbuf);
      lth2 = lth + strlen(message);
      if (lth2 > 118) lth2 = 118;
      j=0;
      for (i=lth; i<lth2; i++, j++) mesbuf[i] = message[j];
      mesbuf[i] = '\n'; 
      mesbuf[i+1] = '\0';

      /* send the message */
      if (envHOSTED == 1){ /*if hosted, send to host to print on standard out */
        lth = strlen(mesbuf)+1;
        ENVSEND(mesbuf, lth, MSG0, HOST0, 0);
        }
        else{ /* else, send to standard out directly */
        ENVPRINTF("%s",mesbuf);
        };

      };

    }
    else envERRMSG("message0: communication channel not open - exiting");

}

