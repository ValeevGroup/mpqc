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
#include "envopen.h"

/* node trace initialization routine */
/* should be first executable instruction in program */
/* and should only be called once */
/* - tracesize is the number of bytes to be allocated for trace data storage */
/* - flush == 1, if trace array fills up, send it back to secondary storage */
/*               and flush */
/*         != 1, if trace array fills up, stop tracing */
/* - sync  == 0, do nothing */
/*         == 1, execute clocksync0 to sync the processors. */

tracenode(tracesize, flush, sync)
int tracesize, flush, sync;
{
  long sendmess, time1, time2;
  double clock0();
  char *malloc();

  /* don't allocate multiple work arrays */
  /* if tracenode called multiple times */
  if (envTOPN != 0){
    free(envTDATA);
    free(envTPTR);
    envTOPN = 0;
    };

  /* initialize communication and tracing statistics variables */
  envOSPC = 0;
  envCTIM1 = 0;
  envCTIM2 = 0;
  envSNT = 0;
  envSVOL = 0;
  envRCV = 0;
  envRVOL = 0;
  envPRB = 0;

  /* allocate the trace array */
  envTSZ = tracesize/sizeof(long);
  if (envTSZ < MINTSZ){

    sprintf(envMBUF,
"tracenode: %d bytes too small for trace array; using %d bytes",
            tracesize,MINTSZ*sizeof(long));
    envTSZ = MINTSZ;
    sendmess = 1;

    }
    else sendmess = 0;

  envTDATA = (long *) malloc(envTSZ*sizeof(long));

  /* calculate how many packets required to send trace data back */
  /* and allocate an array of indices to record beginning indices */
  /* of each packet */
  /* (this algorithm can be wasteful if tmaxlong is "too small"-i.e.<200 ?*/
  /* Each packet will be less than tmaxlong because of need to */
  /* stop at the end of a trace record. Thus, each packet is smaller than */
  /* expected, and more than envSPTR packets could be stuffed into envTDATA.*/
  /* Since it isn't this small, this is not a concern) */
  envSPTR = (envTSZ - MINTSZ)/TMAXLONG + 1;
  envTPTR = (long *) malloc((envSPTR+1)*sizeof(long));

  /* synchronize (if needed) */
  if (sync == 1){
    if (envOPN == 0){
      if (envMOPN == 0){
        ENVOPEN();
        envMOPN = 1;
        };
      envOPN = 1;
      clocksync0();
      envOPN = 0;
      }
      else clocksync0();
    };

  /* set global switches */
  envTOPN = 1;
  envFLU = flush;

  /* send trace messages if they are required. (Done here so that timestamp */
  /* will be consistent) */
  if (sendmess == 1) tracemsg(envMBUF);

  if ((envTDATA == NULL) || (envTPTR == NULL)){

    sprintf(envMBUF,
"tracenode: allocating %d byte trace buffer exhausted memory; tracing disabled",
            tracesize);
    if (envTDATA != NULL) free(envTDATA);
    if (envTPTR != NULL) free(envTPTR);
    tracemsg(envMBUF);

    envTOPN = 0;

    }
    else{

    /* set beginning index of first trace data packet */
    envTPTR[0] = 0;
    envNPTR = 1;

    /* set first limit on amount of tracing data (before checking whether */
    /* need to ship data to host or stop tracing) */
    envMLTH = min(TMAXLONG,envTSZ) - FLUSHTRIGGER;

    };

  if (envTOPN != 0){

    /* record time tracing started */
    ENVCLOCK(&time1, &time2);

    envTDATA[KEY] = TSTART;
    envTDATA[CLOCK1] = time1;
    envTDATA[CLOCK2] = time2;
    envTDATA[CLKSYNC1] = envCLKSTART1;
    envTDATA[CLKSYNC2] = envCLKSTART2;
    envTDATA[OFFSET1] = envCLKOFFSET;
    envTDATA[OFFSET2] = 1000000*(envCLKOFFSET - envTDATA[OFFSET1]);
    envTDATA[DRIFT1] = envCLKDRIFT;
    envTDATA[DRIFT2] = 1000000*(envCLKDRIFT - envTDATA[DRIFT1]);
    envTDATA[CLKREF1] = envCLKREF;
    envTDATA[CLKREF2] = 1000000*(envCLKREF - envTDATA[CLKREF1]);
    envTDATA[TRACE2] = envTRA;
    envTDATA[COMPF2] = envTCP;
    envTDATA[COMMF2] = envTCM;

    envNXT = NTSTARTSIZE;

    };

}

