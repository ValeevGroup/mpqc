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

/* set the types of tracing data collected */
/* if type <= 0, record only open, close, tracestart, tracelevel */
/*                    traceexit calls and summary statisitics) */
/* if type >= 1, also record "mark" calls*/
/* if type >= 2, also record higher level routines, and */
/*               send/recv/sync calls except */
/*               the send/recv calls within the higher level routines */
/* if type >= 3, also record send/recv/sync calls within higher level routines*/

tracelevel(trace,compstats,commstats)
int trace,compstats,commstats;
{
  long time1, time2;

  if (envFUL == 1){

    /* array full right now, but save parameters in case array */
    /* is ever flushed */
    envSTRA = trace;
    envSTCP = compstats;
    envSTCM = commstats;

    }
    else{

    envTRA = trace;
    envTCP = compstats;
    envTCM = commstats;

    };

  /* check whether tracing enabled */
  if (envTOPN == 1){

    ENVCLOCK(&time1, &time2);

    if (envFUL != 1){

      envTDATA[envNXT+KEY] = TLEVEL;
      envTDATA[envNXT+CLOCK1] = time1;
      envTDATA[envNXT+CLOCK2] = time2;
      envTDATA[envNXT+TRACE] = trace;
      envTDATA[envNXT+COMPF] = compstats;
      envTDATA[envNXT+COMMF] = commstats;
      envNXT = envNXT + TLEVELSIZE;

      if (envNXT > envMLTH) envFLUSH(time1, time2);

      };

    if (envTCP > 0){

      envTDATA[envNXT+KEY] = COMPSTATS;
      envTDATA[envNXT+CLOCK1] = time1;
      envTDATA[envNXT+CLOCK2] = time2;
      envTDATA[envNXT+IDLE1] = envCTIM1;
      envTDATA[envNXT+IDLE2] = envCTIM2;
      envNXT = envNXT + COMPSTATSIZE;

      if (envNXT > envMLTH) envFLUSH(time1, time2);

      };

    if (envTCM > 0){

      envTDATA[envNXT+KEY] = COMMSTATS;
      envTDATA[envNXT+CLOCK1] = time1;
      envTDATA[envNXT+CLOCK2] = time2;
      envTDATA[envNXT+MRECV] = envRCV;
      envTDATA[envNXT+MSENT] = envSNT;
      envTDATA[envNXT+MRVOL] = envRVOL;
      envTDATA[envNXT+MSVOL] = envSVOL;
      envTDATA[envNXT+MPROB] = envPRB;
      envNXT = envNXT + COMMSTATSIZE;

      if (envNXT > envMLTH) envFLUSH(time1, time2);

      };

    };

}

