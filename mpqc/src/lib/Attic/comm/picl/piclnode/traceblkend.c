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

/* traceblockend marks the end of a piece of code in */
/* which trace records should not generated when tracelevel == 2 */
/* and envTRACEBLOCK = 1 */

traceblockend(type, location, parameter)
int type, location, parameter;
{
  int timechk;
  long time1, time2;

  /* collect trace info? */
  if (envTOPN == 1){

    /* generate trace records for blocks? */
    if (envTRACEBLOCK == 1){

      /* reset tracing levels */
      envTRA++; envTCP++; envTCM++;

      /* reset clock indicator */
      timechk = 0;

      if (envTRA > 1){

        ENVCLOCK(&time1, &time2);
        timechk = 1;

        /* write out a trace record for traceblockend */
        envTDATA[envNXT+KEY] = TBLOCKEND;
        envTDATA[envNXT+CLOCK1] = time1;
        envTDATA[envNXT+CLOCK2] = time2;
        envTDATA[envNXT+TYPE] = type;
        envTDATA[envNXT+LOCATION] = location;
        envTDATA[envNXT+PARAM] = parameter;
        envNXT = envNXT + TBLOCKSIZE;
    
        if (envNXT > envMLTH) envFLUSH(time1, time2);

        };

      if (envTCP == 2){

        if (timechk == 0){
          ENVCLOCK(&time1, &time2);
          timechk = 1;
          };

        /* write out a compstats record */
        envTDATA[envNXT+KEY] = COMPSTATS;
        envTDATA[envNXT+CLOCK1] = time1;
        envTDATA[envNXT+CLOCK2] = time2;
        envTDATA[envNXT+IDLE1] = envCTIM1;
        envTDATA[envNXT+IDLE2] = envCTIM2;
        envNXT = envNXT + COMPSTATSIZE;

        if (envNXT > envMLTH) envFLUSH(time1, time2);

        };

      if (envTCM == 2){

        if (timechk == 0){
          ENVCLOCK(&time1, &time2);
          };

        /* write out a commstats record */
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

    };

}

