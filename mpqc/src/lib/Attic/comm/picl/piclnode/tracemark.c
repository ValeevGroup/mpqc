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

/* "mark" routine */

tracemark(marktype)
int marktype;
{
  long time1, time2;

  /* check whether tracing enabled */
  if (envTOPN == 1){

    ENVCLOCK(&time1, &time2);

    if (envTRA > 0){

      envTDATA[envNXT+KEY] = TMARK;
      envTDATA[envNXT+CLOCK1] = time1;
      envTDATA[envNXT+CLOCK2] = time2;
      envTDATA[envNXT+TYPE] = marktype;
      envNXT = envNXT + TMARKSIZE;

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

