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

/* write out node trace data and disable node tracing */

traceexit()
{
  long time1, time2;

  /* check whether tracing currently enabled */
  if (envTOPN == 1){

    ENVCLOCK(&time1, &time2);

    envTDATA[envNXT+KEY] = TEXIT;
    envTDATA[envNXT+CLOCK1] = time1;
    envTDATA[envNXT+CLOCK2] = time2;
    envTDATA[envNXT+OSPACE] = envOSPC + (envNXT + TEXITSIZE + 
                              COMPSTATSIZE + COMMSTATSIZE)*sizeof(long);
    envNXT = envNXT + TEXITSIZE; 

    /* write out compstats summary statistics */
    envTDATA[envNXT+KEY] = COMPSTATS;
    envTDATA[envNXT+CLOCK1] = time1;
    envTDATA[envNXT+CLOCK2] = time2;
    envTDATA[envNXT+IDLE1] = envCTIM1;
    envTDATA[envNXT+IDLE2] = envCTIM2;
    envNXT = envNXT + COMPSTATSIZE;
   
    /* write out commstats summary statistics */
    envTDATA[envNXT+KEY] = COMMSTATS;
    envTDATA[envNXT+CLOCK1] = time1;
    envTDATA[envNXT+CLOCK2] = time2;
    envTDATA[envNXT+MRECV] = envRCV;
    envTDATA[envNXT+MSENT] = envSNT;
    envTDATA[envNXT+MRVOL] = envRVOL;
    envTDATA[envNXT+MSVOL] = envSVOL;
    envTDATA[envNXT+MPROB] = envPRB;
    envNXT = envNXT + COMMSTATSIZE;

    /* disable tracing on the node */
    envTOPN = -1;

    };

}
