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

/* either stop tracing, and record this fact, or flush the node trace */
/* array */

envFLUSH(time1, time2)
long time1, time2;
{

  /* can we add another packet? */
  if (envNPTR == envSPTR){

    /* no; therefore, either flush or stop */
    if (envFLU == 1) traceflush();
      else{

      envTDATA[envNXT+KEY] = TFULL;
      envTDATA[envNXT+CLOCK1] = time1;
      envTDATA[envNXT+CLOCK2] = time2;
      envNXT = envNXT + TFULLSIZE;

      /* save tracing parameters */
      envSTRA = envTRA;
      envSTCP = envTCP;
      envSTCM = envTCM;

      /* disable tracing */
      envTRA = 0;
      envTCP = 0;
      envTCM = 0;

      /* mark trace array as full */
      envFUL = 1;

      };

    }
    else{

    /* save beginning of next packet */
    envTPTR[envNPTR] = envNXT;
    envNPTR++;

    /* calculate the new limit on envNXT */
    envMLTH = envNXT + TMAXLONG;
    if (envMLTH > envTSZ) envMLTH = envTSZ;
    envMLTH = envMLTH - FLUSHTRIGGER;

    };

}

