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

#include "globaldefs.h"
#include "envglobals.h"
#include "envclock.h"
#include "envopen.h"
#include "envupdatetime.h"

/* node version of open communication channels.  */

open0(numproc, me, host)
int *numproc, *me, *host;
{
  long oldctim1, oldctim2, start1, start2, stop1, stop2;

  if (envOPN == 0){

    /* record trace information */
    if (envTOPN == 1){

      /* save old idle time values */
      oldctim1 = envCTIM1;
      oldctim2 = envCTIM2;

      /* record start time */
      ENVCLOCK(&start1, &start2);

      /* open communication channel */
      if (envMOPN == 0){
        ENVOPEN();
        envMOPN = 1;
        };
      envOPN = 1;

      /* return requested values and set globals if not already set */
      *numproc = envNPA;
      *me = envME;
      *host = HOST0;

      /* record stop time */
      ENVCLOCK(&stop1, &stop2);

      /* update running total */
      ENVUPDATETIME(&envCTIM1, &envCTIM2, start1, start2, stop1, stop2);

      /* output open trace record */
      envTDATA[envNXT+KEY] = OPEN;
      envTDATA[envNXT+CLOCK1] = start1;
      envTDATA[envNXT+CLOCK2] = start2;
      envNXT = envNXT + OPENSIZE; 

      if (envNXT > envMLTH) envFLUSH(start1, start2);

      if (envTCP > 0){
    
        /* write out info preceding open*/
        envTDATA[envNXT+KEY] = COMPSTATS;
        envTDATA[envNXT+CLOCK1] = start1;
        envTDATA[envNXT+CLOCK2] = start2;
        envTDATA[envNXT+IDLE1] = oldctim1;
        envTDATA[envNXT+IDLE2] = oldctim2;
        envNXT = envNXT + COMPSTATSIZE;
       
        /* write out info following open */
        envTDATA[envNXT+KEY] = COMPSTATS;
        envTDATA[envNXT+CLOCK1] = stop1;
        envTDATA[envNXT+CLOCK2] = stop2;
        envTDATA[envNXT+IDLE1] = envCTIM1;
        envTDATA[envNXT+IDLE2] = envCTIM2;
        envNXT = envNXT + COMPSTATSIZE;
      
        if (envNXT > envMLTH) envFLUSH(stop1, stop2);
     
        };
    
      if (envTCM > 0){
  
        envTDATA[envNXT+KEY] = COMMSTATS;
        envTDATA[envNXT+CLOCK1] = start1;
        envTDATA[envNXT+CLOCK2] = start2;
        envTDATA[envNXT+MRECV] = envRCV;
        envTDATA[envNXT+MSENT] = envSNT;
        envTDATA[envNXT+MRVOL] = envRVOL;
        envTDATA[envNXT+MSVOL] = envSVOL;
        envTDATA[envNXT+MPROB] = envPRB;
        envNXT = envNXT + COMMSTATSIZE;

        if (envNXT > envMLTH) envFLUSH(start1, start2);

        };

      }
      else{

      /* open communication channel */
      if (envMOPN == 0){
        ENVOPEN();
        envMOPN = 1;
        };
      envOPN = 1;

      /* return requested values and set globals if not already set */
      *numproc = envNPA;
      *me = envME;
      *host = HOST0;

      };

    }
    else envERRMSG("open0: communication channel already open - exiting");

}

