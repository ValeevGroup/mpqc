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
#include "envfflush.h"
#include "envfprintf.h"
#include "envopen.h"
#include "envsend.h"
#include "envupdatetime.h"
#include "envclock.h"

/* send a message of type TRACE_MESG back to the trace file. */
/* the message should be a string of length less than 80 characters. */

tracemsg(message)
char *message;
{
  int i, j, lth, lth2, strlen();
  long time1, time2, start1, start2, stop1, stop2;
  double dtime, clock0();
  char mesbuf[120];

  /* check whether tracing enabled */
  if (envTOPN == 1){

    /* record start time for tracemsg */
    dtime = clock0();

    /* open communication (if necessary) */
    if (envMOPN == 0){

      /* record start time */
      ENVCLOCK(&start1, &start2);

      /* open communication channel */
      ENVOPEN();
      envMOPN = 1;

      /* record stop time */
      ENVCLOCK(&stop1, &stop2);

      /* update running total */
      ENVUPDATETIME(&envCTIM1, &envCTIM2, start1, start2, stop1, stop2);

      };

    /* tracing enabled, so prepare a message for the trace file */         
    time1 = dtime;
    time2 = 1000000*(dtime - time1);

    /* (note: if no local trace file, envVER == 0) */
    if (envVER == 0)
      sprintf(mesbuf,"%d %ld %ld %d ", TMSG, time1, time2, envME);
      else
      sprintf(mesbuf,"trace_message clock %ld %ld node %d ",time1,time2, envME);
    lth = strlen(mesbuf);
    lth2 = lth + strlen(message);
    if (lth2 > 118) lth2 = 118;
    j=0;
    for (i=lth; i<lth2; i++, j++) mesbuf[i] = message[j];
    mesbuf[i] = '\n'; 
    mesbuf[i+1] = '\0';

    if (envTFP != NULL){ /* if local trace file, send it there */
      ENVFPRINTF(envTFP,"%s",mesbuf);
      ENVFFLUSH(envTFP);
      }
      else{
      if (envHOSTED == 1){ /* else if hosted, send it to host trace file */
        lth = strlen(mesbuf)+1;
        ENVSEND(mesbuf, lth, TRACE_MESG, HOST0, 0);
        };
      };

    };
}
