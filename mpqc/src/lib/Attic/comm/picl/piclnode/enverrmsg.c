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
#include "envclose.h"
#include "envfflush.h"
#include "envfprintf.h"
#include "envopen.h"
#include "envprintf.h"
#include "envsend.h"

/* Flush and disable tracing, send messages of types TRACE_MESG and MSG0 */
/* to a trace file and standard out, respectively, and exit. If there is */
/* a host, then the MSG0 message goes via the host, else it goes directly */
/* to standard out. */
/* (All this happens whether a communication channel has been opened or not).*/ 
/* The message should be a string of length less than 80 characters. */

envERRMSG(message)
char *message;
{
  int i, j, lth, lth2, strlen();
  long time1, time2;
  double dtime, clock0();
  char mesbuf[120];

  /* open communication (if necessary) */
  if (envMOPN == 0){
    ENVOPEN();
    envMOPN = 1;
    };

  /* record time */
  dtime = clock0();
  time1 = dtime;
  time2 = 1000000*(dtime - time1);

  /* send message to standard out */
  sprintf(mesbuf, "message  node %d clock %ld %ld ", envME, time1, time2);
  lth = strlen(mesbuf);
  lth2 = lth + strlen(message);
  if (lth2 > 118) lth2 = 118;
  j=0;
  for (i=lth; i<lth2; i++, j++) mesbuf[i] = message[j];
  mesbuf[i] = '\n'; 
  mesbuf[i+1] = '\0';

  if (envHOSTED == 1){ /* if hosted, send to host to print on standard out */
    lth = strlen(mesbuf)+1;
    ENVSEND(mesbuf, lth, MSG0, HOST0, 0);
    }
    else{ /* else, send to standard out directly */
    ENVPRINTF("%s",mesbuf);
    };

  /* send message to trace file */
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

  /* flush and close tracing */
  traceexit();
  traceflush();

  /* exit */
  ENVCLOSE();
  exit(1);

}
