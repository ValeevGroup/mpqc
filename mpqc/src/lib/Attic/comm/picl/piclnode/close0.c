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
#include "envclose.h"
#include "envrecv.h"
#include "envsend.h"
#include "envsync.h"

/* node version of close communication channels. */

close0()
{
  int tmp;
  long tracing, time1, time2;

  /* check whether open statement has been executed */
  if (envOPN == 1){

    /* check whether tracing disabled */
    if (envTOPN != 0){

      /* output close trace record and exit tracing */

      ENVCLOCK(&time1, &time2);

      envTDATA[envNXT+KEY] = CLOSE;
      envTDATA[envNXT+CLOCK1] = time1;
      envTDATA[envNXT+CLOCK2] = time2;
      envNXT = envNXT + CLOSESIZE; 

      traceexit();

      /* finish tracing before the communication channel is closed, */
      /* but wait for prompt from host or neighboring node, so that  */
      /* traceflush'ing doesn't effect the behavior of the program */

      if (envHOSTED == 1){

        /* tell host that I have trace data to transmit */
        tracing = -1;
        ENVSEND(&tracing, sizeof(long), CLOSE0, HOST0, 0);

        /* wait until host tells me to send it */
        ENVRECV(&tracing, sizeof(long), TRACE_ARRAY, 0);

        /* process the trace information and disable tracing */
        if (tracing == 1) traceflush();

        /* tell host that done flushing (if it doesn't already know) */
        if (envTFP != NULL){
          tracing = 0;
          ENVSEND(&tracing, sizeof(long), CLOSE0, HOST0, 0);
          };

        }
        else{

        /* wait until everyone done, and take turns flushing trace data */
        ENVSYNC();
        if (envME != 0){
          ENVRECV(&tracing, sizeof(long), TRACE_ARRAY, 0);
          };
        
        /* process the trace information and disable tracing */
        traceflush();

        if (envME != envNPA-1){
          ENVSEND(&tracing, sizeof(long), TRACE_ARRAY, envME+1, 0);
          };

        };

      /* close the communication channel */

      ENVCLOSE();
      envOPN = -1;

      }
      else{

      if (envHOSTED == 1){

        /* tell host that I am closing the communication channel */
        tracing = 0;
        ENVSEND(&tracing, sizeof(long), CLOSE0, HOST0, 0);

        }
        else{

        /* wait until everyone done */
        ENVSYNC();
        if (envME != 0){
          ENVRECV(&tracing, sizeof(long), TRACE_ARRAY, 0);
          };
        if (envME != envNPA-1){
          ENVSEND(&tracing, sizeof(long), TRACE_ARRAY, envME+1, 0);
          };

        };

      /* close the communication channel */

      ENVCLOSE();
      envOPN = -1;

      }

    }
    else envERRMSG("close0: communication channel not open - exiting");

}
