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
#include "envfflush.h"
#include "envfread.h"
#include "envfwrite.h"
#include "envupdatetime.h"
#include "envopen.h"
#include "envrewind.h"
#include "envsend.h"

/* send trace data back to the host and flush the data array */

traceflush()
{
  int size, j, i;
  long begloc, start1, start2, stop1, stop2;

  /* check whether tracing enabled */
  if (envTOPN != 0){

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

      }
      else{

      /* record start time */
      ENVCLOCK(&start1, &start2);

      };

    /* record end of last trace record */
    envTPTR[envNPTR] = envNXT;

    if (envTOPN == 1){ /* not final flush */

      if (envTMPW != NULL){ /* send data to temp file */

        for (i=1; i<=envNPTR; i++){

          begloc = envTPTR[i-1];
          size = (envTPTR[i] - begloc);
          ENVFWRITE(&size, sizeof(long), 1, envTMPW);
          ENVFWRITE(&envTDATA[begloc], sizeof(long), size, envTMPW);
          envTMPBUFS++;

          };

        }
        else{

        if (envTFP != NULL){ /* send trace data to TFP */

          for (i=1; i<=envNPTR; i++){

            begloc = envTPTR[i-1];
            size = (envTPTR[i] - begloc);
            envNFLUSH(&envTDATA[begloc], size);

            };

          }
          else{

          if (envHOSTED == 1){ /* send trace data back to host */

            for (i=1; i<=envNPTR; i++){

              begloc = envTPTR[i-1];
              size = (envTPTR[i] - begloc)*sizeof(long);
              ENVSEND(&envTDATA[begloc], size, TRACE_ARRAY, HOST0, 0);

              };

            };

          };

        };

      /* flush the trace array */
      envTPTR[0] = 0;
      envNPTR = 1;
      envMLTH = min(TMAXLONG,envTSZ) - FLUSHTRIGGER;

      /* record "finished communication" time of traceflush */
      /* and update overhead info */
      ENVCLOCK(&stop1, &stop2);
      envOSPC = envOSPC + envNXT*sizeof(long);

      /* record reinitialization of trace array */
      envTDATA[KEY] = TFLUSH;
      envTDATA[CLOCKSTART1] = start1;
      envTDATA[CLOCKSTART2] = start2;
      envTDATA[CLOCKSTOP1] = stop1;
      envTDATA[CLOCKSTOP2] = stop2;

      envNXT = TFLUSHSIZE;

      }
      else{ /* final flush before disabling tracing */

      if (envTMPBUFS > 0){ /* some trace data in temporary storage */

        for (i=1; i<=envNPTR; i++){ /* first, put all trace data there */

          begloc = envTPTR[i-1];
          size = (envTPTR[i] - begloc);
          ENVFWRITE(&size, sizeof(long), 1, envTMPW);
          ENVFWRITE(&envTDATA[begloc], sizeof(long), size, envTMPW);
          envTMPBUFS++;

          };

        ENVFFLUSH(envTMPW);
        ENVREWIND(envTMPR);

        while (envTMPBUFS > 0){ /* then bring it back one buffer at a time */

          ENVFREAD(&size, sizeof(long), 1, envTMPR);
          ENVFREAD(envTDATA, sizeof(long), size, envTMPR);
          envTMPBUFS--;

          /* and try to send it to a trace file */
          if (envTFP != NULL) envNFLUSH(envTDATA, size); 
            else{ 
            if (envHOSTED == 1){ 
              ENVSEND(envTDATA, size*sizeof(long), TRACE_ARRAY, HOST0, 0);
              };

            };

          };

        }
        else{ /* no trace data in temporary storage, so process normally */

        if (envTFP != NULL){ /* send trace data to TFP */

          for (i=1; i<=envNPTR; i++){

            begloc = envTPTR[i-1];
            size = (envTPTR[i] - begloc);
            envNFLUSH(&envTDATA[begloc], size);

            };

          }
          else{

          if (envHOSTED == 1){ /* send trace data back to host */

            for (i=1; i<=envNPTR; i++){

              begloc = envTPTR[i-1];
              size = (envTPTR[i] - begloc)*sizeof(long);
              ENVSEND(&envTDATA[begloc], size, TRACE_ARRAY, HOST0, 0);

              };

            };

          };

        };

      /* really close tracing now */
      envTOPN = 0;
      free(envTDATA);
      free(envTPTR);

      };

    /* restore the trace variables and mark the array as not being full */
    if (envFUL == 1){

      envTRA = envSTRA;
      envTCP = envSTCP;
      envTCM = envSTCM;
      envFUL = 0;

      };

    };

}

