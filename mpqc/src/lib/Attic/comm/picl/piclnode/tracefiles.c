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
#include "envfclose.h"
#include "envfflush.h"
#include "envfopen.h"
#include "envfread.h"
#include "envfwrite.h"
#include "envprintf.h"
#include "envrewind.h"
#include "envopen.h"
#include "envsend.h"
#include "envunlink.h"
#include "envupdatetime.h"

/* node trace file specification routine */
/* used for specifying temporary disk storage for trace data, */
/*   and for specifying a final trace file distinct from that */
/*   specified by the host */
/* necessary for tracing on a hostless system */
/* useful for tracing on a system with disks attached directly to the */
/*   multiprocessor
/* - tempfile is the prefix (including directory) of the name of the disk */
/*            to be used for temporary storage of trace data. A suffix (the */
/*            node number) is tacked on to make all temporary files unique. */
/* - permfile is the name of the disk where this node's trace data should be */
/*            sent for "permanent" storage. */
/* - verbose == 0, compressed trace records are output */
/*           != 0, normal trace records are output */

tracefiles(tempfile, permfile, verbose)
char *tempfile, *permfile;
int verbose;
{
  int i, lth, times, strlen();
  long begloc, size, index, start1, start2, stop1, stop2;
  char *temptemp, *malloc();

  /* check whether tracing is enabled */
  if (envTOPN == 1){

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

    /* potentially changing tracefile and tempfile, so write out */
    /* everything that has logically been flushed, but not yet sent */
    /* back to a permanent tracefile */

    if (envTMPBUFS > 0){ /* some trace data in temporary storage */

      envTPTR[envNPTR] = envNXT;
      for (i=1; i<=envNPTR; i++){ /* first, put all trace data there */

        begloc = envTPTR[i-1];
        size = (envTPTR[i] - begloc);
        ENVFWRITE(&size, sizeof(long), 1, envTMPW);
        ENVFWRITE(&envTDATA[begloc], sizeof(long), size, envTMPW);
        envTMPBUFS++;

        };

      ENVFFLUSH(envTMPW);
      ENVREWIND(envTMPR);

      while (envTMPBUFS > envNPTR){ /* then bring it back one buffer at a time*/

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

      index = 0;
      while (envTMPBUFS > 0){ /* restore trace buffer to original condition */

        ENVFREAD(&size, sizeof(long), 1, envTMPR);
        ENVFREAD(&envTDATA[index], sizeof(long), size, envTMPR);
        index += size;
        envTMPBUFS--;

        };

      };

    /* close current permanent file */
    if (envTFP != NULL){
      ENVFCLOSE(envTFP);
      };

    /* if requested, set up new permanent file */
    if (permfile != NULL){

      /* try ten times */
      envTFP = NULL;
      times = 10;
      while ((envTFP == NULL) && (times > 0)){
        ENVFOPEN(&envTFP, permfile, "a");
        times--;
        };

      if (envTFP != NULL) envVER = verbose;
        else{
        envVER = 0;
        sprintf(envMBUF,"tracefiles: could not open %-.50s\n",permfile);
        if (envHOSTED == 1){ 
          /* if hosted, send to host to print on standard out */
          lth = strlen(envMBUF)+1;
          ENVSEND(envMBUF, lth, MSG0, HOST0, 0);
          }
          else{ /* else, send to standard out directly */
          ENVPRINTF("%s",envMBUF);
          };
        };

      }
      else envVER = 0;

    /* close current temporary file */
    if (envTMPR != NULL){
      ENVFCLOSE(envTMPR);
      };

    if (envTMPW != NULL){
      ENVFCLOSE(envTMPW);
      };

    /* if requested, set up new temporary file */
    if (tempfile != NULL){

      lth = strlen(tempfile);
      temptemp = malloc0(lth+7);
      sprintf(temptemp, "%s%-.6d", tempfile, envME);
 
      /* try ten times */
      envTMPW = NULL;
      times = 10;
      while ((envTMPW == NULL) && (times > 0)){
        ENVFOPEN(&envTMPW, temptemp, "w");
        times--;
        };

      /* try ten times */
      envTMPR = NULL;
      if (envTMPW != NULL){
        times = 10;
        while ((envTMPR == NULL) && (times > 0)){
          ENVFOPEN(&envTMPR, temptemp, "r");
          times--;
          };
        ENVUNLINK(temptemp);
        };

      if ((envTMPR == NULL) || (envTMPW == NULL)){
        envTMPR = NULL;
        envTMPW = NULL;
        sprintf(envMBUF,"tracefiles: could not open %-.50s\n",temptemp);
        tracemsg(envMBUF);
        };

      free(temptemp);
      };

    /* record stop time */
    ENVCLOCK(&stop1, &stop2);

    /* record changing of trace files */
    envTDATA[envNXT+KEY] = TFILES;
    envTDATA[envNXT+CLOCKSTART1] = start1;
    envTDATA[envNXT+CLOCKSTART2] = start2;
    envTDATA[envNXT+CLOCKSTOP1] = stop1;
    envTDATA[envNXT+CLOCKSTOP2] = stop2;

    envNXT = envNXT + TFILESIZE;

    if (envNXT > envMLTH) envFLUSH(stop1, stop2);
      
    }
    else{

    /* open communication (if necessary) */
    if (envMOPN == 0){

      /* open communication channel */
      ENVOPEN();
      envMOPN = 1;

      };

    /* potentially changing tracefile and tempfile, so write out */
    /* everything that has logically been flushed, but not yet sent */
    /* back to a permanent tracefile */

    if (envTMPBUFS > 0){ /* some trace data in temporary storage */

      envTPTR[envNPTR] = envNXT;
      for (i=1; i<=envNPTR; i++){ /* first, put all trace data there */

        begloc = envTPTR[i-1];
        size = (envTPTR[i] - begloc);
        ENVFWRITE(&size, sizeof(long), 1, envTMPW);
        ENVFWRITE(&envTDATA[begloc], sizeof(long), size, envTMPW);
        envTMPBUFS++;

        };

      ENVFFLUSH(envTMPW);
      ENVREWIND(envTMPR);

      while (envTMPBUFS > envNPTR){ /* then bring it back one buffer at a time*/

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

      index = 0;
      while (envTMPBUFS > 0){ /* restore trace buffer to original condition */

        ENVFREAD(&size, sizeof(long), 1, envTMPR);
        ENVFREAD(&envTDATA[index], sizeof(long), size, envTMPR);
        index += size;
        envTMPBUFS--;

        };

      };

    /* close current permanent file */
    if (envTFP != NULL){
      ENVFCLOSE(envTFP);
      };

    /* if requested, set up new permanent file */
    if (permfile != NULL){

      /* try ten times */
      envTFP = NULL;
      times = 10;
      while ((envTFP == NULL) && (times > 0)){
        ENVFOPEN(&envTFP, permfile, "a");
        times--;
        };

      if (envTFP != NULL) envVER = verbose;
        else{
        sprintf(envMBUF,"tracefiles: could not open %-.50s\n",permfile);
        if (envHOSTED == 1){ 
          /* if hosted, send to host to print on standard out */
          lth = strlen(envMBUF)+1;
          ENVSEND(envMBUF, lth, MSG0, HOST0, 0);
          }
          else{ /* else, send to standard out directly */
          ENVPRINTF("%s",envMBUF);
          };
        envVER = 0;
        };

      }
      else envVER = 0;

    /* close current temporary file */
    if (envTMPR != NULL){
      ENVFCLOSE(envTMPR);
      };

    if (envTMPW != NULL){
      ENVFCLOSE(envTMPW);
      };

    /* if requested, set up new temporary file */
    if (tempfile != NULL){

      lth = strlen(tempfile);
      temptemp = malloc0(lth+7);
      sprintf(temptemp, "%s%-.6d", tempfile, envME);
 
      /* try ten times */
      envTMPW = NULL;
      times = 10;
      while ((envTMPW == NULL) && (times > 0)){
        ENVFOPEN(&envTMPW, temptemp, "w");
        times--;
        };

      /* try ten times */
      envTMPR = NULL;
      if (envTMPW != NULL){
        times = 10;
        while ((envTMPR == NULL) && (times > 0)){
          ENVFOPEN(&envTMPR, temptemp, "r");
          times--;
          };
        ENVUNLINK(temptemp);
        };

      if ((envTMPR == NULL) || (envTMPW == NULL)){
        envTMPR = NULL;
        envTMPW = NULL;
        sprintf(envMBUF,"tracefiles: could not open %-.50s\n",temptemp);
        tracemsg(envMBUF);
        };

      free(temptemp);
      };

    };

}

