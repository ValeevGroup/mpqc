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
#include "envfprintf.h"

/* this routine processes the trace data for a node */

envNFLUSH(tarray, size)
long *tarray, size;
{
  int k, klump, lth, bufsize, strlen();
  long i, arraysize, maxarray;
  long time1, time2, timef1, timef2, idle1, idle2, busy1, busy2;
  double hldclkoffset, hldclkdrift, hldclkref, hldclkstart;
  double dtime, didle, dbusy, dtimef, envNODETIME(), envNODEIDLE();
  static double envTST;
  char mesbuf[200], *tmpbuf;

  /* save current clock normalization values */
  hldclkoffset = envCLKOFFSET;
  hldclkdrift = envCLKDRIFT;
  hldclkref = envCLKREF;
  hldclkstart = envCLKSTART;

  k = 0;
  while (k<size) {

   switch (tarray[k+KEY] & 077) {

    case TSTART:
     envCLKSTART = 0.0;
     envCLKOFFSET = tarray[k+OFFSET1] + .000001*tarray[k+OFFSET2];
     envCLKDRIFT = tarray[k+DRIFT1] + .000001*tarray[k+DRIFT2];
     envCLKREF = tarray[k+CLKREF1] + .000001*tarray[k+CLKREF2];
     envCLKSTART = envNODETIME(tarray[k+CLKSYNC1],tarray[k+CLKSYNC2]);
     envTST = envNODETIME(tarray[k+CLOCK1],tarray[k+CLOCK2]);
     time1 =  envTST;
     time2 =  1000000*(envTST - time1);

     if (envVER == 0)
       sprintf(mesbuf,"%d %ld %ld %d %ld %ld %ld\n",TSTART,
               time1,time2,envME,tarray[k+TRACE2],tarray[k+COMPF2],
               tarray[k+COMMF2]);
       else sprintf(mesbuf,
"trace_start  clock %ld %ld node %d event %ld compstats %ld commstats %ld\n", 
                    time1,time2,envME,tarray[k+TRACE2],tarray[k+COMPF2],
                    tarray[k+COMMF2]);

     ENVFPRINTF(envTFP,"%.199s",mesbuf);
     k = k + NTSTARTSIZE;
     break;

    case OPEN:
     dtime = envNODETIME(tarray[k+CLOCK1],tarray[k+CLOCK2]);
     time1 = dtime;
     time2 = 1000000*(dtime - time1);

     if (envVER == 0) 
       sprintf(mesbuf,"%d %ld %ld %d\n",OPEN,time1,time2,envME);
       else sprintf(mesbuf,"open  clock %ld %ld node %d\n",
                    time1,time2,envME);

     ENVFPRINTF(envTFP,"%.199s",mesbuf);
     k = k + OPENSIZE;
     break;

    case SEND:
     dtime = envNODETIME(tarray[k+CLOCK1],tarray[k+CLOCK2]);
     time1 =  dtime;
     time2 =  1000000*(dtime - time1);

     if (envVER == 0) 
       sprintf(mesbuf,"%d %ld %ld %d %ld %ld %ld\n",SEND,time1,time2,
               envME,tarray[k+TO],tarray[k+TYPE],tarray[k+LENGTH]);
       else sprintf(mesbuf,
"send  clock %ld %ld node %d to %ld type %ld lth %ld\n",
                   time1,time2,envME,tarray[k+TO],tarray[k+TYPE],
                   tarray[k+LENGTH]);

     ENVFPRINTF(envTFP,"%.199s",mesbuf);
     k = k + SENDSIZE;
     break;

    case PROBE:
     dtime = envNODETIME(tarray[k+CLOCK1],tarray[k+CLOCK2]);
     time1 =  dtime;
     time2 =  1000000*(dtime - time1);

     if (envVER == 0) 
       sprintf(mesbuf,"%d %ld %ld %d %ld \n",PROBE,time1,time2,envME,
        tarray[k+TYPE]);
       else sprintf(mesbuf,"probe  clock %ld %ld node %d type %ld\n",
                    time1,time2,envME,tarray[k+TYPE]);

     ENVFPRINTF(envTFP,"%.199s",mesbuf);
     k = k + PROBESIZE;
     break;

    case RECV:
     dtime = envNODETIME(tarray[k+CLOCK1],tarray[k+CLOCK2]);
     time1 =  dtime;
     time2 =  1000000*(dtime - time1);

     if (envVER == 0) 
       sprintf(mesbuf,"%d %ld %ld %d %ld %ld %ld\n",RECV,time1,time2,
        envME,tarray[k+FROM],tarray[k+TYPE],tarray[k+LENGTH]);
       else sprintf(mesbuf,
"recv  clock %ld %ld node %d from %ld type %ld lth %ld\n",
                    time1,time2,envME,tarray[k+FROM],tarray[k+TYPE],
                    tarray[k+LENGTH]);

     ENVFPRINTF(envTFP,"%.199s",mesbuf);
     k = k + RECVSIZE;
     break;

    case BLOCK:
     dtime = envNODETIME(tarray[k+CLOCK1],tarray[k+CLOCK2]);
     time1 =  dtime;
     time2 =  1000000*(dtime - time1);

     if (envVER == 0) 
       sprintf(mesbuf,"%d %ld %ld %d %ld\n",BLOCK,time1,time2,envME,
        tarray[k+TYPE]);
       else 
       sprintf(mesbuf,"recv_blocking  clock %ld %ld node %d type %ld\n",
               time1,time2,envME,tarray[k+TYPE]);

     ENVFPRINTF(envTFP,"%.199s",mesbuf);
     k = k + BLOCKSIZE;
     break;

    case WAKE:
     dtime = envNODETIME(tarray[k+CLOCK1],tarray[k+CLOCK2]);
     time1 =  dtime;
     time2 =  1000000*(dtime - time1);

     if (envVER == 0) 
       sprintf(mesbuf,"%d %ld %ld %d %ld %ld %ld\n",WAKE,time1,time2,
               envME,tarray[k+FROM],tarray[k+TYPE],tarray[k+LENGTH]);
       else 
       sprintf(mesbuf,
"recv_waking  clock %ld %ld node %d from %ld type %ld lth %ld\n",
               time1,time2,envME,tarray[k+FROM],tarray[k+TYPE],
               tarray[k+LENGTH]);

     ENVFPRINTF(envTFP,"%.199s",mesbuf);
     k = k + WAKESIZE;
     break;

    case MESSAGE:  
     dtime = envNODETIME(tarray[k+CLOCK1],tarray[k+CLOCK2]);
     time1 =  dtime;
     time2 =  1000000*(dtime - time1);

     if (envVER == 0) 
       sprintf(mesbuf,"%d %ld %ld %d\n",MESSAGE,time1,time2,envME);
       else sprintf(mesbuf,"message  clock %ld %ld node %d\n",
                    time1,time2,envME);

     ENVFPRINTF(envTFP,"%.199s",mesbuf);
     k = k + MESSAGESIZE;
     break;

    case SYNC:  
     dtime = envNODETIME(tarray[k+CLOCK1],tarray[k+CLOCK2]);
     time1 =  dtime;
     time2 =  1000000*(dtime - time1);

     if (envVER == 0) 
       sprintf(mesbuf,"%d %ld %ld %d\n",SYNC,time1,time2,envME);
       else sprintf(mesbuf,"sync  clock %ld %ld node %d\n",
                    time1,time2,envME);

     ENVFPRINTF(envTFP,"%.199s",mesbuf);
     k = k + SYNCSIZE;
     break;

    case COMPSTATS:
     dtime = envNODETIME(tarray[k+CLOCK1],tarray[k+CLOCK2]);
     time1 = dtime;
     time2 = 1000000*(dtime - time1);

     didle = envNODEIDLE(tarray[k+IDLE1],tarray[k+IDLE2]);
     idle1 = didle;
     idle2 = 1000000*(didle - idle1);

     dbusy = dtime - didle - envTST;
     busy1 = dbusy;
     busy2 = 1000000*(dbusy - busy1);

     if (envVER == 0) 
       sprintf(mesbuf,"%d %ld %ld %d %ld %ld\n",COMPSTATS,time1,
               time2,envME,idle1,idle2);
       else 
       sprintf(mesbuf,
"compstats  clock %ld %ld node %d busy %ld %ld idle %ld %ld\n",
               time1,time2,envME,busy1,busy2,idle1,idle2);

     ENVFPRINTF(envTFP,"%.199s",mesbuf);
     k = k + COMPSTATSIZE;
     break;

    case COMMSTATS:
     dtime = envNODETIME(tarray[k+CLOCK1],tarray[k+CLOCK2]);
     time1 =  dtime;
     time2 =  1000000*(dtime - time1);

     if (envVER == 0) 
       sprintf(mesbuf,"%d %ld %ld %d %ld %ld %ld %ld %ld\n",COMMSTATS,
               time1,time2,envME,tarray[k+MRECV],tarray[k+MRVOL],
               tarray[k+MSENT],tarray[k+MSVOL],tarray[k+MPROB]);
       else 
       sprintf(mesbuf,
"%s  %s %ld %ld %s %d received %ld volume %ld sent %ld volume %ld probed %ld\n",
               "commstats","clock",time1,time2,"node",envME,tarray[k+MRECV],
               tarray[k+MRVOL],tarray[k+MSENT],tarray[k+MSVOL],
               tarray[k+MPROB]);

     ENVFPRINTF(envTFP,"%.199s",mesbuf);
     k = k + COMMSTATSIZE;
     break;

    case CLOSE:
     dtime = envNODETIME(tarray[k+CLOCK1],tarray[k+CLOCK2]);
     time1 =  dtime;
     time2 =  1000000*(dtime - time1);

     if (envVER == 0) 
       sprintf(mesbuf,"%d %ld %ld %d\n",CLOSE,time1,time2,envME);
       else 
       sprintf(mesbuf,"close  clock %ld %ld node %d\n",time1,time2,
               envME);

     ENVFPRINTF(envTFP,"%.199s",mesbuf);
     k = k + CLOSESIZE;
     break;

    case TLEVEL:
     dtime = envNODETIME(tarray[k+CLOCK1],tarray[k+CLOCK2]);
     time1 =  dtime;
     time2 =  1000000*(dtime - time1);

     if (envVER == 0) 
       sprintf(mesbuf,"%d %ld %ld %d %ld %ld %ld\n",TLEVEL,time1,time2,
               envME,tarray[k+TRACE],tarray[k+COMPF],tarray[k+COMMF]);
       else 
       sprintf(mesbuf,
"trace_level  clock %ld %ld node %d event %ld compstats %ld commstats %ld\n",
               time1,time2,envME,tarray[k+TRACE],tarray[k+COMPF],
               tarray[k+COMMF]);

     ENVFPRINTF(envTFP,"%.199s",mesbuf);
     k = k + TLEVELSIZE;
     break;

    case TMARK:  
     dtime = envNODETIME(tarray[k+CLOCK1],tarray[k+CLOCK2]);
     time1 =  dtime;
     time2 =  1000000*(dtime - time1);

     if (envVER == 0) 
       sprintf(mesbuf,"%d %ld %ld %d %ld\n",TMARK,time1,time2,envME,
               tarray[k+TYPE]);
       else 
       sprintf(mesbuf,"trace_mark  clock %ld %ld node %d type %ld\n",
               time1,time2,envME,tarray[k+TYPE]);

     ENVFPRINTF(envTFP,"%.199s",mesbuf);
     k = k + TMARKSIZE;
     break;

    case TMARKS:  
     dtime = envNODETIME(tarray[k+CLOCK1],tarray[k+CLOCK2]);
     time1 =  dtime;
     time2 =  1000000*(dtime - time1);
     arraysize = tarray[k+TYPE];

     if (arraysize < 0){
       maxarray = TMAXLONG - FLUSHTRIGGER - TMARKSIZE;
       if (envVER == 0) 
         sprintf(mesbuf,"%d %ld %ld %d %ld\n",TMARKS,time1,time2,envME,
                 arraysize);
         else 
         sprintf(mesbuf,
"trace_mark_array  clock %ld %ld node %d integer array of size %ld %s %d\n",
                 time1,time2,envME,-arraysize,
                 "too large: can be no larger than",maxarray);
       ENVFPRINTF(envTFP,"%.199s",mesbuf);
       k = k + TMARKSIZE;
       }
       else{
       if (envVER == 0){

         bufsize = 11*arraysize + 60;
         if (bufsize > 195){ /* allocate a buffer large enough for mark array */
           tmpbuf = malloc(bufsize);
           if (tmpbuf == NULL){ /* can't get a large enough array - punt */
             arraysize = -arraysize;
             sprintf(mesbuf,"%d %ld %ld %d %ld\n",TMARKS,time1,time2,envME,
                     arraysize);
             ENVFPRINTF(envTFP,"%.199s",mesbuf);
             }
             else{ /* fill buffer and send it off */
             sprintf(tmpbuf,"%d %ld %ld %d %ld",TMARKS,time1,time2,envME,
                     arraysize);
             lth = strlen(mesbuf)-1;
             for (i=1; i<=arraysize; i++){
               sprintf(&tmpbuf[lth]," %ld",tarray[k+TYPE+i]);
               lth = strlen(mesbuf)-1;
               };
             sprintf(&tmpbuf[lth],"\n");
             ENVFPRINTF(envTFP,"%s",tmpbuf);
             free(tmpbuf);
             };
           }
           else{ /* default buffer is big enough */
           sprintf(mesbuf,"%d %ld %ld %d %ld",TMARKS,time1,time2,envME,
                   arraysize);
           lth = strlen(mesbuf)-1;
           for (i=1; i<=arraysize; i++){
             sprintf(&mesbuf[lth]," %ld",tarray[k+TYPE+i]);
             lth = strlen(mesbuf)-1;
             };
           sprintf(&mesbuf[lth],"\n");
           ENVFPRINTF(envTFP,"%.199s",mesbuf);
           };

         }
         else{

         bufsize = 11*arraysize + 90;
         if (bufsize > 195){ /* allocate a buffer large enough for mark array */
           tmpbuf = malloc(bufsize);
           if (tmpbuf == NULL){ /* can't get a large enough array - punt */
             arraysize = -arraysize;
             sprintf(mesbuf, 
"trace_mark_array  clock %ld %ld node %d %s %ld\n",time1,time2,envME,
"failed trying to allocate buffer space for a tracemarks record of size",
arraysize);
             ENVFPRINTF(envTFP,"%.199s",mesbuf);
             }
             else{ /* fill buffer and send it off */
             sprintf(tmpbuf,
"trace_mark_array  clock %ld %ld node %d size %ld array",
                   time1,time2,envME,arraysize);
             lth = strlen(mesbuf)-1;
             for (i=1; i<=arraysize; i++){
               sprintf(&tmpbuf[lth]," %ld",tarray[k+TYPE+i]);
               lth = strlen(mesbuf)-1;
               };
             sprintf(&tmpbuf[lth],"\n");
             ENVFPRINTF(envTFP,"%s",tmpbuf);
             free(tmpbuf);
             };
           }
           else{ /* default buffer is big enough */
           sprintf(mesbuf,
"trace_mark_array  clock %ld %ld node %d size %ld array",
                   time1,time2,envME,arraysize);
           lth = strlen(mesbuf)-1;
           for (i=1; i<=arraysize; i++){
             sprintf(&mesbuf[lth]," %ld",tarray[k+TYPE+i]);
             lth = strlen(mesbuf)-1;
             };
           sprintf(&mesbuf[lth],"\n");
           ENVFPRINTF(envTFP,"%.199s",mesbuf);
           };

         };

       k = k + TMARKSIZE + arraysize;
       };
     break;

    case TFULL:  
     dtime = envNODETIME(tarray[k+CLOCK1],tarray[k+CLOCK2]);
     time1 =  dtime;
     time2 =  1000000*(dtime - time1);

     if (envVER == 0) 
       sprintf(mesbuf,"%d %ld %ld %d\n",TFULL,time1,time2,envME);
       else 
       sprintf(mesbuf,
"trace_stop  clock %ld %ld node %d : ran out of memory\n",
               time1,time2,envME);

     ENVFPRINTF(envTFP,"%.199s",mesbuf);
     k = k + TFULLSIZE;
     break;

    case TFLUSH:  
     dtime = envNODETIME(tarray[k+CLOCK1],tarray[k+CLOCK2]);
     time1 =  dtime;
     time2 =  1000000*(dtime - time1);

     dtimef = envNODETIME(tarray[k+CLOCKSTOP1],tarray[k+CLOCKSTOP2]);
     timef1 =  dtimef;
     timef2 =  1000000*(dtimef - timef1);

     if (envVER == 0) 
       sprintf(mesbuf,"%d %ld %ld %d %ld %ld\n",TFLUSH,
               time1,time2,envME,timef1,timef2);
       else 
       sprintf(mesbuf,
"trace_flush  clock %ld %ld node %d finished %ld %ld\n",
               time1,time2,envME,timef1,timef2);

     ENVFPRINTF(envTFP,"%.199s",mesbuf);
     k = k + TFLUSHSIZE;
     break;

    case TFILES:  
     dtime = envNODETIME(tarray[k+CLOCK1],tarray[k+CLOCK2]);
     time1 =  dtime;
     time2 =  1000000*(dtime - time1);

     dtimef = envNODETIME(tarray[k+CLOCKSTOP1],tarray[k+CLOCKSTOP2]);
     timef1 =  dtimef;
     timef2 =  1000000*(dtimef - timef1);

     if (envVER == 0) 
       sprintf(mesbuf,"%d %ld %ld %d %ld %ld\n",TFILES,
               time1,time2,envME,timef1,timef2);
       else 
       sprintf(mesbuf,
"trace_files  clock %ld %ld node %d finished %ld %ld\n",
               time1,time2,envME,timef1,timef2);

     ENVFPRINTF(envTFP,"%.199s",mesbuf);
     k = k + TFILESIZE;
     break;

    case TEXIT:
     dtime = envNODETIME(tarray[k+CLOCK1],tarray[k+CLOCK2]);
     time1 =  dtime;
     time2 =  1000000*(dtime - time1);

     if (envVER == 0) 
       sprintf(mesbuf,"%d %ld %ld %d %ld\n",TEXIT,time1,time2,envME,
               tarray[k+OSPACE]);
       else 
       sprintf(mesbuf,"trace_exit  clock %ld %ld node %d space %ld\n",
               time1,time2,envME,tarray[k+OSPACE]);

     ENVFPRINTF(envTFP,"%.199s",mesbuf);
     k = k + TEXITSIZE;
     break;

    case TBLOCKBEGIN:
     dtime = envNODETIME(tarray[k+CLOCK1],tarray[k+CLOCK2]);
     time1 =  dtime;
     time2 =  1000000*(dtime - time1);

     if (envVER == 0) 
       sprintf(mesbuf,"%d %ld %ld %d %ld %ld %ld\n",TBLOCKBEGIN,time1,
               time2,envME,tarray[k+TYPE],tarray[k+LOCATION],
               tarray[k+PARAM]);
       else
       sprintf(mesbuf,
"%s  %s %ld %ld node %d block_type %ld location_type %ld parameter_type %ld\n",
               "block_begin","clock",time1,time2,envME,tarray[k+TYPE],
               tarray[k+LOCATION],tarray[k+PARAM]);

     ENVFPRINTF(envTFP,"%.199s",mesbuf);
     k = k + TBLOCKSIZE;
     break;

    case TBLOCKEND:
     dtime = envNODETIME(tarray[k+CLOCK1],tarray[k+CLOCK2]);
     time1 =  dtime;
     time2 =  1000000*(dtime - time1);

     if (envVER == 0) 
       sprintf(mesbuf,"%d %ld %ld %d %ld %ld %ld\n",TBLOCKEND,time1,
               time2,envME,tarray[k+TYPE],tarray[k+LOCATION],
               tarray[k+PARAM]);
       else
       sprintf(mesbuf,
"%s  %s %ld %ld node %d block_type %ld location_type %ld parameter_type %ld\n",
               "block_end","clock",time1,time2,envME,tarray[k+TYPE],
               tarray[k+LOCATION],tarray[k+PARAM]);

     ENVFPRINTF(envTFP,"%.199s",mesbuf);
     k = k + TBLOCKSIZE;
     break;

    case CLOCKSYNC:  
     envCLKSTART = 0.0;
     envCLKOFFSET = tarray[k+OFFSET1] + .000001*tarray[k+OFFSET2];
     envCLKDRIFT = tarray[k+DRIFT1] + .000001*tarray[k+DRIFT2];
     envCLKREF = tarray[k+CLKREF1] + .000001*tarray[k+CLKREF2];
     envCLKSTART = 
      envNODETIME(tarray[k+CLKSYNC1],tarray[k+CLKSYNC2]);
     dtime = envNODETIME(tarray[k+CLOCK1],tarray[k+CLOCK2]);
     time1 =  dtime;
     time2 =  1000000*(dtime - time1);

     if (envVER == 0) 
       sprintf(mesbuf,"%d %ld %ld %d\n",CLOCKSYNC,
               time1,time2,envME);
       else 
       sprintf(mesbuf,"clock_sync  clock %ld %ld node %d\n",
               time1,time2,envME);

     ENVFPRINTF(envTFP,"%.199s",mesbuf);
     k = k + CLKSYNCSIZE;
     break;

    default:
     sprintf(mesbuf,
"%d is not a legal record type: raw trace data for node %d follows\n",
             tarray[k+KEY],envME);
     ENVFPRINTF(envTFP,"%.199s",mesbuf);
     sprintf(mesbuf,"(error identified at location k=%d)\n",k);
     ENVFPRINTF(envTFP,"%.199s",mesbuf);
     klump = 9*(size/9);
     for (k=0; k<klump; k=k+9){ 
       sprintf(mesbuf,"%ld %ld %ld %ld %ld %ld %ld %ld %ld %ld\n", 
               tarray[k], tarray[k+1], tarray[k+2], tarray[k+3], tarray[k+4], 
               tarray[k+5], tarray[k+6], tarray[k+7], tarray[k+8]);
       ENVFPRINTF(envTFP,"%.199s",mesbuf);
       };
     lth=0;
     for (; k<size; k++){
       sprintf(&mesbuf[lth], "%ld ", tarray[k]);
       lth = strlen(mesbuf);
       };
     sprintf(&mesbuf[lth], "\n");
     ENVFPRINTF(envTFP,"%.199s",mesbuf);

    };

   };

  fflush(envTFP);

  /* restore current clock normalization values */
  envCLKOFFSET = hldclkoffset;
  envCLKDRIFT = hldclkdrift;
  envCLKREF = hldclkref;
  envCLKSTART = hldclkstart;

}

