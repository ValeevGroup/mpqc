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

/* globaldefs.h for nodelib */

/********************** standard definitions and declarations ***************/

#include <stdio.h>
#include <math.h>

#define max(x,y) ((x)>(y)?(x):(y))
#define min(x,y) ((x)<(y)?(x):(y))

char *malloc(), *malloc0() ;

/********************** environmental definitions and variables *************/

         /**************** constant definitions *********************/

                       /* PICL host ID */
#define HOST0  -32768
                       /* number of repetitions to use when sync'ing clocks */
                       /* in clocksync0 */
#define REPS    4

         /***************** variable definitions ******************/

                      /* switch indicating whether open has been called */
int envOPN = 0 ;
                      /* switch indicating whether envopen has been called */
int envMOPN = 0 ;
                      /* node ID */
int envME;                                               
                      /* real host ID */
int envHST;
                      /* switch indicating whether there is a host */
int envHOSTED = -1;
                      /* number of processors allocated */
int envNPA;
                      /* source of last message received when tracing */
int envSRC;
                      /* type of last message received when tracing */
int envTYP;
                      /* number of bytes in last message received when */
                      /* tracing */
int envLTH;
                      /* buffer for error and trace messages */
char envMBUF[80];
                      /* switch indicating whether to normalize */
                      /* clock0 using clocksync0 parameters */
int envCLKSYNC = 0;
                      /* drift of local clock from node 0's clock */
double envCLKDRIFT = 0.0;
                      /* offset of local clock from node 0's clock */
                      /* at reference time */
double envCLKOFFSET = 0.0;
                      /* local reference time when clock drift was */
                      /* calculated */
double envCLKREF = 0.0;
                      /* broadcast plus sync time - used for clocksyncing */
double envWAIT = 0.0;
                      /* time clock started for clock0 and tracing */
double envCLKSTART = 0.0;
long envCLKSTART1 = 0;
long envCLKSTART2 = 0;

/************************ trace definitions and variables ******************/

         /**************** constant definitions *********************/

                       /* key of a "trace start" trace record */
#define TSTART  1
                       /* size of a node "trace start" trace record */
#define NTSTARTSIZE  14
                       /* key of an "open" trace record */
#define OPEN  2
                       /* size of an "open" trace record */
#define OPENSIZE  4
                       /* key of a "load" trace record */
#define LOAD  3
                       /* size of a "load" trace record */
#define LOADSIZE  4
                       /* key of a "send" trace record */
#define SEND  4
                       /* size of a "send" trace record */
#define SENDSIZE  6
                       /* key of a "probe" trace record */
#define PROBE  5
                       /* size of a "probe" trace record */
#define PROBESIZE  4
                       /* key of a "receive" trace record */
#define RECV  6
                       /* size of a "receive" trace record */
#define RECVSIZE  6
                       /* key of a "blocked" trace record */
#define BLOCK  7
                       /* size of a "blocked" trace record */
#define BLOCKSIZE  4
                       /* key of a "waking" trace record */
#define WAKE  8
                       /* size of a "waking" trace record */
#define WAKESIZE  6
                       /* key of a "message" trace record */
#define MESSAGE  9
                       /* size of a "message" trace record */
#define MESSAGESIZE  3
                       /* key of a "synchronization" trace record */
#define SYNC 10
                       /* size of a "synchronization" trace record */
#define SYNCSIZE 3
                       /* key of a "busy/idle time" trace record */
#define COMPSTATS  11
                       /* size of a "busy/idle time" trace record */
#define COMPSTATSIZE  5
                       /* key of a "message passing statistics" trace record */
#define COMMSTATS  12
                       /* size of a "message passing statistics" trace record */
#define COMMSTATSIZE  8
                       /* key of a "close" trace record */
#define CLOSE  13
                       /* size of a "close" trace record */
#define CLOSESIZE  3
                       /* key of a "trace level" trace record */
#define TLEVEL  14
                       /* size of a "trace level" trace record */
#define TLEVELSIZE  6
                       /* key of a "trace mark" trace record */
#define TMARK  15
                       /* size of a "trace mark" trace record */
#define TMARKSIZE  4
                       /* key of a "trace message" trace record */
#define TMSG  16
                       /* number of fields in a "trace message" trace record */
#define TMSGSIZE  5
                       /* key of a "trace full" trace record */
#define TFULL  17
                       /* size of a "trace full" trace record */
#define TFULLSIZE  3
                       /* key of a "trace flush" trace record */
#define TFLUSH  18
                       /* size of a "trace flush" trace record */
#define TFLUSHSIZE  5
                       /* key of a "trace exit" trace record */
#define TEXIT  19
                       /* size of a "trace exit" trace record */
#define TEXITSIZE  4
                       /* key of a "trace block begin" trace record */
#define TBLOCKBEGIN  20
                       /* key of a "trace block end" trace record */
#define TBLOCKEND  21
                       /* size of a "trace block" trace record */
#define TBLOCKSIZE  6
                       /* key of a "trace mark array" trace record */
#define TMARKS  22
                       /* key of a "clock sync" trace record */
#define CLOCKSYNC  23
                       /* size of a "clock sync" trace record */
#define CLKSYNCSIZE  11
                       /* key of a "trace files" trace record */
#define TFILES  24
                       /* size of a "trace files" trace record */
#define TFILESIZE  5

         /***************** variable definitions ******************/

                      /* switch (0/1) indicating whether the parameter */
                      /* checking has been enabled */
int envCHECKING = 1 ; 
                      /* switch indicating whether the tracing facility */
                      /* is "turned on" on the node */
int envTOPN = 0 ; 
                      /* variable indicating if/when standard trace data */
                      /* being collected */
                      /* (description of possible values in tracelevel) */
int envTRA = 0 ;
                      /* variable indicating if/when busy/idle time data */
                      /* being collected */
                      /* (description of possible values in tracelevel) */
int envTCP = 0 ;
                      /* variable indicating if/when message passing */
                      /* statistics being collected */
                      /* (description of possible values in tracelevel) */
int envTCM = 0 ;
                      /* switch indicating whether the trace array should */
                      /* be flushed when it fills up - if so, the application */
                      /* may take a performance hit */
int envFLU = 0 ;
                      /* switch indicating whether the trace array is */
                      /* full or not */
int envFUL = 0 ;
                      /* switch indicating whether the trace data */
                      /* should be expanded before writing to the */
                      /* trace file, or left compressed */
                      /* (compressed: 0 expanded: 1) */
                      /* default is compressed, so that tracemsg and enverrmsg*/
                      /* send things to the host in the correct format */
                      /* (for code simplification) */
int envVER = 0;
                      /* variable where value of envTRA is saved when */
                      /* envFUL = 1 (envTRA set to 0 to disable further */
                      /* tracing) */
int envSTRA;
                      /* variable where value of envTCP is saved when */
                      /* envFUL = 1 (envTCP set to 0 to disable further */
                      /* tracing) */
int envSTCP;
                      /* variable where value of envTCM is saved when */
                      /* envFUL = 1 (envTCM set to 0 to disable further */
                      /* tracing) */
int envSTCM;
                      /* amount of storage allocated for tracing */
long envTSZ;
                      /* indicator as to whether to generate trace records */
                      /* for traceblkbeg and traceblkend */
int envTRACEBLOCK = 1 ;
                      /* indicator as to whether the sampling rating */
                      /* of compstats and commstats is controlled */
int envSAMPLING = 0 ;
                      /* indicator as to whether the next compstat */
                      /* or commstat record(s) should be generated */
int envSAMPLE = 1 ;
                      /* pointer to permanent trace file */
FILE *envTFP = NULL;
                      /* pointer to temporary trace file for reading */
FILE *envTMPR = NULL;
                      /* pointer to temporary trace file for writing */
FILE *envTMPW = NULL;
                      /* counter of number of buffers of trace information */
                      /* written to temp file and not yet reread. */
int envTMPBUFS = 0;
                      
/************************ trace array definitions and variables *************/

         /**************** constant definitions *********************/

                       /* minimum size of trace array (must be able to hold */
                       /* start/tflush, open, 4 largest "regular" records, */
                       /* close, tfull, and texit) */
#define MINTSZ  60
                       /* trigger for trace flushing (minimum amount of buffer*/
                       /* space that must be left in order for tracing to */
                       /* work correctly - four "large" regular records, */
                       /* close, tfull, and texit) */
#define FLUSHTRIGGER  42
                       /* maximum size of a trace record other than start */ 
                       /* and clocksync */
#define MAXTRECORD  8
                       /* offset of "key" field  */
#define KEY  0
                       /* offset of "clock" fields  */
#define CLOCK1  1
#define CLOCK2  2
                       /* offset of "start clock" fields  */
#define CLOCKSTART1  1
#define CLOCKSTART2  2
                       /* offset of "stop clock" fields  */
#define CLOCKSTOP1  3
#define CLOCKSTOP2  4
                       /* offset of "clock start" fields  */
#define CLKSYNC1  3
#define CLKSYNC2  4
                       /* offset of "type" field  */
#define TYPE  3
                       /* offset of "idle time" fields  */
#define IDLE1  3
#define IDLE2  4
                       /* offset of field indicating number of messages */
                       /* received  */
#define MRECV  3
                       /* offset of field indicating number of bytes of */
                       /* trace data collected  */
#define OSPACE  3
                       /* offset of "trace" field  */
#define TRACE  3
                       /* offset of "destination" field  */
#define TO  4
                       /* offset of "source" field  */
#define FROM  4
                       /* offset of field indicating number of messages */
                       /* sent  */
#define MSENT  4
                       /* offset of "compstats" field  */
#define COMPF  4
                       /* offset of "location type" field  */
#define LOCATION  4
                       /* offset of "clock offset" fields  */
#define OFFSET1  5
#define OFFSET2  6
                       /* offset of "message length" field  */
#define LENGTH  5
                       /* offset of field indicating number of bytes */
                       /* received  */
#define MRVOL  5
                       /* offset of "commstats" field  */
#define COMMF  5
                       /* offset of "parameter type" field  */
#define PARAM  5
                       /* offset of field indicating number of bytes */
                       /* sent  */
#define MSVOL  6
                       /* offset of "clock drift" fields  */
#define DRIFT1  7
#define DRIFT2  8
                       /* offset of field indicating number of message */
                       /* queue probes  */
#define MPROB  7
                       /* offset of "clock reference" fields  */
#define CLKREF1  9
#define CLKREF2  10
                       /* offset of "trace" field  */
#define TRACE2  11
                       /* offset of "compstats" field  */
#define COMPF2  12
                       /* offset of "commstats" field  */
#define COMMF2  13

         /***************** variable definitions ***********************/

                       /* number of bytes of trace data collected */
long envOSPC;
                       /* "communicating" time counters */
long envCTIM1;
long envCTIM2;
                       /* number of messages sent */
long envSNT;
                       /* number of bytes sent */
long envSVOL;
                       /* number of messages received */
long envRCV;
                       /* number of bytes received */
long envRVOL;
                       /* number of probe calls */
long envPRB;
                       /* pointer to trace array */
long *envTDATA;
                       /* index to start of next trace record */
long envNXT;
                       /* current limit on envNXT - it represents the */
                       /* "boundary" of the next TMAXLONG length packet */
long envMLTH;  
                       /* pointer to array of indices into envTDATA */
                       /* indicating where the the packets start when sending */
                       /* back trace data */
long *envTPTR;
                       /* size of envTPTR array */
int envSPTR;
                       /* index to next unused location in envTPTR */
                       /* (defined and initialized in tracenode) */
int envNPTR;
