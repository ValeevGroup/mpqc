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

/* envglobals.h for iPSC/860 nodelib */

/********************** environmental definitions and variables *************/

         /**************** data structures **************************/

                      /* cell structure for linked list of outstanding */
                      /* isends and irecvs */
struct cell {
    long msgtype;
    long msgid;
    struct cell *next;
    };

         /**************** constant definitions *********************/

                      /* send/receive to/from any processor */
#define  ANY  -1
                      /* load all processors */
#define  ALL  -1
                      /* (default) process identification number */
#define  DPID  0
                      /* maximum length of a message originating on */
                      /* or arriving at the host (in bytes) */
#define  HMAXLENGTH  256000
                      /* packet size for TRACE_ARRAY messages */
                      /* (in long integers: 64000/sizeof(long) )*/
#define  TMAXLONG  16000
                      /* intrinsic clock second normalization factors: */
                      /* (10^-7 of 10^-7 second time and ) */
                      /*  (2^32)/10^7 of (2^32)/10^7 second time) */
#define  NSFACTOR1  429.4967296
#define  NSFACTOR2  .0000001
                      /* clock synchronization algorithm indicators */
                      /* if CALCOFFSET, calculate offset */
                      /* if CALCDRIFT, calculate drift using an initial */
                      /*  sleep(DSLEEP) */
#define  CALCOFFSET  1
#define  CALCDRIFT   1
#define  DSLEEP     60
                      /* message type to tell the nodes what envNPA should be */
                      /* also used to sync the nodes with the host */
#define HOST_INIT   999999995
                      /* message type used by node 0 to initiate clock */
                      /* synchronization logic */
#define CLOCK_PREP  999999994
                      /* message type used by nodes to calculate clock */
                      /* synchronization parameters */
#define CLOCK_SYNC  999999993
                      /* message type lower limit used (if desired) by */
                      /* porthost/node routines */
#define PORT_INIT   800000000
                      /* smallest "prohibited" type (types reserved for use */
                      /* by portable library and OS messages) */
#define TYPELIMIT1   900000000
                      /* largest "prohibited" type (types reserved for use */
                      /* by portable library and OS messages) */
#define TYPELIMIT2   1073741823
                      /* beginning of next "prohibited" range of types */
                      /* (types reserved for use by OS messages) */
#define TYPELIMIT3   2000000000
                      /* maximum number of outstanding sendbegins */
#define ISNDLIMIT 129
                      /* maximum number of outstanding recvbegins */
#define IRCVLIMIT 129

         /***************** variable definitions ******************/

                      /* array to hold two integer timestamp */
extern unsigned long envTIMESTAMP[2];
                      /* global index variable */
extern long envI;
                      /* array of message IDs of pending isends */
extern struct cell *envISARRAY;
                      /* pointer to head of list of pending isends */
extern struct cell *envISHEAD;
                      /* pointer to tail of list of pending isends */
extern struct cell *envISTAIL;
                      /* pointer to head of empty isend cells */
extern struct cell *envISFREE;
                      /* array of message ID of pending irecvs */
extern struct cell *envIRARRAY;
                      /* pointer to head of list of pending irecvs */
extern struct cell *envIRHEAD;
                      /* pointer to tail of list of pending irecvs */
extern struct cell *envIRTAIL;
                      /* pointer to head of empty irecv cells */
extern struct cell *envIRFREE;

         /**************** function definitions *********************/
long infocount(), infonode(), infotype(), numnodes(), iprobe();
long envPROBE(), envRECVEND(), envSYNC();


/************************ trace definitions and variables ******************/

         /**************** constant definitions *********************/

                       /* message type for trace data */
#define TRACE_ARRAY  999999999
                      /* tracemesg message type */
#define TRACE_MESG   999999998
                      /* msg0 message type */
#define MSG0         999999997
                      /* close0 message type */
#define CLOSE0       999999995

