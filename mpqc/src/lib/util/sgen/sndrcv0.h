
#ifdef ALLOC_GLOBALS
#define EXTERN
#define INITIALIZE(x,y) x = y
#else
#define EXTERN extern
#define INITIALIZE(x,y) (x)
#endif

EXTERN int sgen_sndrcv0_print_send;
EXTERN int sgen_sndrcv0_print_receive;
EXTERN int sgen_sndrcv0_print_host;

/* The message types >= MIN_TYPE and <= MAX_TYPE, as well as
 * SYNC_START and SYNC_END are reserved for use by the sgen generated routines.
 */
#define MIN_TYPE    50000
#define MAX_TYPE    65500
#define SYNC_START  65501
#define SYNC_END    65502

#define USE_INCREMENTED_TYPE

#ifdef USE_INCREMENTED_TYPE
EXTERN int INITIALIZE(sgen_sndrcv0_type,MIN_TYPE);
#define TYPENOINC() sgen_sndrcv0_type

#define TYPE() ((sgen_sndrcv0_type > MAX_TYPE)? \
                  (sgen_sndrcv0_type = MIN_TYPE) \
                 :(sgen_sndrcv0_type++))

#define PRINT(sorr,type,root,size) \
       if (  (sgen_sndrcv0_print_send && (sorr == 's')) \
           ||(sgen_sndrcv0_print_receive && (sorr == 'r'))) \
           { sgen_sndrcv_print(sorr,sgen_sndrcv0_type,root,size); }
#else
#define TYPENOINC() type
#define TYPE() type

#define PRINT(sorr,type,root,size) \
       if (  (sgen_sndrcv0_print_send && (sorr == 's')) \
           ||(sgen_sndrcv0_print_receive && (sorr == 'r'))) \
           { sgen_sndrcv_print(sorr,type,root,size); }
#endif

#define PRINT_DATA(sorr,conv,data) \
       if (  ((sgen_sndrcv0_print_send & 2) && (sorr == 's')) \
           ||((sgen_sndrcv0_print_receive & 2) && (sorr == 'r'))) \
           { if (    (mynode0()==sgen_sndrcv0_print_host) \
                  || (sgen_sndrcv0_print_host == -1)) \
               printf(conv,data); }
