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

/* get tracing information */
/* - remaining: approximate number of trace messages that can be saved in */
/*              the remaining free storage in the trace array */
/* - trace: regular tracing type */
/* - compstats: busy/idle tracing type */
/* - commstats: communication tracing type */
/* if type <= 0, record only open, close, tracestart, tracelevel */
/*                    traceexit calls and summary statisitics) */
/* if type >= 1, also record "mark" calls*/
/* if type >= 2, also record higher level routines, and */
/*               send/recv/sync calls except */
/*               the send/recv calls within the higher level routines */
/* if type >= 3, also record send/recv/sync calls within higher level routines*/

traceinfo (remaining, trace, compstats, commstats)
int *remaining, *trace, *compstats, *commstats;
{

  if (envFUL == 1){

    *trace = envSTRA;
    *compstats = envSTCP;
    *commstats = envSTCM;
    if (envTOPN <= 0) *remaining = -1;
      else *remaining = 0;

    }
    else{

    *trace = envTRA;
    *compstats = envTCP;
    *commstats = envTCM;
    if (envTOPN <= 0) *remaining = -1;
      else *remaining = (envTSZ - envNXT - FLUSHTRIGGER)/MAXTRECORD;

    };

}
