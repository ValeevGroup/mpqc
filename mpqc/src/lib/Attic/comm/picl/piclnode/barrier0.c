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

#include "envglobals.h"
#include "portglobals.h"

barrier0 ( )

/* 
 *  Barrier synchronization using dimensional exchange.
 *  This routine differs from sync0 in that barrier0 is portable.
 */
{
  int p, me, host, top, ord, dir, type, i;
  char s;

  traceblockbegin(BARRIER, 0, 0);

  /**** calculate the barrier ****/
  who0(&p, &me, &host);
  getarc0(&p, &top, &ord, &dir);

  if (me < envHLF) {

    if (envTWN < p) recv0(&s, 1, PORT_INIT);
    for (i=1; i<envHLF; i<<=2){

      type = i+PORT_INIT ;
      send0(&s, 1, type, i^me) ;
      recv0(&s, 1, type) ;

      };
    if (envTWN < p) send0(&s, 1, PORT_INIT, envTWN);

    }
    else{
    send0(&s, 1, PORT_INIT, envTWN) ;
    recv0(&s, 1, PORT_INIT) ;
    };
  /*******************************/

  traceblockend(BARRIER, 0, 0);

}


