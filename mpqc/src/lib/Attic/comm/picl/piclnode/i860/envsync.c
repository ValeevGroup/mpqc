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

/* sync0 - barrier synchronization using dimensional exchange */
/* written in machine dependent routines for efficiency */

long envSYNC()
{
  int i, half, twin, type;
  long cnt;
  char s;

  cnt = 0;
  if (envNPA > 1){

    for (i=1; i<envNPA; i<<=1);
    half = (i == envNPA) ? i : i/2;
    if (envME < half){

      twin = envME + half;
      if (twin < envNPA) crecv(TYPELIMIT1, &s, 1);
      for (i=1; i<half; i<<=1){

        type = i + TYPELIMIT1;
        csend(type, &s, 1, i^envME, DPID);
        crecv(type, &s, 1);
        cnt++;

        };

      if (twin < envNPA){

        csend(TYPELIMIT1, &s, 1, twin, DPID);
        cnt++;

        }

      }
      else{

      csend(TYPELIMIT1, &s, 1, envME-half, DPID);
      crecv(TYPELIMIT1, &s, 1);
      cnt++;

      };

    };

  return(cnt);

}

