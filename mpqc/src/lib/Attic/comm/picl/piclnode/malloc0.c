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

/* node version of dynamic memory allocation */

char *malloc0(bytes)
int bytes;
{
  char *malloc();
  char *p;

  p = malloc(bytes);
  if (p == NULL){

    sprintf(envMBUF, "malloc0: insufficient memory, %d bytes required", bytes);
    tracemsg (envMBUF);

    };

  return(p) ;
}

