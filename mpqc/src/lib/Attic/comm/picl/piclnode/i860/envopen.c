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
#include "envglobaldefs.h"

/* open communication channel at node */

envOPEN()
{
  long mynode(), myhost(), numnodes();
  char *getenv();
  struct cell *i;

  if (getenv("HOSTED") != NULL){
    /* get envNPA value from the host */
    crecv(HOST_INIT, &envNPA, sizeof(int));
    envHOSTED = 1;
    }
    else{
    envNPA = numnodes();
    envHOSTED = 0;
    };

  /* set up linked list for isends and irecvs */
  envISFREE = envISARRAY;
  for (i=envISARRAY; i<envISARRAY+ISNDLIMIT-1; i++) i->next = i+1; 
  envIRFREE = envIRARRAY;
  for (i=envIRARRAY; i<envIRARRAY+IRCVLIMIT-1; i++) i->next = i+1;

  /* initialize node and host IDs */
  envME = mynode();
  envHST = myhost();

}
