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

/* routine to test whether a message of the specified type is waiting */
/* to be read */
/* if message not found, 0 is returned, else 1 is returned */

long envPROBE(type)
int type;
{
  long ftype, found, infotype(), iprobe();

  /* check whether type is legal */
  if (type >= TYPELIMIT1){
    if ((type <= TYPELIMIT2) || (type >= TYPELIMIT3)){
      sprintf(envMBUF, 
"probe0: message type %d is prohibited on iPSC/860 - exiting",type);
      envERRMSG(envMBUF);
      };
    };

  /* test for message */
  found = iprobe(type);

  if (found){

    ftype = infotype();

    if (ftype >= TYPELIMIT1){
      if ((ftype <= TYPELIMIT2) || (ftype >= TYPELIMIT3)){

        sprintf(envMBUF, 
"probe0:  probe type %d matched a prohibited type - exiting",type);
        envERRMSG(envMBUF);

        };

      };

    };

  return(found);

}


