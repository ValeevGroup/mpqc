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
#include "envrecv.h"
#include "envsend.h"

/* clock synchronization routine */
/* (calculates offset, drift, and reference time for normalizing) */
/* (clock with respect to node 0) */
/* (modification of program due to Tom Dunigan) */
envCLOCKSYNC()
{
  int i, j, first;
  double p, s, t, tmin, clk_delay, zeroclk_ref, localclk_ref;
  double clock0();
  
  if (envNPA > 1){

    /* test whether this is the first time clocksync has been called */
    first = ( envCLKSYNC == 0 ? 1 : 0 );

    /* disable clock0 normalization */
    envCLKSYNC = 0;

    switch(envME){
     case 0:
      if (first){

        if (CALCOFFSET){

          /* calculate offset */

          /* send offset time to each node */
          /* (spanning tree would be faster, but error grows) */
          for (i=1; i<envNPA;  i++){ 

            /* tell him to ask */
            ENVSEND(&p, sizeof(double), CLOCK_PREP, i, 0) ;  

            /* do it a few times to get best */
            for (j=0; j<REPS; j++){  

              /* await request */
              ENVRECV(&p, sizeof(double), CLOCK_SYNC, 0);  
              p = clock0();
              /*send him time */
              ENVSEND(&p, sizeof(double), CLOCK_SYNC, i, 0) ;

              };

            };
  
          };

        if (CALCDRIFT){

          /* use DSLEEP to calculate drift */

          /* "sleep" for dsleep seconds */
          p = clock0() + DSLEEP;
          while (clock0() < p) ;

          /* start calculation of drift */
          for (i=1; i<envNPA;  i++){ 

            /* send time to each node */
            /* (tell him to ask) */
            ENVSEND(&p, sizeof(double), CLOCK_PREP, i, 0) ;  

            /* do it a few times to get best */
            for (j=0; j<REPS; j++){  

              /* await request */
              ENVRECV(&p, sizeof(double), CLOCK_SYNC, 0);  
              p = clock0();
              /*send him time */
              ENVSEND(&p, sizeof(double), CLOCK_SYNC, i, 0) ; 

              };
  
            };

          };

        }
        else{

        if (CALCDRIFT){

          /* start calculation of drift */
          for (i=1; i<envNPA;  i++){ 

            /* send time to each node */
            /* (tell him to ask) */
            ENVSEND(&p, sizeof(double), CLOCK_PREP, i, 0) ;  

            /* do it a few times to get best */
            for (j=0; j<REPS; j++){  

              /* await request */
              ENVRECV(&p, sizeof(double), CLOCK_SYNC, 0);  
              p = clock0();
              /*send him time */
              ENVSEND(&p, sizeof(double), CLOCK_SYNC, i, 0) ; 

              };
  
            };

          };

        };

      break;
  
     default:

      if (first){

        if (CALCOFFSET){

          /* await start from 0 */
          ENVRECV(&p, sizeof(double), CLOCK_PREP, 0);  
  
          /* try a few times for min */
          tmin = 1000000.0;
          for (j=0; j<REPS; j++){
  
            s = clock0();
            /* ask 0 for time */
            ENVSEND(&p, sizeof(double), CLOCK_SYNC, 0, 0) ; 
            /* await time from 0 */
            ENVRECV(&p, sizeof(double), CLOCK_SYNC, 0);  
            t = clock0();
            clk_delay = (t - s)/2.0;
            if (clk_delay < tmin){ 
  
              /*save best */
              tmin = clk_delay;
              zeroclk_ref = p + clk_delay;
              localclk_ref = t;
  
              };
  
            };
  
          envCLKREF = localclk_ref;
          /* difference in approximation to node 0's time and local time*/
          envCLKOFFSET = zeroclk_ref - localclk_ref;

          };

        };

      if (CALCDRIFT){
  
        /* await start from 0 */
        ENVRECV(&p, sizeof(double), CLOCK_PREP, 0);  

        /* try a few times for min */
        tmin = 1000000.0;
        for (j=0; j<REPS; j++){
  
          s = clock0();
          /* ask 0 for time */
          ENVSEND(&p, sizeof(double), CLOCK_SYNC, 0, 0) ; 
          /* await time from 0 */
          ENVRECV(&p, sizeof(double), CLOCK_SYNC, 0);  
          t = clock0();
          clk_delay = (t - s)/2.0;
          if (clk_delay < tmin){ 
  
            /*save best */
            tmin = clk_delay;
            zeroclk_ref = p + clk_delay;
            localclk_ref = t;

            };
  
          };
  
        /* difference in new and old offsets divided by */
        /* difference in when new and old were calculated */
        envCLKDRIFT = ((zeroclk_ref - localclk_ref)  - envCLKOFFSET)/
                       (localclk_ref - envCLKREF);
  
        };
  

      }; /*end switch */

    };

    /* enable clock0 normalization */
    envCLKSYNC = 1;

}



