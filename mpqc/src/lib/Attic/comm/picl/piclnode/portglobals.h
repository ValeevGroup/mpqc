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

/* portglobal definitions for high level nodelib routines */

/********************** environmental definitions and variables *************/

         /***************** variable definitions ******************/

                      /* number of processors allocated in open0 */
                      /* (value saved here since setarc0 resets) */
                      /* (envNPA) */
extern int envNPH;
                      /* interconnection topology associated with specified */
                      /* machine model */
                      /* (values and meanings described in setarc0) */
                      /* (default value is fully connected) */
extern int envTOP;
                      /* variable specifying whether natural (0) or */
                      /* Gray ordering (1) of processors should be used */
                      /* (default value is natural ordering) */
extern int envORD;
                      /* direction associated with specified */
                      /* machine model */
                      /* (values and meanings described in setarc0) */
                      /* (default value is "forward") */
extern int envDIR;
                      /* the "half-way" node in the specified topology */
extern int envHLF;
                      /* the "mirror-node" of me (reflected around the */
                      /* half-way node) */
extern int envTWN;

/************************ trace definitions and variables ******************/

         /**************** constant definitions *********************/

                       /* block type of a "barrier" trace record */
#define BARRIER  -1
                       /* block type of a "bcast0" trace record */
#define BCAST0  -2
                       /* block type of a "bcast1" trace record */
#define BCAST1  -3
                       /* block type of a "gcomb" trace record */
#define GCOMB  -4
                       /* block type of a "gather" trace record */
#define GATHER  -5


