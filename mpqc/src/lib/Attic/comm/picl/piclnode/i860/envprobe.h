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

/* macro to test whether a message of the specified type is waiting */
/* to be read */
/* if message not found, 0 is returned, else 1 is returned */
/* argument should be an integer value */

#define ENVPROBE(type, checking) \
 ( checking == 1 ? envPROBE( type ) : iprobe( type ) )


