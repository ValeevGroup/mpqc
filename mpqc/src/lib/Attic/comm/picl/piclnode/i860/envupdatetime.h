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

/* update running total given start and stop times */
/* arguments should be pointers to two long integers */
/* and four long integer values */
/* (essentially calculating T = E - (B - T) twice) */

#define ENVUPDATETIME(t1, t2, b1, b2, e1, e2) \
if ( (unsigned long) *( t2 ) > (unsigned long) b2 ){ \
  *( t2 ) = (4294967295 - (unsigned long) *( t2 )) + (unsigned long) b2 + 1; \
  *( t1 ) = b1 - *( t1 ) - 1; \
  } \
  else { \
  *( t2 ) = (unsigned long) b2 - (unsigned long) *( t2 ); \
  *( t1 ) = b1 - *( t1 ); \
  }; \
if ( (unsigned long) *( t2 ) > (unsigned long) e2 ){ \
  *( t2 ) = (4294967295 - (unsigned long) *( t2 )) + (unsigned long) e2 + 1; \
  *( t1 ) = e1 - *( t1 ) - 1; \
  } \
  else { \
  *( t2 ) = (unsigned long) e2 - (unsigned long) *( t2 ); \
  *( t1 ) = e1 - *( t1 ); \
  }




