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

/* envCLOCK returns time in 200s of nanoseconds */
/* in two long integers */
/* arguments should be pointers to two long integers */

#define ENVCLOCK(t1, t2) \
 _hwclock(envTIMESTAMP); \
 *((unsigned long *) t2 ) = envTIMESTAMP[0]; \
 *((unsigned long *) t1 ) = envTIMESTAMP[1]




