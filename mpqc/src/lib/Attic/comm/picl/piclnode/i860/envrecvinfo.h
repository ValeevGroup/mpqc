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

/* node version of received message info */
/* arguments should be pointers to integers */

#define ENVRECVINFO(bytes, type, node) \
 *( bytes ) = infocount(); \
  *( type ) = infotype(); \
  *( node ) = infonode(); \
  if (*( node ) == envHST) *( node ) = -32768

