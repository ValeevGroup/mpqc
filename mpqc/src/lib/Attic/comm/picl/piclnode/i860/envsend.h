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

/* macro to send a message */
/* arguments should be a pointer and 3 integers, respectively */

#define ENVSEND(buf, bytes, type, node, checking) \
  if ( checking == 1 ) envSEND( buf , bytes , type , node ); \
    else{ \
    if ( node != -1){ \
 \
    if ( node == -32768) _csend( type , buf , bytes , envHST , 0); \
      else _csend( type , buf , bytes , node , 0); \
 \
    } \
    else{ \
 \
    if (numnodes() != envNPA){ \
 \
      for (envI=0; envI<envNPA; envI++) \
        if (envI != envME) _csend( type , buf , bytes , envI , 0); \
 \
      } \
      else _csend( type , buf , bytes , node , 0); \
 \
     } \
   }






