
/* $Log$
 * Revision 1.2  1994/08/26 17:57:42  etseidl
 * get rid of rcs id's and fix a few warnings
 *
 * Revision 1.1.1.1  1993/12/29  12:53:30  etseidl
 * SC source tree 0.1
 *
 * Revision 1.1.1.1  1992/03/17  16:35:08  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:35:07  seidl
 * Initial revision
 *
 * Revision 1.3  1992/01/09  12:06:41  seidl
 * add send0 and recv0 modules
 *
 * Revision 1.2  1992/01/06  11:34:49  seidl
 * add asgn, iseq, and zero modules
 *
 * Revision 1.1  1991/12/20  16:07:47  seidl
 * Initial revision
 *
 * Revision 1.4  91/12/02  17:13:44  cljanss
 * include the generated fwrite and fread headers
 * 
 * Revision 1.3  91/09/28  20:43:22  cljanss
 * switched to new naming convention
 * 
 * Revision 1.2  91/09/28  19:35:56  cljanss
 * added matrix_util
 * 
 * Revision 1.1  1991/06/19  14:44:23  janssen
 * Initial revision
 * */

#ifndef _math_lib_h
#define _math_lib_h

#include <math/array/matrix.h>
#include <math/array/matrix_utl.gbl>
#include <math/array/diag.gbl>
#include <math/array/flin.gbl>
#include <math/array/mxmm.gbl>
#include <math/array/mxmv.gbl>
#include <math/array/print.gbl>

#include <math/array/matrixbc0.h>
#include <math/array/matrixrbc0.h>
#include <math/array/matrixsbc0.h>
#include <math/array/matrixbwr.h>
#include <math/array/matrixbrd.h>
#include <math/array/matrixfree.h>
#include <math/array/matrixfrd.h>
#include <math/array/matrixfwr.h>
#include <math/array/matrixallc.h>
#include <math/array/matrixinit.h>
#include <math/array/matrixprnt.h>
#include <math/array/matrixasgn.h>
#include <math/array/matrixiseq.h>
#include <math/array/matrixzero.h>
#include <math/array/matrixsnd0.h>
#include <math/array/matrixrcv0.h>

#endif /* _math_lib_h */
