
/* $Log$
 * Revision 1.1  1993/12/29 12:53:01  etseidl
 * Initial revision
 *
 * Revision 1.4  1993/04/28  00:32:04  jannsen
 * added some include files and c++ support
 *
 * Revision 1.3  1992/06/17  22:04:48  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.2  1992/05/13  18:29:40  jannsen
 * added bounds checking for derivative integrals
 *
 * Revision 1.1.1.1  1992/03/17  16:33:03  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:33:02  seidl
 * Initial revision
 *
 * Revision 1.1  1992/02/03  15:49:49  seidl
 * Initial revision
 *
 * Revision 1.2  1992/01/10  18:01:24  cljanss
 * added storage.gbl
 *
 * Revision 1.1  1991/12/14  00:19:36  cljanss
 * Initial revision
 * */
/*
 * Revision 1.3  91/09/28  19:26:51  cljanss
 * Switch to new file naming scheme
 * 
 * Revision 1.2  91/08/09  16:56:40  cljanss
 * added some more of the sgen generated included files
 * 
 * Revision 1.1  1991/06/16  16:40:07  janssen
 * Initial revision
 * */

#ifndef _chemistry_qc_intv2_int_libv2_h
#define _chemistry_qc_intv2_int_libv2_h

#ifdef __cplusplus
extern "C" {
#endif

#include <tmpl.h>
#include <chemistry/qc/intv2/atoms.h>
#include <chemistry/qc/intv2/atomsbc0.h>
#include <chemistry/qc/intv2/atomsrbc0.h>
#include <chemistry/qc/intv2/atomssbc0.h>
#include <chemistry/qc/intv2/atomsbrd.h>
#include <chemistry/qc/intv2/atomsbwr.h>
#include <chemistry/qc/intv2/atomsallc.h>
#include <chemistry/qc/intv2/atomsinit.h>
#include <chemistry/qc/intv2/atomsip.h>
#include <chemistry/qc/intv2/atomsprnt.h>
#include <chemistry/qc/intv2/atomsfree.h>
#include <chemistry/qc/intv2/atomsasgn.h>
#include <chemistry/qc/intv2/int_macros.h>
#include <chemistry/qc/intv2/int_flags.h>
#include <chemistry/qc/intv2/int_types.h>
#include <chemistry/qc/intv2/init2e.gbl>
#include <chemistry/qc/intv2/comp_erep.gbl>
#include <chemistry/qc/intv2/atominfo.gbl>
#include <chemistry/qc/intv2/int_print.gbl>
#include <chemistry/qc/intv2/offsets.gbl>
#include <chemistry/qc/intv2/comp_0e.gbl>
#include <chemistry/qc/intv2/comp_1e.gbl>
#include <chemistry/qc/intv2/storage.gbl>
#include <chemistry/qc/intv2/bounds.gbl>
#include <chemistry/qc/intv2/bounds.h>
#include <chemistry/qc/intv2/normalize.gbl>
#include <chemistry/qc/intv2/basis.gbl>

#ifdef __cplusplus
}
#endif


#endif /* _chemistry_qc_intv2_int_libv2_h */
