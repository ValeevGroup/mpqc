/*
 * int_libv2.h
 *
 * Copyright (C) 1996 Limit Point Systems, Inc.
 *
 * Author: Curtis Janssen <cljanss@ca.sandia.gov>
 * Maintainer: LPS
 *
 * This file is part of the SC Toolkit.
 *
 * The SC Toolkit is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Library General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * The SC Toolkit is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public License
 * along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
 * the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * The U.S. Government is granted a limited license as per AL 91-7.
 */

/* $Log$
 * Revision 1.13  1996/10/25 23:28:07  etseidl
 * add copyleft notice
 *
 * Revision 1.12  1996/04/15 15:49:01  etseidl
 * add {One,Two}BodyDerivInt, split c++ stuff into more files
 *
 * Revision 1.11  1996/03/29 01:13:38  etseidl
 * moved SphericalTransformIter to basis.  Replaced it with
 * SphericalTransformIterV2.  Also moved integralv2 from integral to here
 *
 * Revision 1.10  1995/11/16 00:47:35  cljanss
 * Removed normalization for individual basis functions.
 *
 * Revision 1.9  1995/08/21 19:36:21  cljanss
 * 1) New integral storage scheme using AVL trees.
 * 2) Updated bounds routines so the SCF program could use them.
 * 3) Got inttest working again.
 *
 * Revision 1.8  1994/10/31  21:30:10  cljanss
 * Switch KeyVal& args to const RefKeyVal&.
 *
 * Revision 1.7  1994/08/26  22:45:36  etseidl
 * fix a bunch of warnings, get rid of rcs id's, get rid of bread/bwrite and
 * fread/fwrite modules
 *
 * Revision 1.6  1994/08/25  16:49:15  etseidl
 * fix error in prototype for int_read_basis()
 *
 * Revision 1.5  1994/08/24  16:07:02  etseidl
 * no longer include atomsip.h, add prototypes for keyval functions
 *
 * Revision 1.4  1994/08/16  20:21:23  etseidl
 * include utils.gbl
 *
 * Revision 1.3  1994/05/27  23:51:22  cljanss
 * Added support for 2 and 3 center 2 electron integrals.  Added a test porgram.
 *
 * Revision 1.2  1993/12/30  13:32:53  etseidl
 * mostly rcs id stuff
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

#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/storage.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <tmpl.h>
#include <chemistry/qc/intv2/atoms.h>
#include <chemistry/qc/intv2/atomsbc0.h>
#include <chemistry/qc/intv2/atomsrbc0.h>
#include <chemistry/qc/intv2/atomssbc0.h>
#include <chemistry/qc/intv2/atomsallc.h>
#include <chemistry/qc/intv2/atomsinit.h>
#include <chemistry/qc/intv2/atomsprnt.h>
#include <chemistry/qc/intv2/atomsfree.h>
#include <chemistry/qc/intv2/atomsasgn.h>
#include <chemistry/qc/intv2/int_macros.h>
#include <chemistry/qc/intv2/int_flags.h>
#include <chemistry/qc/intv2/int_types.h>
#include <chemistry/qc/intv2/init2e.gbl>
#include <chemistry/qc/intv2/comp_erep.gbl>
#include <chemistry/qc/intv2/comp_erep23.gbl>
#include <chemistry/qc/intv2/atominfo.gbl>
#include <chemistry/qc/intv2/int_print.gbl>
#include <chemistry/qc/intv2/offsets.gbl>
#include <chemistry/qc/intv2/comp_0e.gbl>
#include <chemistry/qc/intv2/comp_1e.gbl>
#include <chemistry/qc/intv2/bounds.gbl>
#include <chemistry/qc/intv2/bounds.h>
#include <chemistry/qc/intv2/normalize.gbl>
#include <chemistry/qc/intv2/basis.gbl>
#include <chemistry/qc/intv2/utils.gbl>

#ifdef __cplusplus
}

#include <chemistry/qc/intv2/int_cplus.h>

#endif


#endif /* _chemistry_qc_intv2_int_libv2_h */
