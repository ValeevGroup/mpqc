/*
 * sym_dmt.h
 *
 * Copyright (C) 1996 Limit Point Systems, Inc.
 *
 * Author: Edward Seidl <seidl@janed.com>
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

#ifndef _sym_dmt_h
#define _sym_dmt_h

#ifdef __cplusplus

extern "C" {
#endif

#include <math/dmt/libdmt.h>

#include <chemistry/qc/intv2/atoms.h>

#include <chemistry/qc/dmtsym/symm.h>
#include <chemistry/qc/dmtsym/symm_mac.h>
#include <chemistry/qc/dmtsym/symminit.h>
#include <chemistry/qc/dmtsym/symmallc.h>
#include <chemistry/qc/dmtsym/symmasgn.h>
#include <chemistry/qc/dmtsym/symmiseq.h>
#include <chemistry/qc/dmtsym/symmzero.h>
#include <chemistry/qc/dmtsym/symmfree.h>
#include <chemistry/qc/dmtsym/symmprnt.h>
#include <chemistry/qc/dmtsym/symmbc0.h>
#include <chemistry/qc/dmtsym/symmrbc0.h>
#include <chemistry/qc/dmtsym/symmsbc0.h>
#include <chemistry/qc/dmtsym/symmrcv0.h>
#include <chemistry/qc/dmtsym/symmsnd0.h>

#include <chemistry/qc/dmtsym/create_r.gbl>
#include <chemistry/qc/dmtsym/mksym.gbl>
#include <chemistry/qc/dmtsym/skeleton.gbl>
#include <chemistry/qc/dmtsym/syminit.gbl>

#ifdef __cplusplus
}

#include <util/keyval/keyval.h>
#include <math/symmetry/pointgrp.h>
#include <chemistry/qc/basis/gaussbas.h>

int sym_struct_from_gbs(const RefGaussianBasisSet&, sym_struct_t&);
int sym_struct_from_pg(const PointGroup&, centers_t&, sym_struct_t&);
int sym_init_centers(const RefKeyVal&, centers_t&, sym_struct_t&);

#endif

#endif
