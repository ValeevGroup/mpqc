/*
 * int_cplus.h
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

#ifndef _chemistry_qc_intv2_int_cplus_h
#define _chemistry_qc_intv2_int_cplus_h

#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/intv2/atoms.h>

class RefKeyVal;
int int_read_basis(const RefKeyVal&, char*, const char*, basis_t&);
int int_read_centers(const RefKeyVal&, centers_t&);

centers_t * int_centers_from_gbs(const RefGaussianBasisSet&);

#endif
