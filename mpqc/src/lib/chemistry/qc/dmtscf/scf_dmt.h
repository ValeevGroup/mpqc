/*
 * scf_dmt.h --- files to include when using the dmtscf library
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

#ifndef _scf_dmt_h
#define _scf_dmt_h

#include <chemistry/qc/dmtsym/sym_dmt.h>

#ifdef __cplusplus
class KeyVal;
extern "C" {
#endif /* __cplusplus */

#include <chemistry/qc/dmtscf/scf.h>
#include <chemistry/qc/dmtscf/scfallc.h>
#include <chemistry/qc/dmtscf/scffree.h>
#include <chemistry/qc/dmtscf/scfzero.h>
#include <chemistry/qc/dmtscf/scfasgn.h>
#include <chemistry/qc/dmtscf/scfinit.h>
#include <chemistry/qc/dmtscf/scfiseq.h>
#include <chemistry/qc/dmtscf/scfprnt.h>
#include <chemistry/qc/dmtscf/scfbc0.h>
#include <chemistry/qc/dmtscf/scfrbc0.h>
#include <chemistry/qc/dmtscf/scfsbc0.h>
#include <chemistry/qc/dmtscf/scfrcv0.h>
#include <chemistry/qc/dmtscf/scfsnd0.h>

#include <chemistry/qc/dmtscf/scf_bnd.gbl>
#include <chemistry/qc/dmtscf/scf_core.gbl>
#include <chemistry/qc/dmtscf/scf_diis.gbl>
#include <chemistry/qc/dmtscf/scf_en.gbl>
#include <chemistry/qc/dmtscf/scf_ex.gbl>
#include <chemistry/qc/dmtscf/scf_gmat.gbl>
#include <chemistry/qc/dmtscf/scf_iter.gbl>
#include <chemistry/qc/dmtscf/scf_loopg.gbl>
#include <chemistry/qc/dmtscf/scf_loopj.gbl>
#include <chemistry/qc/dmtscf/scf_loopk.gbl>
#include <chemistry/qc/dmtscf/scf_mkden.gbl>
#include <chemistry/qc/dmtscf/scf_mkgd.gbl>
#include <chemistry/qc/dmtscf/scf_mkgdlb.gbl>
#include <chemistry/qc/dmtscf/scf_mkj.gbl>
#include <chemistry/qc/dmtscf/scf_mkk.gbl>
#include <chemistry/qc/dmtscf/scf_oeis.gbl>
#include <chemistry/qc/dmtscf/scf_orth.gbl>
#include <chemistry/qc/dmtscf/scf_proj.gbl>
#include <chemistry/qc/dmtscf/scf_vect.gbl>

#ifdef __cplusplus
}

class RefKeyVal;

/* The preferred input reader for C++ programs. */
int scf_init_scf_struct(const RefKeyVal&, centers_t&, scf_struct_t&);
int scf_init_scf(const RefKeyVal&, centers_t&, scf_struct_t&, sym_struct_t&);
int scf_make_old_centers(const RefKeyVal&, centers_t&, centers_t&);

void scf_print_options(FILE*, scf_struct_t&);

class StateOut;
class StateIn;

void put_scf_struct(StateOut&, scf_struct_t&);
void get_scf_struct(StateIn&, scf_struct_t&);

#endif /* __cplusplus */

#endif /* _scf_dmt_h */

/*************************************************************************
 * Local Variables:
 * mode: c
 * eval: (c-set-style "ETS")
 * End:
 */
