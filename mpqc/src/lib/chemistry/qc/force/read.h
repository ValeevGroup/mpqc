/*
 * read.h
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

#ifndef _chemistry_qc_force_read_h
#define _chemistry_qc_force_read_h

#ifdef __cplusplus
void
dmt_force_osscf_keyval_init(KeyVal*,FILE*);
void
dmt_force_csscf_keyval_init(KeyVal*,FILE*);

void dmt_get_csscf_force(StateIn&);
void dmt_put_csscf_force(StateOut&);
void dmt_get_osscf_force(StateIn&);
void dmt_put_osscf_force(StateOut&);

extern "C" {
int
dmt_force_read_and_bcast_boolean(KeyVal*keyval,FILE*,char*name,int*boolval);
int
dmt_force_read_and_bcast_int(KeyVal*keyval,FILE*,char*name,int*intval);
int
dmt_force_read_string(KeyVal*keyval,FILE*,char*name,char**val);
}
#else
int
dmt_force_read_and_bcast_boolean(void*keyval,FILE*,char*name,int*boolval);
int
dmt_force_read_and_bcast_int(void*keyval,FILE*,char*name,int*intval);
int
dmt_force_read_string(void*keyval,FILE*,char*name,char**val);
#endif

#endif
