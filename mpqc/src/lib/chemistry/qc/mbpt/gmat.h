/*
 * gmat.h
 *
 * Copyright (C) 1996 Limit Point Systems, Inc.
 *
 * Author: Ida Nielsen <ibniels@ca.sandia.gov>
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

#ifndef _chemistry_qc_mbpt_gmat_h
#define _chemistry_qc_mbpt_gmat_h
 
int
mbpt_make_gmat(scf_struct_t *scf_info, centers_t *centers,
          RefSymmSCMatrix& Gmat, double *DPmat, FILE *outfile);

int
mbpt_init_gmat(centers_t *centers, scf_struct_t *scf_info, double *intbuf);

void
mbpt_done_gmat(centers_t *centers, scf_struct_t *scf_info);

#endif
