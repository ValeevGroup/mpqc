/*
 * cphf.h
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

#ifndef _mpqcic_cphf_h
#define _mpqcic_cphf_h
 
void
mbpt_cphf(centers_t *centers, scf_struct_t *scf_info, FILE *outfile,
          int nbasis, int nvir, int nocc, double **Scf_Vec, double *L,
          double *eigval, RefSCMatrix& P2);

#endif
