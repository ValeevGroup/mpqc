/*
 * opt2.h
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

#ifndef _chemistry_qc_mbpt_opt2_h
#define _chemistry_qc_mbpt_opt2_h

int
mbpt_opt2(centers_t &centers, scf_struct_t &scf_info, sym_struct_t &sym_info,
          dmt_matrix Scf_Vec, dmt_matrix Fock, dmt_matrix FockO,
          int nfzc, int nfzv, int mem,
          int v1, int v2, int lb,
          FILE* outfile);

int
mbpt_opt2_v1(centers_t *centers, scf_struct_t *scf_info, dmt_matrix Scf_Vec,
             double_vector_t *_evals, int nfzc, int nfzv, int mem,
             FILE* outfile);

int
mbpt_opt2_v2(centers_t *centers, scf_struct_t *scf_info, dmt_matrix Scf_Vec,
             double_vector_t *_evals, int nfzc, int nfzv, int mem,
             FILE* outfile);

int
mbpt_opt2v2lb(centers_t *centers, scf_struct_t *scf_info, dmt_matrix Scf_Vec,
              double_vector_t *_evals, int nfzc, int nfzv, int mem,
              FILE* outfile);

#endif
