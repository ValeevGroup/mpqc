/*
 * matrix.h
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

#ifndef _math_dmt_matrix_h
#define _math_dmt_matrix_h

typedef int dmt_matrix;
typedef int mdescr;
typedef int booln;

#define MAXARGS 50

struct ngl_descr {
  int type;
  mdescr m;
  double *a;
  double *t;
  int readonly;
  };
typedef struct ngl_descr ngl_descr_t;

struct ngl_list {
  ngl_descr_t args[MAXARGS];
  int nargs;
  int rmsgtype;
  int smsgtype;
  void *buf;
  int current_bufsize;
  int total_bufsize;
  int loop_count;
  int readonly;
  int nextblock;
  int innerloopmat;
  struct ngl_list *p;
  };
typedef struct ngl_list loop_t;

#define COLUMNS   1
#define SCATTERED 2

struct dmt_block_info {
  int i;
  int j;
  int magnitude; /* the value of Qvec for this block */
  int ami;
  int amj;
  int nconi;
  int nconj;
  int nprimi;
  int nprimj;
  int dimi;
  int dimj;
};
typedef struct dmt_block_info dmt_cost_t;


#endif
