/*
 * mgd.h --- header for direct G build
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

#ifndef _dmtscf_mgd_h
#define _dmtscf_mgd_h

#define MIN0(a,b) ((a)<(b)) ? (a) : (b)
#define MAX0(a,b) ((a)>(b)) ? (a) : (b)

#ifndef IOFF
#define IOFF(a,b) ((a)>(b))?(a)*((a)+1)/2+(b):(b)*((b)+1)/2+(a)
#endif

struct mgd {
  int si;
  int sj;
  int sk;
  int sl;
  int isz;
  int jsz;
  int ksz;
  int lsz;
  double *glp;
  double *glpo;
  double *plp;
  double *plpo;
  double *gloc;
  double *gloco;
  double *ploc;
  double *ploco;
} ;

#endif

/* Local Variables:
 * mode: c
 * eval: (c-set-style "ETS")
 * End:
 */
