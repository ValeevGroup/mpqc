/*
 * symm_mac.h
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

#ifndef _sym_symm_mac_h
#define _sym_symm_mac_h

enum pgroups {_PG_C1, _PG_CS, _PG_CI, _PG_CN, _PG_CNV, _PG_CNH, _PG_DN,
              _PG_DND, _PG_DNH, _PG_SN, _PG_T, _PG_TH, _PG_TD, _PG_O, _PG_OH,
              _PG_I, _PG_IH};

enum cart_p {y, z, x};
enum cart_d {yy, yz, zz, xy, xz, xx};
enum pure_d {pzz, xmy, pxy, pxz, pyz};
enum cart_f {yyy, yyz, yzz, zzz, xyy, xyz, xzz, xxy, xxz, xxx};
enum cart_g {yyyy, yyyz, yyzz, yzzz, zzzz, xyyy, xyyz, xyzz, xzzz,
             xxyy, xxyz, xxzz, xxxy, xxxz, xxxx};

enum angmom {_AM_S,_AM_P,_AM_D,_AM_F,_AM_G,_AM_H,_AM_I,_AM_J,_AM_K,_AM_L,
            _AM_M,_AM_N,_AM_O,_AM_Q,_AM_R,_AM_T,_AM_U,_AM_V,_AM_X,_AM_Y,_AM_Z};

#endif
