/*
 * buildinter.h
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

/* $Log$
 * Revision 1.3  1996/10/25 23:27:56  etseidl
 * add copyleft notice
 *
 * Revision 1.2  1993/12/30 13:32:45  etseidl
 * mostly rcs id stuff
 *
 * Revision 1.1.1.1  1992/03/17  16:32:32  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:32:31  seidl
 * Initial revision
 *
 * Revision 1.1  1992/02/03  15:49:49  seidl
 * Initial revision
 *
 * Revision 1.1  1991/12/14  00:19:36  cljanss
 * Initial revision
 * */
/*
 * Revision 1.2  91/09/28  18:17:13  cljanss
 * Intermediates are no longer stored, if requested with flags.
 * 
 * Revision 1.1  91/06/16  16:40:07  janssen
 * Initial revision
 *  */

#ifdef ALLOC_BUILDINTER
# define EXTERN
#else
# define EXTERN extern
#endif

EXTERN double int_v_ooze;
EXTERN double int_v_zeta12;
EXTERN double int_v_zeta34;
EXTERN double int_v_oo2zeta12;
EXTERN double int_v_oo2zeta34;
EXTERN double int_v_W0;
EXTERN double int_v_W1;
EXTERN double int_v_W2;
EXTERN double int_v_p120;
EXTERN double int_v_p121;
EXTERN double int_v_p122;
EXTERN double int_v_p340;
EXTERN double int_v_p341;
EXTERN double int_v_p342;
EXTERN double int_v_r10;
EXTERN double int_v_r11;
EXTERN double int_v_r12;
EXTERN double int_v_r20;
EXTERN double int_v_r21;
EXTERN double int_v_r22;
EXTERN double int_v_r30;
EXTERN double int_v_r31;
EXTERN double int_v_r32;
EXTERN double int_v_r40;
EXTERN double int_v_r41;
EXTERN double int_v_r42;
EXTERN double int_v_k12;
EXTERN double int_v_k34;
EXTERN doublep_array3_t int_v_list;

#undef EXTERN

