/*
 * sndrcv0.h
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

#ifdef ALLOC_GLOBALS
#define EXTERN
#define INITIALIZE(x,y) x = y
#else
#define EXTERN extern
#define INITIALIZE(x,y) (x)
#endif

EXTERN int sgen_sndrcv0_print_send;
EXTERN int sgen_sndrcv0_print_receive;
EXTERN int sgen_sndrcv0_print_host;

/* The message types >= MIN_TYPE and <= MAX_TYPE, as well as
 * SYNC_START and SYNC_END are reserved for use by the sgen generated routines.
 */
#define MIN_TYPE    50000
#define MAX_TYPE    65500
#define SYNC_START  65501
#define SYNC_END    65502

#define USE_INCREMENTED_TYPE

#ifdef USE_INCREMENTED_TYPE
EXTERN int INITIALIZE(sgen_sndrcv0_type,MIN_TYPE);
#define TYPENOINC() sgen_sndrcv0_type

#define TYPE() ((sgen_sndrcv0_type > MAX_TYPE)? \
                  (sgen_sndrcv0_type = MIN_TYPE) \
                 :(sgen_sndrcv0_type++))

#define PRINT(sorr,type,root,size) \
       if (  (sgen_sndrcv0_print_send && (sorr == 's')) \
           ||(sgen_sndrcv0_print_receive && (sorr == 'r'))) \
           { sgen_sndrcv_print(sorr,sgen_sndrcv0_type,root,size); }
#else
#define TYPENOINC() type
#define TYPE() type

#define PRINT(sorr,type,root,size) \
       if (  (sgen_sndrcv0_print_send && (sorr == 's')) \
           ||(sgen_sndrcv0_print_receive && (sorr == 'r'))) \
           { sgen_sndrcv_print(sorr,type,root,size); }
#endif

#define PRINT_DATA(sorr,conv,data) \
       if (  ((sgen_sndrcv0_print_send & 2) && (sorr == 's')) \
           ||((sgen_sndrcv0_print_receive & 2) && (sorr == 'r'))) \
           { if (    (mynode0()==sgen_sndrcv0_print_host) \
                  || (sgen_sndrcv0_print_host == -1)) \
               printf(conv,data); }
