/*
 picl.h

 Copyright (C) 1996 Limit Point Systems, Inc.

 Author: Curtis Janssen <cljanss@limitpt.com>
 Maintainer: LPS

 This file is part of the SC Toolkit.

 The SC Toolkit is free software; you can redistribute it and/or modify
 it under the terms of the GNU Library General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 The SC Toolkit is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Library General Public License for more details.

 You should have received a copy of the GNU Library General Public License
 along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
 the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.

 The U.S. Government is granted a limited license as per AL 91-7.
*/

/* All of the official picl routines are declared here.  These
 * declarations apply to all of the picl packages. */

#ifndef _PICL_H
#define _PICL_H

#ifdef __cplusplus
#include <util/group/message.h>
void open0_messagegrp(int*,int*,int*, const RefMessageGrp&);

extern "C" {
#endif

#ifdef __STDC__
#define ARGS0 void
#define ARGS1(a) a
#define ARGS2(a,b) a,b
#define ARGS3(a,b,c) a,b,c
#define ARGS4(a,b,c,d) a,b,c,d
#define ARGS5(a,b,c,d,e) a,b,c,d,e
#define ARGS6(a,b,c,d,e,f) a,b,c,d,e,f
#else
#define ARGS0
#define ARGS1(a)
#define ARGS2(a,b)
#define ARGS3(a,b,c)
#define ARGS4(a,b,c,d)
#define ARGS5(a,b,c,d,e)
#define ARGS6(a,b,c,d,e,f)
#endif

/* Low level routines. */
int host0( ARGS0 );
void check0( ARGS1(int) );
double clock0( ARGS0 );
void close0( ARGS1(int) );
void load0( ARGS2(char*,int) );
void message0( ARGS1(char*) );
void open0( ARGS3(int*,int*,int*) );
int probe0( ARGS1(int) );
void recv0( ARGS3(void*,int,int) );
void recvinfo0( ARGS3(int*,int*,int*) );
void send0( ARGS4(void*,int,int,int) );
void sync0( ARGS0 );
void who0( ARGS3(int*,int*,int*) );

/* High level routines. */
void setarc0( ARGS4(int*,int*,int*,int*) );
void getarc0( ARGS4(int*,int*,int*,int*) );

void bcast0( ARGS4(void*,int,int,int ) );

void gcomb0( ARGS6(void*,int,int,int,int,void(*comb)(void*,void*,int,int)) );
void gand0( ARGS5(void*,int,int,int,int) );
void gmax0( ARGS5(void*,int,int,int,int) );
void gmin0( ARGS5(void*,int,int,int,int) );
void gor0( ARGS5(void*,int,int,int,int) );
void gprod0( ARGS5(void*,int,int,int,int) );
void gsum0( ARGS5(void*,int,int,int,int) );
void gxor0( ARGS5(void*,int,int,int,int) );

void barrier0( ARGS0 );

int gray0( ARGS1(int) );
int ginv0( ARGS1(int) );

/****************************************************************************/
/* These routines are not a part picl, but their use crept into the code.   */

#if !defined(HAVE_NX)
int cubedim( ARGS0 );
int mynode ( ARGS0 );
int numnodes( ARGS0 );
int infocount( ARGS0 );
#endif
int mynode0( ARGS0 );
int cubedim0( ARGS0 );
int numnodes0( ARGS0 );

int loop_out_neigh( ARGS0 );
int my_loop_index( ARGS0 );
int loop_in_neigh( ARGS0 );
int loop_neigh( ARGS1(int) );
void gop0( ARGS4(double*,int,int,int) );
void gop1( ARGS5(double*,int,double*,int,int) );
void gop0_sc( ARGS4(signed char*,int,int,int) );

void picl_prober( ARGS0 );

#ifdef I860
void gdcomb(ARGS5(int,double*,double*,int,int));
#endif

#undef ARGS0
#undef ARGS1
#undef ARGS2
#undef ARGS3
#undef ARGS4
#undef ARGS5
#undef ARGS6

#ifdef __cplusplus
}
#endif

#endif /* _PICL_H */

/*
 Local Variables:
 mode: c++
 c-file-style: "CLJ"
 End:
*/
