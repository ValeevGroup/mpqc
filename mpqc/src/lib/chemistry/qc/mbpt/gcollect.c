/*
 * gcollect.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <util/group/picl.h>

#ifdef HAVE_CUBE_H
#include <cube.h>
#endif
#ifdef HAVE_NX_H
#include <cube.h>
#endif

void
gcollect(double *x, int *lens, double *y)
{
#if defined(HAVE_CUBE_H) || defined(HAVE_NX_H)
  gcolx(x,lens,y);
#else

  int i,leny;
  int me = mynode0();
  int nproc = numnodes0();
  int dim = cubedim0();
  int myoff=0;
  double *temp;

  if (dim==0) {
    memcpy(y,x,lens[0]);
    return;
    }
    
  for (i=0; i < me; i++) myoff += lens[i];
  myoff /= sizeof(double);

  leny = 0;
  for (i=0; i < nproc; i++) leny += lens[i];

  memset(y,0,leny);

  memcpy(&y[myoff],x,lens[me]);

  temp = (double *) malloc(leny);
  if (!temp) {
    fprintf(stderr,"gcollect: could not malloc temp %d\n",leny);
    exit(1);
    }

  gop1(y,leny/sizeof(double),temp,'+',3);

  free(temp);
#endif
  }

/* Local Variables:
 * mode: c++
 * eval: (c-set-style "CLJ-CONDENSED")
 * End:
 */
