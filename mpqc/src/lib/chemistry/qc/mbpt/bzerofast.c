/*
 * bzerofast.c
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

#include <chemistry/qc/mbpt/bzerofast.h>

/* Commented out this version since the compiler
 * cannot handle it */
/*int bzerofast(double *d, int dimension)
  {
    int i;
  
    for (i=dimension; i; i--) {
      *d++ = 0.0;
      }
  
    return(0);
  }                     */
int bzerofast(double *d, int dimension)
{
  int i;

  for (i=0; i<dimension; i++) {
    *d++ = 0.0;
    }

  return(0);
}

/* Local Variables:
 * mode: c++
 * eval: (c-set-style "CLJ-CONDENSED")
 * End:
 */
