//
// bitarray.cc
//
// Modifications are
// Copyright (C) 1998 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

/* bitarray.h -- definition of the BitArray Class
 *
 *      THIS SOFTWARE FITS THE DESCRIPTION IN THE U.S. COPYRIGHT ACT OF A
 *      "UNITED STATES GOVERNMENT WORK".  IT WAS WRITTEN AS A PART OF THE
 *      AUTHOR'S OFFICIAL DUTIES AS A GOVERNMENT EMPLOYEE.  THIS MEANS IT
 *      CANNOT BE COPYRIGHTED.  THIS SOFTWARE IS FREELY AVAILABLE TO THE
 *      PUBLIC FOR USE WITHOUT A COPYRIGHT NOTICE, AND THERE ARE NO
 *      RESTRICTIONS ON ITS USE, NOW OR SUBSEQUENTLY.
 *
 *  Author:
 *      E. T. Seidl
 *      Bldg. 12A, Rm. 2033
 *      Computer Systems Laboratory
 *      Division of Computer Research and Technology
 *      National Institutes of Health
 *      Bethesda, Maryland 20892
 *      Internet: seidl@alw.nih.gov
 *      July, 1993
 */

#include <util/container/bitarray.h>
#include <util/misc/exenv.h>

using namespace std;

using sc::BitArrayLTri;

BitArrayLTri::BitArrayLTri(int i, int j)
  : a(), n(0), nm(0), na(0)
{
  if (i!=j) {
    ExEnv::err0() << indent << "BitArrayLTri(int,int): i != j"
                 << endl;
    abort();
  }
  int sz = i*(i+1)/2;

  if(sz) {
    n = sz/8 + ((sz%8)?1:0);
    a.resize(n, '\0');
    na=sz; nm=i;
  }
}

BitArrayLTri::~BitArrayLTri()
{
  n=0; na = 0; nm = 0;
}

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
