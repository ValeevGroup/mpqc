//
// bitarray.h
//
// Modifications are
// Copyright (C) 1996 Limit Point Systems, Inc.
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

#ifndef _util_container_bitarray_h
#define _util_container_bitarray_h

#include <string.h>
#include <stdlib.h>

#include <util/misc/formio.h>

namespace sc {

//
// class BitArrayLTri is used as the lower triangle of a boolean matrix.
// rather than storing an int or a char, just use one bit for each, so
// instead of n(n+1)/2 bytes of storage you have n(n+1)/16 bytes.  A
// further savings of n bits could be obtained by setting the diagonal to
// always true or always false depending on the application, but this would
// probably be more expensive computationally than it's worth.
//

class BitArrayLTri {
  private:
    unsigned char *a;
    int n;
    int nm;
    int na;

    static int
    ij_offset(int i, int j)
        {
          return (i>j) ? (((i*(i+1)) >> 1) + j) : (((j*(j+1)) >> 1) + i);
        }

  public:
    BitArrayLTri(int =0, int =0);
    ~BitArrayLTri();

    void set(unsigned int i) { a[(i>>3)] |= (1 << (i&7)); }
    void set(unsigned int i, unsigned int j) { set(ij_offset(i,j)); }

    int is_set(unsigned int i, unsigned int j) const
      { int ij = ij_offset(i,j); return (a[(ij>>3)] & (1 << (ij&7))); }
    int is_set(unsigned int i) const
      { return (a[(i>>3)] & (1 << (i&7))); }

    int operator()(unsigned int i, unsigned int j) const
      { int ij = ij_offset(i,j); return (a[(ij>>3)] & (1 << (ij&7))); }
    int operator()(unsigned int i) const
      { return (a[(i>>3)] & (1 << (i&7))); }
    int operator[](unsigned int i) const
      { return (a[(i>>3)] & (1 << (i&7))); }

    int dim() const { return na; }
    int nrow() const { return nm; }
    int ncol() const { return nm; }

    int degree(unsigned int i) const {
      int nedge=0;
      for (int j=0; j < nm; j++) if ((*this)(i,j)) nedge++;
      return nedge;
    }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
