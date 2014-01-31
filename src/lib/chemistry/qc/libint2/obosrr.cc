//
// obosrr.cc
//
// Copyright (C) 2001 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
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

#include <chemistry/qc/libint2/int1e.h>
#include <chemistry/qc/basis/fjt.h>

using namespace std;
using namespace sc;

inline void fail()
{
  ExEnv::errn() << scprintf("failing module:\n%s",__FILE__) << endl;
  abort();
}

void Int1eLibint2::OI_OSrecurs_(double **OIX, double **OIY, double **OIZ, double PA[3], double PB[3],
			      double gamma, int lmaxi, int lmaxj)
{
  int i,j,k;
  double pp = 1/(2*gamma);

  OIX[0][0] = OIY[0][0] = OIZ[0][0] = 1.0;

	/* Upward recursion in j for i=0 */

  OIX[0][1] = PB[0];
  OIY[0][1] = PB[1];
  OIZ[0][1] = PB[2];

  for(j=1;j<lmaxj;j++) {
    OIX[0][j+1] = PB[0]*OIX[0][j];
    OIY[0][j+1] = PB[1]*OIY[0][j];
    OIZ[0][j+1] = PB[2]*OIZ[0][j];
    OIX[0][j+1] += j*pp*OIX[0][j-1];
    OIY[0][j+1] += j*pp*OIY[0][j-1];
    OIZ[0][j+1] += j*pp*OIZ[0][j-1];
  }

	/* Upward recursion in i for all j's */

  OIX[1][0] = PA[0];
  OIY[1][0] = PA[1];
  OIZ[1][0] = PA[2];
  for(j=1;j<=lmaxj;j++) {
    OIX[1][j] = PA[0]*OIX[0][j];
    OIY[1][j] = PA[1]*OIY[0][j];
    OIZ[1][j] = PA[2]*OIZ[0][j];
    OIX[1][j] += j*pp*OIX[0][j-1];
    OIY[1][j] += j*pp*OIY[0][j-1];
    OIZ[1][j] += j*pp*OIZ[0][j-1];
  }
  for(i=1;i<lmaxi;i++) {
    OIX[i+1][0] = PA[0]*OIX[i][0];
    OIY[i+1][0] = PA[1]*OIY[i][0];
    OIZ[i+1][0] = PA[2]*OIZ[i][0];
    OIX[i+1][0] += i*pp*OIX[i-1][0];
    OIY[i+1][0] += i*pp*OIY[i-1][0];
    OIZ[i+1][0] += i*pp*OIZ[i-1][0];
    for(j=1;j<=lmaxj;j++) {
      OIX[i+1][j] = PA[0]*OIX[i][j];
      OIY[i+1][j] = PA[1]*OIY[i][j];
      OIZ[i+1][j] = PA[2]*OIZ[i][j];
      OIX[i+1][j] += i*pp*OIX[i-1][j];
      OIY[i+1][j] += i*pp*OIY[i-1][j];
      OIZ[i+1][j] += i*pp*OIZ[i-1][j];
      OIX[i+1][j] += j*pp*OIX[i][j-1];
      OIY[i+1][j] += j*pp*OIY[i][j-1];
      OIZ[i+1][j] += j*pp*OIZ[i][j-1];
    }
  }

  return;
}
