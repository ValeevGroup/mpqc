//
// comp2e3c.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
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

/* These routines compute two and three center electron repulsion
 * integrals. */

#include <stdexcept>

#include <stdlib.h>
#include <math.h>

#include <chemistry/qc/intv3/flags.h>
#include <chemistry/qc/intv3/macros.h>
#include <chemistry/qc/intv3/types.h>

#include <chemistry/qc/intv3/int2e.h>

using namespace sc;

void
Int2eV3::make_int_unit_shell()
{
  double *exp = new double[1];
  int *am = new int[1];
  int *pure = new int[1];
  double **c = new double*[1];
  *c = new double[1];

  exp[0] = 0.0;
  am[0] = 0;
  pure[0] = 0;
  c[0][0] = 1.0;

  int_unit_shell = new GaussianShell(1,1,exp,am,pure,c,
                                     GaussianShell::Unnormalized,
                                     false);
}

void
Int2eV3::delete_int_unit_shell()
{
  delete int_unit_shell;
  int_unit_shell = 0;
}

/* Compute a 2 center electron repulsion integral.  Electron 1 is in
 * shell psh1 and electron 2 is in psh2, that is (1 | 2).  To avoid
 * confusing the user of these routines, the INT_NOPERM is set.
 */
void
Int2eV3::erep_2center(int &psh1, int &psh2)
{
  if (bs2_ || bs4_) {
      throw std::runtime_error("erep_2center: bs2 or bs4 not null");
    }
  //int shd = 0x11111111; /* a dummy shell that will cause death if used */
  int shd = 0; // shell = 0 is used so intermediate lookup will work
  int oldperm = permute();
  set_permute(0);
  erep(psh1,shd,psh2,shd);
  set_permute(oldperm);
}

/* This is an alternate interface to int_erep2.  It takes
 * as arguments the flags, an integer vector of shell numbers
 * and an integer vector which will be filled in with size
 * information, if it is non-NULL. */
void
Int2eV3::erep_2center(int *shells, int  *sizes)
{
  erep_2center(shells[0],shells[1]);
  if (sizes) {
      sizes[0] = bs1_->shell(shells[0]).nfunction();
      sizes[1] = bs3_->shell(shells[1]).nfunction();
    }
}


/* Computes a 3 center two electron integral.  Electron 1 is in psh1
 * and electron 2 is in psh2 and psh3, that is (1 2 | 3).  To avoid
 * confusing the user of these routines, the INT_NOPERM is set.
 */
void
Int2eV3::erep_3center(int &psh1, int &psh2, int &psh3)
{
  if (bs4_) {
      throw std::runtime_error("erep_3center: bs4 not null");
    }
  //int shd = 0x11111111; /* a dummy shell that will cause death if used */
  int shd = 0; // shell = 0 is used so intermediate lookup will work
  int oldperm = permute();
  set_permute(0);
  erep(psh1,psh2,psh3,shd);
  set_permute(oldperm);
}

/* This is an alternate interface to int_erep3.  It takes
 * as arguments the flags, an integer vector of shell numbers
 * and an integer vector which will be filled in with size
 * information, if it is non-NULL. */
void
Int2eV3::erep_3center(int *shells, int  *sizes)
{
  erep_3center(shells[0],shells[1],shells[2]);
  if (sizes) {
      sizes[0] = bs1_->shell(shells[0]).nfunction();
      sizes[1] = bs2_->shell(shells[1]).nfunction();
      sizes[2] = bs3_->shell(shells[2]).nfunction();
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
