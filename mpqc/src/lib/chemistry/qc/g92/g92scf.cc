//
// g92scf.cc
//
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

#include "g92.h"

////////////////////////////////////////////////////////////////////////////

#define CLASSNAME Gaussian92SCF
#define PARENTS public Gaussian92
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
Gaussian92SCF::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = Gaussian92::_castdown(cd);
  return do_castdowns(casts,cd);
}

char *
Gaussian92SCF::emethod()
{
  static char method[32];
  int conv = (int) -log10(desired_value_accuracy());

  // closed shell
  if (multiplicity_==1)
    sprintf(method,"rhf scf=direct scfcon=%d",conv);

  // open shell
  else
    sprintf(method,"rohf scf=direct scfcon=%d",conv);

  return method;
}

char *
Gaussian92SCF::gmethod()
{
  static char method[36];
  int conv = (int) -log10(desired_value_accuracy());
  
  // closed shell
  if (multiplicity_==1)
    sprintf(method,"rhf force scf=direct scfcon=%d",conv);

  // open shell
  else
    sprintf(method,"rohf force scf=direct scfcon=%d",conv);

  return method;
}

char *
Gaussian92SCF::hmethod()
{
  static char method[48];
  int conv = (int) -log10(desired_value_accuracy());
  int hconv = (int) -log10(desired_hessian_accuracy());
  
  // closed shell
  if (multiplicity_==1)
    sprintf(method,"rhf freq scf=direct scfcon=%d cphf=conv=%d",conv,hconv);
  
  // open shell
  else
    sprintf(method,"rohf freq scf=direct scfcon=%d cphf=conv=%d",conv,hconv);

  return method;
}

Gaussian92SCF::Gaussian92SCF(const RefKeyVal&keyval):
  Gaussian92(keyval)
{
  if (!basis_) {
    fprintf(stderr,"Gaussian92SCF needs a basis\n");
    abort();
  }
}

Gaussian92SCF::~Gaussian92SCF()
{
}

Gaussian92SCF::Gaussian92SCF(StateIn&s) :
  Gaussian92(s)
  maybe_SavableState(s)
{
}

void
Gaussian92SCF::save_data_state(StateOut&s)
{
  Gaussian92::save_data_state(s);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")
// End:
