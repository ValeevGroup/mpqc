//
// g92test.cc
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

#include <string.h>

#include <sys/stat.h>
#include <unistd.h>

#include <util/keyval/keyval.h>

#include "g92.h"

#include <chemistry/molecule/coor.h>
#include <math/optimize/qnewton.h>
#include <math/optimize/gdiis.h>
#include <math/optimize/efc.h>
#include <math/optimize/update.h>

// Force linkages:
#ifndef __PIC__
const ClassDesc &fl0  = Gaussian92SCF::class_desc_;
const ClassDesc &fl0a = Gaussian92UHF::class_desc_;

const ClassDesc &fl1 = IntMolecularCoor::class_desc_;
const ClassDesc &fl2 = QNewtonOpt::class_desc_;
const ClassDesc &fl3 = GDIISOpt::class_desc_;
const ClassDesc &fl4 = EFCOpt::class_desc_;
const ClassDesc &fl5 = BFGSUpdate::class_desc_;
#endif

main(int argc, char**argv)
{
  printf("\n     SC G92 driver program\n\n");

  // the output stream is standard out
  SCostream& o = SCostream::cout;

  char *input =      (argc > 1)? argv[1] : SRCDIR "/mpqc.in";
  char *keyword =    (argc > 2)? argv[2] : "mole";
  char *optkeyword = (argc > 3)? argv[3] : "opt";

  struct stat sb;
  RefMolecularEnergy mole;
  RefOptimize opt;

  if (stat("g92test.ckpt",&sb)==0 && sb.st_size) {
    //StateInText si("g92test.ckpt");
    StateInBin si("g92test.ckpt");
    opt.restore_state(si);
    mole = opt->nlp();
  } else {
    // open keyval input
    RefKeyVal rpkv(new ParsedKeyVal(input));

    mole = rpkv->describedclassvalue(keyword);
    opt = rpkv->describedclassvalue(optkeyword);
    opt->set_checkpoint();
    opt->set_checkpoint_file("g92test.ckpt");
  }

  if (mole->gradient_implemented()) {
    if (opt.nonnull()) {
      opt->optimize();
    } else {
      o << "opt is null\n";
    }
  }
  
  mole->hessian();
  Gaussian92::castdown(mole)->normal_modes()->print("normal modes");
  Gaussian92::castdown(mole)->frequencies()->print("frequencies");

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
