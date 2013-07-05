//
// scftest.cc --- test program
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifdef HAVE_CONFIG_H
#include <mpqc_config.h>
#endif

#include <string.h>

#include <sys/stat.h>
#include <unistd.h>
#include <new>

#include <util/keyval/keyval.h>
#include <util/group/message.h>
#include <util/group/pregtime.h>
#include <util/misc/bug.h>
#include <util/misc/formio.h>
#include <util/state/state_bin.h>

#include <math/optimize/opt.h>

#include <chemistry/molecule/coor.h>
#include <chemistry/molecule/energy.h>

// Force linkages:
#include <chemistry/qc/scf/linkage.h>

using namespace std;
using namespace sc;

Ref<MessageGrp> grp;

static Ref<MessageGrp>
init_mp(const Ref<KeyVal>& keyval)
{
  // if we are on a paragon then use a ParagonMessageGrp
  // otherwise read the message group from the input file
  grp << keyval->describedclassvalue("message");

  if (grp.nonnull()) MessageGrp::set_default_messagegrp(grp);
  else grp = MessageGrp::get_default_messagegrp();

  Ref<Debugger> debugger; debugger << keyval->describedclassvalue(":debug");
  // Let the debugger know the name of the executable and the node
  if (debugger.nonnull()) {
    debugger->set_exec("scftest");
    debugger->set_prefix(grp->me());
    debugger->debug("curt is a hog");
  }
  
  RegionTimer::set_default_regiontimer(
    new ParallelRegionTimer(grp,"scftest",1,0));

  SCFormIO::set_printnode(0);
  //SCFormIO::set_debug(1);

  SCFormIO::setindent(cout, 2);
  SCFormIO::setindent(cerr, 2);
  
  return grp;
}

int
main(int argc, char**argv)
{
  const char *input =      (argc > 1)? argv[1] : SRCDIR "/mpqc.in";
  const char *keyword =    (argc > 2)? argv[2] : "mole";
  const char *optkeyword = (argc > 3)? argv[3] : "opt";

  // open keyval input
  Ref<KeyVal> rpkv(new ParsedKeyVal(input));

  init_mp(rpkv);

  Timer tim;
  tim.enter("input");
  
  if (rpkv->exists("matrixkit")) {
    Ref<SCMatrixKit> kit; kit << rpkv->describedclassvalue("matrixkit");
    SCMatrixKit::set_default_matrixkit(kit);
  }
  
  struct stat sb;
  Ref<MolecularEnergy> mole;
  Ref<Optimize> opt;

  if (stat("scftest.ckpt",&sb)==0 && sb.st_size) {
    StateInBin si("scftest.ckpt");
    opt << SavableState::restore_state(si);
    mole << opt->function();
  } else {
    mole << rpkv->describedclassvalue(keyword);
    opt << rpkv->describedclassvalue(optkeyword);
    if (opt.nonnull()) {
      opt->set_checkpoint();
      opt->set_checkpoint_file("scftest.ckpt");
    }
  }

  tim.exit("input");

  if (mole.nonnull()) {
    if (mole->gradient_implemented()) {
      if (opt.nonnull()) {
        opt->optimize();
      } else {
        mole->gradient().print("gradient");
      }
    } else if (mole->value_implemented()) {
      ExEnv::out0() << indent
           << scprintf("value of mole is %15.10f\n\n", mole->energy());
    }
  }

  mole->print(ExEnv::out0());

  StateOutBin so("scftest.wfn");
  SavableState::save_state(mole.pointer(),so);
  
  tim.print(ExEnv::out0());

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
