//
// etraintest.cc
//
// Copyright (C) 2011 Edward Valeev
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

#include <util/keyval/keyval.h>
#include <util/group/message.h>
#include <util/group/pregtime.h>
#include <util/group/thread.h>
#include <util/group/memory.h>
#include <util/misc/bug.h>
#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <util/state/state_bin.h>
#include <chemistry/qc/etrain/etrain.h>

// Force linkages:
#include <util/group/linkage.h>
#include <chemistry/qc/basis/linkage.h>
#include <chemistry/qc/wfn/linkage.h>
#include <chemistry/qc/scf/linkage.h>
#include <chemistry/qc/dft/linkage.h>
#include <util/state/linkage.h>

using namespace sc;
using namespace std;

// globals and extern
Ref<RegionTimer> tim;

void
do_main(const Ref<KeyVal>& rpkv, const Ref<RegionTimer>& tim)
{
  const char* keyword = "etrain";
  tim->enter("test");
  tim->enter("input");
  Ref<ETraIn> etrain;
  etrain << rpkv->describedclassvalue(keyword);
  tim->exit("input");
  if (etrain.null())
      throw InputError("ETraIn object not found in the input",__FILE__,__LINE__);

  const double tmp = etrain->value();
  etrain->print(ExEnv::out0());

  tim->exit("test");
}

static void
init(int argc, char** argv, const char* input, const char* output, const Ref<KeyVal>& keyval)
{
  // get the thread group.  first try the commandline and environment
  Ref<ThreadGrp> thrgrp = ThreadGrp::initial_threadgrp(argc, argv);

  // if we still don't have a group, try reading the thread group
  // from the input
  if (thrgrp.null()) {
    thrgrp << keyval->describedclassvalue("threadgrp");
  }

  if (thrgrp)
    ThreadGrp::set_default_threadgrp(thrgrp);
  else
    thrgrp = ThreadGrp::get_default_threadgrp();

  Ref<MessageGrp> msggrp = MessageGrp::initial_messagegrp(argc, argv);
  if (msggrp.null()) msggrp << keyval->describedclassvalue("messagegrp");
  if (msggrp)
    MessageGrp::set_default_messagegrp(msggrp);
  else
    msggrp = MessageGrp::get_default_messagegrp();

  if (keyval->exists("matrixkit")) {
    Ref<SCMatrixKit> kit; kit << keyval->describedclassvalue("matrixkit");
    SCMatrixKit::set_default_matrixkit(kit);
  }

  // get the integral factory. first try commandline and environment
  Ref<Integral> integral = Integral::initial_integral(argc, argv);
  // if we still don't have a integral, try reading the integral
  // from the input
  if (integral.null()) {
    integral << keyval->describedclassvalue("integrals");
  }
  if (integral)
    Integral::set_default_integral(integral);
  else
    integral = Integral::get_default_integral();
  ExEnv::out0() << endl << indent
       << "Using " << integral->class_name()
       << " by default for molecular integrals evaluation" << endl << endl;


  Ref<Debugger> debugger; debugger << keyval->describedclassvalue("debug");
  // Let the debugger know the name of the executable and the node
  if (debugger) {
    debugger->set_exec("mpqc-test");
    debugger->set_prefix(msggrp->me());
    debugger->debug("curt is a hog");
  }

  tim = new ParallelRegionTimer(msggrp,"mpqc-test",1,0);
  RegionTimer::set_default_regiontimer(tim);

  SCFormIO::set_printnode(0);
  SCFormIO::init_mp(msggrp->me());
  //SCFormIO::set_debug(1);

  SCFormIO::setindent(ExEnv::out0(), 2);
  SCFormIO::setindent(cerr, 2);

  // get the basename for output files
  const char *basename_source;
  if (output)
    basename_source = output;
  else
    basename_source = input;
  int nfilebase = (int) (::strrchr(basename_source, '.') - basename_source);
  char *basename = new char[nfilebase + 1];
  strncpy(basename, basename_source, nfilebase);
  basename[nfilebase] = '\0';
  SCFormIO::set_default_basename(basename);
}

int main(int argc, char**argv)
{
  char *infile = new char[strlen(SRCDIR) + strlen("/etraintest.in") + 1];
  sprintf(infile,SRCDIR "/etraintest.in");
  if (argc == 2) {
    delete[] infile;
    infile = argv[1];
  }

  // open keyval input
  Ref<KeyVal> rpkv(new ParsedKeyVal(infile));

  init(argc, argv, infile, 0, rpkv);

  ///////////
  //
  // Do the work now
  //
  ///////////

  do_main(rpkv,tim);

  tim->print(ExEnv::out0());

  tim = 0;
  RegionTimer::set_default_regiontimer(0);
  MessageGrp::set_default_messagegrp(0);
  ThreadGrp::set_default_threadgrp(0);
  MemoryGrp::set_default_memorygrp(0);

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
