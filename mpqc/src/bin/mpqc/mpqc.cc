//
// mpqc.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
// Maintainer: LPS
//
// This file is part of MPQC.
//
// MPQC is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// MPQC is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifdef HAVE_CONFIG_H
#include <scconfig.h>
#endif

#include <new.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fstream.h>
#include <strstream.h>

#include <util/options/GetLongOpt.h>
#include <util/misc/newstring.h>
#include <util/keyval/keyval.h>
#include <util/group/message.h>
#include <util/group/thread.h>
#include <util/group/pregtime.h>
#include <util/misc/bug.h>
#include <util/misc/formio.h>
#include <util/misc/exenv.h>

#include <math/optimize/opt.h>

#include <chemistry/molecule/coor.h>
#include <chemistry/molecule/energy.h>
#include <chemistry/molecule/molfreq.h>
#include <chemistry/qc/wfn/wfn.h>

// Force linkages:
#include <chemistry/qc/mbpt/linkage.h>
//#include <chemistry/qc/psi/linkage.h>
#include <util/state/linkage.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "disclaimer.h"

//////////////////////////////////////////////////////////////////////////

static void
clean_up(void)
{
  MessageGrp::set_default_messagegrp(0);
  SCMatrixKit::set_default_matrixkit(0);
  RegionTimer::set_default_regiontimer(0);
}

int
main(int argc, char *argv[])
{
  atexit(clean_up);

  ExEnv::set_args(argc, argv);

#ifdef HAVE_MPI
  // MPI is allowed wait until MPI_Init to fill in argc and argv,
  // so we may have to call MPI_Init before we even know that we
  // want an MPIMessageGrp.  The command name is used to let mpqc
  // know that an early init is needed.
  if (!strcmp(ExEnv::program_name(), "mpqc-mpi")) {
    MPI_Init(&argc, &argv);
  }
#endif

  // parse commandline options
  GetLongOpt options;

  options.enroll("f", GetLongOpt::MandatoryValue,
                 "the name of the input file", "mpqc.in");
  options.enroll("messagegrp", GetLongOpt::MandatoryValue,
                 "which message group to use", 0);
  options.enroll("threadgrp", GetLongOpt::MandatoryValue,
                 "which thread group to use", 0);
  options.enroll("memorygrp", GetLongOpt::MandatoryValue,
                 "which thread group to use", 0);
  options.enroll("l", GetLongOpt::MandatoryValue, "basis set limit", "0");
  options.enroll("W", GetLongOpt::MandatoryValue,
                 "set the working directory", ".");
  options.enroll("c", GetLongOpt::NoValue, "check input then exit", 0);
  options.enroll("v", GetLongOpt::NoValue, "print the version number", 0);
  options.enroll("w", GetLongOpt::NoValue, "print the warranty", 0);
  options.enroll("d", GetLongOpt::NoValue, "debug", 0);
  options.enroll("h", GetLongOpt::NoValue, "print this message", 0);

  options.parse(argc, argv);

  // set the working dir
  if (strcmp(options.retrieve("W"),"."))
    chdir(options.retrieve("W"));

  // get the message group.  first try the commandline and environment
  RefMessageGrp grp = MessageGrp::initial_messagegrp(argc, argv);
  if (grp.nonnull())
    MessageGrp::set_default_messagegrp(grp);
  else
    grp = MessageGrp::get_default_messagegrp();

  // initialize keyval input
  const char *input = options.retrieve("f");
  RefParsedKeyVal parsedkv;
  if (grp->n() == 1) {
    parsedkv = new ParsedKeyVal(input);
  }
  else {
    // read the input file on only node 0
    parsedkv = new ParsedKeyVal();
    char *in_char_array;
    if (grp->me() == 0) {
      ifstream is(input);
      ostrstream ostrs;
      is >> ostrs.rdbuf();
      ostrs << ends;
      in_char_array = ostrs.str();
      int n = ostrs.pcount();
      grp->bcast(n);
      grp->bcast(in_char_array, n);
    }
    else {
      int n;
      grp->bcast(n);
      in_char_array = new char[n];
      grp->bcast(in_char_array, n);
    }
    parsedkv->parse_string(in_char_array);
    delete[] in_char_array;
  }
  RefKeyVal pkv = parsedkv.pointer();
  RefKeyVal ppkv(new AggregateKeyVal(new PrefixKeyVal(":mpqc",pkv),
                                     new PrefixKeyVal(":default",pkv)));
  pkv = new ParsedKeyVal("input",ppkv);
  RefKeyVal keyval = new AggregateKeyVal(ppkv,pkv);

  pkv = ppkv = 0;

  // get the basename for output files
  int nfilebase = (int) (strrchr(input, '.') - input);
  char *basename = new char[nfilebase + 1];
  strncpy(basename, input, nfilebase);
  basename[nfilebase] = '\0';
  SCFormIO::set_default_basename(basename);

  // get the thread group.  first try the commandline and environment
  RefThreadGrp thread = ThreadGrp::initial_threadgrp(argc, argv);
  
  // if we still don't have a group, try reading the thread group
  // from the input
  if (thread.null()) {
    thread = keyval->describedclassvalue("thread");
  }

  if (thread.nonnull())
    ThreadGrp::set_default_threadgrp(thread);
  else
    thread = ThreadGrp::get_default_threadgrp();

  // set up output classes
  SCFormIO::setindent(cout, 2);
  SCFormIO::setindent(cerr, 2);

  SCFormIO::set_printnode(0);
  if (grp->n() > 1)
    SCFormIO::init_mp(grp->me());

  if (options.retrieve("d"))
    SCFormIO::set_debug(1);

  if (options.retrieve("h")) {
    cout << node0 << endl
         << indent << "MPQC version 1.0" << endl << endl;
    options.usage(cout);
    exit(0);
  }
  
  if (options.retrieve("v")) {
    cout << node0 << endl
         << indent << "MPQC version 1.0" << endl << endl;
    exit(0);
  }
  
  if (options.retrieve("w")) {
    cout << node0 << endl
         << indent << "MPQC version 1.0" << endl << endl;
    print_disclaimer(cout);
    exit(0);
  }

  // initialize timing for mpqc
  RefRegionTimer tim = new ParallelRegionTimer(grp,"mpqc",1,1);
  RegionTimer::set_default_regiontimer(tim);

  tim->enter("input");
  
  // now set up the debugger
  RefDebugger debugger = keyval->describedclassvalue("debug");
  if (debugger.nonnull()) {
    debugger->set_exec(argv[0]);
    debugger->set_prefix(grp->me());
    if (options.retrieve("d"))
      debugger->debug("Starting debugger because -d given on command line.");
  }

  // now check to see what matrix kit to use
  if (keyval->exists("matrixkit"))
    SCMatrixKit::set_default_matrixkit(
      keyval->describedclassvalue("matrixkit"));
  
  // announce ourselves
  cout << node0 << endl
       << indent << "       MPQC: Massively Parallel Quantum Chemistry\n\n\n"
       << indent
       << scprintf("Running on a %s with %d nodes.", TARGET_ARCH, grp->n())
       << endl << endl;

  // see if frequencies are wanted
  RefMolecularFrequencies molfreq = keyval->describedclassvalue("freq");
  
  // check for a molecular energy and optimizer
  char * molname = keyval->pcharvalue("filename");
  if (!molname)
    molname = new_string(basename);
  
  char * ckptfile = new char[strlen(molname)+6];
  sprintf(ckptfile,"%s.ckpt",molname);
  
  int restart = 1;
  if (keyval->exists("restart"))
    restart = keyval->booleanvalue("restart");

  int checkpoint = keyval->booleanvalue("checkpoint");
  if (keyval->error() != KeyVal::OK)
    checkpoint=1;

  int savestate = keyval->booleanvalue("savestate");
  if (keyval->error() != KeyVal::OK)
    savestate=1;

  struct stat sb;
  RefMolecularEnergy mole;
  RefOptimize opt;

  if (restart && stat(ckptfile,&sb)==0 && sb.st_size) {
    StateInBinXDR si(ckptfile);
    opt.restore_state(si);
    mole = opt->function();
  } else {
    mole = keyval->describedclassvalue("mole");
    opt = keyval->describedclassvalue("opt");
    if (checkpoint && opt.nonnull()) {
      opt->set_checkpoint();
      opt->set_checkpoint_file(ckptfile);
    }
  }

  delete[] ckptfile;

  int check = (options.retrieve("c") != 0);
  int limit = atoi(options.retrieve("l"));
  if (limit) {
    RefWavefunction wfn(mole);
    if (wfn.nonnull() && wfn->basis_dimension()->n() > limit) {
      cerr << node0 << endl << indent
           << "The limit of " << limit << " basis functions has been exceeded."
           << endl;
      check = 1;
    }
  }

  if (check) {
    cout << node0 << endl << indent
         << "Exiting since the check option is on." << endl;
    exit(0);
  }
  
  tim->change("calc");

  int do_energy = keyval->booleanvalue("do_energy");
  if (keyval->error() != KeyVal::OK)
    do_energy=1;
  
  int do_grad = keyval->booleanvalue("do_gradient");
  if (keyval->error() != KeyVal::OK)
    do_grad=0;

  int do_opt = keyval->booleanvalue("optimize");
  if (keyval->error() != KeyVal::OK)
    do_opt=1;
  
  int do_pdb = keyval->booleanvalue("write_pdb");
  if (keyval->error() != KeyVal::OK)
    do_pdb=0;
  
  int print_mole = keyval->booleanvalue("print_mole");
  if (keyval->error() != KeyVal::OK)
    print_mole=1;
  
  int ready_for_freq = 1;
  if (mole.nonnull()) {
    if ((do_opt || do_grad) && !mole->gradient_implemented()) {
      cout << node0 << indent
           << "WARNING: optimization or gradient requested but the given"
           << endl
           << "         MolecularEnergy object cannot do gradients."
           << endl;
    }

    if (do_opt && opt.nonnull() && mole->gradient_implemented()) {
      int result = opt->optimize();
      if (result) {
        cout << node0 << indent
             << "The optimization has converged." << endl << endl;
      } else {
        cout << node0 << indent
             << "The optimization has NOT converged." << endl << endl;
        ready_for_freq = 0;
      }
    } else if (do_grad && mole->gradient_implemented()) {
      mole->gradient().print("gradient");
    } else if (do_energy && mole->value_implemented()) {
      cout << node0 << indent
           << scprintf("value of mole is %20.15f\n", mole->energy())
           << endl;
    }
  }

  tim->exit("calc");

  if (ready_for_freq && molfreq.nonnull()) {
    tim->enter("frequencies");
    molfreq->compute_displacements();
    cout << node0 << indent
         << "Computing molecular frequencies from "
         << molfreq->ndisplace() << " displacements:" << endl;

    for (int i=0; i<molfreq->ndisplace(); i++) {
      // This produces side-effects in mol and may even change
      // its symmetry.
      cout << node0
           << "Beginning displacement " << i << ":" << endl;
      molfreq->displace(i);

      mole->obsolete();
      RefSCVector gradv = mole->get_cartesian_gradient();
      molfreq->set_gradient(i, gradv);
    }

    molfreq->compute_frequencies_from_gradients();
    //molfreq->thermochemistry(scf_info.nopen+1);
    molfreq->thermochemistry(1);
    tim->exit("frequencies");
  }

  if (mole.nonnull()) {
    if (print_mole)
      mole->print(cout);

    if (do_pdb) {
      ckptfile = new char[strlen(molname)+5];
      sprintf(ckptfile, "%s.pdb", molname);
      ofstream pdbfile(ckptfile);
      mole->molecule()->print_pdb(pdbfile);
      delete[] ckptfile;
    }
    
  }
  else {
    cout << node0 << "mpqc: The molecular energy object is null" << endl
         << " make sure \"mole\" specificies a MolecularEnergy derivative"
         << endl;
  }

  if (savestate) {
    if (opt.nonnull()) {
      ckptfile = new char[strlen(molname)+6];
      sprintf(ckptfile,"%s.ckpt",molname);

      StateOutBinXDR so(ckptfile);
      opt.save_state(so);
      so.close();

      delete[] ckptfile;
    }

    if (mole.nonnull()) {
      ckptfile = new char[strlen(molname)+5];
      sprintf(ckptfile, "%s.wfn",molname);
  
      StateOutBinXDR so(ckptfile);
      mole.save_state(so);
      so.close();

      delete[] ckptfile;
    }
  }
  
  tim->print(cout);

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
