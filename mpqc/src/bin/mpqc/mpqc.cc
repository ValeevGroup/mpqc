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
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with the MPQC; see the file COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifdef HAVE_CONFIG_H
#include <scconfig.h>
#endif

#ifdef HAVE_TIME_H
#include <time.h>
#endif

#include <scdirlist.h>

#include <new.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fstream.h>
#include <strstream.h>

#include <util/options/GetLongOpt.h>
#include <util/misc/newstring.h>
#include <util/keyval/keyval.h>
#include <util/state/state_bin.h>
#include <util/group/message.h>
#include <util/group/mstate.h>
#include <util/group/thread.h>
#include <util/group/pregtime.h>
#include <util/misc/bug.h>
#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <util/render/render.h>

#include <math/optimize/opt.h>

#include <chemistry/molecule/coor.h>
#include <chemistry/molecule/energy.h>
#include <chemistry/molecule/molfreq.h>
#include <chemistry/molecule/fdhess.h>
#include <chemistry/molecule/formula.h>
#include <chemistry/qc/wfn/wfn.h>

// Force linkages:
#include <chemistry/qc/wfn/linkage.h>
#include <chemistry/qc/scf/linkage.h>
#include <chemistry/qc/dft/linkage.h>
#include <chemistry/qc/mbpt/linkage.h>
//#include <chemistry/qc/psi/linkage.h>
#include <util/state/linkage.h>
#ifdef HAVE_SC_SRC_LIB_CHEMISTRY_QC_CC
#  include <chemistry/qc/cc/linkage.h>
#endif
#ifdef HAVE_SC_SRC_LIB_CHEMISTRY_QC_PSI
#  include <chemistry/qc/psi/linkage.h>
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "version.h"

//////////////////////////////////////////////////////////////////////////

static void
trash_stack_b(int &i, char *&ichar)
{
  char stack;
  ichar = &stack;
  ichar -= 10;
  for (i=0; i<1000; i++) {
    *ichar-- = 0xfe;
  }
}

static void
trash_stack()
{
  int i;
  char *ichar;
  trash_stack_b(i,ichar);
}

static void
clean_up(void)
{
  MessageGrp::set_default_messagegrp(0);
  ThreadGrp::set_default_threadgrp(0);
  SCMatrixKit::set_default_matrixkit(0);
  RegionTimer::set_default_regiontimer(0);
}

static void
out_of_memory()
{
  cerr << "ERROR: mpqc: out of memory" << endl;
  cout.flush();
  cerr.flush();
  abort();
}

#if defined(__alpha__)
#include <signal.h>
#include <asm/fpu.h>
extern "C" {
  void ieee_set_fp_control(long);
  long ieee_get_fp_control(void);
};
void sigfpe_handler(int)
{
  long fp_control = ieee_get_fp_control();
  int fatal = 0;
  if (fp_control & IEEE_STATUS_INV) {
      cout << "SGIFPE: invalid operation" << endl;;
      fatal = 1;
    }
  if (fp_control & IEEE_STATUS_DZE) {
      cout << "SGIFPE: divide by zero" << endl;;
      fatal = 1;
    }
  if (fp_control & IEEE_STATUS_OVF) {
      cout << "SGIFPE: overflow" << endl;;
      fatal = 1;
    }
  if (fp_control & IEEE_STATUS_UNF) {
      //cout << "SGIFPE: underflow" << endl;;
    }
  if (fp_control & IEEE_STATUS_INE) {
      //cout << "SGIFPE: inexact" << endl;;
    }
  if (fatal) abort();
}
#endif

int
main(int argc, char *argv[])
{
  //trash_stack();

  KeyValValueboolean truevalue(1), falsevalue(0);
  int i;
  const char *devnull = "/dev/null";
  atexit(clean_up);
  set_new_handler(out_of_memory);

#if defined(__i386__) && defined(__GNUC__)
  // make floating point errors cause an exception (except for denormalized
  // operands, since small numbers are denormalized)
  asm("fldcw %0" : : "o" (0x372));
#endif

#if defined(__alpha__)
  signal(SIGFPE,sigfpe_handler);
#endif

  ExEnv::init(argc, argv);

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
  options.enroll("L", GetLongOpt::NoValue, "print the license", 0);
  options.enroll("k", GetLongOpt::NoValue, "print key/value assignments", 0);
  options.enroll("d", GetLongOpt::NoValue, "debug", 0);
  options.enroll("h", GetLongOpt::NoValue, "print this message", 0);

  options.parse(argc, argv);

  if (options.retrieve("h")) {
    cout << node0
         << indent << "MPQC version " << MPQC_VERSION << endl
         << indent << "compiled for " << TARGET_ARCH << endl
         << SCFormIO::copyright << endl;
    options.usage(cout);
    exit(0);
  }
  
  if (options.retrieve("v")) {
    cout << node0
         << indent << "MPQC version " << MPQC_VERSION << endl
         << indent << "compiled for " << TARGET_ARCH << endl
         << SCFormIO::copyright;
    exit(0);
  }
  
  if (options.retrieve("w")) {
    cout << node0
         << indent << "MPQC version " << MPQC_VERSION << endl
         << indent << "compiled for " << TARGET_ARCH << endl
         << SCFormIO::copyright << endl
         << SCFormIO::warranty;
    exit(0);
  }
  
  if (options.retrieve("L")) {
    cout << node0
         << indent << "MPQC version " << MPQC_VERSION << endl
         << indent << "compiled for " << TARGET_ARCH << endl
         << SCFormIO::copyright << endl
         << SCFormIO::license;
    exit(0);
  }

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

  if (options.retrieve("k")) parsedkv->verbose(1);
  RefKeyVal keyval = new PrefixKeyVal("mpqc",parsedkv.pointer());

  // get the basename for output files
  int nfilebase = (int) (strrchr(input, '.') - input);
  char *basename = new char[nfilebase + 1];
  strncpy(basename, input, nfilebase);
  basename[nfilebase] = '\0';
  SCFormIO::set_default_basename(basename);

  // set up output classes
  SCFormIO::setindent(cout, 2);
  SCFormIO::setindent(cerr, 2);

  SCFormIO::set_printnode(0);
  if (grp->n() > 1)
    SCFormIO::init_mp(grp->me());

  if (options.retrieve("d"))
    SCFormIO::set_debug(1);

  // initialize timing for mpqc
  RefRegionTimer tim;
  if (keyval->exists("timer")) tim = keyval->describedclassvalue("timer");
  else                         tim = new ParallelRegionTimer(grp,"mpqc",1,1);
  RegionTimer::set_default_regiontimer(tim);

  if (tim.nonnull()) tim->enter("input");
  
  // announce ourselves
  const char title1[] = "MPQC: Massively Parallel Quantum Chemistry";
  int ntitle1 = sizeof(title1);
  const char title2[] = "Version " MPQC_VERSION;
  int ntitle2 = sizeof(title2);
  cout << node0 << endl;
  cout << node0 << indent;
  for (i=0; i<(80-ntitle1)/2; i++) cout << node0 << ' ';
  cout << node0 << title1 << endl;
  cout << node0 << indent;
  for (i=0; i<(80-ntitle2)/2; i++) cout << node0 << ' ';
  cout << node0 << title2 << endl << endl;

  const char *tstr = 0;
#if defined(HAVE_TIME) && defined(HAVE_CTIME)
  time_t t;
  time(&t);
  tstr = ctime(&t);
#endif
  if (!tstr) {
    tstr = "UNKNOWN";
  }

  cout << node0
       << indent << scprintf("Architecture is %s.", TARGET_ARCH) << endl
       << indent << scprintf("Hostname is %s.", ExEnv::hostname()) << endl
       << indent << scprintf("Username is %s.", ExEnv::username()) << endl
       << indent << scprintf("Time is %s", tstr);
  cout.flush();

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

  cout << node0 << indent
       << "Using " << grp->class_name()
       << " for message passing (number of nodes = " << grp->n() << ")." << endl
       << indent
       << "Using " << thread->class_name()
       << " for threading (number of threads = " << thread->nthread() << ")." << endl
       << indent
       << "Total number of processors = " << grp->n() * thread->nthread() << endl;

  // now set up the debugger
  RefDebugger debugger = keyval->describedclassvalue("debug");
  if (debugger.nonnull()) {
    Debugger::set_default_debugger(debugger);
    debugger->set_exec(argv[0]);
    debugger->set_prefix(grp->me());
    if (options.retrieve("d"))
      debugger->debug("Starting debugger because -d given on command line.");
  }

  // now check to see what matrix kit to use
  if (keyval->exists("matrixkit"))
    SCMatrixKit::set_default_matrixkit(
      keyval->describedclassvalue("matrixkit"));

  // check for a molecular energy and optimizer
  KeyValValueString molnamedef(basename);
  char * molname = keyval->pcharvalue("filename", molnamedef);
  if (strcmp(molname, basename))
    SCFormIO::set_default_basename(molname);

  char * ckptfile = new char[strlen(molname)+6];
  sprintf(ckptfile,"%s.ckpt",molname);
  
  KeyValValueString restartfiledef(ckptfile);
  char * restartfile = keyval->pcharvalue("restart_file", restartfiledef);
  
  int restart = keyval->booleanvalue("restart",truevalue);

  int checkpoint = keyval->booleanvalue("checkpoint",truevalue);

  int savestate = keyval->booleanvalue("savestate",truevalue);

  struct stat sb;
  RefMolecularEnergy mole;
  RefOptimize opt;

  int statresult, statsize;
  if (restart) {
    if (grp->me() == 0) {
      statresult = stat(restartfile,&sb);
      statsize = (statresult==0) ? sb.st_size : 0;
    }
    grp->bcast(statresult);
    grp->bcast(statsize);
  }
  if (restart && statresult==0 && statsize) {
    BcastStateInBin si(grp,restartfile);
    if (keyval->exists("override")) {
      si.set_override(new PrefixKeyVal(keyval,"override"));
    }
    char *suf = strrchr(restartfile,'.');
    if (!strcmp(suf,".wfn")) {
      mole.key_restore_state(si,"mole");
      cout << node0 << endl << indent
           << "Restored <MolecularEnergy> from " << restartfile << endl;

      opt = keyval->describedclassvalue("opt");
      if (opt.nonnull())
        opt->set_function(mole);
    }
    else {
      opt.key_restore_state(si,"opt");
      if (opt.nonnull()) {
        mole = opt->function();
        cout << node0 << endl << indent
             << "Restored <Optimize> from " << restartfile << endl;
      }
    }
  } else {
    mole = keyval->describedclassvalue("mole");
    opt = keyval->describedclassvalue("opt");
  }

  if (mole.nonnull()) {
    MolecularFormula mf(mole->molecule());
    cout << node0 << endl << indent
         << "Molecular formula " << mf.formula() << endl;
  }

  if (checkpoint && opt.nonnull()) {
    opt->set_checkpoint();
    if (grp->me() == 0) opt->set_checkpoint_file(ckptfile);
    else opt->set_checkpoint_file(devnull);
  }

  // see if frequencies are wanted

  RefMolecularHessian molhess = keyval->describedclassvalue("hess");
  RefMolecularFrequencies molfreq = keyval->describedclassvalue("freq");
  
  int check = (options.retrieve("c") != 0);
  int limit = atoi(options.retrieve("l"));
  if (limit) {
    RefWavefunction wfn(mole);
    if (wfn.nonnull() && wfn->ao_dimension()->n() > limit) {
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
  
  if (tim.nonnull()) tim->change("calc");

  int do_energy = keyval->booleanvalue("do_energy",truevalue);
  
  int do_grad = keyval->booleanvalue("do_gradient",falsevalue);

  int do_opt = keyval->booleanvalue("optimize",truevalue);
  
  int do_pdb = keyval->booleanvalue("write_pdb",falsevalue);
  
  int print_mole = keyval->booleanvalue("print_mole",truevalue);
  
  int print_timings = keyval->booleanvalue("print_timings",truevalue);
  
  // sanity checks for the benefit of reasonable looking output
  if (opt.null())
    do_opt=0;
  else
    do_grad=1;
  
  cout << node0 << endl << indent
       << "MPQC options:" << endl << incindent
       << indent << "matrixkit     = <"
       << SCMatrixKit::default_matrixkit()->class_name() << ">" << endl
       << indent << "filename      = " << molname << endl
       << indent << "restart_file  = " << restartfile << endl
       << indent << "restart       = " << (restart ? "yes" : "no") << endl
       << indent << "checkpoint    = " << (checkpoint ? "yes" : "no") << endl
       << indent << "savestate     = " << (savestate ? "yes" : "no") << endl
       << indent << "do_energy     = " << (do_energy ? "yes" : "no") << endl
       << indent << "do_gradient   = " << (do_grad ? "yes" : "no") << endl
       << indent << "optimize      = " << (do_opt ? "yes" : "no") << endl
       << indent << "write_pdb     = " << (do_pdb ? "yes" : "no") << endl
       << indent << "print_mole    = " << (print_mole ? "yes" : "no") << endl
       << indent << "print_timings = " << (print_timings ? "yes" : "no")
       << endl << decindent;

  delete[] restartfile;
  delete[] ckptfile;
  
  int ready_for_freq = 1;
  if (mole.nonnull()) {
    if (((do_opt && opt.nonnull()) || do_grad)
        && !mole->gradient_implemented()) {
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
        cout << node0 << indent
             << scprintf("Value of the MolecularEnergy: %15.10f",
                         mole->energy())
             << endl << endl;
      } else {
        cout << node0 << indent
             << "The optimization has NOT converged." << endl << endl;
        ready_for_freq = 0;
      }
    } else if (do_grad && mole->gradient_implemented()) {
      mole->do_gradient(1);
      cout << node0 << endl << indent
           << scprintf("Value of the MolecularEnergy: %15.10f",
                       mole->energy())
           << endl;
      if (mole->value_result().actual_accuracy()
          > mole->value_result().desired_accuracy()) {
        cout << node0 << indent
             << "WARNING: desired accuracy not achieved in energy" << endl;
      }
      cout << node0 << endl;
      // Use result_noupdate since the energy might not have converged
      // to the desired accuracy in which case grabbing the result will
      // start up the calculation again.  However the gradient might
      // not have been computed (if we are restarting and the gradient
      // isn't in the save file for example).
      RefSCVector grad;
      if (mole->gradient_result().computed()) {
        grad = mole->gradient_result().result_noupdate();
      }
      else {
        grad = mole->gradient();
      }
      if (grad.nonnull()) {
        grad.print("Gradient of the MolecularEnergy:");
        if (mole->gradient_result().actual_accuracy()
            > mole->gradient_result().desired_accuracy()) {
          cout << node0 << indent
               << "WARNING: desired accuracy not achieved in gradient" << endl;
        }
      }
    } else if (do_energy && mole->value_implemented()) {
      cout << node0 << endl << indent
           << scprintf("Value of the MolecularEnergy: %15.10f",
                       mole->energy())
           << endl << endl;
    }
  }

  if (tim.nonnull()) tim->exit("calc");

  // save this before doing the frequency stuff since that obsoletes the
  // function stuff
  if (savestate) {
    if (opt.nonnull()) {
      if (grp->me() == 0) {
        ckptfile = new char[strlen(molname)+6];
        sprintf(ckptfile,"%s.ckpt",molname);
      }
      else {
        ckptfile = new char[strlen(devnull)+1];
        strcpy(ckptfile, devnull);
      }

      StateOutBin so(ckptfile);
      opt.save_state(so);
      so.close();

      delete[] ckptfile;
    }

    if (mole.nonnull()) {
      if (grp->me() == 0) {
        ckptfile = new char[strlen(molname)+6];
        sprintf(ckptfile,"%s.wfn",molname);
      }
      else {
        ckptfile = new char[strlen(devnull)+1];
        strcpy(ckptfile, devnull);
      }
  
      StateOutBin so(ckptfile);
      mole.save_state(so);
      so.close();

      delete[] ckptfile;
    }
  }

  // Frequency calculation.
  if (ready_for_freq && molfreq.nonnull()) {
    RefSymmSCMatrix xhessian;
    if (molhess.nonnull()) {
      // if "hess" input was given, use it to compute the hessian
      xhessian = molhess->cartesian_hessian();
    }
    else if (mole->hessian_implemented()) {
      // if mole can compute the hessian, use that hessian
      xhessian = mole->get_cartesian_hessian();
    }
    else if (mole->gradient_implemented()) {
      // if mole can compute gradients, use gradients at finite
      // displacements to compute the hessian
      molhess = new FinDispMolecularHessian(mole);
      xhessian = molhess->cartesian_hessian();
    }
    else {
      cout << "mpqc: WARNING: Frequencies cannot be computed" << endl;
    }

    if (xhessian.nonnull()) {
      char *hessfile = SCFormIO::fileext_to_filename(".hess");
      MolecularHessian::write_cartesian_hessian(hessfile,
                                                mole->molecule(), xhessian);
      delete[] hessfile;

      molfreq->compute_frequencies(xhessian);
      // DEGENERACY IS NOT CORRECT FOR NON-SINGLET CASES:
      molfreq->thermochemistry(1);
    }
  }

  // see if any pictures are desired
  RefRender renderer = keyval->describedclassvalue("renderer");
  RefRenderedObject rendered = keyval->describedclassvalue("rendered");
  RefAnimatedObject animated = keyval->describedclassvalue("rendered");
  if (renderer.nonnull() && rendered.nonnull()) {
    if (tim.nonnull()) tim->enter("render");
    if (grp->me() == 0) renderer->render(rendered);
    if (tim.nonnull()) tim->exit("render");
  }
  else if (renderer.nonnull() && animated.nonnull()) {
    if (tim.nonnull()) tim->enter("render");
    if (grp->me() == 0) renderer->animate(animated);
    if (tim.nonnull()) tim->exit("render");
  }
  else if (renderer.nonnull()) {
    if (tim.nonnull()) tim->enter("render");
    int n = keyval->count("rendered");
    for (i=0; i<n; i++) {
      rendered = keyval->describedclassvalue("rendered",i);
      animated = keyval->describedclassvalue("rendered",i);
      if (rendered.nonnull()) {
        // make sure the object has a name so we don't overwrite its file
        if (rendered->name() == 0) {
          char ic[64];
          sprintf(ic,"%02d",i);
          rendered->set_name(ic);
        }
        if (grp->me() == 0) renderer->render(rendered);
      }
      else if (animated.nonnull()) {
        // make sure the object has a name so we don't overwrite its file
        if (animated->name() == 0) {
          char ic[64];
          sprintf(ic,"%02d",i);
          animated->set_name(ic);
        }
        if (grp->me() == 0) renderer->animate(animated);
      }
    }
    if (tim.nonnull()) tim->exit("render");
  }
  RefMolFreqAnimate molfreqanim = keyval->describedclassvalue("animate_modes");
  if (ready_for_freq && molfreq.nonnull()
      && molfreqanim.nonnull() && renderer.nonnull()) {
    if (tim.nonnull()) tim->enter("render");
    molfreq->animate(renderer, molfreqanim);
    if (tim.nonnull()) tim->exit("render");
  }

  if (mole.nonnull()) {
    if (print_mole)
      mole->print(cout);

    if (do_pdb && grp->me() == 0) {
      ckptfile = new char[strlen(molname)+5];
      sprintf(ckptfile, "%s.pdb", molname);
      ofstream pdbfile(ckptfile);
      mole->molecule()->print_pdb(pdbfile);
      delete[] ckptfile;
    }
    
  }
  else {
    cout << node0 << "mpqc: The molecular energy object is null" << endl
         << " make sure \"mole\" specifies a MolecularEnergy derivative"
         << endl;
  }

  if (parsedkv->have_unseen()) {
    cout << node0 << indent
         << "The following keywords in \"" << input << "\" were ignored:"
         << endl;
    cout << incindent;
    parsedkv->print_unseen(cout<<node0);
    cout << decindent;
    cout << node0 << endl;
  }

  if (print_timings)
    if (tim.nonnull()) tim->print(cout);

  delete[] basename;
  delete[] molname;
  SCFormIO::set_default_basename(0);

  molfreqanim = 0;
  animated = 0;
  rendered = 0;
  renderer = 0;
  molfreq = 0;
  molhess = 0;
  opt = 0;
  mole = 0;
  debugger = 0;
  thread = 0;
  tim = 0;
  keyval = 0;
  parsedkv = 0;
  grp = 0;
  clean_up();

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
