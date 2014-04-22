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

// This is needed to make GNU extensions available, such as
// feenableexcept and fedisableexcept.
#ifndef _GNU_SOURCE
#  define _GNU_SOURCE
#endif

#ifdef HAVE_CONFIG_H
#include <mpqc_config.h>
#endif

#ifdef HAVE_TIME_H
#include <time.h>
#endif

#include <new>
#include <stdexcept>
#include <cassert>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fstream>
#include <libgen.h>

#include <mpqc_config.h>
#include <sstream>

#include "mpqcinit.h"

#include <util/options/GetLongOpt.h>
#include <util/misc/scexception.h>
#include <util/misc/newstring.h>
#include <util/keyval/keyval.h>
#include <util/state/state_bin.h>
#include <util/group/message.h>
#include <util/group/memory.h>
#include <util/group/mstate.h>
#include <util/group/thread.h>
#include <util/group/pregtime.h>
#include <util/misc/bug.h>
#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <util/misc/runnable.h>
#include <util/misc/consumableresources.h>
#include <util/render/render.h>

#include <math/optimize/opt.h>

#include <chemistry/molecule/coor.h>
#include <chemistry/molecule/energy.h>
#include <chemistry/molecule/molfreq.h>
#include <chemistry/molecule/findisp.h>
#include <chemistry/molecule/formula.h>
#include <chemistry/qc/wfn/wfn.h>

#include <util/group/linkage.h>
#include <util/state/linkage.h>
#include <chemistry/qc/basis/linkage.h>
#include <chemistry/qc/wfn/linkage.h>
#include <chemistry/qc/scf/linkage.h>
#include <chemistry/qc/dft/linkage.h>
#include <chemistry/qc/etrain/linkage.h>
#include <chemistry/qc/mbpt/linkage.h>
#include <chemistry/qc/lmp2/linkage.h>

#ifdef HAVE_LIBINT2
# include <chemistry/qc/libint2/linkage.h>
#endif // HAVE_LIBINT2

#ifdef MPQC_R12
#include <chemistry/qc/mbptr12/linkage.h>
#include <chemistry/qc/ccr12/linkage.h>
#endif // MPQC_R12

#ifdef HAVE_PSI3
# include <chemistry/qc/psi/linkage.h>
#endif // HAVE_PSI3

#ifdef MPQC_CI
# include <chemistry/qc/ci/linkage.h>
#endif

#ifdef HAVE_MPI
#define OMPI_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#include <mpi.h>
#include <util/group/messmpi.h>
#endif

using namespace std;
using namespace sc;

#include "mpqcin.h"

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

#include <signal.h>

#ifdef HAVE_FENV_H
#  include <fenv.h>
#endif

static void
print_unseen(const Ref<ParsedKeyVal> &parsedkv,
             const char *input)
{
  if (parsedkv->have_unseen()) {
    ExEnv::out0() << endl;
    ExEnv::out0() << indent
         << "The following keywords in \"" << input << "\" were ignored:"
         << endl;
    ExEnv::out0() << incindent;
    parsedkv->print_unseen(ExEnv::out0());
    ExEnv::out0() << decindent;
  }
}

int
try_main(int argc, char *argv[])
{
  KeyValValueboolean truevalue(1), falsevalue(0);
  int i;
  const char *devnull = "/dev/null";

  // parse commandline options
  GetLongOpt options;

  options.usage("[options] [filename]");
  options.enroll("f", GetLongOpt::MandatoryValue,
                 "the name of an object format input file", 0);
  options.enroll("o", GetLongOpt::MandatoryValue,
                 "the name of the output file", 0);
  options.enroll("l", GetLongOpt::MandatoryValue, "basis set limit", "0");
  options.enroll("W", GetLongOpt::MandatoryValue,
                 "set the working directory", ".");
  options.enroll("c", GetLongOpt::NoValue, "check input then exit", 0);
  options.enroll("v", GetLongOpt::NoValue, "print the version number", 0);
  options.enroll("w", GetLongOpt::NoValue, "print the warranty", 0);
  options.enroll("L", GetLongOpt::NoValue, "print the license", 0);
  options.enroll("k", GetLongOpt::NoValue, "print key/value assignments", 0);
  options.enroll("i", GetLongOpt::NoValue, "convert simple to OO input", 0);
  options.enroll("d", GetLongOpt::NoValue, "debug", 0);
  options.enroll("h", GetLongOpt::NoValue, "print this message", 0);

  MPQCInit init(options,argc,argv);

  int optind = options.parse(argc, argv);

  const char *output = options.retrieve("o");
  ostream *outstream = 0;
  if (output != 0) {
    outstream = new ofstream(output);
    ExEnv::set_out(outstream);
  }

  if (options.retrieve("h")) {
    ExEnv::out0()
         << indent << "MPQC version " << MPQC_VERSION << endl
         << indent << "compiled for " << TARGET_ARCH << endl
         << SCFormIO::copyright << endl;
    options.usage(ExEnv::out0());
    exit(0);
  }

  if (options.retrieve("v")) {
    ExEnv::out0()
         << indent << "MPQC version " << MPQC_VERSION << endl
         << indent << "compiled for " << TARGET_ARCH << endl
         << SCFormIO::copyright;
    exit(0);
  }

  if (options.retrieve("w")) {
    ExEnv::out0()
         << indent << "MPQC version " << MPQC_VERSION << endl
         << indent << "compiled for " << TARGET_ARCH << endl
         << SCFormIO::copyright << endl
         << SCFormIO::warranty;
    exit(0);
  }

  if (options.retrieve("L")) {
    ExEnv::out0()
         << indent << "MPQC version " << MPQC_VERSION << endl
         << indent << "compiled for " << TARGET_ARCH << endl
         << SCFormIO::copyright << endl
         << SCFormIO::license;
    exit(0);
  }

  // set the working dir
  if (strcmp(options.retrieve("W"),".")) {
      int err = chdir(options.retrieve("W"));
      MPQC_ASSERT(!err);
  }

  // initialize keyval input
  const char *object_input = options.retrieve("f");
  const char *generic_input;
  if (argc - optind == 0) {
    generic_input = 0;
  }
  else if (argc - optind == 1) {
    generic_input = argv[optind];
  }
  else {
    options.usage();
    throw invalid_argument("extra arguments given");
  }

  ExEnv::init(argc, argv);
  init.init_fp();
  init.init_limits();
  Ref<MessageGrp> grp = init.init_messagegrp();
  init.init_io(grp);

  if (object_input == 0 && generic_input == 0) {
    generic_input = "mpqc.in";
    }
  else if (object_input && generic_input) {
    options.usage();
    throw invalid_argument("only one of -f and a file argument can be given");
    }

  const char *input;
  if (object_input) input = object_input;
  if (generic_input) input = generic_input;

  Ref<ParsedKeyVal> parsedkv;
  // read the input file on only node 0
  char *in_char_array;
  if (grp->me() == 0) {
    ifstream is(input);
    ostringstream ostrs;
    is >> ostrs.rdbuf();
    int n = 1 + strlen(ostrs.str().c_str());
    in_char_array = strcpy(new char[n],ostrs.str().c_str());
    grp->bcast(n);
    grp->bcast(in_char_array, n);
    }
  else {
    int n;
    grp->bcast(n);
    in_char_array = new char[n];
    grp->bcast(in_char_array, n);
    }

  int use_simple_input;
  if (generic_input && grp->me() == 0) {
    MPQCIn mpqcin;
    use_simple_input = mpqcin.check_string(in_char_array);
  }
  else {
    use_simple_input = 0;
  }
  grp->bcast(use_simple_input);

  if (use_simple_input) {
    MPQCIn mpqcin;
    char *simple_input_text = mpqcin.parse_string(in_char_array);
    if (options.retrieve("i")) {
      ExEnv::out0() << "Generated object-oriented input file:" << endl
                   << simple_input_text
                   << endl;
      exit(0);
    }
    parsedkv = new ParsedKeyVal();
    parsedkv->parse_string(simple_input_text);
    delete[] simple_input_text;
    }
  else {
    parsedkv = new ParsedKeyVal();
    parsedkv->parse_string(in_char_array);
    }
  delete[] in_char_array;

  if (options.retrieve("k")) parsedkv->verbose(1);
  Ref<KeyVal> keyval = new PrefixKeyVal(parsedkv.pointer(),"mpqc");

  init.init_basename(input, (output?std::string(output):std::string("")));

  if (options.retrieve("d"))
    SCFormIO::set_debug(1);

  // initialize timing for mpqc
  init.init_timer(grp,keyval);

  Timer timer;

  timer.enter("input");

  // announce ourselves
  const char title1[] = "MPQC: Massively Parallel Quantum Chemistry";
  int ntitle1 = sizeof(title1);
  const char title2[] = "Version " MPQC_VERSION;
  int ntitle2 = sizeof(title2);
  ExEnv::out0() << endl;
  ExEnv::out0() << indent;
  for (i=0; i<(80-ntitle1)/2; i++) ExEnv::out0() << ' ';
  ExEnv::out0() << title1 << endl;
  ExEnv::out0() << indent;
  for (i=0; i<(80-ntitle2)/2; i++) ExEnv::out0() << ' ';
  ExEnv::out0() << title2 << endl << endl;

  const char *tstr = 0;
#if defined(HAVE_TIME) && defined(HAVE_CTIME)
  time_t t;
  time(&t);
  tstr = ctime(&t);
#endif
  if (!tstr) {
    tstr = "UNKNOWN";
  }

  ExEnv::out0()
       << indent << scprintf("Machine:    %s", TARGET_ARCH) << endl
       << indent << scprintf("User:       %s@%s",
                             ExEnv::username(), ExEnv::hostname()) << endl
       << indent << scprintf("Start Time: %s", tstr) << endl;

  Ref<ThreadGrp> thread = init.init_threadgrp(keyval);
  Ref<MemoryGrp> memory = init.init_memorygrp(keyval);

  ExEnv::out0() << indent
       << "Using " << grp->class_name()
       << " for message passing (number of nodes = " << grp->n() << ")." << endl
       << indent
       << "Using " << thread->class_name()
       << " for threading (number of threads = " << thread->nthread() << ")." << endl
       << indent
       << "Using " << memory->class_name()
       << " for distributed shared memory." << endl
       << indent
       << "Total number of processors = " << grp->n() * thread->nthread() << endl;

  // now set up the debugger
  Ref<Debugger> debugger; debugger << keyval->describedclassvalue("debug");
  if (debugger) {
    Debugger::set_default_debugger(debugger);
    debugger->set_exec(argv[0]);
    debugger->set_prefix(grp->me());
    if (options.retrieve("d"))
      debugger->debug("Starting debugger because -d given on command line.");
  }

  // now check to see what matrix kit to use
  if (keyval->exists("matrixkit"))
    SCMatrixKit::set_default_matrixkit(
      dynamic_cast<SCMatrixKit*>(
        keyval->describedclassvalue("matrixkit").pointer()));

  init.init_integrals(keyval);
  Ref<Integral> integral = Integral::get_default_integral();
  ExEnv::out0() << endl << indent
       << "Using " << integral->class_name()
       << " by default for molecular integrals evaluation" << endl;

  init.init_resources(keyval);
  Ref<ConsumableResources> resources = ConsumableResources::get_default_instance();
  ExEnv::out0() << indent
       << "Given resources: " << resources->sprint() << endl << endl;

  // check for a molecular energy and optimizer
  std::string basename = SCFormIO::default_basename();
  KeyValValuestring molnamedef(basename);
  std::string molname = keyval->stringvalue("filename", molnamedef);
  if (molname != basename)
    SCFormIO::set_default_basename(molname.c_str());

  std::string ckptfile = molname;
  ckptfile += ".ckpt";

  std::string restartfile = keyval->stringvalue("restart_file",
                                                KeyValValuestring(ckptfile));

  std::string wfn_file = keyval->stringvalue("wfn_file");
  if (wfn_file.empty()) {
    wfn_file = molname;
    wfn_file += ".wfn";
  }
  std::string mole_ckpt_file;
  mole_ckpt_file = wfn_file;

  int restart = keyval->booleanvalue("restart",truevalue);

  int checkpoint = keyval->booleanvalue("checkpoint",truevalue);
  int checkpoint_freq = keyval->intvalue("checkpoint_freq",KeyValValueint(1));

  int savestate = keyval->booleanvalue("savestate",truevalue);

  struct stat sb;
  Ref<MolecularEnergy> mole;
  Ref<Optimize> opt;

  int statresult, statsize;
  if (restart) {
    if (grp->me() == 0) {
      statresult = stat(restartfile.c_str(),&sb);
      statsize = (statresult==0) ? sb.st_size : 0;
    }
    grp->bcast(statresult);
    grp->bcast(statsize);
  }
  if (restart && statresult==0 && statsize) {
    BcastStateInBin si(grp,restartfile.c_str());
    if (keyval->exists("override")) {
      si.set_override(new PrefixKeyVal(keyval,"override"));
    }
    const char *suf = strrchr(restartfile.c_str(),'.');
    if (!strcmp(suf,".wfn")) {
      mole << SavableState::key_restore_state(si,"mole");
      ExEnv::out0() << endl
                   << indent << "Restored <" << mole->class_name()
                   << "> from " << restartfile << endl;

      opt << keyval->describedclassvalue("opt");
      if (opt)
        opt->set_function(mole.pointer());
    }
    else {
      opt << SavableState::key_restore_state(si,"opt");
      if (opt) {
        mole << opt->function();
        ExEnv::out0() << endl << indent
             << "Restored <Optimize> from " << restartfile << endl;
      }
    }
  } else {
    mole << keyval->describedclassvalue("mole");
    opt << keyval->describedclassvalue("opt");
  }

  if (mole) {
    MolecularFormula mf(mole->molecule());
    ExEnv::out0() << endl << indent
         << "Molecular formula " << mf.formula() << endl;
  }

  int check = (options.retrieve("c") != 0);
  int limit = atoi(options.retrieve("l"));
  if (limit) {
    Ref<Wavefunction> wfn; wfn << mole;
    if (wfn && wfn->ao_dimension()->n() > limit) {
      ExEnv::out0() << endl << indent
           << "The limit of " << limit << " basis functions has been exceeded."
           << endl;
      check = 1;
    }
  }

  if (check) {
    ExEnv::out0() << endl << indent
         << "Exiting since the check option is on." << endl;
    exit(0);
  }

  timer.change("calc");

  const int do_energy = keyval->booleanvalue("do_energy",truevalue);

  const int do_grad = keyval->booleanvalue("do_gradient",falsevalue);

  const int do_freq = keyval->booleanvalue("do_freq",falsevalue);

  const int do_xyz = keyval->booleanvalue("write_xyz",truevalue);

  const int print_mole = keyval->booleanvalue("print_mole",truevalue);

  const int print_resources = keyval->booleanvalue("print_resources",truevalue);

  const int print_timings = keyval->booleanvalue("print_timings",truevalue);

  // default value for optimize is true if opt is given, and false if it is not
  const int do_opt = keyval->booleanvalue("optimize", (!opt ? falsevalue : truevalue));
  if (do_opt && !opt && mole) { // if user requested optimization but no Optimize object given pick a default object
    Ref<AssignedKeyVal> akv = new AssignedKeyVal;
    akv->assign("function", mole.pointer());
    Ref<KeyVal> kv = akv;
    opt = new QNewtonOpt(kv);
  }

  // Read in all of the runnable objects now, so we can get rid of
  // the reference to the input file.
  const int nrun = keyval->count("run");
  std::vector<Ref<Runnable> > run(nrun);
  for (int i=0; i<nrun; i++) {
    run[i] << keyval->describedclassvalue("run",i);
    if (run[i]) {
      ExEnv::out0() << indent
                    << "Read run:" << i
                    << " of type " << run[i]->class_name()
                    << "."
                    << std::endl;
    }
    else {
      ExEnv::out0() << indent
                    << "Problem reading run:" << i
                    << ", it will be ignored."
                    << std::endl;
    }
  }

  // see if the user specified how to compute gradients, hessians, and frequencies
  Ref<MolecularGradient> molgrad;
  molgrad << keyval->describedclassvalue("grad");
  Ref<MolecularHessian> molhess;
  molhess << keyval->describedclassvalue("hess");
  Ref<MolecularFrequencies> molfreq;
  molfreq << keyval->describedclassvalue("freq");

  // see if any pictures are desired
  Ref<Render> renderer;
  renderer << keyval->describedclassvalue("renderer");

  // If we have a renderer, then we will read in some more info
  // below.  Otherwise we can get rid of the keyval's, to eliminate
  // superfluous references to objects that we might otherwise be
  // able to delete.  We cannot read in the remaining rendering
  // objects now, since some of their KeyVal CTOR's are heavyweight,
  // requiring optimized geometries, etc.
  if (!renderer) {
    if (parsedkv) print_unseen(parsedkv, input);
    keyval = 0;
    parsedkv = 0;
  }

  ExEnv::out0() << endl << indent
       << "MPQC options:" << endl << incindent
       << indent << "matrixkit       = <"
       << SCMatrixKit::default_matrixkit()->class_name() << ">" << endl
       << indent << "filename        = " << molname << endl
       << indent << "restart_file    = " << restartfile << endl
       << indent << "restart         = " << (restart ? "yes" : "no") << endl
       << indent << "checkpoint      = " << (checkpoint ? "yes" : "no") << endl
       << indent << "savestate       = " << (savestate ? "yes" : "no") << endl
       << indent << "do_energy       = " << (do_energy ? "yes" : "no") << endl
       << indent << "do_gradient     = " << (do_grad ? "yes" : "no") << endl
       << indent << "do_freq         = " << (do_freq ? "yes" : "no") << endl
       << indent << "optimize        = " << (do_opt ? "yes" : "no") << endl
       << indent << "write_xyz       = " << (do_xyz ? "yes" : "no") << endl
       << indent << "print_mole      = " << (print_mole ? "yes" : "no") << endl
       << indent << "print_timings   = " << (print_timings ? "yes" : "no") << endl
       << indent << "print_resources = " << (print_resources ? "yes" : "no")
       << endl << decindent;

  // Prepare to checkpoint the computation
  if (checkpoint && mole) {
    mole->set_checkpoint();
    if (grp->me() == 0) mole->set_checkpoint_file(mole_ckpt_file.c_str());
    else mole->set_checkpoint_file(devnull);
    mole->set_checkpoint_freq(checkpoint_freq);
  }
  if (checkpoint && opt) {
    opt->set_checkpoint();
    if (grp->me() == 0) opt->set_checkpoint_file(ckptfile.c_str());
    else opt->set_checkpoint_file(devnull);
  }

  int ready_for_freq = 1;
  if (mole) {

    const bool need_gradient = ((do_opt && opt) || do_grad);
    bool have_gradient = mole->gradient_implemented() || molgrad;
    if (need_gradient && !have_gradient) {
      ExEnv::out0() << indent
           << "WARNING: optimization or gradient requested but the given"
           << endl
           << "         MolecularEnergy object cannot do gradients."
           << endl
           << "         Gradient will be computed numerically by finite differences."
           << endl;
      molgrad = new FinDispMolecularGradient(mole);
      have_gradient = true;
    }

    if (do_opt && opt && have_gradient) { // optimize geometry

      const bool use_mpqc_molgrad = !mole->gradient_implemented() && molgrad;
      if (use_mpqc_molgrad) {
        mole->set_molgrad(molgrad);
      }

      ExEnv::out0() << std::endl
                    << indent
                    << "Optimizing using:"
                    << std::endl
                    << incindent;
      opt->print();
      ExEnv::out0() << decindent;
      int result = opt->optimize();
      if (result) {
        ExEnv::out0() << indent
             << "The optimization has converged." << endl << endl;
        ExEnv::out0() << indent
             << scprintf("Value of the MolecularEnergy: %15.10f",
                         mole->energy())
             << endl << endl;
      } else {
        ExEnv::out0() << indent
             << "The optimization has NOT converged." << endl << endl;
        ready_for_freq = 0;
      }

      if (use_mpqc_molgrad) {
        mole->set_molgrad(0);
      }

      // EFV Feb. 6 2011
      // good idea to print optimized coordinates here
      {
        Ref<MolecularCoor> mc = mole->molecularcoor();
        if (mc)
          mc->print(ExEnv::out0());
        else {
          mole->molecule()->print(ExEnv::out0());
          // generate simple internals and print them out just FYI
          Ref<SetIntCoor> cset = new SetIntCoor;
          Ref<IntCoorGen> intcoorgen = new IntCoorGen(mole->molecule());
          intcoorgen->generate(cset);
          cset->update_values(mole->molecule());
          cset->print_details(mole->molecule(), ExEnv::out0());
        }
      }

    } else if (do_grad && have_gradient) { // only compute the gradient

      RefSCVector grad;

      if (mole->gradient_implemented()) {
        mole->do_gradient(1);
        ExEnv::out0() << endl << indent
            << scprintf("Value of the MolecularEnergy: %15.10f", mole->energy())
            << endl;
        if (mole->value_result().actual_accuracy()
            > mole->value_result().desired_accuracy()) {
          ExEnv::out0() << indent
              << "WARNING: desired accuracy not achieved in energy" << endl;
        }
        if (mole->gradient_result().actual_accuracy()
            > mole->gradient_result().desired_accuracy()) {
          ExEnv::out0() << indent
               << "WARNING: desired accuracy not achieved in gradient" << endl;
        }
        ExEnv::out0() << endl;
        // Use result_noupdate since the energy might not have converged
        // to the desired accuracy in which case grabbing the result will
        // start up the calculation again.  However the gradient might
        // not have been computed (if we are restarting and the gradient
        // isn't in the save file for example).
        if (mole->gradient_result().computed()) {
          grad = mole->gradient_result().result_noupdate();
        } else {
          grad = mole->gradient();
        }
      }
      else { // use molgrad
        MPQC_ASSERT(molgrad);
        mole->set_molgrad(molgrad);
        grad = mole->gradient();
        mole->set_molgrad(0);
      }

      if (grad) {
        grad.print("Gradient of the MolecularEnergy:");
      }
    }
    else if (do_energy && mole->value_implemented()) { // only compute the energy
      ExEnv::out0() << endl << indent
           << scprintf("Value of the MolecularEnergy: %15.10f",
                       mole->energy())
           << endl << endl;
    }

  }

  timer.exit("calc");

  // save this before doing the frequency stuff since that obsoletes the
  // function stuff
  if (savestate) {
    if (opt) {
      if (grp->me() == 0) {
        ckptfile = molname;
        ckptfile += ".ckpt";
      }
      else {
        ckptfile = devnull;
      }

      StateOutBin so(ckptfile.c_str());
      SavableState::save_state(opt.pointer(),so);
      so.close();
    }

    if (mole) {
      if (grp->me() != 0) {
        wfn_file = devnull;
      }

      StateOutBin so(wfn_file.c_str());
      SavableState::save_state(mole.pointer(),so);
      so.close();

    }
  }

  // Frequency calculation
  if (do_freq) {

    if ((opt && ready_for_freq) || !opt) {
      RefSymmSCMatrix xhessian;
      if (mole->hessian_implemented()) { // if mole can compute the hessian, use that hessian
        xhessian = mole->get_cartesian_hessian();
      }
      else if (molhess) { // else use user-provided hessian
        const bool molhess_needs_mole = (molhess->energy() == 0);
        if (molhess_needs_mole) molhess->set_energy(mole);
        xhessian = molhess->cartesian_hessian();
        if (molhess_needs_mole) molhess->set_energy(0);
      }
      else if (mole->gradient_implemented() ||
               mole->value_implemented()) {
        // if mole can compute energies or gradients, use finite
        // differences to compute the hessian
        Ref<FinDispMolecularHessian> fdmolhess = new FinDispMolecularHessian(mole);
        fdmolhess->params()->set_restart(restart);
        fdmolhess->params()->set_checkpoint(checkpoint);
        molhess = fdmolhess;
        xhessian = molhess->cartesian_hessian();
        molhess = 0;
      }
      else {
        throw FeatureNotImplemented("mpqc: frequencies cannot be computed for this type of MolecularEnergy",
                                    __FILE__, __LINE__);
      }

      if (xhessian) {
        char *hessfile = SCFormIO::fileext_to_filename(".hess");
        MolecularHessian::write_cartesian_hessian(hessfile,
                                                  mole->molecule(), xhessian);
        delete[] hessfile;

        if (!molfreq)
          molfreq = new MolecularFrequencies(mole->molecule());
        molfreq->compute_frequencies(xhessian);
        // DEGENERACY IS NOT CORRECT FOR NON-SINGLET CASES:
        molfreq->thermochemistry(1);
      }
    }
  }

  if (renderer) {
    Ref<RenderedObject> rendered;
    rendered << keyval->describedclassvalue("rendered");
    Ref<AnimatedObject> animated;
    animated << keyval->describedclassvalue("rendered");
    if (rendered) {
      timer.enter("render");
      if (grp->me() == 0) renderer->render(rendered);
      timer.exit("render");
    }
    else if (animated) {
      timer.enter("render");
      if (grp->me() == 0) renderer->animate(animated);
      timer.exit("render");
    }
    else {
      timer.enter("render");
      int n = keyval->count("rendered");
      for (i=0; i<n; i++) {
        rendered << keyval->describedclassvalue("rendered",i);
        animated << keyval->describedclassvalue("rendered",i);
        if (rendered) {
          // make sure the object has a name so we don't overwrite its file
          if (rendered->name() == 0) {
            char ic[64];
            sprintf(ic,"%02d",i);
            rendered->set_name(ic);
          }
          if (grp->me() == 0) renderer->render(rendered);
        }
        else if (animated) {
          // make sure the object has a name so we don't overwrite its file
          if (animated->name() == 0) {
            char ic[64];
            sprintf(ic,"%02d",i);
            animated->set_name(ic);
          }
          if (grp->me() == 0) renderer->animate(animated);
        }
      }
      timer.exit("render");
    }
    Ref<MolFreqAnimate> molfreqanim;
    molfreqanim << keyval->describedclassvalue("animate_modes");
    if (molfreq && molfreqanim) {
      timer.enter("render");
      molfreq->animate(renderer, molfreqanim);
      timer.exit("render");
    }
  }

  for (int i=0; i<nrun; i++) {
    if (run[i]) {
      ExEnv::out0() << indent << "Running object run:" << i << std::endl;
      ExEnv::out0() << incindent;
      run[i]->run();
      ExEnv::out0() << decindent;
    }
  }
  if (nrun) {
    ExEnv::out0() << std::endl;
    run.resize(0);
  }

  if (mole) {
    if (print_mole)
      mole->print(ExEnv::out0());

    if (do_xyz && grp->me() == 0) {
      ckptfile = molname;
      ckptfile += ".xyz";
      ofstream file(ckptfile.c_str());
      mole->molecule()->print_xyz(file);
    }

  }
  else {
    ExEnv::out0() << "mpqc: The molecular energy object is null" << endl
         << " make sure \"mole\" specifies a MolecularEnergy derivative"
         << endl;
  }

  if (parsedkv) print_unseen(parsedkv, input);

  if (print_resources) {
    const bool print_resources_state = false;
    const bool print_resources_stats = true;
    ExEnv::out0() << endl;
    resources->print_summary(ExEnv::out0(),print_resources_state,print_resources_stats);
  }
  if (print_timings)
    timer.print(ExEnv::out0());

  renderer = 0;
  molfreq = 0;
  molhess = 0;
  molgrad = 0;
  opt = 0;
  mole = 0;
  integral = 0;
  debugger = 0;
  thread = 0;
  keyval = 0;
  parsedkv = 0;
  grp = 0;
  memory = 0;
  resources = 0;

#if defined(HAVE_TIME) && defined(HAVE_CTIME)
  time(&t);
  tstr = ctime(&t);
#endif
  if (!tstr) {
    tstr = "UNKNOWN";
  }
  ExEnv::out0() << endl
               << indent << scprintf("End Time: %s", tstr) << endl;

  if (output != 0) {
    ExEnv::set_out(&cout);
    delete outstream;
  }

  return 0;
}

int
main(int argc, char *argv[])
{
  try {
  try_main(argc, argv);
  }
  catch (SCException &e) {
    cout << argv[0] << ": ERROR: SC EXCEPTION RAISED:" << endl
         << e.what()
         << endl;
    throw;
  }
  catch (bad_alloc &e) {
    cout << argv[0] << ": ERROR: MEMORY ALLOCATION FAILED:" << endl
         << e.what()
         << endl;
    throw;
  }
  catch (exception &e) {
    cout << argv[0] << ": ERROR: EXCEPTION RAISED:" << endl
         << e.what()
         << endl;
    throw;
  }
  catch (...) {
    cout << argv[0] << ": ERROR: UNKNOWN EXCEPTION RAISED" << endl;
    throw;
  }
  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
