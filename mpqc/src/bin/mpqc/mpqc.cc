
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <new.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

#include <util/options/GetLongOpt.h>
#include <util/misc/newstring.h>
#include <util/keyval/keyval.h>
#include <util/group/message.h>
#include <util/group/pregtime.h>
#include <util/misc/bug.h>
#include <util/misc/formio.h>

#include <math/optimize/opt.h>

#include <chemistry/molecule/coor.h>
#include <chemistry/molecule/energy.h>
#include <chemistry/molecule/molfreq.h>
#include <chemistry/qc/wfn/wfn.h>

// Force linkages:
#include <chemistry/qc/mbpt/linkage.h>

#include "disclaimer.h"

//////////////////////////////////////////////////////////////////////////

static void
clean_up(void)
{
  MessageGrp::set_default_messagegrp(0);
}

int
main(int argc, char *argv[])
{
  atexit(clean_up);
  
  // first thing, set up output classes
  SCFormIO::setindent(cout, 2);
  SCFormIO::setindent(cerr, 2);

  // parse commandline options
  GetLongOpt options;

  options.enroll("f", GetLongOpt::MandatoryValue,
                 "the name of the input file", "mpqc.in");
  options.enroll("messagegrp", GetLongOpt::MandatoryValue,
                 "which message group to use", 0);
  options.enroll("l", GetLongOpt::MandatoryValue, "basis set limit", "0");
  options.enroll("c", GetLongOpt::NoValue, "check input then exit", 0);
  options.enroll("v", GetLongOpt::NoValue, "print the version number", 0);
  options.enroll("w", GetLongOpt::NoValue, "print the warranty", 0);
  options.enroll("d", GetLongOpt::NoValue, "debug", 0);
  options.enroll("h", GetLongOpt::NoValue, "print this message", 0);

  options.parse(argc, argv);

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

  // open keyval input
  const char *input = options.retrieve("f");
  RefKeyVal pkv(new ParsedKeyVal(input));
  RefKeyVal ppkv(new AggregateKeyVal(new PrefixKeyVal(":mpqc",pkv),
                                     new PrefixKeyVal(":default",pkv)));
  pkv = new ParsedKeyVal("input",ppkv);
  RefKeyVal keyval = new AggregateKeyVal(ppkv,pkv);

  pkv = ppkv = 0;

  // get the message group.  the commandline take precedence, then what is
  // in the input file.
  //RefMessageGrp grp = MessageGrp::initial_messagegrp(argc, argv);
  RefMessageGrp grp;
  
  // if we are on a paragon then use a ParagonMessageGrp
  // otherwise read the message group from the input file
  if (grp.null()) {
#ifdef HAVE_NX_H
    grp = new ParagonMessageGrp;
#else
    grp = keyval->describedclassvalue("message");
#endif
  }

  if (grp.nonnull())
    MessageGrp::set_default_messagegrp(grp);
  else
    grp = MessageGrp::get_default_messagegrp();

  SCFormIO::set_printnode(0);
  SCFormIO::set_messagegrp(grp);
  if (options.retrieve("d"))
    SCFormIO::set_debug(1);

  // initialize timing for mpqc
  RefRegionTimer tim = new ParallelRegionTimer(grp,"mpqc",1,0);
  RegionTimer::set_default_regiontimer(tim);

  tim->enter("input");
  
  // now set up the debugger
  RefDebugger debugger = keyval->describedclassvalue(":debug");
  if (debugger.nonnull()) {
    debugger->set_exec(argv[0]);
    debugger->set_prefix(grp->me());
    if (options.retrieve("d"))
      debugger->debug("curt is a hog");
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
    molname = new_string("mpqc");
  
  char * ckptfile = new char[strlen(molname)+6];
  sprintf(ckptfile,"%s.ckpt",molname);
  
  int restart = 1;
  if (keyval->exists("restart"))
    restart = keyval->booleanvalue("restart");
  
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
    if (opt.nonnull()) {
      opt->set_checkpoint();
      opt->set_checkpoint_file(ckptfile);
    }
  }

  delete[] ckptfile;

  int check = (int) options.retrieve("c");
  int limit = atoi(options.retrieve("l"));
  if (limit) {
    Wavefunction *wfn = Wavefunction::castdown(mole);
    if (wfn && wfn->basis_dimension()->n() > limit) {
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
    do_grad=1;

  int do_opt = keyval->booleanvalue("optimize");
  if (keyval->error() != KeyVal::OK)
    do_opt=1;
  
  if (mole.nonnull()) {
    if (do_opt && opt.nonnull() && mole->gradient_implemented()) {
      opt->optimize();
    } else if (do_grad && mole->gradient_implemented()) {
      mole->gradient().print("gradient");
    } else if (do_energy && mole->value_implemented()) {
      cout << node0 << indent
           << scprintf("value of mole is %20.15f\n\n", mole->energy());
    }
  }

  tim->exit("calc");

  if (molfreq.nonnull()) {
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
      RefSCVector gradv = mole->gradient();
      molfreq->set_gradient(i, gradv);
    }

    molfreq->compute_frequencies_from_gradients();
    //molfreq->thermochemistry(scf_info.nopen+1);
    molfreq->thermochemistry(1);
    tim->exit("frequencies");
  }

  if (mole.nonnull()) {
    mole->print(cout);
  }
  else {
    cout << node0 << "mpqc: The molecular energy object is null" << endl
         << " make sure \"mole\" specificies a MolecularEnergy derivative"
         << endl;
  }

  ckptfile = new char[strlen(molname+5)];
  sprintf(ckptfile, "%s.wfn",molname);
  
  StateOutBinXDR so(ckptfile);
  mole.save_state(so);
  
  tim->print(cout);

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")
// End:
