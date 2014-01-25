
#include <util/misc/bug.h>
#include <util/ref/ref.h>
#include <util/misc/bug.h>
#include <util/group/pregtime.h>
#include <util/group/message.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/scf/clhf.h>
#include <chemistry/qc/scf/fbclhf.h>
#include <chemistry/qc/lcao/fockdist.h>

#include <chemistry/qc/scf/linkage.h>

#undef DEBUG
#define DEBUG 0

using namespace sc;
using namespace std;

Ref<MessageGrp> grp;

static Ref<ThreadGrp>
init_thread(const Ref<KeyVal>& keyval, int &argc, char **&argv)
{
  ///////////////////////////////////////////////////////////////
  // Initialize the thread group

  // get the thread group.  first try the commandline and environment
  Ref<ThreadGrp> thread = ThreadGrp::initial_threadgrp(argc, argv);
  
  // if we still don't have a group, try reading the thread group
  // from the input
  if (thread == 0) {
    thread << keyval->describedclassvalue("thread");
  }

  if (thread)
    ThreadGrp::set_default_threadgrp(thread);
  else
    thread = ThreadGrp::get_default_threadgrp();

  return thread;
}

static Ref<MessageGrp>
init_message(const Ref<KeyVal>& keyval,int &argc, char **&argv)
{
  Ref<MessageGrp> grp;

  grp << keyval->describedclassvalue("message");

  if (grp == 0) grp = MessageGrp::initial_messagegrp(argc, argv);

  if (grp == 0) grp = MessageGrp::get_default_messagegrp();

  MessageGrp::set_default_messagegrp(grp);

  RegionTimer::set_default_regiontimer(
      new ParallelRegionTimer(grp,"focktest",1,0));

  SCFormIO::set_printnode(0);

  SCFormIO::setindent(std::cout, 2);
  SCFormIO::setindent(std::cerr, 2);
  
  return grp;
}

static void
clean_up(void)
{
  MemoryGrp::set_default_memorygrp(0);
  MessageGrp::set_default_messagegrp(0);
  ThreadGrp::set_default_threadgrp(0);
  SCMatrixKit::set_default_matrixkit(0);
  Integral::set_default_integral(0);
  RegionTimer::set_default_regiontimer(0);
}

int
try_main(int argc, char **argv)
{
  const char *input =      (argc > 1)? argv[1] : SRCDIR "/focktest.in";
  const char *keyword =    (argc > 2)? argv[2] : "mole";
  const char *optkeyword = (argc > 3)? argv[3] : "opt";

  // open keyval input
  Ref<KeyVal> rpkv(new ParsedKeyVal(input));

  Ref<KeyVal> mole(new PrefixKeyVal(rpkv,"mole"));

  ///////////////////////////////////////////////////////////////

  Ref<ThreadGrp> thr = init_thread(rpkv,argc,argv);
  Ref<MessageGrp> msg = init_message(rpkv,argc,argv);

  ExEnv::out0() << indent
                << "nthread = " << thr->nthread() << std::endl;
  ExEnv::out0() << indent
                << "nnode   = " << msg->n() << std::endl;

  ///////////////////////////////////////////////////////////////

  Ref<RegionTimer> regtim;
  if (rpkv->exists("timer")) regtim << rpkv->describedclassvalue("timer");
  else                       regtim = new ParallelRegionTimer(msg,"focktest",1,1);
  RegionTimer::set_default_regiontimer(regtim);

  ///////////////////////////////////////////////////////////////

  Timer tim;
  tim.enter("input");

  Ref<Debugger> debugger; debugger << rpkv->describedclassvalue("debug");
  if (debugger) {
      Debugger::set_default_debugger(debugger);
      debugger->set_exec(argv[0]);
      debugger->set_prefix(msg->me());
    }
  
  if (rpkv->exists("matrixkit")) {
    Ref<SCMatrixKit> kit; kit << rpkv->describedclassvalue("matrixkit");
    SCMatrixKit::set_default_matrixkit(kit);
  }

  // get the integral factory. first try commandline and environment
  Ref<Integral> integral = Integral::initial_integral(argc, argv);
  
  // if we still don't have a integral, try reading the integral
  // from the input
  if (integral == 0) {
    integral << rpkv->describedclassvalue("integrals");
  }

  if (integral)
    Integral::set_default_integral(integral);
  else
    integral = Integral::get_default_integral();

  tim.exit("input");

  Ref<CLHF> clhf;
  //clhf << rpkv->describedclassvalue(keyword);

  std::string fockbuildmatrixtype = rpkv->stringvalue("fockbuildmatrixtype");

  bool prefetch_blocks = false;
  if (fockbuildmatrixtype == "prefetched_distributed") {
      prefetch_blocks = true;
  }

#if !DEBUG

  KeyValValueboolean def_compute_clhf_energy(1);
  bool compute_clhf_energy = rpkv->booleanvalue("compute_clhf_energy",
                                                def_compute_clhf_energy);

  KeyValValueboolean def_compute_fockbuild_energy(1);
  bool compute_fockbuild_energy = rpkv->booleanvalue("compute_fockbuild_energy",
                                                     def_compute_fockbuild_energy);

  double eclhf;
  if (compute_clhf_energy) {
      ExEnv::out0() << indent << "computing CLHF energy" << std::endl;
      tim.enter("CLHF energy/density");
      clhf = new CLHF(mole);
      eclhf = clhf->energy();
      tim.exit();
    }

  double enewclhf;
  Ref<CLHF> newclhf;
  if (compute_fockbuild_energy) {
      ExEnv::out0() << indent << "computing FockBuildCLHF energy" << std::endl;
      tim.enter("New CLHF energy/density");
      newclhf = new FockBuildCLHF(mole);
      enewclhf = newclhf->energy();
      tim.exit();
    }

  if (compute_clhf_energy) {
      ExEnv::out0() << indent << scprintf("E(CLHF)    = %12.8f", eclhf) << std::endl;
    }
  if (compute_fockbuild_energy) {
      ExEnv::out0() << indent << scprintf("E(NewCLHF) = %12.8f", enewclhf) << std::endl;
    }
  if (compute_clhf_energy && compute_fockbuild_energy) {
      ExEnv::out0() << indent << scprintf("DeltaE     = %12.8f", enewclhf-eclhf) << std::endl;
    }
  if (compute_fockbuild_energy) {
      newclhf->print();
    }

  Ref<FockDistribution> fockdist = new FockDistribution;

  ExEnv::out0() << indent << "building multi-basis Fock matrix" << std::endl;
  Ref<GaussianBasisSet> basis2; basis2 << rpkv->describedclassvalue("basis2");
  if (basis2) {
      RefSymmSCMatrix density = clhf->ao_density();
      Ref<GaussianBasisSet> basis1 = clhf->basis();

      ExEnv::out0() << indent << "basis1:" << std::endl;
      basis1->print();

      ExEnv::out0() << indent << "basis2:" << std::endl;
      basis2->print();

      // compute the Fock matrix in the original basis
      RefSymmSCMatrix f11(density.dim(),density.kit());
      f11.assign(0.0);
      Ref<FockContribution> fc11
          = new CLHFContribution(basis1,basis1,basis1,fockbuildmatrixtype);
      fc11->set_fmat(0, f11);
      fc11->set_pmat(0, density);
      Ref<FockBuild> fb11 = new FockBuild(fockdist, fc11, prefetch_blocks,
                                          basis1,basis1,basis1);
      fb11->set_accuracy(1e-12);
      fb11->build();
      f11.print("Fock matrix in original basis");

      // compute the Fock matrix in the new basis
      RefSymmSCMatrix f22(basis2->basisdim(),basis2->matrixkit());
      f22.assign(0.0);
      Ref<FockContribution> fc22
          = new CLHFContribution(basis2,basis2,basis1,fockbuildmatrixtype);
      fc22->set_fmat(0, f22);
      fc22->set_pmat(0, density);
      Ref<FockBuild> fb22 = new FockBuild(fockdist, fc22, prefetch_blocks,
                                          basis2,basis2,basis1);
      fb22->set_accuracy(1e-12);
      fb22->build();
      f22.print("Fock matrix in basis2");
    }

#else
  ExEnv::out0() << indent << "computing CLHF energy" << std::endl;
  tim.enter("CLHF energy/density");
  clhf = new CLHF(mole);
  double eclhf = clhf->energy();
  tim.exit();

  RefSymmSCMatrix density = clhf->ao_density();

  tim.enter("CLHF Fock build");
  Ref<GaussianBasisSet> gbs = clhf->basis();
  gbs->print();
  RefSymmSCMatrix f(density.dim(),density.kit());
  f.assign(0.0);
  Ref<FockContribution> fc
      = new CLHFContribution(gbs,gbs,gbs,fockbuildmatrixtype);
  fc->set_fmat(0, f);
  fc->set_pmat(0, density);
  Ref<FockBuild> fb = new FockBuild(fockdist,fc,prefetch_blocks,gbs);
  fb->set_accuracy(1e-12);
  fb->build();
  tim.exit("CLHF Fock build");

  f.print("FockBuilder Skeleton G Matrix AO Basis (CLHF Density)");
  // now symmetrize the skeleton G matrix in f
  Ref<PetiteList> pl = clhf->integral()->petite_list();
  f.scale(1.0/(double)pl->order());
  f.print("Scaled FockBuilder Skeleton G Matrix AO Basis (CLHF Density)");
  RefSymmSCMatrix fsymm = clhf->fock(0)->clone();
  pl->symmetrize(f,fsymm);
  fsymm.print("FockBuilder G Matrix SO Basis (CLHF Density)");

  density.print("Density");
  (clhf->fock(0)-clhf->core_hamiltonian())
      ->print("CLHF Fock Matrix Electronic Contribution");
  
  (clhf->fock(0)-clhf->core_hamiltonian()-fsymm)
      ->print("Difference");
#endif

  tim.print(ExEnv::out0());

  MessageGrp::set_default_messagegrp(0);
  ThreadGrp::set_default_threadgrp(0);

  clean_up();

  return 0;
}

// extern "C" begin_vmon();
// extern "C" end_vmon();

int
main(int argc, char *argv[])
{
//   begin_vmon();
  try {
      try_main(argc, argv);
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
//   end_vmon();
  return 0;
}
