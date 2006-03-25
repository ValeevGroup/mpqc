
#include <util/ref/ref.h>
#include <util/misc/bug.h>
#include <util/group/pregtime.h>
#include <util/group/message.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/scf/fockbuild.h>
#include <chemistry/qc/scf/clhfcontrib.h>
#include <chemistry/qc/scf/clhf.h>
#include <chemistry/qc/scf/fbclhf.h>

#include <chemistry/qc/scf/linkage.h>

#undef DEBUG
#define DEBUG 1

using namespace sc;
using namespace std;

Ref<RegionTimer> tim;
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
  
  tim = new ParallelRegionTimer(grp,"scftest",1,0);
  RegionTimer::set_default_regiontimer(tim);

  SCFormIO::set_printnode(0);
  //SCFormIO::set_debug(1);

  SCFormIO::setindent(cout, 2);
  SCFormIO::setindent(cerr, 2);
  
  return grp;
}

int
try_main(int argc, char **argv)
{
  const char *input =      (argc > 1)? argv[1] : SRCDIR "/mpqc.in";
  const char *keyword =    (argc > 2)? argv[2] : "mole";
  const char *optkeyword = (argc > 3)? argv[3] : "opt";

  // open keyval input
  Ref<KeyVal> rpkv(new ParsedKeyVal(input));

  Ref<KeyVal> mole(new PrefixKeyVal(rpkv,"mole"));

  init_mp(rpkv);

  tim->enter("input");
  
  if (rpkv->exists("matrixkit")) {
    Ref<SCMatrixKit> kit; kit << rpkv->describedclassvalue("matrixkit");
    SCMatrixKit::set_default_matrixkit(kit);
  }

  tim->exit("input");

  Ref<CLHF> clhf;
  //clhf << rpkv->describedclassvalue(keyword);

#if !DEBUG
  std::cout << "computing CLHF energy" << std::endl;
  tim->enter("CLHF energy/density");
  clhf = new CLHF(mole);
  double eclhf = clhf->energy();
  tim->exit();

  std::cout << "computing NewCLHF energy" << std::endl;
  tim->enter("New CLHF energy/density");
  Ref<CLHF> newclhf;
  newclhf = new FockBuildCLHF(mole);
  double enewclhf = newclhf->energy();
  tim->exit();

  std::cout << scprintf("E(CLHF)    = %12.8f", eclhf) << std::endl;
  std::cout << scprintf("E(NewCLHF) = %12.8f", enewclhf) << std::endl;
  std::cout << scprintf("DeltaE     = %12.8f", enewclhf-eclhf) << std::endl;

#else
  std::cout << "computing CLHF energy" << std::endl;
  tim->enter("CLHF energy/density");
  clhf = new CLHF(mole);
  double eclhf = clhf->energy();
  tim->exit();

  RefSymmSCMatrix density = clhf->ao_density();

  tim->enter("CLHF Fock build");
  Ref<GaussianBasisSet> gbs = clhf->basis();
  gbs->print();
  RefSymmSCMatrix f(density.dim(),density.kit());
  f.assign(0.0);
  Ref<FockContribution> fc
      = new CLHFContribution(gbs,gbs,gbs,gbs);
  fc->set_fmat(0, f);
  fc->set_pmat(0, density);
  Ref<TwoBodyInt> eri = clhf->integral()->electron_repulsion();
  Ref<FockBuild> fb = new FockBuild(fc,1e-12,&eri,gbs);
  fb->build();
  tim->exit("CLHF Fock build");

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

  tim->print(ExEnv::out0());

  return 0;
}


int
main(int argc, char *argv[])
{
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
  return 0;
}
