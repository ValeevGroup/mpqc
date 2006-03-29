
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
#define DEBUG 0

using namespace sc;
using namespace std;

Ref<RegionTimer> tim;
Ref<MessageGrp> grp;

static Ref<MessageGrp>
init_mp(const Ref<KeyVal>& keyval)
{
  grp << keyval->describedclassvalue("message");

  if (grp.nonnull()) MessageGrp::set_default_messagegrp(grp);
  else grp = MessageGrp::get_default_messagegrp();
  
  tim = new ParallelRegionTimer(grp,"scftest",1,0);
  RegionTimer::set_default_regiontimer(tim);

  SCFormIO::set_printnode(0);

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

  // get the integral factory. first try commandline and environment
  Ref<Integral> integral = Integral::initial_integral(argc, argv);
  
  // if we still don't have a integral, try reading the integral
  // from the input
  if (integral.null()) {
    integral << rpkv->describedclassvalue("integrals");
  }

  if (integral.nonnull())
    Integral::set_default_integral(integral);
  else
    integral = Integral::get_default_integral();

  tim->exit("input");

  Ref<CLHF> clhf;
  //clhf << rpkv->describedclassvalue(keyword);

#if !DEBUG
  std::cout << "computing CLHF energy" << std::endl;
  tim->enter("CLHF energy/density");
  clhf = new CLHF(mole);
  double eclhf = clhf->energy();
  tim->exit();

  std::cout << "computing FockBuildCLHF energy" << std::endl;
  tim->enter("New CLHF energy/density");
  Ref<CLHF> newclhf;
  newclhf = new FockBuildCLHF(mole);
  double enewclhf = newclhf->energy();
  tim->exit();

  std::cout << scprintf("E(CLHF)    = %12.8f", eclhf) << std::endl;
  std::cout << scprintf("E(NewCLHF) = %12.8f", enewclhf) << std::endl;
  std::cout << scprintf("DeltaE     = %12.8f", enewclhf-eclhf) << std::endl;

  std::cout << "building multi-basis Fock matrix" << std::endl;
  Ref<GaussianBasisSet> basis2; basis2 << rpkv->describedclassvalue("basis2");
  if (basis2.nonnull()) {
      RefSymmSCMatrix density = clhf->ao_density();
      Ref<GaussianBasisSet> basis1 = clhf->basis();

      std::cout << "basis1:" << std::endl;
      basis1->print();

      std::cout << "basis2:" << std::endl;
      basis2->print();

      // compute the Fock matrix in the original basis
      RefSymmSCMatrix f11(density.dim(),density.kit());
      f11.assign(0.0);
      Ref<FockContribution> fc11
          = new CLHFContribution(basis1,basis1,basis1,basis1);
      fc11->set_fmat(0, f11);
      fc11->set_pmat(0, density);
      clhf->integral()->set_basis(basis1,basis1,basis1,basis1);
      Ref<FockBuild> fb11 = new FockBuild(fc11,1e-12,
                                          basis1,basis1,basis1,basis1);
      fb11->build();
      f11.print("Fock matrix in original basis");

      // compute the Fock matrix in the new basis
      RefSymmSCMatrix f22(basis2->basisdim(),basis2->matrixkit());
      f22.assign(0.0);
      Ref<FockContribution> fc22
          = new CLHFContribution(basis2,basis2,basis1,basis1);
      fc22->set_fmat(0, f22);
      fc22->set_pmat(0, density);
      clhf->integral()->set_basis(basis2,basis2,basis1,basis1);
      Ref<FockBuild> fb22 = new FockBuild(fc22,1e-12,
                                          basis2,basis2,basis1,basis1);
      fb22->build();
      f22.print("Fock matrix in basis2");
    }

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
  Ref<FockBuild> fb = new FockBuild(fc,1e-12,gbs);
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
