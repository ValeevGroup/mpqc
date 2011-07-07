
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/symmint.h>
#include <chemistry/qc/basis/orthog.h>
#include <chemistry/qc/basis/files.h>
#include <chemistry/qc/intv3/intv3.h>
#include <chemistry/qc/lcao/fockbuild.h>
#include <chemistry/qc/lcao/fockdist.h>
#include <chemistry/qc/lcao/clhfcontrib.h>
#include <math/optimize/scextrap.h>
#include <math/optimize/diis.h>
#include <math/optimize/scextrapmat.h>
#include <util/group/mstate.h>
#include <util/group/pregtime.h>
#include <util/group/memory.h>
#include <util/group/thread.h>
#include <exception>
#include <mpqcinit.h>

using namespace sc;
using std::cout;
using std::endl;

int main_gamess(int argc, char **argv)
{
  GetLongOpt opt;
  opt.usage("[options] [filename]");
  opt.enroll("verbose", GetLongOpt::NoValue, "enable extra printing", 0);

  MPQCInit init(opt,argc,argv);

  int optind = opt.parse(argc, argv);
  std::string inputfile;
  if (argc - optind == 0) {
    inputfile = "hf.in";
  }
  else if (argc - optind == 1) {
    inputfile = argv[optind];
  }
  else {
    opt.usage(std::cout);
    return 1;
  }

  // initialize molecule given description
  Ref<Molecule> molecule = new Molecule;
  {
    molecule->add_atom(8, 0.0, 0.0, 0.0);
  }

  // initialize basis set given name
  Ref<GaussianBasisSet> basis;
  {
    Ref<AssignedKeyVal> tmpkv = new AssignedKeyVal;
    tmpkv->assign("name", "cc-pVDZ");
    tmpkv->assign("molecule", molecule.pointer());
    Ref<KeyVal> kv = tmpkv;
    basis = new GaussianBasisSet(kv);
  }


  Ref<sc::KeyVal> keyval = init.init(inputfile);

  Ref<sc::ThreadGrp> thr = sc::ThreadGrp::get_default_threadgrp();
  Ref<sc::MessageGrp> msg = sc::MessageGrp::get_default_messagegrp();

  if (opt.retrieve("verbose")) {
    sc::ExEnv::out0() << indent
                      << "nthread = " << thr->nthread() << std::endl;
    sc::ExEnv::out0() << indent
                      << "nnode   = " << msg->n() << std::endl;
  }

  //Ref<GaussianBasisSet> basis;
  //basis << keyval->describedclassvalue("basis");

  // make an integral factory and get the petite list
  Ref<Integral> integral = new IntegralV3(basis);
  Ref<PetiteList> pl = integral->petite_list();
  integral->set_storage(32000000);
  Integral::set_default_integral(integral);

  // construct the overlap matrix, s
  RefSymmSCMatrix s_skel(basis->basisdim(), basis->matrixkit());
  Ref<SCElementOp> s_op
      = new OneBodyIntOp(new SymmOneBodyIntIter(integral->overlap(), pl));

  s_skel.assign(0.0);
  s_skel.element_op(s_op);
  s_op=0;

  RefSymmSCMatrix s_SO(pl->SO_basisdim(), basis->so_matrixkit());
  pl->symmetrize(s_skel,s_SO);

  // construct the basis set orthogonalizer
  Ref<OverlapOrthog> orthog
      = new OverlapOrthog(OverlapOrthog::Symmetric,
                          s_SO,
                          basis->so_matrixkit(),
                          1.0e-6 /*lindep tolerance*/,
                          0 /* debug */
          );

  // construct the core hamiltonian (nuclear attraction + kinetic energy)
  RefSymmSCMatrix hcore_skel(basis->basisdim(), basis->matrixkit());
  Ref<SCElementOp> hcore_op
      = new OneBodyIntOp(new SymmOneBodyIntIter(integral->hcore(), pl));

  hcore_skel.assign(0.0);
  hcore_skel.element_op(hcore_op);
  hcore_op=0;

  RefSymmSCMatrix hcore_SO(pl->SO_basisdim(), basis->so_matrixkit());
  pl->symmetrize(hcore_skel,hcore_SO);

  RefSymmSCMatrix density_AO(pl->AO_basisdim(), basis->so_matrixkit());
  density_AO.assign(0.0);

  RefSCMatrix vector_OSO(pl->SO_basisdim(), orthog->orthog_dim(),
                         basis->so_matrixkit());
  vector_OSO.assign(0.0);

  RefDiagSCMatrix evals(orthog->orthog_dim(), basis->so_matrixkit());

  RefSymmSCMatrix density_SO(pl->SO_basisdim(), basis->so_matrixkit());

  Ref<SelfConsistentExtrapolation> extrap = new DIIS;
  extrap->set_tolerance(1.0e-6);

  // read the occupation numbers
  RefSymmSCMatrix occ(orthog->orthog_dim(), basis->so_matrixkit());
  occ.assign(0.0);
  for (int i=0; i<orthog->orthog_dim(); i++) {
      occ(i,i) = keyval->doublevalue("occ",i);
    }

  // construct the initial guess to vector_OSO and density_AO
  RefSymmSCMatrix hcore_OSO(orthog->orthog_dim(), basis->so_matrixkit());
  hcore_OSO.assign(0.0);
  hcore_OSO.accumulate_transform(orthog->basis_to_orthog_basis(),hcore_SO);
  hcore_OSO.diagonalize(evals, vector_OSO);

  density_SO.assign(0.0);
  RefSCMatrix vector_SO = orthog->basis_to_orthog_basis().t() * vector_OSO;
  density_SO.accumulate_transform(vector_SO, occ);
  density_AO = pl->to_AO_basis(density_SO);

  while (!extrap->converged()) {
      //while (!extrap->converged()) {
      // construct the G matrix
      RefSymmSCMatrix g_skel(basis->basisdim(), basis->matrixkit());
      g_skel.assign(0.0);
      Ref<FockContribution> g_contrib
          = new CLHFContribution(basis,basis,basis,"replicated");
      g_contrib->set_fmat(0, g_skel);
      g_contrib->set_pmat(0, density_AO);
      Ref<FockDistribution> fockdist = new FockDistribution;
      Ref<FockBuild> fb = new FockBuild(fockdist, g_contrib, false,
                                        basis, basis, basis);
      fb->set_accuracy(1e-12);
      fb->build();
  
      g_skel.scale(1.0/(double)pl->order());

      RefSymmSCMatrix g_SO(pl->SO_basisdim(), basis->so_matrixkit());
      pl->symmetrize(g_skel,g_SO);

      // compute the exchange/correlation potential
      RefSymmSCMatrix vxc_SO;
      double exc = 0.0;

      // construct the fock matrix
      RefSymmSCMatrix f_SO = hcore_SO + g_SO;

      // transform the fock matrix into the orthogonal basis
      RefSymmSCMatrix f_OSO(orthog->orthog_dim(), basis->so_matrixkit());
      f_OSO.assign(0.0);
      f_OSO.accumulate_transform(orthog->basis_to_orthog_basis(),f_SO);

      // Extrapolate the fock matrix
      Ref<SCExtrapData> data = new SymmSCMatrixSCExtrapData(f_OSO);
      RefSymmSCMatrix error_MO(orthog->orthog_dim(), basis->so_matrixkit());
      error_MO.assign(0.0);
      error_MO.accumulate_transform(vector_OSO.t(),f_OSO);
      error_MO->scale_diagonal(0.0);
  
      RefSymmSCMatrix error_SO(pl->SO_basisdim(), basis->so_matrixkit());
      error_SO.assign(0.0);
      error_SO.accumulate_transform(orthog->basis_to_orthog_basis().t()
                                    *vector_OSO, error_MO);
      Ref<SCExtrapError> error = new SymmSCMatrixSCExtrapError(error_SO);
      extrap->extrapolate(data,error);

      // Diagonalize the fock matrix
      f_OSO.diagonalize(evals, vector_OSO);

      // Compute the density matrix
      density_SO.assign(0.0);
      vector_SO = orthog->basis_to_orthog_basis().t() * vector_OSO;
      density_SO.accumulate_transform(vector_SO, occ);
      density_AO = pl->to_AO_basis(density_SO);

      double energy =exc + basis->molecule()->nuclear_repulsion_energy()
                     +(density_SO * 0.5 * (hcore_SO + f_SO)).trace();

      ExEnv::out0() << indent
                    << scprintf("energy = %14.8f, error = %10.8f",
                                energy, extrap->error()) << std::endl;
    }

  return 0;
}
