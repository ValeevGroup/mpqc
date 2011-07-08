
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/symmint.h>
#include <chemistry/qc/basis/orthog.h>
#include <chemistry/qc/basis/files.h>
#include <chemistry/qc/intv3/intv3.h>
#include <chemistry/qc/lcao/wfnworld.h>
#include <chemistry/qc/mbptr12/ref.h>
#include <chemistry/qc/scf/clhf.h>
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
    inputfile = "pt2r12.in";
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
  molecule->set_point_group(new PointGroup("d2h"));

  // initialize basis set given name
  Ref<GaussianBasisSet> basis;
  {
    Ref<AssignedKeyVal> tmpkv = new AssignedKeyVal;
    tmpkv->assign("name", "cc-pVDZ");
    tmpkv->assign("molecule", molecule.pointer());
    Ref<KeyVal> kv = tmpkv;
    basis = new GaussianBasisSet(kv);
  }

  // print environment
  Ref<sc::ThreadGrp> thr = sc::ThreadGrp::get_default_threadgrp();
  Ref<sc::MessageGrp> msg = sc::MessageGrp::get_default_messagegrp();
  Ref<sc::Integral> integral = sc::Integral::get_default_integral()->clone();
  if (opt.retrieve("verbose")) {
    sc::ExEnv::out0() << indent
                      << "nthread = " << thr->nthread() << std::endl;
    sc::ExEnv::out0() << indent
                      << "nnode   = " << msg->n() << std::endl;
  }

#define TEST_R12_INFRASTRUCTURE 1
#if TEST_R12_INFRASTRUCTURE
  {
    // create CLHF object
    Ref<CLHF> clhf;
    {
      Ref<AssignedKeyVal> akv = new AssignedKeyVal;
      akv->assign("molecule", molecule.pointer());
      akv->assign("basis", basis.pointer());
      Ref<KeyVal> kv = akv;
      clhf = new CLHF(kv);
    }

    // compute CLHF object
    sc::ExEnv::out0() << "Energy = " << clhf->energy() << std::endl;

    // grab orbitals and their symmetries
    RefSCMatrix C_so = clhf->eigenvectors();
    integral->set_basis(basis);
    Ref<PetiteList> pl = integral->petite_list();
    RefSCMatrix C_ao = pl->evecs_to_AO_basis(C_so);
    std::vector<unsigned int> orbsym(C_ao.coldim().n(), 0);
    {
      Ref<SCBlockInfo> blocks = C_ao.coldim()->blocks();
      if (blocks.nonnull()) {
        const unsigned int nirreps = blocks->nblock();
        for(unsigned int irrep=0; irrep<nirreps; ++irrep) {
          const unsigned int first_orb_in_irrep = blocks->start(irrep);
          const unsigned int norbs_in_irrep = blocks->size(irrep);
          const unsigned int plast_orb_in_irrep = first_orb_in_irrep + norbs_in_irrep;
          for(unsigned int o=first_orb_in_irrep; o<plast_orb_in_irrep; ++o) {
            orbsym[o] = irrep;
          }
        }
      }
    }

    // grab densities
    RefSymmSCMatrix P_so = clhf->density(); P_so.scale(0.5);
    RefSymmSCMatrix P_mo = P_so.kit()->symmmatrix(C_ao.coldim());  P_mo.assign(0.0);
    P_mo.accumulate_transform(clhf->mo_to_so(), P_so, SCMatrix::TransposeTransform);
    P_mo.print("MO density");

    // create World in which we will compute
    // use defaults for all params
    Ref<WavefunctionWorld> world = new WavefunctionWorld(Ref<KeyVal>(new AssignedKeyVal), clhf);

    // use its orbitals to initialize Extern_RefWavefunction
    Ref<RefWavefunction> ref_wfn = new Extern_RefWavefunction(world, basis, integral,
                                                              C_ao, orbsym,
                                                              P_mo, P_mo,
                                                              clhf->nelectron()/2, 3, 1);
  }
#endif
  
  return 0;
}
