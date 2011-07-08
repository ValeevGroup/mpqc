
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
#include <moinfo.h>

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
  Ref<KeyVal> kv = init.init_keyval(MessageGrp::get_default_messagegrp(), inputfile);
  init.init_integrals(kv);

  // print environment
  Ref<sc::ThreadGrp> thr = sc::ThreadGrp::get_default_threadgrp();
  Ref<sc::MessageGrp> msg = sc::MessageGrp::get_default_messagegrp();
  Ref<sc::Integral> integral = sc::Integral::get_default_integral()->clone();
  ExEnv::out0() << indent << "Using " << integral->class_name() << " for integrals by default" << std::endl;
  if (opt.retrieve("verbose")) {
    sc::ExEnv::out0() << indent
                      << "nthread = " << thr->nthread() << std::endl;
    sc::ExEnv::out0() << indent
                      << "nnode   = " << msg->n() << std::endl;
  }

  //
  // Read molecule, basis, and orbitals
  //
  ExternReadMOInfo rdorbs("GAMESS_MOINFO.TXT");
  Ref<GaussianBasisSet> basis = rdorbs.basis();
  RefSCMatrix C_ao = rdorbs.coefs();
  std::vector<unsigned int> orbsym = rdorbs.orbsym();
  const unsigned int nocc = rdorbs.nocc();
  const unsigned int nfzc = rdorbs.nfzc();
  const unsigned int nfzv = rdorbs.nfzv();

  // make an OrbitalSpace object that represents these orbitals
  Ref<OrbitalSpace> orbs = new OrbitalSpace("p", "MOInfo orbitals", C_ao, basis, integral);

  //
  // Read 1-RDM
  //
  Ref<ExternReadRDMOne> rdrdm1 = new ExternReadRDMOne("GAMESS_RDM1.TXT",
                                                      orbs);
  RefSymmSCMatrix P_mo = rdrdm1->scmat(AnySpinCase1);

  //
  // Test orbs and 1-rdm
  //
  // create CLHF object
  Ref<CLHF> clhf;
  {
    Ref<AssignedKeyVal> akv = new AssignedKeyVal;
    akv->assign("molecule", basis->molecule().pointer());
    akv->assign("basis", basis.pointer());
    Ref<KeyVal> kv = akv;
    clhf = new CLHF(kv);
  }

  // compute CLHF object
  sc::ExEnv::out0() << "Energy = " << clhf->energy() << std::endl;
  clhf->eigenvectors().print("MO coefficients from CLHF");

  // create World in which we will compute
  // use defaults for all params
  Ref<WavefunctionWorld> world = new WavefunctionWorld(Ref<KeyVal>(new AssignedKeyVal), clhf);

  // use its orbitals to initialize Extern_RefWavefunction
  P_mo.scale(0.5);
  integral->set_basis(basis);
  Ref<RefWavefunction> ref_wfn = new Extern_RefWavefunction(world, basis, integral,
                                                            C_ao, orbsym,
                                                            P_mo, P_mo,
                                                            nocc, nfzc, nfzv);


#define TEST_R12_INFRASTRUCTURE 0
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
