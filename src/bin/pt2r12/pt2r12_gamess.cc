
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/symmint.h>
#include <chemistry/qc/basis/orthog.h>
#include <chemistry/qc/basis/files.h>
#include <chemistry/qc/basis/cart.h>
#include <chemistry/qc/intv3/intv3.h>
#include <chemistry/qc/lcao/wfnworld.h>
#include <chemistry/qc/lcao/fockbuilder.h>
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

  // GAMESS reports MO coefficients in *cartesian* basis
  // need to transform to the "real" basis?
  if (C_ao.rowdim().n() != basis->nbasis()) { // yes
    integral->set_basis(basis);
    Ref<CartesianBasisSet> cbasis = new CartesianBasisSet(basis, integral);
    assert(C_ao.rowdim().n() == cbasis->nbasis()); // if number of AOs different from computed here, give up

    // now compute overlap between spherical and cartesian basis
    Ref<Integral> localints = integral->clone();
    localints->set_basis(basis,cbasis);
    RefSCMatrix S_sph_cart = detail::onebodyint_ao<&Integral::overlap>(basis, cbasis, localints);

    // SO basis is always blocked, so first make sure S_sph_cart is blocked
    {
      Ref<Integral> braints = integral->clone();  braints->set_basis(basis);
      Ref<PetiteList> brapl = braints->petite_list();
      Ref<Integral> ketints = integral->clone();  ketints->set_basis(cbasis);
      Ref<PetiteList> ketpl = ketints->petite_list();

      RefSCMatrix S_sph_cart_blk = basis->so_matrixkit()->matrix(brapl->AO_basisdim(),ketpl->AO_basisdim());
      S_sph_cart_blk->convert(S_sph_cart);
      S_sph_cart = S_sph_cart_blk;
    }

    // compute projector from cart to spherical basis
    // P(cart->sph) = S^-1 (sph/sph) * S (sph/cart)
    RefSymmSCMatrix S_sph_sph;
    {
      localints->set_basis(basis,basis);
      S_sph_sph = detail::onebodyint<&Integral::overlap>(basis, localints);
    }
    Ref<OverlapOrthog> orthog = new OverlapOrthog(OverlapOrthog::Symmetric, S_sph_sph,
                                                  S_sph_sph.kit(), 1e-8, 0);
    C_ao = orthog->overlap_inverse() * S_sph_cart * C_ao;
  }

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
//  clhf->eigenvectors().print("MO coefficients from CLHF");
//  C_ao.print("MO coefficients from GAMESS");

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

  return 0;
}
