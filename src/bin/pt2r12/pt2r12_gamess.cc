
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
#include <chemistry/qc/mbptr12/pt2r12.h>
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
  opt.usage("[options]");
  opt.enroll("verbose", GetLongOpt::NoValue, "enable extra printing", 0);

  MPQCInit init(opt,argc,argv);

  int optind = opt.parse(argc, argv);
  if (argc - optind != 0) {
    opt.usage(std::cout);
    return 1;
  }
  init.init_integrals();

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
  ExternReadMOInfo rdorbs("R12_COOR.TXT");
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
  { // test the metric
    Ref<Integral> localints = integral->clone();
    RefSymmSCMatrix S_ao = sc::detail::overlap(basis, localints);
    RefSymmSCMatrix S_mo = S_ao.kit()->symmmatrix(C_ao.coldim());
    S_mo.assign(0.0);
    S_mo.accumulate_transform(C_ao, S_ao, SCMatrix::TransposeTransform);
    S_ao.print("AO overlap matrix");
    S_mo.print("MO overlap matrix");
  }

  //
  // Read 1-RDM
  //
  Ref<ExternReadRDMOne> rdrdm1 = new ExternReadRDMOne("R12_RDM1.TXT",
                                                      orbs);
  RefSymmSCMatrix P_mo = rdrdm1->scmat(AnySpinCase1);
  P_mo.eigvals().print("Eigenvalues of 1-RDM read from disk");

#define TEST_GAMESS_DENSITY_AGAINST_PSI 0
#if TEST_GAMESS_DENSITY_AGAINST_PSI
  // Read Psi density and figure out how to transform Psi orbitals to GAMESS orbitals
  // will use this to rotate Psi 2-RDM into GAMESS basis for comparison
  RefSCMatrix U_Psi_to_GAMESS;  // this matrix rotates Psi orbitals to Gamess orbitals
  {
    Ref<ExternReadRDMOne> rdrdm1_psi = new ExternReadRDMOne("R12_RDM1_PSI.TXT",
                                                        orbs);
    RefSymmSCMatrix P_mo_psi = rdrdm1_psi->scmat(AnySpinCase1);
    P_mo_psi.eigvals().print("Eigenvalues of 1-RDM (Psi) read from disk");
    U_Psi_to_GAMESS = P_mo.eigvecs() * P_mo_psi.eigvecs().t();
    U_Psi_to_GAMESS.print("this rotates Psi orbitals to GAMESS orbitals");

    RefSymmSCMatrix P_mo_gpsi = P_mo_psi.clone(); P_mo_gpsi.assign(0.0);
    P_mo_gpsi.accumulate_transform(U_Psi_to_GAMESS, P_mo_psi);
    P_mo_gpsi.print("1-RDM (Psi) transformed to GAMESS orbitals");
    {
      std::ofstream os("R12_RDM1_PSI_GAMESS.TXT");
      const unsigned int nocc = 10;
      for (unsigned int b1 = 0; b1 < nocc; ++b1) {
        for (unsigned int k1 = 0; k1 <= b1; ++k1) {
          const double value = P_mo_gpsi.get_element(b1, k1);
          os << scprintf("%8d%8d%23.17lf",
                         b1 + 1, k1 + 1, value) << std::endl;
        }
      }
      os << scprintf("%8d%8d%23.17lf",
                 -1,-1,-1.0) << std::endl;
      os.close();
    }
  }
#endif

  //
  // Read 2-RDM
  //
  Ref<ExternReadRDMTwo> rdrdm2 = new ExternReadRDMTwo("R12_RDM2.TXT",
                                                      orbs);
  Ref<RDM<One> > rdm1_recomputed = rdrdm2->rdm_m_1();
  rdm1_recomputed->scmat(AnySpinCase1).eigvals().print("Eigenvalues of 1-RDM recomputed from 2-RDM");
  RefSymmSCMatrix P2_mo = rdrdm2->scmat(AnySpinCase2);

#if TEST_GAMESS_DENSITY_AGAINST_PSI
  {
    Ref<ExternReadRDMTwo> rdrdm2_psi = new ExternReadRDMTwo("R12_RDM2_PSI.TXT",
                                                            orbs);
    RefSymmSCMatrix P2_mo_psi = rdrdm2_psi->scmat(AnySpinCase2);
    RefSymmSCMatrix P2_mo_gpsi = P2_mo_psi.clone();
    P2_mo_gpsi.assign(0.0);
    // transform Psi density to GAMESS orbitals
    const unsigned int nmo = U_Psi_to_GAMESS.ncol();
    const unsigned int nocc = 10;
    for(unsigned int b1=0; b1<nocc; ++b1) {
      for(unsigned int b2=0; b2<nocc; ++b2) {
        for(unsigned int k1=0; k1<nocc; ++k1) {
          for(unsigned int k2=0; k2<nocc; ++k2) {
            const unsigned int b12 = b1*nmo + b2;
            const unsigned int k12 = k1*nmo + k2;
            double gpsi_value = 0.0;
            for(unsigned int B1=0; B1<nocc; ++B1) {
              for(unsigned int B2=0; B2<nocc; ++B2) {
                for(unsigned int K1=0; K1<nocc; ++K1) {
                  for(unsigned int K2=0; K2<nocc; ++K2) {
                    const unsigned int B12 = B1*nmo + B2;
                    const unsigned int K12 = K1*nmo + K2;
                    gpsi_value +=
                        U_Psi_to_GAMESS.get_element(b1,B1) *
                        U_Psi_to_GAMESS.get_element(b2,B2) *
                        U_Psi_to_GAMESS.get_element(k1,K1) *
                        U_Psi_to_GAMESS.get_element(k2,K2) *
                        P2_mo_psi.get_element(B12,K12)
                        ;
                  }
                }
              }
            }
            P2_mo_gpsi.set_element(b12,k12,gpsi_value);
          }
        }
      }
    }
    // dump to a file for comparison
    {
      std::ofstream os("R12_RDM2_PSI_GAMESS.TXT");
      for(unsigned int b1=0; b1<nocc; ++b1) {
        for(unsigned int k1=0; k1<=b1; ++k1) {
          for(unsigned int k2=0; k2<=b1; ++k2) {
            const unsigned int b2_max = (k2 == b1) ? k1 : k2;
            for(unsigned int b2=0; b2<=b2_max; ++b2) {
              const unsigned int b12 = b1*nmo + b2;
              const unsigned int k12 = k1*nmo + k2;
              const double value = P2_mo_gpsi.get_element(b12, k12);
              os << scprintf("%8d%8d%8d%8d%23.17lf",
                             b1+1, b2+1, k1+1, k2+1,
                             value) << std::endl;
            }
          }
        }
      }
      os << scprintf("%8d%8d%8d%8d%23.17lf",
                     -1,-1,-1,-1,-1.0) << std::endl;
      os.close();
    }
    // reread and compute 1-rdm
    {
      Ref<ExternReadRDMTwo> rdrdm2 = new ExternReadRDMTwo("R12_RDM2_PSI_GAMESS.TXT",
                                                        orbs);
      Ref<RDM<One> > rdm1_recomputed = rdrdm2->rdm_m_1();
      rdm1_recomputed->scmat(AnySpinCase1).eigvals().print("Eigenvalues of 1-RDM recomputed from 2-RDM transformed from Psi to GAMESS orbitals");
    }

  }
#endif

  //
  // Test orbs and 1-rdm
  //
  // create CLHF object
  Ref<CLHF> clhf;
  {
    Ref<AssignedKeyVal> akv = new AssignedKeyVal;
    akv->assign("molecule", basis->molecule().pointer());
    akv->assign("basis", basis.pointer());
    akv->assign("value_accuracy", 1e-8);
    Ref<KeyVal> kv = akv;
    clhf = new CLHF(kv);
  }

  // compute CLHF object
  sc::ExEnv::out0() << "Energy = " << clhf->energy() << std::endl;
  {
    clhf->eigenvectors().print("MO coefficients (in AO basis) from CLHF");
    C_ao.print("MO coefficients (in AO basis) from GAMESS");
  }
  { // test the metric
    Ref<Integral> localints = integral->clone(); localints->set_basis(basis);
    Ref<PetiteList> plist = localints->petite_list();
    RefSymmSCMatrix S_ao = sc::detail::overlap(basis, localints);
    RefSCMatrix C_ao_mpqc = plist->evecs_to_AO_basis(clhf->eigenvectors());
    RefSymmSCMatrix S_mo = S_ao.kit()->symmmatrix(C_ao.coldim());
    S_mo.assign(0.0);
    S_mo.accumulate_transform(C_ao_mpqc, S_ao, SCMatrix::TransposeTransform);
    S_ao.print("AO overlap matrix");
    S_mo.print("MO overlap matrix");

    RefSCMatrix S12 = C_ao_mpqc.t() * S_ao * C_ao;
    S12.print("Overlap between MPQC and GAMESS MOs");
  }

  // create World in which we will compute
  // use defaults for all params
  Ref<WavefunctionWorld> world = new WavefunctionWorld(Ref<KeyVal>(new AssignedKeyVal));
  world->set_wfn(clhf.pointer());

  // use its orbitals to initialize Extern_RefWavefunction
  P_mo.scale(0.5);
  integral->set_basis(basis);
  Ref<RefWavefunction> ref_wfn = new Extern_RefWavefunction(world, basis, integral,
                                                            C_ao, orbsym,
                                                            P_mo, P_mo,
                                                            nocc, nfzc, nfzv);

  // create PT2R12 object
  Ref<PT2R12> pt2r12;
  {
    Ref<AssignedKeyVal> kva = new AssignedKeyVal;
    kva->assign("molecule", basis->molecule().pointer());
    kva->assign("basis", basis.pointer());
    kva->assign("refwfn", ref_wfn.pointer());
    kva->assign("world", world.pointer());
    kva->assign("rdm2", rdrdm2.pointer());
    kva->assign("corr_factor", "stg-6g");
    kva->assign("corr_param", "1.0");
    {
      Ref<AssignedKeyVal> tmpkv = new AssignedKeyVal;
      tmpkv->assign("name", "cc-pVDZ-F12-CABS");
      tmpkv->assign("puream", "true");
      tmpkv->assign("molecule", basis->molecule().pointer());
      Ref<KeyVal> kv = tmpkv;
      Ref<GaussianBasisSet> aux_basis = new GaussianBasisSet(kv);
      kva->assign("aux_basis", aux_basis.pointer());
    }
    Ref<KeyVal> kv = kva;

    pt2r12 = new PT2R12(kv);
    const double ept2r12 = pt2r12->value();
  }

  return 0;
}
