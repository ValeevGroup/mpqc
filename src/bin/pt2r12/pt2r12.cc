
#include <iostream>
#include <exception>
#include <numeric>

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
#include <mpqcinit.h>
#include <moinfo.h>

// Force linkages:
#include <scdirlist.h>
#include <util/group/linkage.h>
#include <chemistry/qc/basis/linkage.h>
#ifdef HAVE_SC_SRC_LIB_CHEMISTRY_QC_PSI
#include <chemistry/qc/psi/linkage.h>
#endif
#include <chemistry/qc/wfn/linkage.h>
#include <chemistry/qc/scf/linkage.h>
#ifdef HAVE_SC_SRC_LIB_CHEMISTRY_QC_MBPTR12
#  include <chemistry/qc/mbptr12/linkage.h>
#endif
#ifdef HAVE_SC_SRC_LIB_CHEMISTRY_QC_LIBINT2
#  include <chemistry/qc/libint2/linkage.h>
#endif
#include <util/state/linkage.h>

#ifdef HAVE_MPI
#define MPICH_SKIP_MPICXX
#include <mpi.h>
#include <util/group/messmpi.h>
#endif

using std::cout;
using std::endl;
using namespace sc;

extern int try_main(int argc, char *argv[]);

int
main(int argc, char *argv[])
{
  try {
    try_main(argc, argv);
  }
  catch (std::bad_alloc &e) {
      cout << argv[0] << ": ERROR: MEMORY ALLOCATION FAILED:" << endl
           << e.what()
           << endl;
      exit(1);
    }
  catch (std::exception &e) {
      cout << argv[0] << ": ERROR: EXCEPTION RAISED:" << endl
           << e.what()
           << endl;
      exit(1);
    }
  catch (...) {
      cout << argv[0] << ": ERROR: UNKNOWN EXCEPTION RAISED" << endl;
      exit(1);
    }
  return 0;
}

int try_main(int argc, char **argv)
{
  const bool debug_print = false;
  GetLongOpt opt;
  opt.usage("[options]");
  opt.enroll("prefix", GetLongOpt::MandatoryValue, "filename prefix for input data", 0);
  opt.enroll("verbose", GetLongOpt::NoValue, "enable extra printing", 0);

  MPQCInit init(opt,argc,argv);

  int optind = opt.parse(argc, argv);
  // no misc command-line options allowed
  if (argc - optind != 0) {
    opt.usage(std::cout);
    return 1;
  }
  // must receive prefix
  const char* filename_prefix_cstr = opt.retrieve("prefix");
  if (filename_prefix_cstr == 0) {
    opt.usage(std::cout);
    return 1;
  }
  const std::string filename_prefix(filename_prefix_cstr);
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
  ExternMOInfo rdorbs(filename_prefix + ".pt2r12.dat", integral); // all MO info is contained in rdorbs
  Ref<OrbitalSpace> orbs = rdorbs.orbs();

  Ref<GaussianBasisSet> basis = orbs->basis();
  RefSCMatrix C_ao = orbs->coefs();
  const std::vector<unsigned int>& fzcpi = rdorbs.fzcpi();
  const std::vector<unsigned int>& inactpi = rdorbs.inactpi();
  const std::vector<unsigned int>& actpi = rdorbs.actpi();
  const std::vector<unsigned int>& fzvpi = rdorbs.fzvpi();
  const std::vector<unsigned int>& mopi = rdorbs.mopi();
  std::vector<unsigned int> occpi;
  const unsigned int nirrep = fzcpi.size();
  for (int i = 0; i < nirrep; ++i)
  {
    occpi.push_back(fzcpi[i] + inactpi[i] + actpi[i]);
  }
  const unsigned int nfzc = std::accumulate(fzcpi.begin(), fzcpi.end(), 0.0);
  const unsigned int ninact = std::accumulate(inactpi.begin(), inactpi.end(), 0.0);
  const unsigned int nact = std::accumulate(actpi.begin(), actpi.end(), 0.0);
  const unsigned int nfzv = std::accumulate(fzvpi.begin(), fzvpi.end(), 0.0);
  const unsigned int nmo = orbs->rank();
  const unsigned int nuocc = nmo - nfzc - ninact - nact;

  { // test the metric
    Ref<Integral> localints = integral->clone();
    RefSymmSCMatrix S_so = sc::detail::overlap(basis, localints);
    localints->set_basis(basis);
    RefSymmSCMatrix S_ao = localints->petite_list()->to_AO_basis(S_so);
    S_ao.print("AO overlap matrix");
    RefSCMatrix C_ao_bsdim = C_ao.kit()->matrix(S_ao.dim(), C_ao.coldim());
    C_ao_bsdim->convert(C_ao);
    RefSymmSCMatrix S_mo = C_ao.kit()->symmmatrix(C_ao.coldim());
    S_mo.assign(0.0);
    S_mo.accumulate_transform(C_ao_bsdim, S_ao, SCMatrix::TransposeTransform);
    S_ao.print("AO overlap matrix");
    S_mo.print("MO overlap matrix");
  }

  basis = orbs->basis();
  C_ao = orbs->coefs();

  /////////////////////////////////////////////
  // Read 2-RDM
  /////////////////////////////////////////////

  // molcas reports 2-RDM in terms of active occupied orbitals only, indexed occording to molcas convention
  // thus use the map from molcas active occupied orbitals to MPQC occupied range
  // first make an OrbitalSpace for MPQC occupied orbitals
  Ref<OrbitalSpace> occ_orbs = new OrbitalSpace(std::string("z(sym)"),
                                                std::string("symmetry-ordered occupied MOInfo orbitals"),
                                                orbs->coefs(),
                                                orbs->basis(),
                                                orbs->integral(), orbs->evals(),
                                                0, nuocc + nfzv,
                                                OrbitalSpace::symmetry);
#if 1
  sc::ExEnv::out0() << "debug: print occ and orbs" << std::endl;
  occ_orbs->print_detail();
  orbs->print_detail();
#endif

  // Currently we support two choices for reporting 2-RDM:
  // in active space only (MOLCAS) or in the entire occupied space (GAMESS)
  Ref<ExternSpinFreeRDMTwo> rdrdm2;
#if PT2R12GAMESS
  // GAMESS reports the density in occupied orbitals
  rdrdm2 = new ExternSpinFreeRDMTwo(filename_prefix + ".pt2r12.rdm2.dat",
                                    rdorbs.occindexmap_occ(),
                                    occ_orbs);
#else
  // MOLCAS reports the density in active orbitals
  rdrdm2 = new ExternSpinFreeRDMTwo(filename_prefix + ".pt2r12.rdm2.dat",
                                    rdorbs.actindexmap_occ(),
                                    occ_orbs);
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
  if (1) { // test the metric
    Ref<Integral> localints = integral->clone(); localints->set_basis(basis);
    Ref<PetiteList> plist = localints->petite_list();

    RefSCMatrix C_ao_clhf = plist->evecs_to_AO_basis(clhf->eigenvectors());
    C_ao_clhf.print("MO coefficients (in AO basis) from CLHF");
    C_ao.print("MO coefficients (in AO basis) from host program");

    RefSymmSCMatrix S_so = sc::detail::overlap(basis, localints);
    RefSymmSCMatrix S_ao = localints->petite_list()->to_AO_basis(S_so);
    S_ao.print("AO overlap matrix");
    RefSymmSCMatrix S_mo = C_ao_clhf.kit()->symmmatrix(C_ao_clhf.coldim());
    S_mo.assign(0.0);
    S_mo.accumulate_transform(C_ao_clhf, S_ao, SCMatrix::TransposeTransform);
    S_mo.print("MO overlap matrix");

    RefSCMatrix C_ao_bsdim = C_ao.kit()->matrix(S_ao.dim(), C_ao.coldim());
    C_ao_bsdim->convert(C_ao);
    RefSCMatrix S12 = C_ao_clhf.t() * S_ao * C_ao_bsdim;
    S12.print("Overlap between CLHF and host program MOs");
  }

  // create World in which we will compute
  // use defaults for all params
  Ref<WavefunctionWorld> world = new WavefunctionWorld(Ref<KeyVal>(new AssignedKeyVal));
  world->set_wfn(clhf.pointer());

  // to construct Extern_RefWavefunction we need an energy/correlation-ordered orbitals
  // make 1-RDM to the full MO space
  RefSymmSCMatrix P1_mo;
  {
    Ref<SpinFreeRDM<One> > rdrdm1 = rdrdm2->rdm_m_1();
    RefSymmSCMatrix P1_mo_occ = rdrdm1->scmat();
    std::vector<unsigned int> occ_to_orbs_indexmap;
    occ_to_orbs_indexmap = (*orbs) << *(rdrdm1->orbs());
    P1_mo = P1_mo_occ.kit()->symmmatrix(orbs->dim());
    P1_mo.assign(0.0);
    const unsigned int nocc = P1_mo_occ.n();
    for(unsigned int i1=0; i1<nocc; ++i1) {
      const unsigned int ii1 = occ_to_orbs_indexmap[i1];
      for(unsigned int i2=0; i2<=i1; ++i2) {
        const unsigned int ii2 = occ_to_orbs_indexmap[i2];
        P1_mo.set_element(ii1, ii2, P1_mo_occ.get_element(i1, i2));
      }
    }
    P1_mo.scale(0.5);
  }

  if(debug_print)
  {
    sc::ExEnv::out0() << "debug:print before refwfn" << std::endl;
    orbs->print_detail();
    P1_mo->print("P1_mo");
    sc::ExEnv::out0() << std::endl<<orbs->orbsym()[0] << " " << orbs->orbsym()[1] << " " << orbs->orbsym()[2] <<" "<< orbs->orbsym()[3] << std::endl;
  }

  // use its orbitals to initialize Extern_RefWavefunction
  integral->set_basis(basis);
// the following constructor doesn't work correctly
//  Ref<RefWavefunction> ref_wfn = new Extern_RefWavefunction(world, basis, integral,
//                                                            orbs->coefs(), orbs->orbsym(),
//                                                            P1_mo, P1_mo,
//                                                            nfzc+ninact+nact,
//                                                            nfzc,
//                                                            nfzv);
  Ref<RefWavefunction> ref_wfn = new Extern_RefWavefunction(world, basis, integral,
                                                            orbs->coefs(), orbs->orbsym(),
                                                            P1_mo, P1_mo,
                                                            mopi,
                                                            occpi,
                                                            fzcpi,
                                                            fzvpi);
  if(debug_print)
  {
    sc::ExEnv::out0() << "debug:print refwfn orbs " << std::endl;
    ref_wfn->occ_act(Alpha)->print_detail();
    ref_wfn->orbs(Alpha)->print_detail();
  }

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
    kva->assignboolean("spinadapted", 1);
    kva->assignboolean("pt2_correction", 1);
    Ref<KeyVal> kv = kva;

    pt2r12 = new PT2R12(kv);
    const double ept2r12 = pt2r12->value();
  }

  return 0;
}
