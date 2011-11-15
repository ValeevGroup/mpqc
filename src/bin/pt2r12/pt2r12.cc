
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
#include <extern_pt2r12.h>

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
///////////////////////

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
  opt.enroll("prefix", GetLongOpt::MandatoryValue, "mandatory filename prefix, will look for files:\n\
                        $val.pt2r12.dat\n                        $val.pt2r12.rdm2.dat\n                       ", 0);
  opt.enroll("cabs", GetLongOpt::MandatoryValue, "name for CABS; default: construct CABS automatically", 0);
  opt.enroll("dfbs", GetLongOpt::MandatoryValue, "name for DFBS; default: no density fitting", 0);
  opt.enroll("f12exp", GetLongOpt::MandatoryValue, "f12 exponent; default: 1.0", "1.0");
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
  // may receive CABS basis set name
  const char* cabs_name_cstr = opt.retrieve("cabs");
  const std::string cabs_name(cabs_name_cstr ? cabs_name_cstr : "");
  // may receive DFBS basis set name
  const char* dfbs_name_cstr = opt.retrieve("dfbs");
  const std::string dfbs_name(dfbs_name_cstr ? dfbs_name_cstr : "");
  // may receive F12 exponent
  const char* f12exp_cstr = opt.retrieve("f12exp");
  const std::string f12exp_str(f12exp_cstr);

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
  Ref<ExternMOInfo> rdorbs = new ExternMOInfo(filename_prefix + ".pt2r12.dat", integral); // all MO info is contained in rdorbs
  Ref<OrbitalSpace> orbs = rdorbs->orbs();
  Ref<GaussianBasisSet> basis = orbs->basis();
  RefSCMatrix C_ao = orbs->coefs();
  const std::vector<unsigned int>& fzcpi = rdorbs->fzcpi();
  const std::vector<unsigned int>& inactpi = rdorbs->inactpi();
  const std::vector<unsigned int>& actpi = rdorbs->actpi();
  const std::vector<unsigned int>& fzvpi = rdorbs->fzvpi();
  const unsigned int nfzc = std::accumulate(fzcpi.begin(), fzcpi.end(), 0.0);
  const unsigned int ninact = std::accumulate(inactpi.begin(), inactpi.end(), 0.0);
  const unsigned int nact = std::accumulate(actpi.begin(), actpi.end(), 0.0);
  const unsigned int nfzv = std::accumulate(fzvpi.begin(), fzvpi.end(), 0.0);
  const unsigned int nmo = orbs->rank();
  const unsigned int nuocc = nmo - nfzc - ninact - nact;

  if (0) { // test the metric
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
                                    rdorbs->occindexmap_occ(),
                                    occ_orbs);
#else
  // MOLCAS reports the density in active orbitals
  rdrdm2 = new ExternSpinFreeRDMTwo(filename_prefix + ".pt2r12.rdm2.dat",
                                    rdorbs->actindexmap_occ(),
                                    occ_orbs);
#endif

  // create World in which we will compute
  // use defaults for all params
  Ref<WavefunctionWorld> world;
  {
    Ref<AssignedKeyVal> kva = new AssignedKeyVal;
    if (dfbs_name.empty() == false) {
      Ref<AssignedKeyVal> tmpkv = new AssignedKeyVal;
      tmpkv->assign("name", dfbs_name.c_str());
      tmpkv->assign("molecule", basis->molecule().pointer());
      Ref<KeyVal> kv = tmpkv;
      Ref<GaussianBasisSet> df_basis = new GaussianBasisSet(kv);
      kva->assign("df_basis", df_basis.pointer());
    }
    Ref<KeyVal> kv = kva;
    world = new WavefunctionWorld(kv);
  }

  Ref<ExternPT2R12> extern_pt2r12;
  {
    Ref<AssignedKeyVal> kva = new AssignedKeyVal;
    kva->assign("orbs_info", rdorbs.pointer());
    kva->assign("rdm2", rdrdm2.pointer());
    kva->assign("world", world.pointer());
    kva->assign("f12exp", f12exp_str);
    kva->assign("cabs", cabs_name);
    kva->assign("basis", orbs->basis().pointer());
    kva->assign("molecule", orbs->basis()->molecule().pointer());
    Ref<KeyVal> kv = kva;
    extern_pt2r12 = new ExternPT2R12(kv);
  }

  extern_pt2r12->compute();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
