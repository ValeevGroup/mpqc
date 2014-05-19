
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
#include <chemistry/qc/nbody/ref.h>
#include <chemistry/qc/mbptr12/pt2r12.h>
#include <chemistry/qc/scf/clhf.h>
#include <math/optimize/scextrap.h>
#include <math/optimize/diis.h>
#include <math/optimize/scextrapmat.h>
#include <util/group/mstate.h>
#include <util/group/pregtime.h>
#include <util/group/memory.h>
#include <util/group/thread.h>
#include <util/misc/consumableresources.h>
#include "../mpqc/mpqcinit.h"

#ifdef HAVE_MADNESS
# include <util/madness/init.h>
#endif

#include "moinfo.h"
#include "extern_pt2r12.h"

// Force linkages:
#include <mpqc_config.h>
#include <util/group/linkage.h>
#include <chemistry/qc/basis/linkage.h>
#ifdef HAVE_PSI3
#include <chemistry/qc/psi/linkage.h>
#endif
#include <chemistry/qc/wfn/linkage.h>
#include <chemistry/qc/scf/linkage.h>
#ifdef HAVE_LIBINT2
#  include <chemistry/qc/libint2/linkage.h>
#  include <chemistry/qc/mbptr12/linkage.h>
#endif
#include <util/state/linkage.h>

#ifdef HAVE_MPI
#define OMPI_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#include <mpi.h>
#include <util/group/messmpi.h>
#endif

// MUST HAVE LIBINT2
#ifdef HAVE_LIBINT2
#  include <chemistry/qc/libint2/libint2.h>
#else
#  error "this copy of MPQC does not include the Libint2 library (see libint.valeyev.net) -- cannot use F12 methods"
#endif

using std::cout;
using std::endl;
using namespace sc;
///////////////////////

extern int try_main(int argc, char *argv[]);

#ifndef SKIP_MAIN
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
#else

static std::vector<std::string> argv_;
extern "C" void add_arg_(const char * arg, long len){
    std::string s(arg,len);
    s.erase(s.find_last_not_of(' ')+1);
    argv_.push_back(s);
}

extern int try_main(int argc, char *argv[]);

char* arg(const char * value){
    char *p = new char[strlen(value)];
    strcpy(p,value);
    return p;
}

extern "C" int pt2r12_main_() {
    
    int argc = argv_.size() +1;
    char* argv[argc];
    
    
    argv[0] = arg("");
    for (int i=0; i<argv_.size(); i++){
	argv[i+1]= arg(argv_[i].c_str());
    }
    
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

#endif
int try_main(int argc, char **argv)
{
  const bool debug_print = false;

  // pt2r12 receives user instructions via command line
  GetLongOpt opt;
  opt.usage("[options]");
  opt.enroll("prefix", GetLongOpt::MandatoryValue, "mandatory filename prefix, will look for files:\n\
                        $val.pt2r12.dat\n                        $val.pt2r12.rdm2.dat\n                       ", 0);
  opt.enroll("obs", GetLongOpt::MandatoryValue, "name for the orbital basis set; optional, if given will be used to set the defaults for cabs, dfbs, and f12exp", 0);
  opt.enroll("cabs", GetLongOpt::MandatoryValue, "name for CABS; default: construct CABS automatically", 0);
  opt.enroll("dfbs", GetLongOpt::MandatoryValue, "name for DFBS; default: no density fitting; use \"none\" to override the default for the obs", 0);
  opt.enroll("f12exp", GetLongOpt::MandatoryValue, "f12 exponent; default: 1.0", "1.0");
  opt.enroll("r12", GetLongOpt::MandatoryValue, "compute [2]_R12 correction; default: true", 0);
  opt.enroll("verbose", GetLongOpt::NoValue, "enable extra printing", 0);
#if defined(HAVE_MPQC3_RUNTIME)
  opt.enroll("singles", GetLongOpt::MandatoryValue, "compute [2]_s correction; default: false", 0);
  opt.enroll("partitionH", GetLongOpt::MandatoryValue, "How to partition Hamiltonian in [2]_s: fock, dyall_1, dyall_2; default: fock", 0);
  opt.enroll("mpqc3", GetLongOpt::MandatoryValue, "enable MPQC3 runtime features; default: true", 0);
#endif

  // initialize the environment
  MPQCInit init(opt,argc,argv);
  ExEnv::init(argc, argv);
  init.init_fp();
  init.init_limits();
  Ref<MessageGrp> grp = init.init_messagegrp();
  init.init_io(grp);
  init.init_timer(grp,0);

#ifdef HAVE_MADNESS
  MADNESSRuntime::initialize();
#endif

  Timer timer;

  const char *tstr = 0;
#if defined(HAVE_TIME) && defined(HAVE_CTIME)
  time_t t;
  time(&t);
  tstr = ctime(&t);
#endif
  if (!tstr) {
    tstr = "UNKNOWN";
  }
  ExEnv::out0()
       << indent << scprintf("Machine:    %s", TARGET_ARCH) << endl
       << indent << scprintf("User:       %s@%s",
                             ExEnv::username(), ExEnv::hostname()) << endl
       << indent << scprintf("Start Time: %s", tstr) << endl;

  Ref<ThreadGrp> thread = init.init_threadgrp();
  Ref<MemoryGrp> memory = init.init_memorygrp();
#ifdef HAVE_LIBINT2
  Integral::set_default_integral(new IntegralLibint2);
#endif
  init.init_integrals();
  Ref<Integral> integral = Integral::get_default_integral()->clone();
  init.init_resources();
  Ref<ConsumableResources> resources =
      ConsumableResources::get_default_instance();

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
  // may receive OBS basis set name
  const char* obs_name_cstr = opt.retrieve("obs");
  const std::string obs_name(obs_name_cstr ? obs_name_cstr : "");
  // may receive CABS basis set name
  const char* cabs_name_cstr = opt.retrieve("cabs");
  std::string cabs_name(cabs_name_cstr ? cabs_name_cstr : "");
  // if OBS given but CABS basis is not, look up a default value
  if (cabs_name.empty() && not obs_name.empty()) {
    cabs_name = R12Technology::default_cabs_name(obs_name);
  }
  // may receive DFBS basis set name
  const char* dfbs_name_cstr = opt.retrieve("dfbs");
  std::string dfbs_name(dfbs_name_cstr ? dfbs_name_cstr : "");
  // if OBS given but DFBS is not, look up a default DFBS
  if (dfbs_name.empty() && not obs_name.empty())
    dfbs_name = DensityFittingRuntime::default_dfbs_name(obs_name, 1); // for OBS with cardinal number X use DFBS with cardinal number X+1
  if (dfbs_name == "none") dfbs_name = "";
  // may receive F12 exponent
  const char* f12exp_cstr = opt.retrieve("f12exp");
  std::string f12exp_str(f12exp_cstr);
  // if OBS given but F12 exponent is not, look up a default value
  if (f12exp_str.empty() && not obs_name.empty()) {
    const double f12exp_default = R12Technology::default_stg_exponent(obs_name);
    if (f12exp_default != 0.0) {
      std::ostringstream oss;
      oss << f12exp_default;
      f12exp_str = oss.str();
    }
  }

  const char* r12_cstr = opt.retrieve("r12");
  const std::string r12_str = r12_cstr?r12_cstr:"";

#if defined(HAVE_MPQC3_RUNTIME)
  const char* singles_cstr = opt.retrieve("singles");
  const std::string singles_str = singles_cstr?singles_cstr:"";
  const char* partition_cstr = opt.retrieve("partitionH");
  const std::string partition_str = partition_cstr?partition_cstr:"";
  const char* mpqc3_cstr = opt.retrieve("mpqc3");
  const std::string mpqc3_str = mpqc3_cstr?mpqc3_cstr:"";
#endif

  ExEnv::out0() << indent << "Given resources: " << resources->sprint() << endl
      << endl;
  if (opt.retrieve("verbose")) {
    ExEnv::out0() << indent << "Using " << grp->class_name()
        << " for message passing (number of nodes = " << grp->n() << ")."
        << endl << indent << "Using " << thread->class_name()
        << " for threading (number of threads = " << thread->nthread() << ")."
        << endl << indent << "Using " << memory->class_name()
        << " for distributed shared memory." << endl << indent
        << "Total number of processors = " << grp->n() * thread->nthread()
        << endl;
    ExEnv::out0() << indent << "Using " << integral->class_name()
        << " for integrals by default" << std::endl;
  }

  //
  // Read molecule, basis, and orbitals
  //
  Ref<ExternMOInfo> rdorbs = new ExternMOInfo(filename_prefix + ".pt2r12.dat",
                                              integral,
                                              obs_name); // all MO info is contained in rdorbs
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
  const unsigned int nuocc = nmo - nfzc - ninact - nact - nfzv;

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
#if 0
  sc::ExEnv::out0() << "debug: print occ and orbs" << std::endl;
  occ_orbs->print_detail();
  orbs->print_detail();
#endif

//   // Currently we support two choices for reporting 2-RDM:
//   // in active space only (MOLCAS) or in the entire occupied space (GAMESS)
//   Ref<ExternSpinFreeRDMTwo> rdrdm2;
// #if PT2R12GAMESS
//   // GAMESS reports the density in occupied orbitals
//   rdrdm2 = new ExternSpinFreeRDMTwo(filename_prefix + ".pt2r12.rdm2.dat",
//                                     rdorbs->occindexmap_occ(),
//                                     occ_orbs);
// #else
//   // MOLCAS reports the density in active orbitals
//   rdrdm2 = new ExternSpinFreeRDMTwo(filename_prefix + ".pt2r12.rdm2.dat",
//                                     rdorbs->actindexmap_occ(),
//                                     occ_orbs);
// #endif

// the 2-RDM is reported in the active space only (for both MOLCAS and GAMESS))
  Ref<ExternSpinFreeRDMTwo> rdrdm2;
  rdrdm2 = new ExternSpinFreeRDMTwo(filename_prefix + ".pt2r12.rdm2.dat",
                                    rdorbs->actindexmap_occ(),
                                    occ_orbs);


  // create World in which we will compute
  // use defaults for all params
  Ref<WavefunctionWorld> world;
  {
    Ref<AssignedKeyVal> kva = new AssignedKeyVal;
    if (dfbs_name.empty() == false) {
      Ref<AssignedKeyVal> tmpkv = new AssignedKeyVal;
      tmpkv->assign("name", dfbs_name.c_str());
      if (dfbs_name.find("aug-cc-pV") != std::string::npos &&
          dfbs_name.find("Z-RI") != std::string::npos) { // if aug-cc-pVXZ-RI, make one as a union of
                                                         // cc-pVXZ-RI and augmentation-cc-pVXZ-RI
        std::string ccpvxzri_name(dfbs_name, 4, dfbs_name.size()-4);

        Ref<AssignedKeyVal> tmpkv1 = new AssignedKeyVal;
        tmpkv1->assign("name", ccpvxzri_name);
        tmpkv1->assign("molecule", basis->molecule().pointer());

        Ref<GaussianBasisSet> ccpvxzri = new GaussianBasisSet(tmpkv1);
        //Ref<GaussianBasisSet> ccpvxzri = new UncontractedBasisSet(tmpkv1);

        Ref<AssignedKeyVal> tmpkv2 = new AssignedKeyVal;
        tmpkv2->assign("name", std::string("augmentation-") + ccpvxzri_name);
        tmpkv2->assign("molecule", basis->molecule().pointer());
        Ref<GaussianBasisSet> augmentationccpvxzri = new GaussianBasisSet(tmpkv2);

        Ref<GaussianBasisSet> df_basis = new UnionBasisSet(ccpvxzri, augmentationccpvxzri);
        kva->assign("df_basis", df_basis.pointer());
      }
      else { // otherwise assume the basis exists in the library
        tmpkv->assign("molecule", basis->molecule().pointer());
        Ref<KeyVal> kv = tmpkv;
        Ref<GaussianBasisSet> df_basis = new GaussianBasisSet(kv);
        kva->assign("df_basis", df_basis.pointer());
      }
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
    if(not r12_str.empty())
      kva->assign("pt2_correction", r12_str);

#if defined(HAVE_MPQC3_RUNTIME)
    if(not singles_str.empty())
      kva->assign("cabs_singles", singles_str);
    if(not partition_str.empty())
      kva->assign("cabs_singles_h0", partition_str);
    if(not mpqc3_str.empty())
      kva->assign("use_mpqc3", mpqc3_str);
#endif
    Ref<KeyVal> kv = kva;
    extern_pt2r12 = new ExternPT2R12(kv);
  }

  extern_pt2r12->compute();

  extern_pt2r12->print();

  timer.print(ExEnv::out0());

#if defined(HAVE_TIME) && defined(HAVE_CTIME)
  time(&t);
  tstr = ctime(&t);
#endif
  if (!tstr) {
    tstr = "UNKNOWN";
  }
  ExEnv::out0() << std::endl
                << indent << scprintf("End Time: %s", tstr) << std::endl;

#ifdef HAVE_MADNESS
  MADNESSRuntime::finalize();
#endif

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
