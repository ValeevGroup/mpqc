
#include <cmath>
#include <cassert>

#include <util/misc/scexception.h>
#include <util/misc/consumableresources.h>
#include <util/misc/regtime.h>
#include <util/misc/formio.h>
#include <util/state/stateio.h>

#include <chemistry/qc/basis/petite.h>

#include <chemistry/qc/scf/fbclhf.h>
#include <chemistry/qc/lcao/fockbuild.h>
#include <chemistry/qc/lcao/fockbuild_runtime.h>
#include <chemistry/qc/lcao/clhfcontrib.h>

#ifdef MPQC_NEW_FEATURES
#  include <util/misc/xmlwriter.h>
static constexpr bool xml_debug = false;
#endif // MPQC_NEW_FEATURES

using namespace std;
using namespace sc;

static ClassDesc FockBuildCLHF_cd(
  typeid(FockBuildCLHF),"FockBuildCLHF",1,"public CLHF",
  0, create<FockBuildCLHF>, create<FockBuildCLHF>);

FockBuildCLHF::FockBuildCLHF(StateIn& s) :
  SavableState(s),
  CLHF(s)
{
  fockdist_ << SavableState::restore_state(s);
  s.get(fockbuildmatrixtype_);
  s.get(prefetch_blocks_);
}

FockBuildCLHF::FockBuildCLHF(const Ref<KeyVal>& keyval) :
  CLHF(keyval)
{
  fockdist_ << keyval->describedclassvalue("fockdist");
  if (fockdist_.null()) {
      fockdist_ = new FockDistribution;
    }
  KeyValValuestring deffbm("replicated");
  fockbuildmatrixtype_ = keyval->stringvalue("fockbuildmatrixtype",deffbm);
  if (fockbuildmatrixtype_ != "replicated"
      && fockbuildmatrixtype_ != "distributed"
      && fockbuildmatrixtype_ != "prefetched_distributed") {
      throw InputError("fockbuildmatrixtype must be \"replicated\","
                           "\"distributed\", or \"prefetched_distributed\".",
                           __FILE__,
                           __LINE__,
                           "fockbuildmatrixtype",
                           fockbuildmatrixtype_.c_str(),
                           class_desc());
    }
  if (fockbuildmatrixtype_ == "prefetched_distributed") {
      prefetch_blocks_ = true;
  }
  else {
      prefetch_blocks_ = false;
    }
}

FockBuildCLHF::~FockBuildCLHF()
{
}

void
FockBuildCLHF::save_data_state(StateOut& s)
{
  CLHF::save_data_state(s);
  SavableState::save_state(fockdist_.pointer(),s);
  s.put(fockbuildmatrixtype_);
  s.put(prefetch_blocks_);
}

void
FockBuildCLHF::init_threads()
{
  Ref<GaussianBasisSet> gbs(basis());
  Ref<FockContribution> fc
      = new CLHFContribution(gbs,gbs,gbs,fockbuildmatrixtype_);
  fb_ = new FockBuild(fockdist_,
                      fc,
                      prefetch_blocks_,
                      gbs, gbs, gbs,
                      scf_grp_, threadgrp_, integral());
}

void
FockBuildCLHF::done_threads()
{
  fb_ = 0;
}

void
FockBuildCLHF::ao_fock(double accuracy)
{
  Timer routine_tim("ao_fock");
  Timer step_tim("misc");
  int nthread = threadgrp_->nthread();

  Ref<GaussianBasisSet> gbs = basis();
  Ref<PetiteList> pl = integral()->petite_list(gbs);
  
  // transform the density difference to the AO basis
  RefSymmSCMatrix dd = cl_dens_diff_;
  cl_dens_diff_ = pl->to_AO_basis(dd);

  // These two lines are needed to get the same cl_dens_diff_
  // as the old CLHF class has.
  //cl_dens_diff_->scale(2.0);
  //cl_dens_diff_->scale_diagonal(0.5);

  double gmat_accuracy = accuracy;
  if (min_orthog_res() < 1.0) { gmat_accuracy *= min_orthog_res(); }

  if (debug_>1) {
    cl_dens_diff_.print("cl_dens_diff before set_pmat");
  }

  fb_->contrib()->set_fmat(0, cl_gmat_);
  fb_->contrib()->set_pmat(0, cl_dens_diff_);
  fb_->set_accuracy(gmat_accuracy);

  if (debug_>1) {
    cl_gmat_.print("cl_gmat before build");
    cl_dens_diff_.print("cl_dens_diff before build");
  }

  step_tim.change("build");
  fb_->build();

  ExEnv::out0() << indent << scprintf("%20.0f integrals\n",
                                                  fb_->contrib()->nint());

  step_tim.change("misc");

  // get rid of the AO basis density difference
  cl_dens_diff_ = dd;

  // now symmetrize the skeleton G matrix, placing the result in dd
  RefSymmSCMatrix skel_gmat = cl_gmat_.copy();
  skel_gmat.scale(1.0/(double)pl->order());
  if (debug_>1) {
    skel_gmat.print("skel_gmat before symmetrize");
  }
  dd = cl_dens_diff_.clone();
  pl->symmetrize(skel_gmat,dd);
  if (debug_>1) {
    dd.print("dd after symmetrize");
  }

  // F = H+G
  cl_fock_.result_noupdate().assign(hcore_);
  cl_fock_.result_noupdate().accumulate(dd);
  accumddh_->accum(cl_fock_.result_noupdate());
  cl_fock_.computed()=1;
}

void
FockBuildCLHF::print(std::ostream&o) const
{
  CLHF::print(o);
  ExEnv::out0() << indent
                << "fockbuildmatrixtype = " << fockbuildmatrixtype_
                << std::endl;
  fockdist_->print(o);
}

////////////////////////////////////////////////////////

ClassDesc DFCLHF::cd_(
  typeid(DFCLHF),"DFCLHF",1,"public CLHF",
  0, create<DFCLHF>, create<DFCLHF>);

DFCLHF::DFCLHF(StateIn& s) :
  SavableState(s),
  CLHF(s)
{
  world_ << SavableState::restore_state(s);
}

DFCLHF::DFCLHF(const Ref<KeyVal>& keyval) :
  CLHF(keyval)
{
  dens_reset_freq_ = 1;

  // if world not given, make this the center of a new World
  world_ << keyval->describedclassvalue("world", KeyValValueRefDescribedClass(0));
  if (world_.null())
    world_ = new WavefunctionWorld(keyval);
  if (world_.null())
    throw InputError("DFCLHF requires a WavefunctionWorld; input did not specify it, neither could it be constructed",
                     __FILE__, __LINE__, "world");
  if (world_->wfn() == 0) world_->set_wfn(this);

#ifdef MPQC_NEW_FEATURES
  xml_debug_ = keyval->booleanvalue("xml_debug", KeyValValueboolean(false));
#endif // MPQC_NEW_FEATURES

  // need a nonblocked cl_gmat_ in this method
  Ref<PetiteList> pl = integral()->petite_list();
  gmat_ = basis()->so_matrixkit()->symmmatrix(pl->SO_basisdim());
  gmat_.assign(0.0);
}

DFCLHF::~DFCLHF()
{
}

void
DFCLHF::save_data_state(StateOut& s)
{
  CLHF::save_data_state(s);
  SavableState::save_state(world_.pointer(),s);
}

void
DFCLHF::ao_fock(double accuracy)
{
  Timer routine_tim("ao_fock");
  Timer step_tim("misc");
  int nthread = threadgrp_->nthread();

#ifdef MPQC_NEW_FEATURES
  if(xml_debug_) begin_xml_context("compute_fock", "compute_fock.xml");
#endif // MPQC_NEW_FEATURES

  // transform the density difference to the AO basis
  RefSymmSCMatrix dd = cl_dens_diff_;
#ifdef MPQC_NEW_FEATURES
  if(xml_debug_) write_as_xml("cl_dens_diff_", cl_dens_diff_);
#endif // MPQC_NEW_FEATURES
  Ref<PetiteList> pl = integral()->petite_list();
  cl_dens_diff_ = pl->to_AO_basis(dd);

  double gmat_accuracy = accuracy;
  if (min_orthog_res() < 1.0) { gmat_accuracy *= min_orthog_res(); }

  if (debug_>1) {
    cl_dens_diff_.print("cl_dens_diff before set_pmat");
  }

  Ref<OrbitalSpaceRegistry> oreg = world_->moints_runtime()->factory()->orbital_registry();
  Ref<AOSpaceRegistry> aoreg = world_->moints_runtime()->factory()->ao_registry();
  std::string aospace_id = new_unique_key(oreg);
  if (aoreg->key_exists(basis()) == false) {
    Ref<OrbitalSpace> aospace = new AtomicOrbitalSpace(aospace_id, "DFCLHF AO basis set", basis(), integral());
    aoreg->add(basis(), aospace);
    MPQC_ASSERT(oreg->key_exists(aospace_id) == false); // should be ensured by using new_unique_key
    oreg->add(make_keyspace_pair(aospace));
  }
  // feed the spin densities to the builder, cl_dens_diff_ includes total density right now, so halve it
  RefSymmSCMatrix Pa = cl_dens_diff_.copy(); Pa.scale(0.5);
  RefSymmSCMatrix Pb = Pa;

  Ref<FockBuildRuntime> fb_rtime = world_->fockbuild_runtime();
  fb_rtime->set_densities(Pa, Pb);

  step_tim.change("build");
  Ref<OrbitalSpace> aospace = aoreg->value(basis());
  RefSCMatrix G;
#ifdef MPQC_NEW_FEATURES
  if(xml_debug_) write_as_xml("D", Pa);
#endif // MPQC_NEW_FEATURES
  {
    const std::string jkey = ParsedOneBodyIntKey::key(aospace->id(),aospace->id(),std::string("J"));
#ifdef MPQC_NEW_FEATURES
    if(xml_debug_) begin_xml_context("compute_J");
#endif // MPQC_NEW_FEATURES
    RefSCMatrix J = fb_rtime->get(jkey);
#ifdef MPQC_NEW_FEATURES
    if(xml_debug_) write_as_xml("J", J), end_xml_context("compute_J");
#endif // MPQC_NEW_FEATURES
    G = J.copy();
  }
  {
    const std::string kkey = ParsedOneBodyIntKey::key(aospace->id(),aospace->id(),std::string("K"),AnySpinCase1);
#ifdef MPQC_NEW_FEATURES
    if(xml_debug_) begin_xml_context("compute_K");
#endif // MPQC_NEW_FEATURES
    RefSCMatrix K = fb_rtime->get(kkey);
#ifdef MPQC_NEW_FEATURES
    if(xml_debug_) write_as_xml("K", K), end_xml_context("compute_K");
#endif // MPQC_NEW_FEATURES
    G.accumulate( -1.0 * K);
  }
#ifdef MPQC_NEW_FEATURES
  if(xml_debug_) end_xml_context("compute_fock"), assert(false);
#endif // MPQC_NEW_FEATURES
  Ref<SCElementOp> accum_G_op = new SCElementAccumulateSCMatrix(G.pointer());
  RefSymmSCMatrix G_symm = G.kit()->symmmatrix(G.coldim()); G_symm.assign(0.0);
  G_symm.element_op(accum_G_op); G = 0;
  G_symm = pl->to_SO_basis(G_symm);
  gmat_.accumulate(G_symm); G_symm = 0;
  step_tim.change("misc");

  // get rid of the AO basis density difference
  cl_dens_diff_ = dd;

  // F = H+G
  cl_fock_.result_noupdate().assign(hcore_);
  cl_fock_.result_noupdate().accumulate(gmat_);
  accumddh_->accum(cl_fock_.result_noupdate());
  cl_fock_.computed()=1;
  //this->reset_density();
}

void
DFCLHF::reset_density() {
  CLHF::reset_density();
  gmat_.assign(0.0);
}

void
DFCLHF::print(std::ostream&o) const
{
  CLHF::print(o);
  world_->print(o);
}

Ref<DensityFittingInfo>
DFCLHF::dfinfo() const {
  Ref<DensityFittingInfo> result = const_cast<DensityFittingInfo*>(world_->tfactory()->df_info());
  return result;
}


#ifdef MPQC_NEW_FEATURES
boost::property_tree::ptree&
DFCLHF::write_xml(
    boost::property_tree::ptree& parent,
    const XMLWriter& writer
)
{
  using boost::property_tree::ptree;
  ptree& child = get_my_ptree(parent);
  world()->write_xml(child, writer);
  return CLHF::write_xml(parent, writer);
}
#endif // MPQC_NEW_FEATURES
