//
// cadfclhf.cc
//
// Copyright (C) 2013 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: Nov 13, 2013
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

// Standard library includes
#include <memory>
#include <math.h>

// Boost includes
#include <boost/tuple/tuple_io.hpp>
#include <boost/math/special_functions/erf.hpp>

// MPQC includes
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/basis/gaussbas.h>
#include <util/group/messmpi.h>
#include <util/misc/regtime.h>
#include <util/misc/scexception.h>
#include <util/misc/xmlwriter.h>
#include <math/scmat/blas.h>
#include <util/container/conc_cache.h>

#include "cadfclhf.h"

#define EIGEN_NO_AUTOMATIC_RESIZING 1

#define CADF_USE_BLAS 0


using namespace sc;
using namespace std;
using boost::thread;
using boost::thread_group;
using namespace sc::parameter;

typedef std::pair<int, int> IntPair;
typedef CADFCLHF::CoefContainer CoefContainer;
typedef CADFCLHF::Decomposition Decomposition;
typedef std::pair<CoefContainer, CoefContainer> CoefPair;

static boost::mutex debug_print_mutex;

////////////////////////////////////////////////////////////////////////////////

ClassDesc CADFCLHF::cd_(
    typeid(CADFCLHF),
    "CADFCLHF",
    1, // verion number
    "public CLHF",
    0, // default constructor pointer
    create<CADFCLHF>, // KeyVal constructor pointer
    create<CADFCLHF> // StateIn constructor pointer
);

//////////////////////////////////////////////////////////////////////////////////

CADFCLHF::CADFCLHF(StateIn& s) :
    SavableState(s),
    CLHF(s),
    stats_(),
    iter_stats_(0)
{
  //----------------------------------------------------------------------------//
  // Nothing to do yet
  throw FeatureNotImplemented("SavableState construction of CADFCLHF",
      __FILE__, __LINE__, class_desc());
  //----------------------------------------------------------------------------//
}

//////////////////////////////////////////////////////////////////////////////////

CADFCLHF::CADFCLHF(const Ref<KeyVal>& keyval) :
    CLHF(keyval),
#if USE_INTEGRAL_CACHE
    // TODO Initialize the number of bins to a more reasonable size, check for enough memory, etc
    ints3_(make_shared<ThreeCenterIntCache>(
        gbs_->nshell() * gbs_->nshell() * gbs_->nshell() / scf_grp_->n()
    )),
    // TODO Initialize the number of bins to a more reasonable size, check for enough memory, etc
    ints2_(make_shared<TwoCenterIntCache>(
        gbs_->nshell() * gbs_->nshell() / scf_grp_->n()
    )),
#endif
    L_3(gbs_->nshell() * gbs_->nshell() * 3),
    L_DC(gbs_->nshell()),
    decomps_(make_shared<DecompositionCache>(
        molecule()->natom() * molecule()->natom()
    )),
    ints4maxes_(make_shared<FourCenterMaxIntCache>(
        gbs_->nshell() * gbs_->nshell()
    )),
    local_pairs_spot_(0),
    memory_used_(0)
{
  //----------------------------------------------------------------------------//
  // Get the auxiliary basis set
  dfbs_ << keyval->describedclassvalue("df_basis", KeyValValueRefDescribedClass(0));
  if(dfbs_.null()){
    throw InputError("CADFCLHF requires a density fitting basis set",
        __FILE__, __LINE__, "df_basis");
  }
  //----------------------------------------------------------------------------//
  nthread_ = keyval->intvalue("nthread", KeyValValueint(threadgrp_->nthread()));
  //----------------------------------------------------------------------------//
  // get the bound for filtering pairs.  This should be smaller than the general
  //   Schwarz bound.
  // TODO Figure out a reasonable default value for this
  pair_screening_thresh_ = keyval->doublevalue("pair_screening_thresh", KeyValValuedouble(1e-8));
  density_screening_thresh_ = keyval->doublevalue("density_screening_thresh", KeyValValuedouble(1e-8));
  full_screening_thresh_ = keyval->doublevalue("full_screening_thresh", KeyValValuedouble(1e-8));
  distance_screening_thresh_ = keyval->doublevalue("distance_screening_thresh", KeyValValuedouble(full_screening_thresh_));
  B_screening_thresh_ = keyval->doublevalue("B_screening_thresh", KeyValValuedouble(full_screening_thresh_));
  coef_screening_thresh_ = keyval->doublevalue("coef_screening_thresh", KeyValValuedouble(1e-8));
  full_screening_expon_ = keyval->doublevalue("full_screening_expon", KeyValValuedouble(1.0));
  full_screening_thresh_min_ = keyval->doublevalue("full_screening_thresh_min", KeyValValuedouble(full_screening_thresh_min_));
  distance_damping_factor_ = keyval->doublevalue("distance_damping_factor", KeyValValuedouble(1.0));
  //----------------------------------------------------------------------------//
  // For now, use coulomb metric.  We can easily make this a keyword later
  metric_oper_type_ = TwoBodyOper::eri;
  //----------------------------------------------------------------------------//
  do_linK_ = keyval->booleanvalue("do_linK", KeyValValueboolean(false));
  linK_block_rho_ = keyval->booleanvalue("linK_block_rho", KeyValValueboolean(false));
  B_use_buffer_ = keyval->booleanvalue("B_use_buffer", KeyValValueboolean(B_use_buffer_));
  linK_sorted_B_contraction_ = keyval->booleanvalue("linK_sorted_B_contraction", KeyValValueboolean(false));
  linK_use_distance_ = keyval->booleanvalue("linK_use_distance", KeyValValueboolean(false));
  use_extents_ = keyval->booleanvalue("use_extents", KeyValValueboolean(use_extents_));
  use_max_extents_ = keyval->booleanvalue("use_max_extents", KeyValValueboolean(use_max_extents_));
  subtract_extents_ = keyval->booleanvalue("subtract_extents", KeyValValueboolean(subtract_extents_));
  well_separated_thresh_ = keyval->doublevalue("well_separated_thresh", KeyValValuedouble(well_separated_thresh_));
  //----------------------------------------------------------------------------//
  print_screening_stats_ = keyval->intvalue("print_screening_stats", KeyValValueint(0));
  //----------------------------------------------------------------------------//
  print_iteration_timings_ = keyval->booleanvalue("print_iteration_timings", KeyValValueboolean(print_iteration_timings_));
  //----------------------------------------------------------------------------//
  exact_diagonal_J_ = keyval->booleanvalue("exact_diagonal_J", KeyValValueboolean(exact_diagonal_J_));
  exact_diagonal_K_ = keyval->booleanvalue("exact_diagonal_K", KeyValValueboolean(exact_diagonal_K_));
  //----------------------------------------------------------------------------//
  use_sparse_ = keyval->booleanvalue("use_sparse", KeyValValueboolean(use_sparse_));
  use_norms_nu_ = keyval->booleanvalue("use_norms_nu", KeyValValueboolean(use_norms_nu_));
  use_norms_sigma_ = keyval->booleanvalue("use_norms_sigma", KeyValValueboolean(use_norms_sigma_));
  xml_screening_data_ = keyval->booleanvalue("xml_screening_data", KeyValValueboolean(xml_screening_data_));
  all_to_all_L_3_ = keyval->booleanvalue("all_to_all_L_3", KeyValValueboolean(all_to_all_L_3_));
  sig_pairs_J_ = keyval->booleanvalue("sig_pairs_J", KeyValValueboolean(sig_pairs_J_));
  screen_B_ = keyval->booleanvalue("screen_B", KeyValValueboolean(screen_B_));
  screen_B_use_distance_ = keyval->booleanvalue("screen_B_use_distance", KeyValValueboolean(screen_B_use_distance_));
  scale_screening_thresh_ = keyval->booleanvalue("scale_screening_thresh", KeyValValueboolean(scale_screening_thresh_));
  distribute_coefficients_ = keyval->booleanvalue("distribute_coefficients", KeyValValueboolean(distribute_coefficients_));
  store_coefs_transpose_ = keyval->booleanvalue("store_coefs_transpose", KeyValValueboolean(store_coefs_transpose_));
  //----------------------------------------------------------------------------//
  stats_.print_level = print_screening_stats_;
  stats_.xml_stats_saved = xml_screening_data_;
  //----------------------------------------------------------------------------//
  xml_debug_ = keyval->booleanvalue("xml_debug", KeyValValueboolean(false));
  //----------------------------------------------------------------------------//
  for(auto&& ish : shell_range(gbs_)) {
    max_fxn_obs_ = std::max(ish.nbf, max_fxn_obs_);
    max_fxn_atom_obs_ = std::max(ish.atom_nbf, max_fxn_atom_obs_);
  }
  for(auto&& Xsh : shell_range(dfbs_)) {
    max_fxn_dfbs_ = std::max(Xsh.nbf, max_fxn_dfbs_);
    max_fxn_atom_dfbs_ = std::max(Xsh.atom_nbf, max_fxn_atom_dfbs_);
  }
  if(do_linK_) {
    B_buffer_size_ = std::min((unsigned int)DEFAULT_TARGET_BLOCK_SIZE, gbs_->nbasis())
      * sizeof(double) * max_fxn_obs_ * max_fxn_dfbs_;
  }
  else {
    B_buffer_size_ = std::min((unsigned int)DEFAULT_TARGET_BLOCK_SIZE, gbs_->nbasis())
      * sizeof(double) * max_fxn_obs_
      * std::min(DEFAULT_TARGET_BLOCK_SIZE, max_fxn_atom_dfbs_);
  }
  const decltype(B_buffer_size_) buff_tmp = B_buffer_size_;
  B_buffer_size_ = keyval->sizevalue("B_buffer_size",
      KeyValValuesize(B_buffer_size_)
  );
  if(B_buffer_size_ < buff_tmp){
    throw InputError("B_buffer_size is too small",
        __FILE__, __LINE__, "B_buffer_size");
  }
  initialize();
  //----------------------------------------------------------------------------//
}

//////////////////////////////////////////////////////////////////////////////////

CADFCLHF::~CADFCLHF()
{
  // Clean up the coefficient data
  if(coefficients_data_) delete[] coefficients_data_;
  if(dist_coefs_data_) delete[] dist_coefs_data_;
}


void
CADFCLHF::print(ostream&o) const
{
  o << indent << "Closed Shell Hartree-Fock (CLHF):" << endl;
  o << incindent;
  auto bool_str = [](bool val) -> std::string { return val ? "yes" : "no"; };
  std::string fmt = "%3.1e";
  auto double_str = [&fmt](double val) -> std::string {
     return std::string(scprintf(fmt.c_str(), val).str());
  };
  o << indent << "all_to_all_L_3 = " << bool_str(all_to_all_L_3_) << endl;
  o << indent << "B_screening_thresh = " << double_str(B_screening_thresh_) << endl;
  o << indent << "B_use_buffer = " << bool_str(B_use_buffer_) << endl;
  o << indent << "basis name = " << gbs_->label() << endl;
  o << indent << "coef_screening_thresh = " << double_str(coef_screening_thresh_) << endl;
  o << indent << "density_screening_thresh = " << double_str(density_screening_thresh_) << endl;
  o << indent << "dfbasis name = " << dfbs_->label() << endl;
  o << indent << "distance_damping_factor = " << double_str(distance_damping_factor_) << endl;
  o << indent << "distance_screening_thresh = " << double_str(distance_screening_thresh_) << endl;
  o << indent << "distribute_coefficients = " << bool_str(distribute_coefficients_) << endl;
  o << indent << "do_linK = " << bool_str(do_linK_) << endl;
  o << indent << "exact_diagonal_J = " << bool_str(exact_diagonal_J_) << endl;
  o << indent << "exact_diagonal_K = " << bool_str(exact_diagonal_K_) << endl;
  o << indent << "full_screening_expon = " << double_str(full_screening_expon_) << endl;
  o << indent << "full_screening_thresh = " << double_str(full_screening_thresh_) << endl;
  o << indent << "full_screening_thresh_min = " << double_str(full_screening_thresh_min_) << endl;
  o << indent << "linK_block_rho = " << bool_str(linK_block_rho_) << endl;
  o << indent << "linK_sorted_B_contraction = " << bool_str(linK_sorted_B_contraction_) << endl;
  o << indent << "linK_use_distance = " << bool_str(linK_use_distance_) << endl;
  o << indent << "pair_screening_thresh = " << double_str(pair_screening_thresh_) << endl;
  o << indent << "scale_screening_thresh = " << bool_str(scale_screening_thresh_) << endl;
  o << indent << "screen_B = " << bool_str(screen_B_) << endl;
  o << indent << "screen_B_use_distance = " << bool_str(screen_B_use_distance_) << endl;
  o << indent << "store_coefs_transpose = " << bool_str(store_coefs_transpose_) << endl;
  o << indent << "subtract_extents = " << bool_str(subtract_extents_) << endl;
  o << indent << "use_extents = " << bool_str(use_extents_) << endl;
  o << indent << "use_max_extents = " << bool_str(use_max_extents_) << endl;
  o << indent << "use_norms_nu = " << bool_str(use_norms_nu_) << endl;
  o << indent << "use_norms_sigma = " << bool_str(use_norms_sigma_) << endl;
  o << indent << "use_sparse = " << bool_str(use_sparse_) << endl;
  o << indent << "well_separated_thresh = " << double_str(well_separated_thresh_) << endl;
  o << indent << "xml_screening_data = " << bool_str(xml_screening_data_) << endl;


  CLHF::print(o);
  if(print_screening_stats_) {
    // TODO global sum stats at some point.  These numbers will be wrong otherwise
    stats_.global_sum(scf_grp_);
    stats_.print_summary(o, gbs_, dfbs_, print_screening_stats_);
  }
  o << decindent;

  if(xml_screening_data_) {
    begin_xml_context("cadfclhf_screening", "screening_stats.xml");

    write_as_xml("statistics", stats_,
      std::map<std::string, const std::string>({
        {"basis_name", gbs_->label()},
        {"dfbasis_name", dfbs_->label()},
        {"well_separated_thresh", std::string(scprintf("%3.1e", well_separated_thresh_).str())},
        {"distance_damping_factor", std::string(scprintf("%3.1e", distance_damping_factor_).str())},
        {"use_extents", (use_extents_ ? "true" : "false")},
        {"use_max_extents", (use_max_extents_ ? "true" : "false")},
        {"subtract_extents", (subtract_extents_ ? "true" : "false")}
      })
    );
    end_xml_context("cadfclhf_screening");
  }

}

//////////////////////////////////////////////////////////////////////////////////

void
CADFCLHF::done_threads(){
  CLHF::done_threads();
}

//////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////

void
CADFCLHF::save_data_state(StateOut& so)
{
  //----------------------------------------------------------------------------//
  // Nothing to do yet
  throw FeatureNotImplemented("SavableState writing in CADFCLHF",
      __FILE__, __LINE__, class_desc());
  //----------------------------------------------------------------------------//
}

//////////////////////////////////////////////////////////////////////////////////

void
CADFCLHF::ao_fock(double accuracy)
{
  /*=======================================================================================*/
  /* Setup                                                 		                        {{{1 */ #if 1 // begin fold
  //---------------------------------------------------------------------------------------//
  Timer timer("ao_fock");
  //---------------------------------------------------------------------------------------//
  if(xml_debug_){
    begin_xml_context("compute_fock", "compute_fock.xml");
  }
  if(not have_coefficients_) {
    ints_computed_locally_ = 0;
    compute_coefficients();
    ints_computed_ = ints_computed_locally_;
    scf_grp_->sum(&ints_computed_, 1);
    if(scf_grp_->me() == 0) {
      ExEnv::out0() << "  Computed " << ints_computed_ << " integrals to determine coefficients." << endl;
    }
  }
  //---------------------------------------------------------------------------------------//
  timer.enter("misc");
  iter_stats_ = &(stats_.next_iteration());
  if(xml_screening_data_) {
    iter_stats_->set_nthread(nthread_);
  }
  //----------------------------------------//
  // transform the density difference to the AO basis
  RefSymmSCMatrix dd = cl_dens_diff_;
  if(xml_debug_) write_as_xml("cl_dens_diff_", cl_dens_diff_);
  Ref<PetiteList> pl = integral()->petite_list();
  cl_dens_diff_ = pl->to_AO_basis(dd);
  //----------------------------------------//
  double gmat_accuracy = accuracy;
  if (min_orthog_res() < 1.0) { gmat_accuracy *= min_orthog_res(); }
  //----------------------------------------//
  // copy over the density
  D_ = cl_dens_diff_.copy().convert2RefSCMat();
  D_.scale(0.5);
  if(xml_debug_) write_as_xml("D", D_);
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Form G                                                		                        {{{1 */ #if 1 // begin fold
  //---------------------------------------------------------------------------------------//
  // compute J and K
  timer.change("build");
  RefSCMatrix G;
  {
    ints_computed_locally_ = 0;
    if(xml_debug_) begin_xml_context("compute_J");
    RefSCMatrix J = compute_J();
    if(xml_debug_) write_as_xml("J", J), end_xml_context("compute_J");
    G = J.copy();
    ints_computed_ = ints_computed_locally_;
    scf_grp_->sum(&ints_computed_, 1);
    if(scf_grp_->me() == 0) {
      ExEnv::out0() << "        Computed " << ints_computed_ << " integrals for J part" << endl;
    }
  }
  {
    ints_computed_locally_ = 0;
    if(xml_debug_) begin_xml_context("compute_K");
    RefSCMatrix K = compute_K();
    if(xml_debug_) write_as_xml("K", K), end_xml_context("compute_K");
    G.accumulate( -1.0 * K);
    ints_computed_ = ints_computed_locally_;
    scf_grp_->sum(&ints_computed_, 1);
    if(scf_grp_->me() == 0) {
      ExEnv::out0() << "        Computed " << ints_computed_ << " integrals for K part" << endl;
    }
  }
  if(xml_debug_) end_xml_context("compute_fock"), assert(false);
  density_reset_ = false;
  //---------------------------------------------------------------------------------------//
  // Move data back to a RefSymmSCMatrix, transform back to the SO basis
  Ref<SCElementOp> accum_G_op = new SCElementAccumulateSCMatrix(G.pointer());
  RefSymmSCMatrix G_symm = G.kit()->symmmatrix(G.coldim()); G_symm.assign(0.0);
  G_symm.element_op(accum_G_op); G = 0;
  G_symm = pl->to_SO_basis(G_symm);
  //---------------------------------------------------------------------------------------//
  // Accumulate difference back into gmat_
  gmat_.accumulate(G_symm); G_symm = 0;
  //---------------------------------------------------------------------------------------//
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Clean up                                             		                        {{{1 */ #if 1 // begin fold
  //---------------------------------------------------------------------------------------//
  timer.change("misc");
  //----------------------------------------//
  // restore the SO version of the density difference
  cl_dens_diff_ = dd;
  //----------------------------------------//
  if(print_screening_stats_) {
    iter_stats_->global_sum(scf_grp_);
  }
  //----------------------------------------//
  // F = H+G
  cl_fock_.result_noupdate().assign(hcore_);
  cl_fock_.result_noupdate().accumulate(gmat_);
  accumddh_->accum(cl_fock_.result_noupdate());
  cl_fock_.computed()=1;
  //---------------------------------------------------------------------------------------//
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
}

//////////////////////////////////////////////////////////////////////////////////

void
CADFCLHF::reset_density()
{
  CLHF::reset_density();
  gmat_.assign(0.0);
  density_reset_ = true;
}

//////////////////////////////////////////////////////////////////////////////////

void
CADFCLHF::loop_shell_pairs_threaded(
    PairSet pset,
    const std::function<void(int, const ShellData&, const ShellData&)>& f
)
{
  local_pairs_spot_ = 0;
  boost::thread_group compute_threads;
  // Loop over number of threads
  for(int ithr = 0; ithr < nthread_; ++ithr) {
    // ...and create each thread that computes pairs
    compute_threads.create_thread([&,ithr](){
      ShellData ish, jsh;
      //----------------------------------------//
      while(get_shell_pair(ish, jsh, pset)){
        f(ithr, ish, jsh);
      }

    });
  }
  compute_threads.join_all();
}

//////////////////////////////////////////////////////////////////////////////////

bool
CADFCLHF::get_shell_pair(ShellData& mu, ShellData& nu, PairSet pset)
{
  // Atomicly access and increment
  int spot = local_pairs_spot_++;
  const auto& plist = pset == SignificantPairs ? local_pairs_sig_ : local_pairs_all_;
  if(spot < plist.size()) {
    const auto& next_pair = plist[spot];
    //----------------------------------------//
    if(dynamic_) {
      // Here's where we'd need to check if we're running low on pairs and prefetch some more
      // When implemented, this should use a std::async or something like that
      throw FeatureNotImplemented("dynamic load balancing", __FILE__, __LINE__, class_desc());
    }
    //----------------------------------------//
    mu = ShellData(next_pair.first, gbs_.pointer(), dfbs_.pointer());
    nu = ShellData(next_pair.second, gbs_.pointer(), dfbs_.pointer());
  }
  else{
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////////


