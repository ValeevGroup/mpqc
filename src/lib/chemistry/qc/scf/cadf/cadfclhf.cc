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
#include <chemistry/qc/scf/iter_logger.h>

#include "cadfclhf.h"
#include "approx_pairs.h"

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
    L_3(gbs_->nshell() * gbs_->nshell() / scf_grp_->n()),
    L_B(gbs_->nshell() * gbs_->nshell() / scf_grp_->n()),
    L_3_star(gbs_->nshell() * gbs_->nshell() / scf_grp_->n()),
    L_d_over(gbs_->nshell() * gbs_->nshell() / scf_grp_->n()),
    L_d_under_ranges(gbs_->nshell() * gbs_->nshell() / scf_grp_->n()),
    L_DC(gbs_->nshell()),
    L_C_under(gbs_->nshell()),
    L_schwarz(gbs_->nshell()),
    decomps_(make_shared<DecompositionCache>(
        molecule()->natom() * molecule()->natom()
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
  pair_writer_ << keyval->describedclassvalue("pair_writer", KeyValValueRefDescribedClass(0));
  if(pair_writer_.nonnull()) {
    pair_writer_->wfn_ = this;
  }
  //----------------------------------------------------------------------------//
  nthread_ = keyval->intvalue("nthread", KeyValValueint(threadgrp_->nthread()));
  //----------------------------------------------------------------------------//
  screening_thresh_ = keyval->doublevalue("screening_thresh", KeyValValuedouble(screening_thresh_));
  pair_screening_thresh_ = keyval->doublevalue("pair_screening_thresh", KeyValValuedouble(screening_thresh_));
  full_screening_thresh_ = keyval->doublevalue("full_screening_thresh", KeyValValuedouble(screening_thresh_));
  distance_screening_thresh_ = keyval->doublevalue("distance_screening_thresh", KeyValValuedouble(full_screening_thresh_));
  B_screening_thresh_ = keyval->doublevalue("B_screening_thresh", KeyValValuedouble(full_screening_thresh_));
  d_over_screening_thresh_ = keyval->doublevalue("d_over_screening_thresh", KeyValValuedouble(B_screening_thresh_));
  d_under_screening_thresh_ = keyval->doublevalue("d_under_screening_thresh", KeyValValuedouble(B_screening_thresh_));
  // Unused?
  coef_screening_thresh_ = keyval->doublevalue("coef_screening_thresh", KeyValValuedouble(screening_thresh_));
  full_screening_thresh_min_ = keyval->doublevalue("full_screening_thresh_min", KeyValValuedouble(full_screening_thresh_*1e-5));
  // Obscure
  full_screening_expon_ = keyval->doublevalue("full_screening_expon", KeyValValuedouble(1.0));
  distance_damping_factor_ = keyval->doublevalue("distance_damping_factor", KeyValValuedouble(1.0));
  //----------------------------------------------------------------------------//
  // For now, use coulomb metric.  We can easily make this a keyword later
  metric_oper_type_ = TwoBodyOper::eri;
  //----------------------------------------------------------------------------//
  do_linK_ = keyval->booleanvalue("do_linK", KeyValValueboolean(do_linK_));
  screen_B_ = keyval->booleanvalue("screen_B", KeyValValueboolean(screen_B_));
  safe_extents_ = keyval->booleanvalue("safe_extents", KeyValValueboolean(safe_extents_));
  linK_block_rho_ = keyval->booleanvalue("linK_block_rho", KeyValValueboolean(screen_B_));
  B_use_buffer_ = keyval->booleanvalue("B_use_buffer", KeyValValueboolean(B_use_buffer_));
  linK_use_distance_ = keyval->booleanvalue("linK_use_distance", KeyValValueboolean(false));
  use_extents_ = keyval->booleanvalue("use_extents", KeyValValueboolean(use_extents_));
  use_max_extents_ = keyval->booleanvalue("use_max_extents", KeyValValueboolean(use_max_extents_));
  subtract_extents_ = keyval->booleanvalue("subtract_extents", KeyValValueboolean(subtract_extents_));
  well_separated_thresh_ = keyval->doublevalue("well_separated_thresh", KeyValValuedouble(well_separated_thresh_));
  //----------------------------------------------------------------------------//
  print_screening_stats_ = keyval->intvalue("print_screening_stats", KeyValValueint(0));
  min_atoms_per_node_ = keyval->intvalue("min_atoms_per_node", KeyValValueint(min_atoms_per_node_));
  //----------------------------------------------------------------------------//
  count_ints_only_ = keyval->booleanvalue("count_ints_only", KeyValValueboolean(count_ints_only_));
  count_ints_hist_clip_edges_ = keyval->booleanvalue("count_ints_hist_clip_edges", KeyValValueboolean(count_ints_hist_clip_edges_));
  count_ints_hist_min_ratio_ = keyval->doublevalue("count_ints_hist_min_ratio", KeyValValuedouble(count_ints_hist_min_ratio_));
  count_ints_hist_max_ratio_ = keyval->doublevalue("count_ints_hist_max_ratio", KeyValValuedouble(count_ints_hist_max_ratio_));
  count_ints_hist_max_int_ = keyval->doublevalue("count_ints_hist_max_int", KeyValValuedouble(count_ints_hist_max_int_));
  count_ints_hist_bins_ = keyval->intvalue("count_ints_hist_bins", KeyValValueint(count_ints_hist_bins_));
  count_ints_hist_ratio_bins_ = keyval->intvalue("count_ints_hist_ratio_bins", KeyValValueint(count_ints_hist_ratio_bins_));
  count_ints_hist_distance_bins_ = keyval->intvalue("count_ints_hist_distance_bins", KeyValValueint(count_ints_hist_distance_bins_));
  count_ints_histogram_ = keyval->booleanvalue("count_ints_histogram", KeyValValueboolean(count_ints_histogram_));
  count_ints_use_norms_ = keyval->booleanvalue("count_ints_use_norms", KeyValValueboolean(count_ints_use_norms_));
  count_ints_exclude_thresh_ = keyval->doublevalue("count_ints_exclude_thresh", KeyValValuedouble(count_ints_exclude_thresh_));
  //----------------------------------------------------------------------------//
  print_iteration_timings_ = keyval->booleanvalue("print_iteration_timings", KeyValValueboolean(print_iteration_timings_));
  //----------------------------------------------------------------------------//
  exact_diagonal_J_ = keyval->booleanvalue("exact_diagonal_J", KeyValValueboolean(exact_diagonal_J_));
  exact_diagonal_K_ = keyval->booleanvalue("exact_diagonal_K", KeyValValueboolean(exact_diagonal_K_));
  //----------------------------------------------------------------------------//
  debug_coulomb_energy_ = keyval->booleanvalue("debug_coulomb_energy", KeyValValueboolean(debug_coulomb_energy_));
  debug_exchange_energy_ = keyval->booleanvalue("debug_exchange_energy", KeyValValueboolean(debug_exchange_energy_));
  //----------------------------------------------------------------------------//
  use_norms_nu_ = keyval->booleanvalue("use_norms_nu", KeyValValueboolean(use_norms_nu_));
  use_norms_B_ = keyval->booleanvalue("use_norms_B", KeyValValueboolean(use_norms_B_));
  match_orbitals_ = keyval->booleanvalue("match_orbitals", KeyValValueboolean(match_orbitals_));
  match_orbitals_use_svd_ = keyval->booleanvalue("match_orbitals_use_svd", KeyValValueboolean(match_orbitals_use_svd_));
  match_orbitals_max_diff_ = keyval->doublevalue("match_orbitals_max_diff", KeyValValuedouble(match_orbitals_max_diff_));
  match_orbitals_max_homo_offset_ = keyval->doublevalue("match_orbitals_max_homo_offset", KeyValValuedouble(match_orbitals_max_homo_offset_));
  use_norms_sigma_ = keyval->booleanvalue("use_norms_sigma", KeyValValueboolean(use_norms_sigma_));
  sigma_norms_chunk_by_atoms_ = keyval->booleanvalue("sigma_norms_chunk_by_atoms", KeyValValueboolean(sigma_norms_chunk_by_atoms_));
  xml_screening_data_ = keyval->booleanvalue("xml_screening_data", KeyValValueboolean(xml_screening_data_));
  screen_B_use_distance_ = keyval->booleanvalue("screen_B_use_distance", KeyValValueboolean(screen_B_use_distance_));
  scale_screening_thresh_ = keyval->booleanvalue("scale_screening_thresh", KeyValValueboolean(scale_screening_thresh_));
  distribute_coefficients_ = keyval->booleanvalue("distribute_coefficients", KeyValValueboolean(distribute_coefficients_));
  store_coefs_transpose_ = keyval->booleanvalue("store_coefs_transpose", KeyValValueboolean(store_coefs_transpose_));
  if(distribute_coefficients_) store_coefs_transpose_ = false;
  shuffle_L_3_keys_ = keyval->booleanvalue("shuffle_L_3_keys", KeyValValueboolean(shuffle_L_3_keys_));
  shuffle_J_assignments_ = keyval->booleanvalue("shuffle_J_assignments", KeyValValueboolean(shuffle_J_assignments_));
  new_exchange_algorithm_ = keyval->booleanvalue("new_exchange_algorithm", KeyValValueboolean(new_exchange_algorithm_));
  //----------------------------------------------------------------------------//
  if(print_screening_stats_) {
    stats_ = std::make_shared<ScreeningStatistics>();
    stats_->print_level = print_screening_stats_;
    stats_->xml_stats_saved = xml_screening_data_;
  }
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
  if(new_exchange_algorithm_) {
    B_buffer_size_ =
      std::max(
          std::min((unsigned int)DEFAULT_TARGET_BLOCK_SIZE, gbs_->nbasis()) * max_fxn_dfbs_,
          (unsigned int)(max_fxn_atom_dfbs_ * 2 * max_fxn_obs_)
      ) * sizeof(double) * max_fxn_obs_;
  }
  else if(do_linK_) {
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
  if(pair_writer_.nonnull()) pair_writer_->wfn_ = 0;
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
  std::string ifmt = "%5d";
  auto double_str = [&fmt](double val) -> std::string {
     return std::string(scprintf(fmt.c_str(), val).str());
  };
  auto int_str = [&ifmt](int val) -> std::string {
     return std::string(scprintf(ifmt.c_str(), val).str());
  };
  o << indent << "B_screening_thresh = " << double_str(B_screening_thresh_) << endl;
  o << indent << "B_use_buffer = " << bool_str(B_use_buffer_) << endl;
  o << indent << "basis name = " << gbs_->label() << endl;
  o << indent << "coef_screening_thresh = " << double_str(coef_screening_thresh_) << endl;
  o << indent << "count_ints_histogram = " << bool_str(count_ints_histogram_) << endl;
  o << indent << "count_ints_only = " << bool_str(count_ints_only_) << endl;
  o << indent << "count_ints_hist_clip_edges = " << bool_str(count_ints_hist_clip_edges_) << endl;
  o << indent << "count_ints_hist_bins = " << int_str(count_ints_hist_bins_) << endl;
  o << indent << "count_ints_hist_distance_bins = " << int_str(count_ints_hist_distance_bins_) << endl;
  o << indent << "count_ints_hist_ratio_bins = " << int_str(count_ints_hist_ratio_bins_) << endl;
  o << indent << "count_ints_hist_min_ratio = " << double_str(count_ints_hist_min_ratio_) << endl;
  o << indent << "count_ints_hist_max_ratio = " << double_str(count_ints_hist_max_ratio_) << endl;
  o << indent << "count_ints_hist_max_int = " << double_str(count_ints_hist_max_int_) << endl;
  o << indent << "count_ints_use_norms = " << bool_str(count_ints_use_norms_) << endl;
  o << indent << "count_ints_exclude_thresh = " << double_str(count_ints_exclude_thresh_) << endl;
  o << indent << "d_over_screening_thresh = " << double_str(d_over_screening_thresh_) << endl;
  o << indent << "d_under_screening_thresh = " << double_str(d_under_screening_thresh_) << endl;
  o << indent << "debug_coulomb_energy = " << bool_str(debug_coulomb_energy_) << endl;
  o << indent << "debug_exchange_energy = " << bool_str(debug_exchange_energy_) << endl;
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
  o << indent << "linK_use_distance = " << bool_str(linK_use_distance_) << endl;
  o << indent << "match_orbitals = " << bool_str(match_orbitals_) << endl;
  o << indent << "match_orbitals_max_diff = " << double_str(match_orbitals_max_diff_) << endl;
  o << indent << "match_orbitals_max_homo_offset = " << double_str(match_orbitals_max_homo_offset_) << endl;
  o << indent << "match_orbitals_use_svd = " << bool_str(match_orbitals_use_svd_) << endl;
  o << indent << "min_atoms_per_node = " << int_str(min_atoms_per_node_) << endl;
  o << indent << "pair_screening_thresh = " << double_str(pair_screening_thresh_) << endl;
  o << indent << "safe_extents = " << bool_str(safe_extents_) << endl;
  o << indent << "scale_screening_thresh = " << bool_str(scale_screening_thresh_) << endl;
  o << indent << "screen_B = " << bool_str(screen_B_) << endl;
  o << indent << "screen_B_use_distance = " << bool_str(screen_B_use_distance_) << endl;
  o << indent << "screening_thresh = " << double_str(screening_thresh_) << endl;
  o << indent << "shuffle_J_assignments = " << bool_str(shuffle_J_assignments_) << endl;
  o << indent << "shuffle_L_3_keys = " << bool_str(shuffle_L_3_keys_) << endl;
  o << indent << "sigma_norms_chunk_by_atoms = " << bool_str(sigma_norms_chunk_by_atoms_) << endl;
  o << indent << "store_coefs_transpose = " << bool_str(store_coefs_transpose_) << endl;
  o << indent << "subtract_extents = " << bool_str(subtract_extents_) << endl;
  o << indent << "thread_4c_ints = " << bool_str(thread_4c_ints_) << endl;
  o << indent << "use_extents = " << bool_str(use_extents_) << endl;
  o << indent << "use_max_extents = " << bool_str(use_max_extents_) << endl;
  o << indent << "new_exchange_algorithm = " << bool_str(new_exchange_algorithm_) << endl;
  o << indent << "use_norms_nu = " << bool_str(use_norms_nu_) << endl;
  o << indent << "use_norms_B = " << bool_str(use_norms_B_) << endl;
  o << indent << "use_norms_sigma = " << bool_str(use_norms_sigma_) << endl;
  o << indent << "well_separated_thresh = " << double_str(well_separated_thresh_) << endl;
  o << indent << "xml_screening_data = " << bool_str(xml_screening_data_) << endl;


  CLHF::print(o);
  if(print_screening_stats_) {
    // TODO global sum stats at some point.  These numbers will be wrong otherwise
    stats_->global_sum(scf_grp_);
    stats_->print_summary(o, gbs_, dfbs_, print_screening_stats_);
  }
  o << decindent;

  if(xml_screening_data_ or count_ints_only_) {
    begin_xml_context("cadfclhf_screening", "screening_stats.xml");

    write_as_xml("statistics", *stats_,
      std::map<std::string, const std::string>({
        {"basis_name", gbs_->label()},
        {"dfbasis_name", dfbs_->label()},
        {"well_separated_thresh", std::string(scprintf("%3.1e", well_separated_thresh_).str())},
        {"distance_damping_factor", std::string(scprintf("%3.1e", distance_damping_factor_).str())},
        {"use_extents", (use_extents_ ? "true" : "false")},
        {"use_max_extents", (use_max_extents_ ? "true" : "false")},
        {"subtract_extents", (subtract_extents_ ? "true" : "false")},
        {"histogram_mode", (stats_->histogram_mode ? "true" : "false")}
      })
    );
    end_xml_context("cadfclhf_screening");
  }

}

//////////////////////////////////////////////////////////////////////////////////

void
CADFCLHF::done_threads(){
  //for (int i=0; i < threadgrp_->nthread(); i++) tbis_[i] = 0;
  //delete[] tbis_;
  //tbis_ = 0;
}

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

  /*-----------------------------------------------------*/
  /* Count ints only mode                           {{{2 */ #if 2 // begin fold
  if(count_ints_only_) {
    iter_stats_ = &(stats_->next_iteration());
    iter_stats_->set_nthread(nthread_);
    if(count_ints_histogram_) stats_->histogram_mode = true;
    double max_distance = 0.0;
    if(count_ints_histogram_) {
      for(auto&& ish : shell_range(gbs_)) {
        for(auto&& jsh : shell_range(gbs_)) {
          for(auto&& Xsh : shell_range(dfbs_)) {
            const double r = get_R(ish, jsh, Xsh, true);
            if(r > max_distance) max_distance = r;
          }
        }
      }
    }
    double min_pair_exponent = std::numeric_limits<double>::infinity();
    double max_pair_exponent = 0;
    for(auto exp : effective_df_exponents_) {
      if(exp > max_pair_exponent) max_pair_exponent = exp;
      if(exp < min_pair_exponent) min_pair_exponent = exp;
    }
    const int max_am = dfbs_->max_angular_momentum();
    for(int l = 0; l <= max_am; ++l) {
      iter_stats_->int_am_counts.push_back(0.0);
      iter_stats_->int_am_ratio_sums.push_back(0.0);
      iter_stats_->int_am_ratio_log_sums.push_back(0.0);
    }
    for(int am = 0; am < dfbs_->max_angular_momentum()+1; ++am) {
      iter_stats_->distance_hists.emplace_back(
          count_ints_hist_distance_bins_, 0.0, max_distance,
          count_ints_hist_ratio_bins_, count_ints_hist_min_ratio_, count_ints_hist_max_ratio_,
          false, true, count_ints_hist_clip_edges_
      );
      iter_stats_->distance_noschwarz_hists.emplace_back(
          count_ints_hist_distance_bins_, 0.0, max_distance,
          count_ints_hist_ratio_bins_, count_ints_hist_min_ratio_, count_ints_hist_max_ratio_,
          false, true, count_ints_hist_clip_edges_
      );
      iter_stats_->values_hists.emplace_back(
          count_ints_hist_bins_, count_ints_exclude_thresh_, count_ints_hist_max_int_,
          count_ints_hist_bins_, count_ints_exclude_thresh_, count_ints_hist_max_int_,
          true, true, count_ints_hist_clip_edges_
      );
      iter_stats_->exponent_ratio_hists.emplace_back(
          //count_ints_hist_bins_, min_max_pair_exponent.first->second, min_max_pair_exponent.second->second,
          count_ints_hist_bins_, min_pair_exponent, max_pair_exponent,
          count_ints_hist_ratio_bins_, count_ints_hist_min_ratio_, count_ints_hist_max_ratio_,
          true, true, count_ints_hist_clip_edges_
      );
    }
    count_ints();
    // Just give it something to keep it happy
    cl_fock_.result_noupdate().assign(hcore_);
    return;
  }
  /*******************************************************/ #endif //2}}}
  /*-----------------------------------------------------*/

  if(xml_debug_){
    begin_xml_context("compute_fock", "compute_fock.xml");
  }

  if(not have_coefficients_) {

    /*-----------------------------------------------------*/
    /* Precompute the coefficients                    {{{2 */ #if 2 // begin fold

    ints_computed_locally_ = 0;
    compute_coefficients();
    decltype(ints_computed_locally_.load()) ints_computed = ints_computed_locally_;
    scf_grp_->sum(&ints_computed, 1);
    if(scf_grp_->me() == 0) {
      ExEnv::out0() << indent << "Computed " << ints_computed << " integrals to determine coefficients." << endl;
    }
    if(pair_writer_.nonnull()) {
      if(distribute_coefficients_) {
        throw FeatureNotImplemented("Can't write pairs with distributed coefficients", __FILE__, __LINE__, class_desc());
      }
      pair_writer_->write_pairs();
    }
    if(iter_log_.nonnull()) {
      iter_log_->log_global_misc([](ptree& parent, const XMLWriter& writer, RefSymmSCMatrix hcore) {
        ptree& child = writer.insert_child(parent, hcore, "hcore");
      }, hcore_);
    }

    /*******************************************************/ #endif //2}}}
    /*-----------------------------------------------------*/

  }

  //---------------------------------------------------------------------------------------//

  timer.enter("misc");
  if(print_screening_stats_) {
    iter_stats_ = &(stats_->next_iteration());
    if(xml_screening_data_) {
      iter_stats_->set_nthread(nthread_);
    }
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


  /*=======================================================================================*/
  /* Form G                                                		                        {{{1 */ #if 1 // begin fold
  //---------------------------------------------------------------------------------------//

  // compute J and K
  timer.change("build");
  RefSCMatrix G;
  {
    ints_computed_locally_ = 0;
    if(xml_debug_) begin_xml_context("compute_J");

    // Do the actual computation of J
    RefSCMatrix J = compute_J();

    // Log stuff for debugging and profiling purposes
    if(xml_debug_) write_as_xml("J", J), end_xml_context("compute_J");
    if(iter_log_.nonnull()) {
      iter_log_->log_iter_misc([J](ptree& parent, const XMLWriter& writer) {
        ptree& child = writer.insert_child(parent, J, "J");
        child.put("<xmlattr>.name", "Coulomb matrix");
      });
    }
    decltype(ints_computed_locally_.load()) ints_computed = ints_computed_locally_;
    scf_grp_->sum(&ints_computed, 1);
    if(scf_grp_->me() == 0) {
      ExEnv::out0() << indent << "Computed " << ints_computed << " integrals for J part" << endl;
    }

    // Copy J into G
    G = J.copy();

  }
  {
    ints_computed_locally_ = 0;
    if(xml_debug_) begin_xml_context("compute_K");

    // Do the actual computation of K
    RefSCMatrix K = compute_K();

    // Log stuff for debugging and profiling
    if(xml_debug_) write_as_xml("K", K), end_xml_context("compute_K");
    if(iter_log_.nonnull()) {
      iter_log_->log_iter_misc([K](ptree& parent, const XMLWriter& writer) {
        ptree& child = writer.insert_child(parent, K, "K");
        child.put("<xmlattr>.name", "exchange matrix");
      });
    }
    decltype(ints_computed_locally_.load()) ints_computed = ints_computed_locally_;
    scf_grp_->sum(&ints_computed, 1);
    if(scf_grp_->me() == 0) {
      ExEnv::out0() << indent << "Computed " << ints_computed << " integrals for K part" << endl;
    }

    // Accumulate K into G
    G.accumulate( -1.0 * K);

  }

  if(xml_debug_) end_xml_context("compute_fock"), assert(false);

  // Reset the density_reset_ flag; it will get set to true if the density is reset between now
  //   and the next fock build
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
  if(iter_log_.nonnull()) {
    iter_log_->log_iter_misc([](ptree& parent, const XMLWriter& writer, RefSymmSCMatrix F) {
      ptree& child = writer.insert_child(parent, F, "F");
      child.put("<xmlattr>.name", "ao fock matrix");
    }, cl_fock_.result_noupdate());
  }
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

bool
CADFCLHF::get_shell_pair(ShellData& mu, ShellData& nu, PairSet pset)
{
  // Atomicly access and increment
  int spot = local_pairs_spot_++;
  const auto& plist = pset == SignificantPairs ? local_pairs_sig_ : local_pairs_all_;
  if(spot < plist.size()) {
    const auto& next_pair = plist[spot];
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

void
sc::get_split_range_part(
    const sc::Ref<sc::MessageGrp>& msg,
    int full_begin, int full_end,
    int& out_begin, int& out_size
)
{
  sc::get_split_range_part(msg->me(), msg->n(), full_begin, full_end, out_begin, out_size);
}

void
sc::get_split_range_part(
    int me, int n,
    int full_begin, int full_end,
    int& out_begin, int& out_size
)
{
  const int size = full_end - full_begin;
  int n_per = size / n;
  const int remain = size % n;
  if(remain == 0) {
    out_begin = full_begin + n_per * me;
    out_size = n_per;
  }
  else if(size < n) {
    if(me < size) {
      out_begin = full_begin + me;
      out_size = 1;
    }
    else {
      out_begin = full_begin + size - 1;
      if(out_begin < full_begin) out_begin = full_begin;
      out_size = 0;
    }
  }
  else {
    if(me < remain) {
      out_begin = full_begin + me * (n_per + 1);
      out_size = n_per + 1;
    }
    else {
      out_begin = full_begin + remain * (n_per + 1) + (me - remain) * n_per;
      out_size = n_per;
    }
  }
}

