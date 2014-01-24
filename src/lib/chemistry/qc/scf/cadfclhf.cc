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

#include <chemistry/qc/basis/petite.h>
#include <util/group/messmpi.h>
#include <util/misc/regtime.h>
#include <util/misc/scexception.h>
#include <chemistry/qc/scf/cadfclhf.h>
#include <util/misc/xmlwriter.h>
#include <boost/tuple/tuple_io.hpp>
#include <chemistry/qc/basis/gaussbas.h>
#include <memory>
#include <math.h>

#define USE_INTEGRAL_CACHE 1
#define EIGEN_NO_AUTOMATIC_RESIZING 1

#ifdef ECLIPSE_PARSER_ONLY
#include <util/misc/sharedptr.h>
#endif

using namespace sc;
using namespace std;
using boost::thread;
using boost::thread_group;
using namespace sc::parameter;
using std::make_shared;
static constexpr bool xml_debug = false;

typedef std::pair<int, int> IntPair;
typedef CADFCLHF::CoefContainer CoefContainer;
typedef CADFCLHF::Decomposition Decomposition;
typedef std::pair<CoefContainer, CoefContainer> CoefPair;

#if LINK_TIME_LIST_INSERTIONS
std::atomic_uint_fast64_t OrderedShellList::insertion_time_nanos = ATOMIC_VAR_INIT(0);
std::atomic_uint_fast64_t OrderedShellList::n_insertions = ATOMIC_VAR_INIT(0);
#endif

static boost::mutex debug_print_mutex;

////////////////////////////////////////////////////////////////////////////////
// Debugging asserts and outputs

#define M_DUMP(M) std::cout << #M << " is " << M.rows() << " x " << M.cols() << std::endl;

#define M_ROW_ASSERT(M1, M2) \
  if(M1.rows() != M2.rows()) { \
    boost::lock_guard<boost::mutex> _tmp(debug_print_mutex); \
    cout << "assertion failed.  Rows not equal:  M1 => " << M1.rows() << " x " << M1.cols() << ", M2 => " << M2.rows() << " x " << M2.cols() << endl; \
    assert(false); \
  }
#define M_COL_ASSERT(M1, M2) \
  if(M1.cols() != M2.cols()) { \
    boost::lock_guard<boost::mutex> _tmp(debug_print_mutex); \
    cout << "assertion failed. Cols not equal:  M1 => " << M1.rows() << " x " << M1.cols() << ", M2 => " << M2.rows() << " x " << M2.cols() << endl; \
    assert(false); \
  }

#define M_PROD_CHECK(R, M1, M2) \
  if(R.rows() != M1.rows() || R.cols() != M2.cols() || M1.cols() != M2.rows()) { \
    boost::lock_guard<boost::mutex> _tmp(debug_print_mutex); \
    cout << "can't perform multiplication: (" << R.rows() << " x " << R.cols() << ") = (" << M1.rows() << " x " << M1.cols() << ") * (" << M2.rows() << " x " << M2.cols() << ")" << endl; \
    assert(false); \
  }

#define M_DOT_CHECK(M1, M2) \
  if(1 != M1.rows() || 1 != M2.cols() || M1.cols() != M2.rows()) { \
    boost::lock_guard<boost::mutex> _tmp(debug_print_mutex); \
    cout << "can't perform multiplication: (" << 1 << " x " << 1 << ") = (" << M1.rows() << " x " << M1.cols() << ") * (" << M2.rows() << " x " << M2.cols() << ")" << endl; \
    assert(false); \
  }

#define DECOMP_PRINT(D, M) \
  {\
    boost::lock_guard<boost::mutex> _tmp(debug_print_mutex); \
    cout << "Decomposition:  " << D.matrixQR().rows() << " x " << D.matrixQR().cols() << ", M => " << M.rows() << " x " << M.cols() << endl; \
  }

#define M_EQ_ASSERT(M1, M2) M_ROW_ASSERT(M1, M2); M_COL_ASSERT(M1, M2);

#define M_BLOCK_ASSERT(M, b1a, b1b, b2a, b2b) \
  if(b1a < 0) { \
    boost::lock_guard<boost::mutex> _tmp(debug_print_mutex); \
    cout << "assertion 1 failed.  data: " << "(" << M.rows() << ", " << M.cols() << ", " << b1a << ", " << b1b << ", " << b2a << ", " << b2b << ")" << endl; \
  } \
  else if(b1b < 0) { \
    boost::lock_guard<boost::mutex> _tmp(debug_print_mutex); \
    cout << "assertion 2 failed.  data: " << "(" << M.rows() << ", " << M.cols() << ", " << b1a << ", " << b1b << ", " << b2a << ", " << b2b << ")" << endl; \
  } \
  else if(b1a > M.rows() - b2a) { \
    boost::lock_guard<boost::mutex> _tmp(debug_print_mutex); \
    cout << "assertion 3 failed.  data: " << "(" << M.rows() << ", " << M.cols() << ", " << b1a << ", " << b1b << ", " << b2a << ", " << b2b << ")" << endl; \
  } \
  else if(b1b > M.cols() - b2b) { \
    boost::lock_guard<boost::mutex> _tmp(debug_print_mutex); \
    cout << "assertion 4 failed.  data: " << "(" << M.rows() << ", " << M.cols() << ", " << b1a << ", " << b1b << ", " << b2a << ", " << b2b << ")" << endl; \
  }

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
    CLHF(s)
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
    // TODO Initialize the number of bins to a more reasonable size, check for enough memory, etc
    ints3_(basis()->nshell() * basis()->nshell() * basis()->nshell() / scf_grp_->n() / std::min((unsigned int)50, basis()->nbasis())),
    ints2_(basis()->nshell() * basis()->nshell() / scf_grp_->n() / std::min((unsigned int)50, basis()->nbasis())),
    local_pairs_spot_(0)
{
  //----------------------------------------------------------------------------//
  // Get the auxiliary basis set
  dfbs_ << keyval->describedclassvalue("df_basis", KeyValValueRefDescribedClass(0));
  if(dfbs_.null()){
    throw InputError("CADFCLHF requires a density fitting basis set",
        __FILE__, __LINE__, "df_basis");
  }
  //----------------------------------------------------------------------------//
  // get the bound for filtering pairs.  This should be smaller than the general
  //   Schwarz bound.
  // TODO Figure out a reasonable default value for this
  pair_screening_thresh_ = keyval->doublevalue("pair_screening_thresh", KeyValValuedouble(1e-8));
  density_screening_thresh_ = keyval->doublevalue("density_screening_thresh", KeyValValuedouble(1e-8));
  full_screening_thresh_ = keyval->doublevalue("full_screening_thresh", KeyValValuedouble(1e-8));
  coef_screening_thresh_ = keyval->doublevalue("coef_screening_thresh", KeyValValuedouble(1e-8));
  //----------------------------------------------------------------------------//
  // For now, use coulomb metric.  We can easily make this a keyword later
  metric_oper_type_ = TwoBodyOper::eri;
  //----------------------------------------------------------------------------//
  do_linK_ = keyval->booleanvalue("do_linK", KeyValValueboolean(false));
  //----------------------------------------------------------------------------//
  print_screening_stats_ = keyval->booleanvalue("print_screening_stats", KeyValValueboolean(false));
  //----------------------------------------------------------------------------//
  initialize();
  //----------------------------------------------------------------------------//
}

//////////////////////////////////////////////////////////////////////////////////

CADFCLHF::~CADFCLHF()
{
  if(have_coefficients_){
    // Clean up the coefficient data
    deallocate(coefficients_data_);
  }
}

//////////////////////////////////////////////////////////////////////////////////

void
CADFCLHF::initialize()
{
  //----------------------------------------------------------------------------//
  // Check that the density is local
  if(!local_dens_){
    throw FeatureNotImplemented("Can't handle density matrices that don't fit on one node",
        __FILE__, __LINE__, class_desc());
  }
  //----------------------------------------------------------------------------//
  // need a nonblocked cl_gmat_ in this method
  Ref<PetiteList> pl = integral()->petite_list();
  gmat_ = basis()->so_matrixkit()->symmmatrix(pl->SO_basisdim());
  gmat_.assign(0.0);
  //----------------------------------------------------------------------------//
  // Determine if the message group is an instance of MPIMessageGrp
  using_mpi_ = dynamic_cast<MPIMessageGrp*>(scf_grp_.pointer()) ? true : false;
  //----------------------------------------------------------------------------//
  have_coefficients_ = false;
}

//////////////////////////////////////////////////////////////////////////////////

void
CADFCLHF::init_threads()
{
  Timer timer("init threads");
  assert(not threads_initialized_);
  //----------------------------------------------------------------------------//
  // TODO fix this so that deallocating the tbis_ array doesn't cause a seg fault when this isn't called (we don't need it)
  SCF::init_threads();
  //----------------------------------------------------------------------------//
  const int nthread = threadgrp_->nthread();
  const int n_node = scf_grp_->n();
  //----------------------------------------------------------------------------//
  // initialize the two electron integral classes
  // ThreeCenter versions
  integral()->set_basis(basis(), basis(), dfbs_);
  eris_3c_.resize(nthread);
  do_threaded(nthread, [&](int ithr){
    eris_3c_[ithr] = (integral()->coulomb<3>());
  });
  for (int i=0; i < nthread; i++) {
    // TODO different fitting metrics
    if(metric_oper_type_ == coulomb_oper_type_){
      metric_ints_3c_.push_back(eris_3c_[i]);
    }
    else{
      throw FeatureNotImplemented("non-coulomb metrics in CADFCLHF", __FILE__, __LINE__, class_desc());
    }
  }
  // TwoCenter versions
  integral()->set_basis(dfbs_, dfbs_);
  eris_2c_.resize(nthread);
  do_threaded(nthread, [&](int ithr){
    eris_2c_[ithr] = integral()->coulomb<2>();
  });
  for (int i=0; i < nthread; i++) {
    if(metric_oper_type_ == coulomb_oper_type_){
      metric_ints_2c_.push_back(eris_2c_[i]);
    }
    else{
      throw FeatureNotImplemented("non-coulomb metrics in CADFCLHF", __FILE__, __LINE__, class_desc());
    }
  }
  // Reset to normal setup
  integral()->set_basis(basis(), basis(), basis(), basis());
  //----------------------------------------------------------------------------//
  // Set up the all pairs vector, needed to prescreen Schwarz bounds
  if(not dynamic_){
    const int nshell = basis()->nshell();
    for(int ish=0, inode=0; ish < nshell; ++ish){
      for(int jsh=0; jsh <= ish; ++jsh, ++inode){
        IntPair ij(ish, jsh);
        pair_assignments_[AllPairs][ij] = inode % n_node;
      }
    }
    // Make the backwards mapping for the current node
    const int me = scf_grp_->me();
    for(auto it : pair_assignments_[AllPairs]){
      if(it.second == me){
        local_pairs_[AllPairs].push_back(it.first);
      }
    }
  }
  //----------------------------------------------------------------------------//
  init_significant_pairs();
  //----------------------------------------------------------------------------//
  if(not dynamic_){
    const int me = scf_grp_->me();
    int inode = 0;
    for(auto sig_pair : sig_pairs_) {
      const int assignment = inode % n_node;
      pair_assignments_[SignificantPairs][sig_pair] = assignment;
      if(assignment == me){
        local_pairs_[SignificantPairs].push_back(sig_pair);
      }
      ++inode;
    }
    // Make the assignments for the mu, X pairs in K
    for(int mu_set = 0; mu_set < sig_blocks_.size(); ++mu_set) {
      for(auto sig_block : sig_blocks_[mu_set]) {
        const int assignment = inode % n_node; ++inode;
        auto pair = std::make_pair(mu_set, sig_block);
        pair_assignments_k_[SignificantPairs][pair] = assignment;
        if(assignment == me) {
          local_pairs_k_[SignificantPairs].push_back(pair);
        }
      } // end loop over blocks for mu
    } // in loop over mu sets
    for(auto ish : shell_range(basis())) {
      for(auto Yblk : shell_block_range(dfbs_, 0, 0, NoLastIndex, NoRestrictions)) {
        const int assignment = inode % n_node; ++inode;
        auto pair = std::make_pair((int)ish, Yblk);
        pair_assignments_k_[AllPairs][pair] = assignment;
        if(assignment == me) {
          local_pairs_k_[AllPairs].push_back(pair);
        }
      }
    }
  }
  //----------------------------------------------------------------------------//
  threads_initialized_ = true;
}

//////////////////////////////////////////////////////////////////////////////////

void
CADFCLHF::init_significant_pairs()
{
  Timer timer("init significant pairs");
  ExEnv::out0() << "  Computing significant shell pairs" << endl;
  std::atomic_int n_significant_pairs(0);
  const int nthread = threadgrp_->nthread();
  boost::mutex pair_mutex;
  std::vector<std::pair<double, IntPair>> pair_values;
  //----------------------------------------//
  local_pairs_spot_ = 0;
  do_threaded(nthread, [&](int ithr){

    ShellData ish, jsh;
    std::vector<std::pair<double, IntPair>> my_pair_vals;

    while(get_shell_pair(ish, jsh, AllPairs)){
      const double max_val = ints4maxes_.get(ish, jsh, ish, jsh, coulomb_oper_type_, [&]() -> double {
        tbis_[ithr]->compute_shell(ish, jsh, ish, jsh);
        const double* buffer = tbis_[ithr]->buffer(coulomb_oper_type_);
        /*
        return *std::max_element(buffer, buffer + ish.nbf * jsh.nbf * ish.nbf * jsh.nbf,
            [](double a, double b){ return fabs(a) < fabs(b); }
        );*/
        double frob_val = 0.0;
        for(int i = 0; i < ish.nbf*jsh.nbf*ish.nbf*jsh.nbf; ++i) {
          frob_val += fabs(buffer[i]);
        }
        return sqrt(frob_val);
      });
      my_pair_vals.push_back(make_pair(sqrt(fabs(max_val)), IntPair(ish, jsh)));
    } // end while get shell pair
    //----------------------------------------//
    // put our values on the node-level vector
    boost::lock_guard<boost::mutex> lg(pair_mutex);
    for(auto item : my_pair_vals){
      pair_values.push_back(item);
    }

  });
  //----------------------------------------//
  double max_schwarz = std::max_element(pair_values.begin(), pair_values.end())->first;
  // Take the max across all nodes
  scf_grp_->max(max_schwarz);
  //----------------------------------------//
  // Now go through the list and figure out which ones are significant
  shell_to_sig_shells_.resize(basis()->nshell());
  std::vector<double> sig_values;
  schwarz_frob_.resize(basis()->nshell(), basis()->nshell());
  schwarz_frob_ = Eigen::MatrixXd::Zero(basis()->nshell(), basis()->nshell());
  do_threaded(nthread, [&](int ithr){
    std::vector<IntPair> my_sig_pairs;
    std::vector<double> my_sig_values;
    for(int i = ithr; i < pair_values.size(); i += nthread){
      auto item = pair_values[i];
      schwarz_frob_(item.second.first, item.second.second) = item.first;
      schwarz_frob_(item.second.second, item.second.first) = item.first;
      if(item.first * max_schwarz > pair_screening_thresh_){
        my_sig_pairs.push_back(item.second);
        my_sig_values.push_back(item.first);
        ++n_significant_pairs;
      }
    }
    //----------------------------------------//
    // put our values on the node-level vector
    boost::lock_guard<boost::mutex> lg(pair_mutex);
    using boost_ext::zip;
    for(auto item : my_sig_pairs){
      sig_pairs_.push_back(item);
    }
    for(auto item : my_sig_values){
      sig_values.push_back(item);
    }
  });
  //----------------------------------------//
  // Now all-to-all the significant pairs
  // This should be done with an MPI_alltoall_v or something like that
  scf_grp_->sum(schwarz_frob_.data(), basis()->nshell() * basis()->nshell());
  int n_sig = n_significant_pairs;
  int my_n_sig = n_significant_pairs;
  int n_sig_pairs[scf_grp_->n()];
  for(int inode = 0; inode < scf_grp_->n(); n_sig_pairs[inode++] = 0);
  n_sig_pairs[scf_grp_->me()] = n_sig;
  scf_grp_->sum(n_sig_pairs, scf_grp_->n());
  scf_grp_->sum(n_sig);
  {
    int sig_data[n_sig*2];
    double sig_vals[n_sig];
    for(int inode = 0; inode < n_sig*2; sig_vals[inode] = 0.0, sig_data[inode++] = 0);
    int data_offset = 0;
    for(int inode = 0; inode < scf_grp_->me(); inode++){
      data_offset += 2 * n_sig_pairs[inode];
    }
    for(int ipair = 0; ipair < my_n_sig; ++ipair) {
      sig_data[data_offset + 2 * ipair] = sig_pairs_[ipair].first;
      sig_data[data_offset + 2 * ipair + 1] = sig_pairs_[ipair].second;
      sig_vals[data_offset/2 + ipair] = sig_values[ipair];
    }
    scf_grp_->sum(sig_data, n_sig*2);
    scf_grp_->sum(sig_vals, n_sig);
    sig_pairs_.clear();
    sig_values.clear();
    for(int ipair = 0; ipair < n_sig; ++ipair){
      sig_pairs_.push_back(IntPair(sig_data[2*ipair], sig_data[2*ipair + 1]));
      sig_values.push_back(sig_vals[ipair]);
    }
  } // get rid of sig_data and sig_vals
  //----------------------------------------//
  // Now compute the significant pairs for the outer loop of the exchange and the sig_partners_ array
  sig_partners_.resize(basis()->nbasis());
  sig_blocks_.resize(basis()->nbasis());
  int pair_index = 0;
  for(auto pair : sig_pairs_) {
    ShellData ish(pair.first, basis()), jsh(pair.second, basis());
    sig_partners_[ish].insert(jsh);
    L_schwarz[ish].insert(jsh, sig_values[pair_index]);
    if(ish != jsh) {
      sig_partners_[jsh].insert(ish);
      L_schwarz[jsh].insert(ish, sig_values[pair_index]);
    }
    for(auto Xblk : iter_shell_blocks_on_center(dfbs_, ish.center)) {
      sig_blocks_[ish].insert(Xblk);
      sig_blocks_[jsh].insert(Xblk);
    }
    for(auto Xblk : iter_shell_blocks_on_center(dfbs_, jsh.center)) {
      sig_blocks_[ish].insert(Xblk);
      sig_blocks_[jsh].insert(Xblk);
    }
    //----------------------------------------//
    ++pair_index;
  }
  out_assert(L_schwarz.size(), >, 0);
  //----------------------------------------//
  #if !LINK_SORTED_INSERTION
  do_threaded(nthread, [&](int ithr){
    auto L_schwarz_iter = L_schwarz.begin();
    const auto& L_schwarz_end = L_schwarz.end();
    L_schwarz_iter.advance(ithr);
    while(L_schwarz_iter != L_schwarz_end) {
      L_schwarz_iter->second.sort();
      L_schwarz_iter.advance(nthread);
    }
  });
  out_assert(L_schwarz[0].size(), >, 0);
  #endif
  //----------------------------------------//
  max_schwarz_.resize(basis()->nshell());
  for(auto ish : shell_range(basis())) {
    double max_val = 0.0;
    for(auto ish : L_schwarz[ish]) {
      if(ish.value > max_val) {
        max_val = ish.value;
      }
    }
    max_schwarz_[ish] = max_val;
  }
  //----------------------------------------//
  ExEnv::out0() << "  Number of significant shell pairs:  " << n_sig << endl;
  ExEnv::out0() << "  Number of total basis pairs:  " << (basis()->nshell() * (basis()->nshell() + 1) / 2) << endl;
  ExEnv::out0() << "  max_schwarz = " << max_schwarz << endl;

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
  //---------------------------------------------------------------------------------------//
  if(not have_coefficients_) {
    ints_computed_locally_ = 0;
    compute_coefficients();
    ints_computed_ = ints_computed_locally_;
    scf_grp_->sum(ints_computed_);
    if(scf_grp_->me() == 0) {
      ExEnv::out0() << "  Computed " << ints_computed_ << " integrals to determine coefficients." << endl;
    }
  }
  //---------------------------------------------------------------------------------------//
  timer.enter("misc");
  int nthread = threadgrp_->nthread();
  //----------------------------------------//
  // transform the density difference to the AO basis
  RefSymmSCMatrix dd = cl_dens_diff_;
  Ref<PetiteList> pl = integral()->petite_list();
  cl_dens_diff_ = pl->to_AO_basis(dd);
  //----------------------------------------//
  double gmat_accuracy = accuracy;
  if (min_orthog_res() < 1.0) { gmat_accuracy *= min_orthog_res(); }
  //----------------------------------------//
  // copy over the density
  // cl_dens_diff_ includes total density right now, so halve it
  D_ = cl_dens_diff_.copy().convert2RefSCMat(); D_.scale(0.5);
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Form G                                                		                        {{{1 */ #if 1 // begin fold
  //---------------------------------------------------------------------------------------//
  // compute J and K
  timer.change("build");
  if(xml_debug){
    begin_xml_context("compute_fock", "compute_fock.xml");
  }
  RefSCMatrix G;
  {
    ints_computed_locally_ = 0;
    if(xml_debug) begin_xml_context("compute_J");
    RefSCMatrix J = compute_J();
    if(xml_debug) write_as_xml("J", J), end_xml_context("compute_J");
    G = J.copy();
    ints_computed_ = ints_computed_locally_;
    scf_grp_->sum(ints_computed_);
    if(scf_grp_->me() == 0) {
      ExEnv::out0() << "        Computed " << ints_computed_ << " integrals for J part" << endl;
    }
  }
  {
    ints_computed_locally_ = 0;
    if(xml_debug) begin_xml_context("compute_K");
    RefSCMatrix K = compute_K();
    if(xml_debug) write_as_xml("K", K), end_xml_context("compute_K");
    G.accumulate( -1.0 * K);
    ints_computed_ = ints_computed_locally_;
    scf_grp_->sum(ints_computed_);
    if(scf_grp_->me() == 0) {
      ExEnv::out0() << "        Computed " << ints_computed_ << " integrals for K part" << endl;
    }
  }
  if(xml_debug) end_xml_context("compute_fock"), assert(false);
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
}

//////////////////////////////////////////////////////////////////////////////////

RefSCMatrix
CADFCLHF::compute_J()
{
  /*=======================================================================================*/
  /* Setup                                                 		                        {{{1 */ #if 1 // begin fold
  //----------------------------------------//
  // Convenience variables
  Timer timer("compute J");
  MessageGrp& msg = *scf_grp_;
  const int nthread = threadgrp_->nthread();
  const int me = msg.me();
  const int n_node = msg.n();
  const Ref<GaussianBasisSet>& obs = basis();
  const int nbf = obs->nbasis();
  const int dfnbf = dfbs_->nbasis();
  //----------------------------------------//
  // Get the density in an Eigen::Map form
  double *D_ptr = allocate<double>(nbf*nbf);
  D_.convert(D_ptr);
  typedef Eigen::Map<Eigen::VectorXd> VectorMap;
  typedef Eigen::Map<Eigen::MatrixXd> MatrixMap;
  // Matrix and vector wrappers, for convenience
  VectorMap d(D_ptr, nbf*nbf);
  MatrixMap D(D_ptr, nbf, nbf);
  //----------------------------------------//
  Eigen::MatrixXd J(nbf, nbf);
  J = Eigen::MatrixXd::Zero(nbf, nbf);
  //----------------------------------------//
  // reset the iteration over local pairs
  local_pairs_spot_ = 0;
  //----------------------------------------//
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Form C_tilde and d_tilde                              		                        {{{1 */ #if 1 // begin fold
  //----------------------------------------//
  timer.enter("compute C_tilde");
  Eigen::VectorXd C_tilde(dfnbf);
  C_tilde = Eigen::VectorXd::Zero(dfnbf);
  {
    boost::mutex C_tilde_mutex;
    boost::thread_group compute_threads;
    // Loop over number of threads
    for(int ithr = 0; ithr < nthread; ++ithr) {
      // ...and create each thread that computes pairs
      compute_threads.create_thread([&,ithr](){
        /*-----------------------------------------------------*/
        /* C_tilde compute thread                         {{{2 */ #if 2 // begin fold
        Eigen::VectorXd Ct(dfnbf);
        Ct = Eigen::VectorXd::Zero(dfnbf);
        ShellData ish, jsh;
        //----------------------------------------//
        while(get_shell_pair(ish, jsh, SignificantPairs)){
          //----------------------------------------//
          // Permutation prefactor
          double pf = (ish == jsh) ? 2.0 : 4.0;
          //----------------------------------------//
          for(int mu = ish.bfoff; mu <= ish.last_function; ++mu){
            for(int nu = jsh.bfoff; nu <= jsh.last_function; ++nu){
              IntPair ij(mu, nu);
              auto& cpair = coefs_[ij];
              VectorMap& Ca = *(cpair.first);
              VectorMap& Cb = *(cpair.second);
              Ct.segment(ish.atom_dfbfoff, ish.atom_dfnbf) += pf * D(mu, nu) * Ca;
              if(ish.center != jsh.center) {
                Ct.segment(jsh.atom_dfbfoff, jsh.atom_dfnbf) += pf * D(mu, nu) * Cb;
              }
            } // end loop over mu
          } // end loop over nu
        } // end while get shell pair
        //----------------------------------------//
        // add our contribution to the node level C_tilde
        boost::lock_guard<boost::mutex> Ctlock(C_tilde_mutex);
        C_tilde += Ct;
        //----------------------------------------//
        /********************************************************/ #endif //2}}}
        /*-----------------------------------------------------*/
      }); // end create_thread
    } // end enumeration of threads
    compute_threads.join_all();
  } // compute_threads is destroyed here
  //----------------------------------------//
  // Global sum C_tilde
  scf_grp_->sum(C_tilde.data(), dfnbf);
  if(xml_debug) {
    write_as_xml("C_tilde", C_tilde);
  }
  timer.exit("compute C_tilde");
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Form g_tilde                                          		                        {{{1 */ #if 1 // begin fold
  //----------------------------------------//
  // TODO Thread this
  timer.enter("compute g_tilde");
  auto X_g_Y = ints_to_eigen(
    ShellBlockData<>(dfbs_), ShellBlockData<>(dfbs_),
    eris_2c_[0],
    coulomb_oper_type_
  );
  Eigen::VectorXd g_tilde(dfnbf);
  g_tilde = (*X_g_Y) * C_tilde;
  if(xml_debug) {
    write_as_xml("g_tilde", g_tilde);
  }
  timer.exit("compute g_tilde");
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Form dtilde and add in Ctilde contribution to J       		                        {{{1 */ #if 1 // begin fold
  timer.enter("compute d_tilde");
  Eigen::VectorXd d_tilde(dfnbf);
  d_tilde = Eigen::VectorXd::Zero(dfnbf);
  {
    boost::mutex tmp_mutex;
    boost::thread_group compute_threads;
    //----------------------------------------//
    // reset the iteration over local pairs
    local_pairs_spot_ = 0;
    // Loop over number of threads
    for(int ithr = 0; ithr < nthread; ++ithr) {
      // ...and create each thread that computes pairs
      compute_threads.create_thread([&,ithr](){
        /*-----------------------------------------------------*/
        /* d_tilde compute and C_tilde contract thread    {{{2 */ #if 2 // begin fold
        Eigen::VectorXd dt(dfnbf);
        dt = Eigen::VectorXd::Zero(dfnbf);
        Eigen::MatrixXd jpart(nbf, nbf);
        jpart = Eigen::MatrixXd::Zero(nbf, nbf);
        //----------------------------------------//
        ShellData ish, jsh;
        while(get_shell_pair(ish, jsh, SignificantPairs)){
          //----------------------------------------//
          double perm_fact = (ish == jsh) ? 2.0 : 4.0;
          //----------------------------------------//
          for(auto Xblk : shell_block_range(dfbs_)){
            auto g_part = ints_to_eigen(
                ish, jsh, Xblk,
                eris_3c_[ithr],
                coulomb_oper_type_
            );
            for(auto mu : function_range(ish)) {
              dt.segment(Xblk.bfoff, Xblk.nbf).transpose() += perm_fact * D.row(mu).segment(jsh.bfoff, jsh.nbf)
                  * g_part->middleRows(mu.bfoff_in_shell*jsh.nbf, jsh.nbf);
              //----------------------------------------//
              jpart.row(mu).segment(jsh.bfoff, jsh.nbf).transpose() +=
                  g_part->middleRows(mu.bfoff_in_shell*jsh.nbf, jsh.nbf)
                  * C_tilde.segment(Xblk.bfoff, Xblk.nbf);
            }
          } // end loop over kshbf
          // only constructing the lower triangle of the J matrix, so zero the strictly upper part
          jpart.triangularView<Eigen::StrictlyUpper>().setZero();
        } // end while get_shell_pair
        //----------------------------------------//
        // add our contribution to the node level d_tilde
        boost::lock_guard<boost::mutex> tmp_lock(tmp_mutex);
        d_tilde += dt;
        J += jpart;
        /*******************************************************/ #endif //2}}}
        /*-----------------------------------------------------*/
      }); // end create_thread
    } // end enumeration of threads
    compute_threads.join_all();
  } // compute_threads is destroyed here
  //----------------------------------------//
  // Global sum d_tilde
  msg.sum(d_tilde.data(), dfnbf);
  if(xml_debug) {
    write_as_xml("d_tilde", d_tilde);
  }
  timer.exit("compute d_tilde");
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Add in first and third term contributions to J       		                        {{{1 */ #if 1 // begin fold
  {
    timer.enter("J contributions");
    boost::mutex tmp_mutex;
    boost::thread_group compute_threads;
    //----------------------------------------//
    // reset the iteration over local pairs
    local_pairs_spot_ = 0;
    // Loop over number of threads
    for(int ithr = 0; ithr < nthread; ++ithr) {
      // ...and create each thread that computes pairs
      compute_threads.create_thread([&,ithr](){
        /*-----------------------------------------------------*/
        /* first and third terms compute thread           {{{2 */ #if 2 // begin fold
        Eigen::MatrixXd jpart(nbf, nbf);
        jpart = Eigen::MatrixXd::Zero(nbf, nbf);
        //----------------------------------------//
        ShellData ish, jsh;
        while(get_shell_pair(ish, jsh, SignificantPairs)){
          //----------------------------------------//
          for(auto mu : function_range(ish)){
            for(auto nu : function_range(jsh)){
              //----------------------------------------//
              IntPair ij(mu, nu);
              auto cpair = coefs_[ij];
              VectorMap& Ca = *(cpair.first);
              VectorMap& Cb = *(cpair.second);
              //----------------------------------------//
              // First term contribution
              jpart(mu, nu) += d_tilde.segment(ish.atom_dfbfoff, ish.atom_dfnbf).transpose() * Ca;
              if(ish.center != jsh.center){
                jpart(mu, nu) += d_tilde.segment(jsh.atom_dfbfoff, jsh.atom_dfnbf).transpose() * Cb;
              }
              //----------------------------------------//
              // Third term contribution
              jpart(mu, nu) -= g_tilde.segment(ish.atom_dfbfoff, ish.atom_dfnbf).transpose() * Ca;
              if(ish.center != jsh.center){
                jpart(mu, nu) -= g_tilde.segment(jsh.atom_dfbfoff, jsh.atom_dfnbf).transpose() * Cb;
              }
              //----------------------------------------//
            } // end loop over nu
          } // end loop over mu
        } // end while get_shell_pair
        // Sum the thread's contributions to the node-level J
        boost::lock_guard<boost::mutex> tmp_lock(tmp_mutex);
        J += jpart;
        /*******************************************************/ #endif //2}}}
        /*-----------------------------------------------------*/
      }); // end create_thread
    } // end enumeration of threads
    compute_threads.join_all();
    timer.exit("J contributions");
  } // compute_threads is destroyed here
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Global sum J                                         		                        {{{1 */ #if 1 // begin fold
  //----------------------------------------//
  msg.sum(J.data(), nbf*nbf);
  //----------------------------------------//
  // Fill in the upper triangle of J
  for(int mu = 0; mu < nbf; ++mu) {
    for(int nu = 0; nu < mu; ++nu) {
      J(nu, mu) = J(mu, nu);
    }
  }
  //----------------------------------------//
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Transfer J to a RefSCMatrix                           		                        {{{1 */ #if 1 // begin fold
  Ref<Integral> localints = integral()->clone();
  localints->set_basis(obs);
  Ref<PetiteList> pl = localints->petite_list();
  RefSCDimension obsdim = pl->AO_basisdim();
  RefSCMatrix result(
      obsdim,
      obsdim,
      obs->so_matrixkit()
  );
  result.assign(J.data());
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Clean up                                             		                        {{{1 */ #if 1 // begin fold
  //----------------------------------------//
  deallocate(D_ptr);
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  return result;
}

//////////////////////////////////////////////////////////////////////////////////

RefSCMatrix
CADFCLHF::compute_K()
{
  /*=======================================================================================*/
  /* Setup                                                 		                        {{{1 */ #if 1 // begin fold
  //----------------------------------------//
  // Convenience variables
  Timer timer("compute K");
  MessageGrp& msg = *scf_grp_;
  const int nthread = threadgrp_->nthread();
  const int me = msg.me();
  const int n_node = msg.n();
  const Ref<GaussianBasisSet>& obs = basis();
  GaussianBasisSet* obsptr = basis();
  GaussianBasisSet* dfbsptr = dfbs_;
  const int nbf = obs->nbasis();
  const int dfnbf = dfbs_->nbasis();
  //----------------------------------------//
  // Get the density in an Eigen::Map form
  double *D_ptr = allocate<double>(nbf*nbf);
  D_.convert(D_ptr);
  typedef Eigen::Map<Eigen::VectorXd> VectorMap;
  typedef Eigen::Map<Eigen::MatrixXd> MatrixMap;
  // Matrix and vector wrappers, for convenience
  VectorMap d(D_ptr, nbf*nbf);
  MatrixMap D(D_ptr, nbf, nbf);
  // Match density scaling in old code:
  D *= 0.5;
  //----------------------------------------//
  Eigen::MatrixXd Kt(nbf, nbf);
  Eigen::MatrixXd Kt2(nbf, nbf);
  Eigen::MatrixXd K2(nbf, nbf);
  Kt = Eigen::MatrixXd::Zero(nbf, nbf);
  Kt2 = Eigen::MatrixXd::Zero(nbf, nbf);
  K2 = Eigen::MatrixXd::Zero(nbf, nbf);
  //----------------------------------------//
  // reset the iteration over local pairs
  local_pairs_spot_ = 0;
  //----------------------------------------//
  auto get_ish_Xblk_pair = [&](ShellData& ish, ShellBlockData<>& Xblk, PairSet ps) -> bool {
    int spot = local_pairs_spot_++;
    if(spot < local_pairs_k_[ps].size()){
      auto& sig_pair = local_pairs_k_[ps][spot];
      ish = ShellData(sig_pair.first, basis(), dfbs_);
      Xblk = sig_pair.second;
      return true;
    }
    else{
      return false;
    }
  };
  //----------------------------------------//
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Make the CADF-LinK lists                                                         {{{1 */ #if 1 // begin fold
  if(do_linK_){
    timer.enter("LinK lists");
    // First clear all of the lists
    L_D.clear();
    L_3.clear();
    //for(auto L4_key : L_4_keys) {
    //  L_4(L4_key.first, L4_key.second).clear();
    //}
    //L_4_keys.clear();
    //============================================================================//
    std::vector<Eigen::Vector3d> positions;
    auto& atoms = molecule()->atoms();
    positions.resize(atoms.size());
    int i = 0;
    for(auto& atom : atoms) {
      positions[i++] << atom.xyz(0), atom.xyz(1), atom.xyz(2);
    }
    Eigen::MatrixXd D_frob(obs->nshell(), obs->nshell());
    // TODO Distribute over MPI processes as well
    #if LINK_TIME_LIST_INSERTIONS
    std::atomic_uint_fast64_t run_time{ 0 };
    OrderedShellList::insertion_time_nanos = 0;
    OrderedShellList::n_insertions = 0;
    //int n_insertions = 0;
    auto start_link = std::chrono::high_resolution_clock::now();
    #endif
    do_threaded(nthread, [&](int ithr){
      #if LINK_TIME_LIST_INSERTIONS
      auto inner_timer = make_auto_timer(run_time);
      #endif
      for(int lsh_index = ithr; lsh_index < obs->nshell(); lsh_index += nthread) {
        ShellData lsh(lsh_index, obs, dfbs_);
        for(auto jsh : shell_range(obs)) {
          double dnorm = D.block(lsh.bfoff, jsh.bfoff, lsh.nbf, jsh.nbf).norm();
          D_frob(lsh, jsh) = dnorm;
          double screen_val = fabs(dnorm * max_schwarz_[lsh] * max_schwarz_[jsh]);
          // NOTE!!! This is not an exact condition, so it may cause differences in
          //   very small decimal places.
          if(screen_val > density_screening_thresh_) {
            L_D[lsh].insert(jsh, screen_val);
          }
        } // end loop over jsh
        #if !LINK_SORTED_INSERTION
        L_D[lsh].sort();
        #endif
        //----------------------------------------//
        for(auto ksh : L_coefs[lsh]) {
          //----------------------------------------//
          bool lists_changed_outer = false;
          for(auto Xsh : iter_shells_on_centers(dfbs_, lsh.center, ksh.center)) {
            const double CXnorm = lsh.center == Xsh.center ?
                C_trans_frob_[Xsh](lsh.shoff_in_atom, ksh)
                : C_trans_frob_[Xsh](ksh.shoff_in_atom, lsh);
            bool lists_changed_inner = false;
            //------------------------------------------//
            for(auto jsh : L_D[lsh]) {
              bool mu_added = false;
              //----------------------------------------//
              for(auto ish : L_schwarz[jsh]) {
                const double screen_val = fabs(
                    CXnorm
                    * schwarz_df_[Xsh]
                    * D_frob(lsh, jsh)
                    * schwarz_frob_(ish, jsh)
                );
                if(screen_val > full_screening_thresh_) {
                  mu_added = true; lists_changed_outer = true; lists_changed_inner = true;
                  L_3[{ish, Xsh}].insert(jsh, screen_val);
                  //L_4((int)ish, jsh).insert(lsh, screen_val);
                }
              } // end loop over ish_w_val (a.k.a. mu)
              //----------------------------------------//
              if(not mu_added) break;
            //----------------------------------------//
            } // end loop over jsh (a.k.a. rho)
            //if(lists_changed_inner) {
            //  L_K(lsh, Xsh).insert(ksh, 0);
            //  L_K_keys.emplace(lsh, Xsh);
            //}
          } // end loop over Xsh
          //----------------------------------------//
          if(not lists_changed_outer) break;
        } // end loop over ksh (a.k.a. nu)
        //----------------------------------------//
      } // end loop over lsh
    });
    #if not LINK_SORTED_INSERTION
    do_threaded(nthread, [&](int ithr){
      auto L_3_iter = L_3.begin();
      const auto& L_3_end = L_3.end();
      L_3_iter.advance(ithr);
      while(L_3_iter != L_3_end) {
        L_3_iter->second.sort();
        L_3_iter.advance(nthread);
      }
    });
    #endif
    #if LINK_TIME_LIST_INSERTIONS
    auto end_link = std::chrono::high_resolution_clock::now();
    static std::locale loc("");
    auto old_loc = std::cout.getloc();
    std::cout.imbue(loc);
    cout << "        LinK list insertions took "
        << double(OrderedShellList::insertion_time_nanos)/1e6
        << " milliseconds to do " << OrderedShellList::n_insertions
        << " insertions"
        << endl
        << "            ( = "
        << double(OrderedShellList::insertion_time_nanos) / double(OrderedShellList::n_insertions)
        << " ns per insertion)"
        << endl;
    cout << "        LinK list build took "
        << double(run_time)/1e6
        << " milliseconds ("
        << double(std::chrono::duration_cast<std::chrono::nanoseconds>(end_link - start_link).count())/1e6 * double(nthread)
        << " ms ideal)"
        << endl;
    std::cout.imbue(old_loc);
    #endif
    timer.exit("LinK lists");
    if(print_screening_stats_ and scf_grp_->me() == 0) {
      int n_LD = 0, n_L3 = 0, n_L4 = 0;
      const int n_LD_max = obs->nshell() * obs->nshell();
      int n_L3_max = 0;
      for(auto L_D_item : L_D) {
        n_LD += L_D_item.second.size();
      }
      ExEnv::out0() << "        Density pairs skipped: " << (n_LD_max - n_LD)
          << "/" << n_LD_max << " (= " << setprecision(2) << (100 * double(n_LD_max - n_LD)/double(n_LD_max)) << "%)"
          << endl;
      const int max_L3_keys = obs->nshell() * dfbs_->nshell();
      const int n_L3_keys = L_3.size();
      ExEnv::out0() << "        obs/dfbs pairs skipped in K: " << (max_L3_keys - n_L3_keys)
          << "/" << max_L3_keys << " (= "
          << setprecision(2) << (100.0 * double(max_L3_keys - n_L3_keys)/double(max_L3_keys)) << "%)"
          << endl;
    }
  } // end if do_linK_
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Loop over local shell pairs for three body contributions                         {{{1 */ #if 1 // begin fold
  {
    Timer timer("three body contributions");
    boost::mutex tmp_mutex, L_3_mutex;
    auto L_3_key_iter = L_3.begin();
    for(int i = 0; i < scf_grp_->me(); ++i, ++L_3_key_iter);
    if(scf_grp_->me() >= L_3.size()) L_3_key_iter = L_3.end();
    const int n_node = scf_grp_->n();
    boost::thread_group compute_threads;
    //----------------------------------------//
    // reset the iteration over local pairs
    local_pairs_spot_ = 0;
    // Loop over number of threads
    for(int ithr = 0; ithr < nthread; ++ithr) {
      // ...and create each thread that computes pairs
      compute_threads.create_thread([&,ithr](){
        Eigen::MatrixXd Kt_part(nbf, nbf);
        Kt_part = Eigen::MatrixXd::Zero(nbf, nbf);
        //----------------------------------------//
        ShellData ish;
        ShellBlockData<> Xblk;
        //----------------------------------------//
        auto get_ish_Xblk_3 = [&]() -> bool {
          if(do_linK_) {
            boost::lock_guard<boost::mutex> lg(L_3_mutex);
            if(L_3_key_iter == L_3.end()) {
              return false;
            }
            else {
              auto& pair = L_3_key_iter->first;
              L_3_key_iter.advance(n_node);
              ish = ShellData(pair.first, basis(), dfbs_);
              Xblk = ShellBlockData<>(ShellData(pair.second, dfbs_, basis()));
              return true;
            }
          }
          else {
            int spot = local_pairs_spot_++;
            if(spot < local_pairs_k_[SignificantPairs].size()){
              auto& sig_pair = local_pairs_k_[SignificantPairs][spot];
              ish = ShellData(sig_pair.first, basis(), dfbs_);
              Xblk = sig_pair.second;
              return true;
            }
            else{
              return false;
            }
          }
        };
        //----------------------------------------//
        while(get_ish_Xblk_3()) {
          /*-----------------------------------------------------*/
          /* Compute B intermediate                         {{{2 */ #if 2 // begin fold
          std::vector<Eigen::MatrixXd> B_mus(ish.nbf);
          for(auto mu : function_range(ish)) {
            B_mus[mu.bfoff_in_shell].resize(nbf, Xblk.nbf);
            B_mus[mu.bfoff_in_shell] = Eigen::MatrixXd::Zero(nbf, Xblk.nbf);
          }
          if(do_linK_){
            for(auto Xsh : shell_range(Xblk)) {
              for(ShellData jsh : L_3[{ish, Xsh}]){
                assert(Xblk.nshell == 1);
                auto g3_ptr = ints_to_eigen(
                    ish, jsh, Xsh,
                    eris_3c_[ithr], coulomb_oper_type_
                );
                auto& g3 = *g3_ptr;
                for(auto mu : function_range(ish)) {
                  B_mus[mu.bfoff_in_shell] += 2.0 *
                      D.middleRows(jsh.bfoff, jsh.nbf).transpose() *
                      g3.middleRows(mu.bfoff_in_shell*jsh.nbf, jsh.nbf);
                } // end loop over mu
              } // end loop over jsh
            }
          }
          else {
            for(ShellData jsh : iter_significant_partners(ish)){
              auto g3_ptr = ints_to_eigen(
                  ish, jsh, Xblk,
                  eris_3c_[ithr], coulomb_oper_type_
              );
              auto& g3 = *g3_ptr;
              for(auto mu : function_range(ish)) {
                B_mus[mu.bfoff_in_shell] += 2.0 *
                    D.middleRows(jsh.bfoff, jsh.nbf).transpose() *
                    g3.middleRows(mu.bfoff_in_shell*jsh.nbf, jsh.nbf);
              } // end loop over mu
            } // end loop over jsh
          }
          if(xml_debug) {
            for(auto mu : function_range(ish)){
              write_as_xml("B_part", B_mus[mu.bfoff_in_shell], std::map<std::string, int>{
                {"mu", mu}, {"Xbfoff", Xblk.bfoff},
                {"Xnbf", Xblk.nbf}, {"dfnbf", dfbs_->nbasis()}
              });
            }
          }
          /*******************************************************/ #endif //2}}}
          /*-----------------------------------------------------*/
          /* Compute K contributions                        {{{2 */ #if 2 // begin fold
          const int obs_atom_bfoff = obs->shell_to_function(obs->shell_on_center(Xblk.center, 0));
          const int obs_atom_nbf = obs->nbasis_on_center(Xblk.center);
          for(auto X : function_range(Xblk)) {
            Eigen::MatrixXd& C_X = coefs_transpose_[X];
            for(auto mu : function_range(ish)) {
              // B_mus[mu.bfoff_in_shell] is (nbf x Ysh.nbf)
              // C_Y is (Y.{obs_}atom_nbf x nbf)
              // result should be (Y.{obs_}atom_nbf x 1)
              Kt_part.row(mu).segment(obs_atom_bfoff, obs_atom_nbf).transpose() +=
                  C_X * B_mus[mu.bfoff_in_shell].col(X.bfoff_in_block);
              // The sigma <-> nu term
              Kt_part.row(mu).transpose() += C_X.transpose()
                  * B_mus[mu.bfoff_in_shell].col(X.bfoff_in_block).segment(obs_atom_bfoff, obs_atom_nbf);
              // Add back in the nu.center == Y.center part
              Kt_part.row(mu).segment(obs_atom_bfoff, obs_atom_nbf).transpose() -=
                  C_X.middleCols(obs_atom_bfoff, obs_atom_nbf).transpose()
                  * B_mus[mu.bfoff_in_shell].col(X.bfoff_in_block).segment(obs_atom_bfoff, obs_atom_nbf);
              //----------------------------------------//
            }
          }
          /*******************************************************/ #endif //2}}}
          /*-----------------------------------------------------*/
        } // end while get ish Xblk pair
        //============================================================================//
        //----------------------------------------//
        // Sum Kt parts within node
        boost::lock_guard<boost::mutex> lg(tmp_mutex);
        Kt += Kt_part;
      }); // end create_thread
    } // end enumeration of threads
    compute_threads.join_all();
  } // compute_threads is destroyed here
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Loop over local shell pairs and compute two body part 		                        {{{1 */ #if 1 // begin fold
  {
    Timer timer("two body contributions");
    boost::mutex tmp_mutex;
    boost::thread_group compute_threads;
    //----------------------------------------//
    // reset the iteration over local pairs
    local_pairs_spot_ = 0;
    // Loop over number of threads
    auto g2_full_ptr = ints_to_eigen(
        ShellBlockData<>(dfbs_),
        ShellBlockData<>(dfbs_),
        eris_2c_[0], coulomb_oper_type_
    );
    for(int ithr = 0; ithr < nthread; ++ithr) {
      // ...and create each thread that computes pairs
      compute_threads.create_thread([&,ithr](){
        Eigen::MatrixXd Kt_part(nbf, nbf);
        Kt_part = Eigen::MatrixXd::Zero(nbf, nbf);
        ShellData ish;
        ShellBlockData<> Yblk;
        while(get_ish_Xblk_pair(ish, Yblk, SignificantPairs)) {
          std::vector<Eigen::MatrixXd> Ct_mus(ish.nbf);
          for(auto mu : function_range(ish)) {
            Ct_mus[mu.bfoff_in_shell].resize(nbf, Yblk.nbf);
            Ct_mus[mu.bfoff_in_shell] = Eigen::MatrixXd::Zero(nbf, Yblk.nbf);
          }
          for(auto Xblk : shell_block_range(dfbs_, 0, 0, NoLastIndex, SameCenter)){
            //auto g2_ptr = ints_to_eigen(
            //    Yblk, Xblk,
            //    eris_2c_[ithr], coulomb_oper_type_
            //);
            const auto& g2 = g2_full_ptr->block(Yblk.bfoff, Xblk.bfoff, Yblk.nbf, Xblk.nbf);
            for(auto jsh : iter_significant_partners(ish)) {
              for(auto mu : function_range(ish)) {
                for(auto rho : function_range(jsh)) {
                  BasisFunctionData first = mu, second = rho;
                  if(jsh > ish) first = rho, second = mu;
                  IntPair mu_rho(first, second);
                  assert(coefs_.find(mu_rho) != coefs_.end());
                  auto& cpair = coefs_[mu_rho];
                  auto& Cmu = (jsh <= ish || ish.center == jsh.center) ? *cpair.first : *cpair.second;
                  auto& Crho = jsh <= ish ? *cpair.second : *cpair.first;
                  if(Xblk.center == mu.center) {
                    Ct_mus[mu.bfoff_in_shell].row(rho) += 2.0 *
                        g2 * Cmu.segment(Xblk.bfoff_in_atom, Xblk.nbf);
                  }
                  else if(Xblk.center == rho.center){
                    Ct_mus[mu.bfoff_in_shell].row(rho) += 2.0 *
                        g2 * Crho.segment(Xblk.bfoff_in_atom, Xblk.nbf);
                  } // end loop over Y
                } // end loop over rho in jsh
              } // end loop over mu in ish
            } // end loop over jsh
          } // end loop over Xblk
          //----------------------------------------//
          std::vector<Eigen::MatrixXd> dt_mus(ish.nbf);
          for(auto mu : function_range(ish)) {
            dt_mus[mu.bfoff_in_shell].resize(nbf, Yblk.nbf);
            dt_mus[mu.bfoff_in_shell] = D * Ct_mus[mu.bfoff_in_shell];
          }
          if(xml_debug) {
            for(auto mu : function_range(ish)){
              write_as_xml("new_dt_part", dt_mus[mu.bfoff_in_shell], std::map<std::string, int>{
                {"mu", mu}, {"Xbfoff", Yblk.bfoff},
                {"Xnbf", Yblk.nbf}, {"dfnbf", dfbs_->nbasis()}
              });
            }
          }
          //----------------------------------------//
          for(auto Y : function_range(Yblk)) {
            Eigen::MatrixXd& C_Y = coefs_transpose_[Y];
            const int obs_atom_bfoff = obs->shell_to_function(obs->shell_on_center(Y.center, 0));
            const int obs_atom_nbf = obs->nbasis_on_center(Y.center);
            for(auto mu : function_range(ish)) {
              // dt_mus[mu.bfoff_in_shell] is (nbf x Ysh.nbf)
              // C_Y is (Y.{obs_}atom_nbf x nbf)
              // result should be (Y.{obs_}atom_nbf x 1)
              Kt_part.row(mu).segment(obs_atom_bfoff, obs_atom_nbf).transpose() -=
                  0.5 * C_Y * dt_mus[mu.bfoff_in_shell].col(Y.bfoff_in_block);
              // The sigma <-> nu term
              Kt_part.row(mu).transpose() -= 0.5
                  * C_Y.transpose()
                  * dt_mus[mu.bfoff_in_shell].col(Y.bfoff_in_block).segment(obs_atom_bfoff, obs_atom_nbf);
              // Add back in the nu.center == Y.center part
              Kt_part.row(mu).segment(obs_atom_bfoff, obs_atom_nbf).transpose() += 0.5
                  * C_Y.middleCols(obs_atom_bfoff, obs_atom_nbf).transpose()
                  * dt_mus[mu.bfoff_in_shell].col(Y.bfoff_in_block).segment(obs_atom_bfoff, obs_atom_nbf);
              //----------------------------------------//
            }
          }
          //----------------------------------------//
        } // end while get pair
        //----------------------------------------//
        // Sum Kt parts within node
        boost::lock_guard<boost::mutex> lg(tmp_mutex);
        Kt += Kt_part;
        Kt2 += Kt_part;
      }); // end create_thread
    } // end enumeration of threads
    compute_threads.join_all();
  } // compute_threads is destroyed here
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Global sum K                                         		                        {{{1 */ #if 1 // begin fold
  //----------------------------------------//
  msg.sum(Kt.data(), nbf*nbf);
  msg.sum(Kt2.data(), nbf*nbf);
  //----------------------------------------//
  // Symmetrize K
  Eigen::MatrixXd K(nbf, nbf);
  K = Kt + Kt.transpose();
  K2 = Kt2 + Kt2.transpose();
  if(xml_debug) write_as_xml("K2", K2);
  //----------------------------------------//
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Transfer K to a RefSCMatrix                           		                        {{{1 */ #if 1 // begin fold
  Ref<Integral> localints = integral()->clone();
  localints->set_basis(obs);
  Ref<PetiteList> pl = localints->petite_list();
  RefSCDimension obsdim = pl->AO_basisdim();
  RefSCMatrix result(
      obsdim,
      obsdim,
      obs->so_matrixkit()
  );
  result.assign(K.data());
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Clean up                                             		                        {{{1 */ #if 1 // begin fold
  //----------------------------------------//
  deallocate(D_ptr);
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  return result;
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
  const int nthread = threadgrp_->nthread();
  // Loop over number of threads
  for(int ithr = 0; ithr < nthread; ++ithr) {
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

/*
void
CADFCLHF::loop_shell_pairs_threaded(
    const std::function<void(int, const ShellData&, const ShellData&)>& f
)
{
  loop_shell_pairs_threaded(AllPairs, f);
}
*/

bool
CADFCLHF::get_shell_pair(ShellData& mu, ShellData& nu, PairSet pset)
{
  // Atomicly access and increment
  int spot = local_pairs_spot_++;
  if(spot < local_pairs_[pset].size()) {
    IntPair& next_pair = local_pairs_[pset][spot];
    //----------------------------------------//
    if(dynamic_) {
      // Here's where we'd need to check if we're running low on pairs and prefetch some more
      // When implemented, this should use a std::async or something like that
      throw FeatureNotImplemented("dynamic load balancing", __FILE__, __LINE__, class_desc());
    }
    //----------------------------------------//
    mu = ShellData(next_pair.first, basis().pointer(), dfbs_.pointer());
    nu = ShellData(next_pair.second, basis().pointer(), dfbs_.pointer());
  }
  else{
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////////

void
CADFCLHF::compute_coefficients()
{
  /*=======================================================================================*/
  /* Setup                                                		                        {{{1 */ #if 1 // begin fold
  //---------------------------------------------------------------------------------------//
  // References for speed
  Timer timer("compute coefficients");
  const Ref<GaussianBasisSet>& obs = basis();
  // Constants for convenience
  const int nbf = obs->nbasis();
  const int dfnbf = dfbs_->nbasis();
  const int natom = obs->ncenter();
  /*-----------------------------------------------------*/
  /* Initialize coefficient memory                  {{{2 */ #if 2 // begin fold
  // Now initialize the coefficient memory.
  // First, compute the amount of memory needed
  // Coefficients will be stored jbf <= ibf
  int ncoefs = 0;
  for(auto ibf : function_range(obs, dfbs_)){
    for(auto jbf : function_range(obs, dfbs_, 0, ibf)){
      ncoefs += ibf.atom_dfnbf;
      if(ibf.center != jbf.center){
        ncoefs += jbf.atom_dfnbf;
      }
    }
  }
  coefficients_data_ = allocate<double>(ncoefs);
  memset(coefficients_data_, 0, ncoefs * sizeof(double));
  double *spot = coefficients_data_;
  for(auto ibf : function_range(obs, dfbs_)){
    for(auto jbf : function_range(obs, dfbs_, 0, ibf)){
      double *data_spot_a = spot;
      spot += ibf.atom_dfnbf;
      double *data_spot_b = spot;
      if(ibf.center != jbf.center){
        spot += jbf.atom_dfnbf;
      }
      IntPair ij(ibf, jbf);
      CoefContainer coefs_a = make_shared<Eigen::Map<Eigen::VectorXd>>(data_spot_a, ibf.atom_dfnbf);
      CoefContainer coefs_b = make_shared<Eigen::Map<Eigen::VectorXd>>(
          data_spot_b, ibf.center == jbf.center ? 0 : int(jbf.atom_dfnbf));
      coefs_.emplace(ij, std::make_pair(coefs_a, coefs_b));
    }
  }
  //----------------------------------------//
  // We can save a lot of mess by storing references
  //   to the jbf > ibf pairs for the ish == jsh cases only
  for(auto ish : shell_range(obs)){
    for(auto ibf : function_range(ish)){
      for(auto jbf : function_range(obs, ibf + 1, ish.last_function)){
        IntPair ij(ibf, jbf);
        IntPair ji(jbf, ibf);
        // This will do a copy, but not of the data, just of the map
        //   to the data, which is okay
        coefs_.emplace(ij, coefs_[ji]);
      }
    }
  }
  //----------------------------------------//
  // Now make the transposes, for more efficient use later
  coefs_transpose_.resize(dfnbf);
  for(auto Y : function_range(dfbs_)){
    coefs_transpose_[Y].resize(obs->nbasis_on_center(Y.center), nbf);
  }
  /********************************************************/ #endif //2}}}
  /*-----------------------------------------------------*/
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Compute the coefficients in threads                   		                        {{{1 */ #if 1 // begin fold
  //---------------------------------------------------------------------------------------//
  // reset the iteration over local pairs
  local_pairs_spot_ = 0;
  const int nthread = threadgrp_->nthread();
  {
    boost::thread_group compute_threads;
    // Loop over number of threads
    for(int ithr = 0; ithr < nthread; ++ithr) {
      // ...and create each thread that computes pairs
      compute_threads.create_thread([&,ithr](){
        /*-----------------------------------------------------*/
        /* Coefficient compute thread                     {{{2 */ #if 2 // begin fold
        ShellData ish, jsh;
        while(get_shell_pair(ish, jsh)){
          const int dfnbfAB = ish.center == jsh.center ? ish.atom_dfnbf : ish.atom_dfnbf + jsh.atom_dfnbf;
          //----------------------------------------//
          std::shared_ptr<Decomposition> decomp = get_decomposition(
              ish, jsh, metric_ints_2c_[ithr]
          );
          //----------------------------------------//
          Eigen::MatrixXd ij_M_X(ish.nbf*jsh.nbf, dfnbfAB);
          for(auto ksh : iter_shells_on_center(dfbs_, ish.center)){
            auto ij_M_k = ints_to_eigen(
                ish, jsh, ksh,
                metric_ints_3c_[ithr],
                metric_oper_type_
            );
            for(auto ibf : function_range(ish)){
              for(auto jbf : function_range(jsh)){
                const int ijbf = ibf.bfoff_in_shell * jsh.nbf + jbf.bfoff_in_shell;
                ij_M_X.row(ijbf).segment(ksh.bfoff_in_atom, ksh.nbf) = ij_M_k->row(ijbf);
              } // end loop over functions in jsh
            } // end loop over functions in ish
          } // end loop over shells on ish.center
          if(ish.center != jsh.center){
            for(auto ksh : iter_shells_on_center(dfbs_, jsh.center)){
              auto ij_M_k = ints_to_eigen(
                  ish, jsh, ksh,
                  metric_ints_3c_[ithr],
                  metric_oper_type_
              );
              for(auto ibf : function_range(ish)){
                for(auto jbf : function_range(jsh)){
                  const int ijbf = ibf.bfoff_in_shell * jsh.nbf + jbf.bfoff_in_shell;
                  const int dfbfoff = ish.atom_dfnbf + ksh.bfoff_in_atom;
                  ij_M_X.row(ijbf).segment(dfbfoff, ksh.nbf) = ij_M_k->row(ijbf);
                } // end loop over functions in jsh
              } // end loop over functions in ish
            } // end loop over shells on jsh.center
          } // end if ish.center != jsh.center
          //----------------------------------------//
          for(auto mu : function_range(ish)){
            // Since coefficients are only stored jbf <= ibf,
            //   we need to figure out how many jbf's to do
            for(auto nu : function_range(jsh)){
              if(ish == jsh and nu > mu) continue;
              const int ijbf = mu.bfoff_in_shell*jsh.nbf + nu.bfoff_in_shell;
              IntPair mu_nu(mu, nu);
              Eigen::VectorXd Ctmp(dfnbfAB);
              Ctmp = decomp->solve(ij_M_X.row(ijbf).transpose());
              assert((coefs_.find(mu_nu) != coefs_.end())
                  || ((cout << "couldn't find coefficients for " << mu << ", " << nu << endl), false));
              //    || ((cout << "couldn't find coefficients for " << mu << ", " << nu << endl ), false);
              *(coefs_[mu_nu].first) = Ctmp.head(ish.atom_dfnbf);
              if(ish.center != jsh.center){
                *(coefs_[mu_nu].second) = Ctmp.tail(jsh.atom_dfnbf);
              }
            } // end jbf loop
          } // end ibf loop
        } // end while get_shell_pair
        /********************************************************/ #endif //2}}}
        /*-----------------------------------------------------*/
      });  // End threaded compute function and compute_threads.create_thread() call
    } // end loop over number of threads
    compute_threads.join_all();
  } // thread_group compute_threads is destroyed and goes out of scope, threads are destroyed (?)
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Global sum coefficient memory                        		                        {{{1 */ #if 1 // begin fold
  //---------------------------------------------------------------------------------------//
  scf_grp_->sum(coefficients_data_, ncoefs);
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Store the transpose coefficients                     		                        {{{1 */ #if 1 // begin fold
  //---------------------------------------------------------------------------------------//
  for(auto Y : function_range(dfbs_)){
    for(auto ish : iter_shells_on_center(obs, Y.center)){
      for(auto mu : function_range(ish)){
        for(auto nu : function_range(obs)){
          //----------------------------------------//
          if(nu <= mu){
            auto& C = *(coefs_[IntPair(mu, nu)].first);
            coefs_transpose_[Y](mu.bfoff_in_atom, nu) = C[Y.bfoff_in_atom];
          }
          else{
            auto cpair = coefs_[IntPair(nu, mu)];
            if(mu.center != nu.center){
              auto& C = *(cpair.second);
              coefs_transpose_[Y](mu.bfoff_in_atom, nu) = C[Y.bfoff_in_atom];
            }
            else{
              auto& C = *(cpair.first);
              coefs_transpose_[Y](mu.bfoff_in_atom, nu) = C[Y.bfoff_in_atom];
            }
          }
          //----------------------------------------//
        }
      }
    }
  }
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Debugging output                                     		                        {{{1 */ #if 1 // begin fold
  if(xml_debug) {
    begin_xml_context(
        "df_coefficients",
        "coefficients.xml"
    );
    for(auto mu : function_range(obs, dfbs_)){
      for(auto nu : function_range(obs, dfbs_, mu)) {
        IntPair mn(mu, nu);
        auto cpair = coefs_[mn];
        auto& Ca = *(cpair.first);
        auto& Cb = *(cpair.second);
        //----------------------------------------//
        // Fill in the zeros for easy comparison
        Eigen::VectorXd coefs(dfnbf);
        coefs = Eigen::VectorXd::Zero(dfnbf);
        coefs.segment(mu.atom_dfbfoff, mu.atom_dfnbf) = Ca;
        if(mu.center != nu.center){
          coefs.segment(nu.atom_dfbfoff, nu.atom_dfnbf) = Cb;
        }
        write_as_xml(
            "coefficient_vector",
            coefs,
            std::map<std::string, int>{
              { "ao_index1", mu },
              { "ao_index2", nu }
            }
        );
      }
    }
    end_xml_context("df_coefficients");
  }
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Make the CADF-LinK lists                                                         {{{1 */ #if 1 // begin fold
  schwarz_df_.resize(dfbs_->nshell());
  for(auto Xsh : shell_range(dfbs_)){
    auto g_XX_ptr = ints_to_eigen(
        Xsh, Xsh, eris_2c_[0], coulomb_oper_type_
    );
    auto& g_XX = *g_XX_ptr;
    double frob_norm = 0.0;
    for(auto X : function_range(Xsh)) {
      frob_norm += g_XX(X.bfoff_in_shell, X.bfoff_in_shell);
    }
    schwarz_df_[Xsh] = sqrt(frob_norm);
  }
  if(do_linK_) {
    for(auto lsh : shell_range(obs)) {
      out_assert(L_schwarz[lsh].size(), >, 0);
      for(auto ksh : L_schwarz[lsh]) {
        double Ct = 0.0;
        //----------------------------------------//
        for(auto sigma : function_range(lsh)) {
          for(auto nu : function_range(ksh)) {
            BasisFunctionData first, second;
            if(ksh <= lsh) {
              first = sigma;
              second = nu;
            }
            else {
              first = nu;
              second = sigma;
            }
            IntPair mu_sigma(first, second);
            assert(coefs_.find(mu_sigma) != coefs_.end());
            auto& cpair = coefs_[mu_sigma];
            for(auto Xsh : iter_shells_on_center(dfbs_, first.center)){
              auto& C = *cpair.first;
              auto g_XX_ptr = ints_to_eigen(
                  Xsh, Xsh, eris_2c_[0], coulomb_oper_type_
              );
              auto& g_XX = *g_XX_ptr;
              for(auto X : function_range(Xsh)) {

                Ct += C[X.bfoff_in_atom] * C[X.bfoff_in_atom] * g_XX(X.bfoff_in_shell, X.bfoff_in_shell);
              }
            }
            if(first.center != second.center) {
              for(auto Xsh : iter_shells_on_center(dfbs_, second.center)){
                auto& C = *cpair.second;
                auto g_XX_ptr = ints_to_eigen(
                    Xsh, Xsh, eris_2c_[0], coulomb_oper_type_
                );
                auto& g_XX = *g_XX_ptr;
                for(auto X : function_range(Xsh)) {
                  Ct += C[X.bfoff_in_atom] * C[X.bfoff_in_atom] * g_XX(X.bfoff_in_shell, X.bfoff_in_shell);
                }
              }
            }
          } // end loop over nu
        } // end loop over sigma
        //----------------------------------------//
        L_coefs[lsh].insert(ksh, sqrt(Ct));
        //----------------------------------------//
      } // end loop over ksh
    } // end loop over lsh
    #if !LINK_SORTED_INSERTION
    do_threaded(nthread, [&](int ithr){
      auto L_coefs_iter = L_coefs.begin();
      const auto& L_coefs_end = L_coefs.end();
      L_coefs_iter.advance(ithr);
      while(L_coefs_iter != L_coefs_end) {
        L_coefs_iter->second.sort();
        L_coefs_iter.advance(nthread);
      }
    });
    #endif
    //----------------------------------------//
    // Compute the Frobenius norm of C_transpose_ blocks
    C_trans_frob_.resize(dfbs_->nshell());
    for(auto Xsh : shell_range(dfbs_, obs)) {
      // obs is being used as the aux basis, so atom_dfnsh is the number
      //   of shells on the same center in the orbital basis
      C_trans_frob_[Xsh].resize(Xsh.atom_dfnsh, obs->nshell());
      for(auto ish : iter_shells_on_center(obs, Xsh.center)) {
        for(auto jsh : shell_range(obs)) {
          double frob_val = 0.0;
          for(auto X : function_range(Xsh)) {
            for(auto mu : function_range(ish)) {
              for(auto rho : function_range(jsh)) {
                const double cval = coefs_transpose_[X](mu.bfoff_in_atom, rho);
                frob_val += cval * cval;
              }
            }
          } // end loop over X in Xsh
          C_trans_frob_[Xsh](ish.shoff_in_atom, jsh) = sqrt(frob_val);
        } // end loop over jsh
      } // end loop over ish
    } // end loop over Xsh
  } // end if do_linK
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Clean up                                              		                        {{{1 */ #if 1 // begin fold
  //---------------------------------------------------------------------------------------//
  have_coefficients_ = true;
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
}

//////////////////////////////////////////////////////////////////////////////////

std::shared_ptr<Decomposition>
CADFCLHF::get_decomposition(int ish, int jsh, Ref<TwoBodyTwoCenterInt> ints)
{
  // TODO atom swap symmetry
  const int atomA = basis()->shell_to_center(ish);
  const int atomB = basis()->shell_to_center(jsh);
  //----------------------------------------//
  return decomps_.get(atomA, atomB, [&](){
    // Make the decomposition
    std::shared_ptr<Decomposition> decompAB;
    //----------------------------------------//
    const int dfnshA = dfbs_->nshell_on_center(atomA);
    const int dfnbfA = dfbs_->nbasis_on_center(atomA);
    const int dfshoffA = dfbs_->shell_on_center(atomA, 0);
    const int dfbfoffA = dfbs_->shell_to_function(dfshoffA);
    //----------------------------------------//
    const int dfnshB = dfbs_->nshell_on_center(atomB);
    const int dfnbfB = dfbs_->nbasis_on_center(atomB);
    const int dfshoffB = dfbs_->shell_on_center(atomB, 0);
    const int dfbfoffB = dfbs_->shell_to_function(dfshoffB);
    //----------------------------------------//
    // Compute the integrals we need
    const int nrows = atomA == atomB ? dfnbfA : dfnbfA + dfnbfB;
    Eigen::MatrixXd g2AB(nrows, nrows);
    // AA integrals
    for(int ishA = dfshoffA; ishA < dfshoffA + dfnshA; ++ishA){
      const int dfbfoffiA = dfbs_->shell_to_function(ishA);
      const int dfnbfiA = dfbs_->shell(ishA).nfunction();
      for(int jshA = dfshoffA; jshA < dfshoffA + dfnshA; ++jshA){
        const int dfbfoffjA = dfbs_->shell_to_function(jshA);
        const int dfnbfjA = dfbs_->shell(jshA).nfunction();
        auto shell_ints = ints_to_eigen(
            ishA, jshA, ints,
            metric_oper_type_
        );
        g2AB.block(
            dfbfoffiA - dfbfoffA, dfbfoffjA - dfbfoffA,
            dfnbfiA, dfnbfjA
        ) = *shell_ints;
        g2AB.block(
            dfbfoffjA - dfbfoffA, dfbfoffiA - dfbfoffA,
            dfnbfjA, dfnbfiA
        ) = shell_ints->transpose();
      }
    }
    if(atomA != atomB) {
      // AB integrals
      for(int ishA = dfshoffA; ishA < dfshoffA + dfnshA; ++ishA){
        const int dfbfoffiA = dfbs_->shell_to_function(ishA);
        const int dfnbfiA = dfbs_->shell(ishA).nfunction();
        for(int jshB = dfshoffB; jshB < dfshoffB + dfnshB; ++jshB){
          const int dfbfoffjB = dfbs_->shell_to_function(jshB);
          const int dfnbfjB = dfbs_->shell(jshB).nfunction();
          auto shell_ints = ints_to_eigen(
              ishA, jshB, ints,
              metric_oper_type_
          );
          g2AB.block(
              dfbfoffiA - dfbfoffA, dfbfoffjB - dfbfoffB + dfnbfA,
              dfnbfiA, dfnbfjB
          ) = *shell_ints;
          g2AB.block(
              dfbfoffjB - dfbfoffB + dfnbfA, dfbfoffiA - dfbfoffA,
              dfnbfjB, dfnbfiA
          ) = shell_ints->transpose();
        }
      }
      // BB integrals
      for(int ishB = dfshoffB; ishB < dfshoffB + dfnshB; ++ishB){
        const int dfbfoffiB = dfbs_->shell_to_function(ishB);
        const int dfnbfiB = dfbs_->shell(ishB).nfunction();
        for(int jshB = dfshoffB; jshB < dfshoffB + dfnshB; ++jshB){
          const int dfbfoffjB = dfbs_->shell_to_function(jshB);
          const int dfnbfjB = dfbs_->shell(jshB).nfunction();
          auto shell_ints = ints_to_eigen(
              ishB, jshB, ints,
              metric_oper_type_
          );
          g2AB.block(
              dfbfoffiB - dfbfoffB + dfnbfA, dfbfoffjB - dfbfoffB + dfnbfA,
              dfnbfiB, dfnbfjB
          ) = *shell_ints;
          g2AB.block(
              dfbfoffjB - dfbfoffB + dfnbfA, dfbfoffiB - dfbfoffB + dfnbfA,
              dfnbfjB, dfnbfiB
          ) = shell_ints->transpose();
        }
      }
    }
    return make_shared<Decomposition>(g2AB);
  });
}

//////////////////////////////////////////////////////////////////////////////////

CADFCLHF::TwoCenterIntContainerPtr
CADFCLHF::ints_to_eigen(int ish, int jsh, Ref<TwoBodyTwoCenterInt>& ints, TwoBodyOper::type int_type){
#if USE_INTEGRAL_CACHE
  return ints2_.get(ish, jsh, int_type, [ish, jsh, int_type, &ints, this](){
#endif
    GaussianBasisSet& ibs = *(ints->basis1());
    GaussianBasisSet& jbs = *(ints->basis2());
    const int nbfi = ibs.shell(ish).nfunction();
    const int nbfj = jbs.shell(jsh).nfunction();
    auto rv = make_shared<TwoCenterIntContainer>(nbfi, nbfj);
    //----------------------------------------//
    ints->compute_shell(ish, jsh);
    const double* buffer = ints->buffer(int_type);
    //::memcpy(rv->data(), buffer, nbfi*nbfj * sizeof(double));
    std::copy(buffer, buffer + nbfi*nbfj, rv->data());
    //std::move(buffer, buffer + nbfi*nbfj, rv->data());
    //----------------------------------------//
    ints_computed_locally_ += nbfi * nbfj;
    return rv;
#if USE_INTEGRAL_CACHE
  });
#endif
}

//////////////////////////////////////////////////////////////////////////////////

CADFCLHF::ThreeCenterIntContainerPtr
CADFCLHF::ints_to_eigen(int ish, int jsh, int ksh, Ref<TwoBodyThreeCenterInt>& ints, TwoBodyOper::type int_type){
  #if USE_INTEGRAL_CACHE
  return ints3_.get(ish, jsh, ksh, int_type, [ish, jsh, ksh, int_type, &ints, this](){
  #endif
    const int nbfi = ints->basis1()->shell(ish).nfunction();
    const int nbfj = ints->basis2()->shell(jsh).nfunction();
    const int nbfk = ints->basis3()->shell(ksh).nfunction();
    auto rv = make_shared<ThreeCenterIntContainer>(nbfi * nbfj, nbfk);
    //----------------------------------------//
    ints->compute_shell(ish, jsh, ksh);
    const double* buffer = ints->buffer(int_type);
    //::memcpy(rv->data(), buffer, nbfi*nbfj*nbfk * sizeof(double));
    std::copy(buffer, buffer + nbfi*nbfj*nbfk, rv->data());
    //std::move(buffer, buffer + nbfi*nbfj*nbfk, rv->data());
    //----------------------------------------//
    ints_computed_locally_ += nbfi * nbfj * nbfk;
    return rv;
#if USE_INTEGRAL_CACHE
  });
#endif
}

