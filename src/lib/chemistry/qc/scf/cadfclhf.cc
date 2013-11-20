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
#include <chemistry/qc/scf/cadfclhf.h>
#include <util/misc/scexception.h>

using namespace sc;
using namespace std;

typedef std::pair<int, int> IntPair;
typedef CADFCLHF::CoefContainer CoefContainer;
typedef CADFCLHF::IntContainer2 IntContainer2;
typedef CADFCLHF::IntContainer3 IntContainer3;
typedef CADFCLHF::Decomposition Decomposition;
typedef std::pair<CoefContainer, CoefContainer> CoefPair;

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
    local_pairs_spot_(0),
    is_fetching_pairs_(false)
{
  //----------------------------------------------------------------------------//
  // Get the auxiliary basis set
  dfbs_ << keyval->describedclassvalue("df_basis", KeyValValueRefDescribedClass(0));
  if(dfbs_.null()){
    throw InputError("CADFCLHF requires a density fitting basis set",
        __FILE__, __LINE__, "df_basis");
  }
  //----------------------------------------------------------------------------//
  initialize();
  //----------------------------------------------------------------------------//
}

//////////////////////////////////////////////////////////////////////////////////

CADFCLHF::~CADFCLHF()
{
  if(have_coefficients_){
    // Cleanup from init_threads(), but only if we called it
    SCF::done_threads();
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
}

//////////////////////////////////////////////////////////////////////////////////

void
CADFCLHF::init_threads()
{
  const int nthread = threadgrp_->nthread();
  //----------------------------------------------------------------------------//
  // If we're doing static load balancing, set up pair assignments here
  if(not dynamic_){
    const int n_node = scf_grp_->n();
    const int nshell = basis()->nshell();
    for(int ish=0, inode=0; ish < nshell; ++ish){
      for(int jsh=0; jsh <= ish; ++jsh, ++inode){
        IntPair ij(ish, jsh);
        pair_assignments_[ij] = inode % n_node;
      }
    }
    // Make the backwards mapping for the current node
    int me = scf_grp_->me();
    for(auto it : pair_assignments_){
      if(it.second == me) local_pairs_.push_back(it.first);
    }
  }
  //----------------------------------------------------------------------------//
  // initialize the two electron integral classes
  // ThreeCenter versions
  integral()->set_basis(basis(), basis(), dfbs_);
  for (int i=0; i < nthread; i++) {
    eris_3c_.push_back(integral()->coulomb<3>());
    // TODO different fitting metrics
    metric_ints_3c_.push_back(eris_3c_.back());
  }
  // TwoCenter versions
  integral()->set_basis(dfbs_, dfbs_);
  for (int i=0; i < nthread; i++) {
    eris_2c_.push_back(integral()->coulomb<2>());
    metric_ints_2c_.push_back(eris_2c_.back());
  }
  // Reset to normal setup
  integral()->set_basis(basis(), basis(), basis(), basis());
  //----------------------------------------------------------------------------//
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
  if(not have_coefficients_) {
    Timer coeff_tim("compute coefficients");
    compute_coefficients();
  } // Timer coeff_tim destroyed, causing it to stop
  //---------------------------------------------------------------------------------------//
  Timer routine_tim("ao_fock");
  Timer step_tim("misc");
  //---------------------------------------------------------------------------------------//
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
  RefSymmSCMatrix D = cl_dens_diff_.copy(); D.scale(0.5);
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Form G                                                		                        {{{1 */ #if 1 // begin fold
  //---------------------------------------------------------------------------------------//
  // compute J and K
  step_tim.change("build");
  RefSCMatrix G;
  {
    RefSCMatrix J = compute_J();
    G = J.copy();
  }
  {
    RefSCMatrix K = compute_K();
    G.accumulate( -1.0 * K);
  }
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
  step_tim.change("misc");
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
  int nthread = threadgrp_->nthread();
  // reset the iteration over local pairs
  local_pairs_spot_ = 0;

}

//////////////////////////////////////////////////////////////////////////////////

RefSCMatrix
CADFCLHF::compute_K()
{
  int nthread = threadgrp_->nthread();
  // reset the iteration over local pairs
  local_pairs_spot_ = 0;
}

//////////////////////////////////////////////////////////////////////////////////

bool
CADFCLHF::get_shell_pair(int& mu, int& nu)
{
  // Atomicly access and increment
  int spot = local_pairs_spot_++;
  if(spot < local_pairs_.size()) {
    IntPair& next_pair = local_pairs_[spot];
    //----------------------------------------//
    if(dynamic_) {
      // Here's where we'd need to check if we're running low on pairs and prefetch some more
      // When implemented, this should use a std::async or something like that
      throw FeatureNotImplemented("dynamic load balancing", __FILE__, __LINE__, class_desc());
    }
    //----------------------------------------//
    mu = next_pair.first;
    nu = next_pair.second;
  }
  else{
    mu = NoMorePairs;
    nu = NoMorePairs;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////////

void
CADFCLHF::compute_coefficients()
{
  //----------------------------------------//
  // setup threads here, since this only
  //   gets called once per SCF computation
  init_threads();
  //----------------------------------------//
  // References for speed
  GaussianBasisSet& obs = *(basis());
  GaussianBasisSet& dfbs = *(dfbs_);
  // Constants for convenience
  const int nbf = obs.nbasis();
  const int dfnbf = dfbs.nbasis();
  const int natom = obs.ncenter();
  //----------------------------------------//
  // Now initialize the coefficient memory.
  // First, compute the amount of memory needed
  // Coefficients will be stored jbf <= ibf
  int ncoefs = 0;
  for(int ibf = 0; ibf < nbf; ++ibf){
    const int ishA = obs.function_to_shell(ibf);
    const int atomA = obs.shell_to_center(ishA);
    const int dfnbfA = dfbs.nbasis_on_center(atomA);
    for(int jbf = 0; jbf <= ibf; ++jbf){
      const int jshB = obs.function_to_shell(jbf);
      const int atomB = obs.shell_to_center(jshB);
      const int dfnbfB = dfbs.nbasis_on_center(atomB);
      ncoefs += dfnbfA;
      if(atomA != atomB){
        ncoefs += dfnbfB;
      }
    }
  }
  coefficients_data_ = allocate<double>(ncoefs);
  double *spot = coefficients_data_;
  for(int ibf = 0; ibf < nbf; ++ibf){
    const int ishA = obs.function_to_shell(ibf);
    const int atomA = obs.shell_to_center(ishA);
    const int dfnbfA = dfbs.nbasis_on_center(atomA);
    for(int jbf = 0; jbf <= ibf; ++jbf){
      const int jshB = obs.function_to_shell(jbf);
      const int atomB = obs.shell_to_center(jshB);
      const int dfnbfB = dfbs.nbasis_on_center(atomB);
      double *data_spot_a = spot;
      IntPair ij(ibf, jbf);
      spot += dfnbfA;
      if(atomA != atomB){
        double *data_spot_b = spot;
        spot += dfnbfB;
        // This is unreadable...sorry
        coefs_.emplace(ij, std::pair<Eigen::Map<Eigen::VectorXd>, Eigen::Map<Eigen::VectorXd>>(
            Eigen::Map<Eigen::VectorXd>(data_spot_a, dfnbfA),
            Eigen::Map<Eigen::VectorXd>(data_spot_b, dfnbfB)
        ));
      }
      else{
        double *data_spot_b = 0;
        coefs_.emplace(ij, std::pair<Eigen::Map<Eigen::VectorXd>, Eigen::Map<Eigen::VectorXd>>(
            Eigen::Map<Eigen::VectorXd>(data_spot_a, dfnbfA),
            Eigen::Map<Eigen::VectorXd>(data_spot_b, 0)
        ));
      }
    }
  }
  //----------------------------------------//
  // reset the iteration over local pairs
  local_pairs_spot_ = 0;
  //----------------------------------------//
  have_coefficients_ = true;
}

void
detail::compute_coefs_task::operator()(){
  //----------------------------------------//
  // References for speed
  GaussianBasisSet& obs = *(wfn_->basis());
  GaussianBasisSet& dfbs = *(wfn_->dfbs_);
  // Constants for convenience
  const int nbf = obs.nbasis();
  const int dfnbf = dfbs.nbasis();
  const int natom = obs.ncenter();
  int ish = 0;
  int jsh = 0;
  //----------------------------------------//
  while(wfn_->get_shell_pair(ish, jsh)){
    const IntPair ij(ish, jsh);
    //----------------------------------------//
    const int atomA = obs.shell_to_center(ish);
    const int nbfA = obs.nbasis_on_center(atomA);
    const int dfnbfA = dfbs.nbasis_on_center(atomA);
    //----------------------------------------//
    const int atomB = obs.shell_to_center(jsh);
    const int nbfB = obs.nbasis_on_center(atomB);
    const int dfnbfB = dfbs.nbasis_on_center(atomB);
    //----------------------------------------//
    std::shared_ptr<Decomposition> decomp =
        wfn_->get_decomposition(ish, jsh, wfn_->metric_ints_2c_[ithr_]);
  }
}

std::shared_ptr<Decomposition>
CADFCLHF::get_decomposition(int ish, int jsh, Ref<TwoBodyTwoCenterInt> ints)
{
  const int atomA = basis()->shell_to_center(ish);
  const int atomB = basis()->shell_to_center(jsh);
  assert(atomA <= atomB); // for now; could also handle the opposite case and rearrange
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
    Eigen::MatrixXd g2AB;
    for(int ishA = dfshoffA; ishA < dfshoffA + dfnshA; ++ishA){
      const int dfbfoffiA = dfbs_->shell_to_function(ishA);
      const int dfnbfiA = dfbs_->shell(ishA).nfunction();
      for(int jshB = dfshoffA; jshB <= dfshoffA + dfnshB; ++jshB){
        const int dfbfoffjB = dfbs_->shell_to_function(jshB);
        const int dfnbfjB = dfbs_->shell(jshB).nfunction();
        std::shared_ptr<Eigen::MatrixXd> shell_ints = ints_to_eigen(ishA, jshB, ints);
        g2AB.block(
            dfbfoffiA - dfbfoffA, dfbfoffjB - dfbfoffB,
            dfnbfiA, dfnbfjB
        ) = *shell_ints;
      }
    }
    return std::shared_ptr<Decomposition>(new Decomposition(g2AB));
  });
}

std::shared_ptr<Eigen::MatrixXd>
CADFCLHF::ints_to_eigen(int ish, int jsh, Ref<TwoBodyTwoCenterInt> ints){
  // TODO permutational symmetry; keep in mind that a transpose is needed
  return ints2_.get(ish, jsh, [&]{
    const int dfnbfi = dfbs_->shell(ish).nfunction();
    const int dfnbfj = dfbs_->shell(jsh).nfunction();
    std::shared_ptr<Eigen::MatrixXd> rv(
        new Eigen::MatrixXd(dfnbfi, dfnbfj)
    );
    ints->compute_shell(dfnbfi, dfnbfj);
    const double* buffer = ints->buffer();
    // TODO vectorized copy
    int buffer_spot = 0;
    for(int ibf = 0; ibf < dfnbfi; ++ibf){
      for(int jbf = 0; jbf < dfnbfj; ++jbf){
        (*rv)(ibf, jbf) = buffer[buffer_spot++];
      }
    }
    return rv;
  });
}


