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

#ifdef ECLIPSE_PARSER_ONLY
#include <util/misc/sharedptr.h>
#endif

using namespace sc;
using namespace std;
using boost::thread;
using boost::thread_group;
using namespace sc::parameter;
using std::make_shared;
static constexpr bool xml_debug = true;

typedef std::pair<int, int> IntPair;
typedef CADFCLHF::CoefContainer CoefContainer;
typedef CADFCLHF::Decomposition Decomposition;
typedef std::pair<CoefContainer, CoefContainer> CoefPair;

static boost::mutex debug_print_mutex;

////////////////////////////////////////////////////////////////////////////////
// Debugging asserts and outputs

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

#define DECOMP_PRINT(D, M) \
  {\
    boost::lock_guard<boost::mutex> _tmp(debug_print_mutex); \
    cout << "Decomposition:  " << D.matrixQR().rows() << " x " << D.matrixQR().cols() << ", M => " << M.rows() << " x " << M.cols() << endl; \
  }

#define M_EQ_ASSERT(M1, M2) M_ROW_ASSERT(M1, M2); M_COL_ASSERT(M1, M2);

#define M_BLOCK_ASSERT(M, b1a, b1b, b2a, b2b) \
  if(b1a < 0) { \
    boost::lock_guard<boost::mutex> _tmp(debug_print_mutex); \
    cout << "assertion 1 failed.  data: " << boost::make_tuple(M.rows(), M.cols(), b1a, b1b, b2a, b2b) << endl; \
  } \
  else if(b1b < 0) { \
    boost::lock_guard<boost::mutex> _tmp(debug_print_mutex); \
    cout << "assertion 2 failed.  data: " << boost::make_tuple(M.rows(), M.cols(), b1a, b1b, b2a, b2b) << endl; \
  } \
  else if(b1a > M.rows() - b2a) { \
    boost::lock_guard<boost::mutex> _tmp(debug_print_mutex); \
    cout << "assertion 3 failed.  data: " << boost::make_tuple(M.rows(), M.cols(), b1a, b1b, b2a, b2b) << endl; \
  } \
  else if(b1b > M.cols() - b2b) { \
    boost::lock_guard<boost::mutex> _tmp(debug_print_mutex); \
    cout << "assertion 4 failed.  data: " << boost::make_tuple(M.rows(), M.cols(), b1a, b1b, b2a, b2b) << endl; \
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
  // For now, use coulomb metric.  We can easily make this a keyword later
  metric_oper_type_ = TwoBodyOper::eri;
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
  have_coefficients_ = false;
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
    const int me = scf_grp_->me();
    const int dfnsh = dfbs_->nshell();
    for(auto it : pair_assignments_){
      if(it.second == me){
        local_pairs_.push_back(it.first);
        //----------------------------------------//
        const int ish = it.first.first;
        const int atomA = basis()->shell_to_center(ish);
        const int shoffA = basis()->shell_on_center(atomA, 0);
        const int dfshoffA = dfbs_->shell_on_center(atomA, 0);
        const int dfnshA = dfbs_->nshell_on_center(atomA);
        const int jsh = it.first.second;
        const int atomB = basis()->shell_to_center(jsh);
        const int shoffB = basis()->shell_on_center(atomB, 0);
        const int dfshoffB = dfbs_->shell_on_center(atomB, 0);
        const int dfnshB = dfbs_->nshell_on_center(atomB);
        //----------------------------------------//
        // Create mutexes corresponding to the pair
        // ints3_ mutexes
        // This avoids dynamic mutex creation which requires
        //   a global thread lock
        /*
        for(int kshdf = 0; kshdf < dfbs_->nshell(); ++kshdf){
          ints3_.init_mutex(it.first.first, it.first.second, kshdf, coulomb_oper_type_);
          if(coulomb_oper_type_ != metric_oper_type_){
            ints3_.init_mutex(it.first.first, it.first.second, kshdf, metric_oper_type_);
          }
        }
        //----------------------------------------//
        if(coulomb_oper_type_ != metric_oper_type_){
          // Here we only need locks for the pairs that will be involved in the
          //   computation of the coefficients
          for(int ishdf = dfshoffA; ishdf < dfshoffA + dfnshA; ++ishdf){
            for(int jshdf = dfshoffB; jshdf < dfshoffB + dfnshB; ++jshdf){
              ints2_.init_mutex(ishdf, jshdf, metric_oper_type_);
            }
          }
        }
        */
        //----------------------------------------//
        // Finally, one mutex for the decomposition
        //decomps_.init_mutex(atomA, atomB);
      }
    }
    // We need to be able to lock all of the coulomb ints2_ (for now)
    //for(int ishdf = 0; ishdf < dfnsh; ++ishdf){
    //  for(int jshdf = 0; jshdf < dfnsh; ++jshdf){
    //    ints2_.init_mutex(ishdf, jshdf, coulomb_oper_type_);
    //  }
    //}
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
  D_ = cl_dens_diff_.copy().convert2RefSCMat(); D_.scale(0.5);
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Form G                                                		                        {{{1 */ #if 1 // begin fold
  //---------------------------------------------------------------------------------------//
  // compute J and K
  step_tim.change("build");
  RefSCMatrix G;
  {
    if(xml_debug){
      begin_xml_context(
          "compute_J",
          "compute_J.xml"
      );
    }
    RefSCMatrix J = compute_J();
    if(xml_debug){
      write_as_xml("J", J);
      end_xml_context("compute_J");
    }
    G = J.copy();
    assert(false);
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
  /*=======================================================================================*/
  /* Setup                                                 		                        {{{1 */ #if 1 // begin fold
  //----------------------------------------//
  // Convenience variables
  MessageGrp& msg = *scf_grp_;
  const int nthread = threadgrp_->nthread();
  const int me = msg.me();
  const int n_node = msg.n();
  GaussianBasisSet& obs = *basis();
  GaussianBasisSet& dfbs = *dfbs_;
  const int nbf = obs.nbasis();
  const int dfnbf = dfbs.nbasis();
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
        while(get_shell_pair(ish, jsh)){
          //----------------------------------------//
          const int nbfi = obs.shell(ish).nfunction();
          const int bfoffi = obs.shell_to_function(ish);
          //----------------------------------------//
          const int nbfj = obs.shell(jsh).nfunction();
          const int bfoffj = obs.shell_to_function(jsh);
          //----------------------------------------//
          double perm_fact = ish == jsh ? 1.0 : 2.0;
          //----------------------------------------//
          for(int ibf = 0; ibf < nbfi; ++ibf){
            const int jmax = ish == jsh ? ibf : nbfj-1;
            double bf_perm_fact = ish == jsh ? 2.0 : 1.0;
            const int mu = bfoffi + ibf;
            for(int jbf = 0; jbf <= jmax; ++jbf){
              const int nu = bfoffj + jbf;
              if(ibf == jbf) bf_perm_fact = 1.0;
              IntPair ij(mu, nu);
              auto cpair = coefs_[ij];
              VectorMap& Ca = *(cpair.first);
              VectorMap& Cb = *(cpair.second);
              const double pf = perm_fact * bf_perm_fact;
              Ct.segment(ish.atom_dfbfoff, ish.atom_dfnbf) += pf * D(mu, nu) * Ca;
              if(ish.center != jsh.center) {
                Ct.segment(jsh.atom_dfbfoff, jsh.atom_dfnbf) += pf * D(mu, nu) * Cb;
              }
            } // end loop over jbf
          } // end loop over ibf
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
  msg.sum(C_tilde.data(), dfnbf);
  if(xml_debug) {
    write_as_xml("C_tilde", C_tilde);
  }
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Form g_tilde                                          		                        {{{1 */ #if 1 // begin fold
  //----------------------------------------//
  // TODO Thread this
  shared_ptr<Eigen::MatrixXd> X_g_Y = ints_to_eigen(
    boost::irange(0, (int)dfbs.nshell()),
    boost::irange(0, (int)dfbs.nshell()),
    eris_2c_[0],
    coulomb_oper_type_
  );
  Eigen::VectorXd g_tilde(dfnbf);
  g_tilde = (*X_g_Y) * C_tilde;
  if(xml_debug) {
    write_as_xml("g_tilde", g_tilde);
  }
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Form dtilde and add in Ctilde contribution to J       		                        {{{1 */ #if 1 // begin fold
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
        while(get_shell_pair(ish, jsh)){
          //----------------------------------------//
          const int bfoffi = obs.shell_to_function(ish);
          //----------------------------------------//
          const int bfoffj = obs.shell_to_function(jsh);
          //----------------------------------------//
          double perm_fact = ish == jsh ? 1.0 : 2.0;
          //----------------------------------------//
          //for(int kshdf = 0; kshdf < dfbs.nshell(); ++kshdf){
          for(auto kshdf : shell_iterator(&dfbs)){
            std::shared_ptr<Eigen::MatrixXd> g_part = ints_to_eigen(
                ish, jsh, kshdf,
                eris_3c_[ithr],
                coulomb_oper_type_
            );
            //for(int ibf = 0; ibf < ish.nbf; ++ibf){
            for(auto ibf : iter_functions_in_shell(&obs, ish)){
              // Don't need to do permutational symmetry
              //   within shells for this part
              //   since nothing here is stored with
              //   symmetric indices
              for(auto jbf : iter_functions_in_shell(&obs, jsh)){
              //for(int jbf = 0; jbf < jsh.nbf; ++jbf){
                const int ijbf = ibf.bfoff_in_shell*jsh.nbf + jbf.bfoff_in_shell;
                //----------------------------------------//
                // compute dtilde contribution
                dt.segment(kshdf.bfoff, kshdf.nbf) += perm_fact * D(ibf, jbf) * g_part->row(ijbf);
                //----------------------------------------//
                // add C_tilde * g_part to J
                if(jbf <= ibf){
                  // only constructing the lower triangle of the J matrix
                  jpart(ibf, jbf) += 0.5 * g_part->row(ijbf) * C_tilde.segment(kshdf.bfoff, kshdf.nbf);
                }
                //----------------------------------------//
              } // end loop over jbf
            } // end loop over ibf
          } // end loop over kshbf
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
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Add in first and third term contributions to J       		                        {{{1 */ #if 1 // begin fold
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
        /* first and third terms compute thread           {{{2 */ #if 2 // begin fold
        Eigen::MatrixXd jpart(nbf, nbf);
        jpart = Eigen::MatrixXd::Zero(nbf, nbf);
        //----------------------------------------//
        ShellData ish, jsh;
        while(get_shell_pair(ish, jsh)){
          //----------------------------------------//
          const int nbfi = obs.shell(ish).nfunction();
          const int bfoffi = obs.shell_to_function(ish);
          //----------------------------------------//
          const int nbfj = obs.shell(jsh).nfunction();
          const int bfoffj = obs.shell_to_function(jsh);
          //----------------------------------------//
          double perm_fact = (ish == jsh ? 0.5 : 1.0);
          //----------------------------------------//
          for(int ibf = 0; ibf < nbfi; ++ibf){
            const int jmax = ish == jsh ? ibf : nbfj-1;
            double bf_perm_fact = (ish == jsh ? 1.0 : 0.5);
            const int mu = bfoffi + ibf;
            for(int jbf = 0; jbf <= jmax; ++jbf){
              const int nu = bfoffj + jbf;
              if(mu == nu) bf_perm_fact = 1.0;
              const double pf = perm_fact * bf_perm_fact;
              const int ijbf = ibf*nbfj + jbf;
              IntPair ij(mu, nu);
              //----------------------------------------//
              auto cpair = coefs_[ij];
              VectorMap& Ca = *(cpair.first);
              VectorMap& Cb = *(cpair.second);
              //----------------------------------------//
              // First term contribution
              jpart(mu, nu) += pf * d_tilde.segment(ish.atom_dfbfoff, ish.atom_dfnbf).transpose() * Ca;
              if(ish.center != jsh.center){
                jpart(mu, nu) += pf * d_tilde.segment(jsh.atom_dfbfoff, jsh.atom_dfnbf).transpose() * Cb;
              }
              //----------------------------------------//
              // Third term contribution
              jpart(mu, nu) -= pf * g_tilde.segment(ish.atom_dfbfoff, ish.atom_dfnbf).transpose() * Ca;
              if(ish.center != jsh.center){
                jpart(mu, nu) -= pf * g_tilde.segment(jsh.atom_dfbfoff, jsh.atom_dfnbf).transpose() * Cb;
              }
              //----------------------------------------//

            } // end loop over jbf
          } // end loop over ibf
        } // end while get_shell_pair
        // Sum the thread's contributions to the node-level J
        boost::lock_guard<boost::mutex> tmp_lock(tmp_mutex);
        J += jpart;
        /*******************************************************/ #endif //2}}}
        /*-----------------------------------------------------*/
      }); // end create_thread
    } // end enumeration of threads
    compute_threads.join_all();
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
  localints->set_basis(&obs);
  Ref<PetiteList> pl = localints->petite_list();
  RefSCDimension obsdim = pl->AO_basisdim();
  RefSCMatrix result(
      obsdim,
      obsdim,
      obs.so_matrixkit()
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
  int nthread = threadgrp_->nthread();
  // reset the iteration over local pairs
  local_pairs_spot_ = 0;
  assert(false);
}

//////////////////////////////////////////////////////////////////////////////////

bool
CADFCLHF::get_shell_pair(ShellData& mu, ShellData& nu)
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
    mu = ShellData(next_pair.first, basis(), dfbs_);
    nu = ShellData(next_pair.second, basis(), dfbs_);
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
  // setup threads here, since this only
  //   gets called once per SCF computation
  init_threads();
  /*-----------------------------------------------------*/
  // References for speed
  GaussianBasisSet& obs = *(basis());
  GaussianBasisSet& dfbs = *(dfbs_);
  // Constants for convenience
  const int nbf = obs.nbasis();
  const int dfnbf = dfbs.nbasis();
  const int natom = obs.ncenter();
  /*-----------------------------------------------------*/
  /* Initialize coefficient memory                  {{{2 */ #if 2 // begin fold
  // Now initialize the coefficient memory.
  // First, compute the amount of memory needed
  // Coefficients will be stored jbf <= ibf
  int ncoefs = 0;
  for(auto ibf : function_iterator(&obs, &dfbs)){
    for(auto jbf : function_iterator(&obs, &dfbs, 0, ibf)){
      ncoefs += ibf.atom_dfnbf;
      if(ibf.center != jbf.center){
        ncoefs += jbf.atom_dfnbf;
      }
    }
  }
  coefficients_data_ = allocate<double>(ncoefs);
  memset(coefficients_data_, 0, ncoefs * sizeof(double));
  double *spot = coefficients_data_;
  for(auto ibf : function_iterator(&obs, &dfbs)){
    for(auto jbf : function_iterator(&obs, &dfbs, 0, ibf)){
      double *data_spot_a = spot;
      spot += ibf.atom_dfnbf;
      double *data_spot_b = spot;
      if(ibf.center != jbf.center){
        spot += jbf.atom_dfnbf;
      }
      IntPair ij(ibf, jbf);
      CoefContainer coefs_a = make_shared<Eigen::Map<Eigen::VectorXd>>(data_spot_a, ibf.atom_dfnbf);
      CoefContainer coefs_b = make_shared<Eigen::Map<Eigen::VectorXd>>(
          data_spot_b, ibf.center == jbf.center ? 0 : jbf.atom_dfnbf);
      coefs_.emplace(ij, std::make_pair(coefs_a, coefs_b));
    }
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
          for(auto ksh : iter_shells_on_center(&dfbs, ish.center)){
            std::shared_ptr<Eigen::MatrixXd> ij_M_k = ints_to_eigen(
                ish, jsh, ksh,
                metric_ints_3c_[ithr],
                metric_oper_type_
            );
            for(auto ibf : iter_functions_in_shell(&obs, ish)){
              for(auto jbf : iter_functions_in_shell(&obs, jsh)){
                const int ijbf = ibf.bfoff_in_shell * jsh.nbf + jbf.bfoff_in_shell;
                ij_M_X.row(ijbf).segment(ksh.bfoff_in_atom, ksh.nbf) = ij_M_k->row(ijbf);
              } // end loop over functions in jsh
            } // end loop over functions in ish
          } // end loop over shells on ish.center
          if(ish.center != jsh.center){
            for(auto ksh : iter_shells_on_center(&dfbs, jsh.center)){
              std::shared_ptr<Eigen::MatrixXd> ij_M_k = ints_to_eigen(
                  ish, jsh, ksh,
                  metric_ints_3c_[ithr],
                  metric_oper_type_
              );
              for(auto ibf : iter_functions_in_shell(&obs, ish)){
                for(auto jbf : iter_functions_in_shell(&obs, jsh)){
                  const int ijbf = ibf.bfoff_in_shell * jsh.nbf + jbf.bfoff_in_shell;
                  const int dfbfoff = ish.atom_dfnbf + ksh.bfoff_in_atom;
                  ij_M_X.row(ijbf).segment(dfbfoff, ksh.nbf) = ij_M_k->row(ijbf);
                } // end loop over functions in jsh
              } // end loop over functions in ish
            } // end loop over shells on jsh.center
          } // end if ish.center != jsh.center
          //----------------------------------------//
          for(int ibf = 0; ibf < ish.nbf; ++ibf){
            // Since coefficients are only stored jbf <= ibf,
            //   we need to figure out how many jbf's to do
            const int jmax = ish == jsh ? ibf : jsh.nbf-1;
            for(int jbf = 0; jbf <= jmax; ++jbf){
              const int ijbf = ibf*jsh.nbf + jbf;
              IntPair ij(ibf + ish.bfoff, jbf + jsh.bfoff);
              Eigen::VectorXd Ctmp(dfnbfAB);
              Ctmp = decomp->solve(ij_M_X.row(ijbf).transpose());
              *(coefs_[ij].first) = Ctmp.head(ish.atom_dfnbf);
              if(ish.center != jsh.center){
                *(coefs_[ij].second) = Ctmp.tail(jsh.atom_dfnbf);
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
  /* Debugging output                                     		                        {{{1 */ #if 1 // begin fold
  if(xml_debug) {
    begin_xml_context(
        "df_coefficients",
        "coefficients.xml"
    );
    for(auto mu : function_iterator(&obs, &dfbs)){
      for(auto nu : function_iterator(&obs, &dfbs, mu)) {
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
        std::shared_ptr<Eigen::MatrixXd> shell_ints = ints_to_eigen(
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
          std::shared_ptr<Eigen::MatrixXd> shell_ints = ints_to_eigen(
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
          std::shared_ptr<Eigen::MatrixXd> shell_ints = ints_to_eigen(
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

std::shared_ptr<Eigen::MatrixXd>
CADFCLHF::ints_to_eigen(int ish, int jsh, Ref<TwoBodyTwoCenterInt> ints, TwoBodyOper::type int_type){
  // TODO permutational symmetry; keep in mind that a transpose is needed
  return ints2_.get(ish, jsh, int_type, [ish, jsh, int_type, &ints](){
    GaussianBasisSet& ibs = *(ints->basis1());
    GaussianBasisSet& jbs = *(ints->basis2());
    const int nbfi = ibs.shell(ish).nfunction();
    const int nbfj = jbs.shell(jsh).nfunction();
    auto rv = make_shared<Eigen::MatrixXd>(nbfi, nbfj);
    ints->compute_shell(ish, jsh);
    const double* buffer = ints->buffer(int_type);
    // TODO vectorized copy
    int buffer_spot = 0;
    for(int ibf = 0; ibf < nbfi; ++ibf){
      for(int jbf = 0; jbf < nbfj; ++jbf){
        (*rv)(ibf, jbf) = buffer[buffer_spot++];
      }
    }
    return rv;
  });
}

//////////////////////////////////////////////////////////////////////////////////

std::shared_ptr<Eigen::MatrixXd>
CADFCLHF::ints_to_eigen(
    int_range&& ishs,
    int_range&& jshs,
    Ref<TwoBodyTwoCenterInt> ints,
    TwoBodyOper::type int_type
)
{
  GaussianBasisSet& ibs = *(ints->basis1());
  GaussianBasisSet& jbs = *(ints->basis2());
  //----------------------------------------//
  // First figure out the sizes
  int nbfis = 0;
  for(auto ish : ishs)
    nbfis += ibs.shell(ish).nfunction();
  int nbfjs = 0;
  for(auto jsh : jshs)
    nbfjs += jbs.shell(jsh).nfunction();
  //----------------------------------------//
  auto rv = make_shared<Eigen::MatrixXd>(nbfis, nbfjs);
  //----------------------------------------//
  int ioffset = 0, joffset = 0;
  for(int ish : ishs) {
    const int nbfi = ibs.shell(ish).nfunction();
    joffset = 0;
    for(int jsh : jshs) {
      const int nbfj = jbs.shell(jsh).nfunction();
      rv->block(
          ioffset, joffset,
          nbfi, nbfj
      ) = *(ints_to_eigen(ish, jsh, ints, int_type));
      joffset += nbfj;
    }
    ioffset += nbfi;
  }
  return rv;
}

//////////////////////////////////////////////////////////////////////////////////

std::shared_ptr<Eigen::MatrixXd>
CADFCLHF::ints_to_eigen(int ish, int jsh, int ksh, Ref<TwoBodyThreeCenterInt> ints, TwoBodyOper::type int_type){
  // TODO permutational symmetry; keep in mind that a transpose is needed
  return ints3_.get(ish, jsh, ksh, int_type, [ish, jsh, ksh, int_type, &ints](){
    const int nbfi = ints->basis1()->shell(ish).nfunction();
    const int nbfj = ints->basis2()->shell(jsh).nfunction();
    const int nbfk = ints->basis3()->shell(ksh).nfunction();
    auto rv = make_shared<Eigen::MatrixXd>(nbfi * nbfj, nbfk);
    ints->compute_shell(ish, jsh, ksh);
    const double* buffer = ints->buffer(int_type);
    // TODO vectorized copy
    int buffer_spot = 0;
    for(int ibf = 0; ibf < nbfi; ++ibf){
      for(int jbf = 0; jbf < nbfj; ++jbf){
        for(int kbf = 0; kbf < nbfk; ++kbf){
          (*rv)(ibf*nbfj + jbf, kbf) = buffer[buffer_spot++];
        }
      }
    }
    return rv;
  });
}

/////////////////////////////////////////////////////////////////////////////
// Iterator stuff

//============================================================================//
// ShellIter

ShellIter::ShellIter(
    const Ref<GaussianBasisSet>& basis,
    int position
) : pos_(position),
    basis_(basis),
    dfbasis_(0)
{

}

ShellIter::ShellIter(
    const Ref<GaussianBasisSet>& basis,
    const Ref<GaussianBasisSet>& dfbasis,
    int position
) : pos_(position),
    basis_(basis),
    dfbasis_(dfbasis)
{

}

bool
ShellIter::operator!=(const ShellIter& other) const
{
  return pos_ != other.pos_;
  // This doesn't seem to work. Maybe there's a mistake?  It's not that important
  //    or not basis_->equiv(other.basis_)
  //    or (dfbasis_.nonnull() and not other.basis_.nonnull())
  //    or (dfbasis_.null() and not other.basis_.null())
  //    or (dfbasis_.nonnull() and not dfbasis_->equiv(other.dfbasis_));
}

const ShellIter&
ShellIter::operator++()
{
  ++pos_;
  return *this;
}

const ShellData
ShellIter::operator*() const
{
  if(dfbasis_.nonnull()) {
    ShellData rv(pos_, basis_, dfbasis_);
    return rv;
  }
  else {
    ShellData rv(pos_, basis_);
    return rv;
  }

}

//============================================================================//
// BasisFunctionIter

BasisFunctionIter::BasisFunctionIter(
    const Ref<GaussianBasisSet>& basis,
    int position
) : pos_(position),
    basis_(basis),
    dfbasis_(0)
{

}

BasisFunctionIter::BasisFunctionIter(
    const Ref<GaussianBasisSet>& basis,
    const Ref<GaussianBasisSet>& dfbasis,
    int position
) : pos_(position),
    basis_(basis),
    dfbasis_(dfbasis)
{

}

bool
BasisFunctionIter::operator!=(const BasisFunctionIter& other) const
{
  return pos_ != other.pos_;
  // This doesn't seem to work. Maybe there's a mistake?  It's not that important
  //    or not basis_->equiv(other.basis_)
  //    or (dfbasis_.nonnull() and not other.basis_.nonnull())
  //    or (dfbasis_.null() and not other.basis_.null())
  //    or (dfbasis_.nonnull() and not dfbasis_->equiv(other.dfbasis_));
}

const BasisFunctionIter&
BasisFunctionIter::operator++()
{
  ++pos_;
  return *this;
}

const BasisFunctionData
BasisFunctionIter::operator*() const
{
  if(dfbasis_.nonnull()) {
    BasisFunctionData rv(pos_, basis_, dfbasis_);
    return rv;
  }
  else {
    BasisFunctionData rv(pos_, basis_);
    return rv;
  }
}

//============================================================================//
// ShellData

ShellData::ShellData()
  : index(-1),
    basis_(0),
    dfbasis_(0),
    bfoff(this, &ShellData::assert_not_initialized),
    nbf(this, &ShellData::assert_not_initialized),
    center(this, &ShellData::assert_not_initialized),
    atom_bfoff(this, &ShellData::assert_not_initialized),
    atom_shoff(this, &ShellData::assert_not_initialized),
    atom_nsh(this, &ShellData::assert_not_initialized),
    atom_nbf(this, &ShellData::assert_not_initialized),
    bfoff_in_atom(this, &ShellData::assert_not_initialized),
    shoff_in_atom(this, &ShellData::assert_not_initialized),
    atom_last_function(this, &ShellData::assert_not_initialized),
    atom_last_shell(this, &ShellData::assert_not_initialized),
    last_function(this, &ShellData::assert_not_initialized),
    atom_dfshoff(this, &ShellData::assert_not_initialized),
    atom_dfbfoff(this, &ShellData::assert_not_initialized),
    atom_dfnbf(this, &ShellData::assert_not_initialized),
    atom_dfnsh(this, &ShellData::assert_not_initialized),
    atom_df_last_function(this, &ShellData::assert_not_initialized),
    atom_df_last_shell(this, &ShellData::assert_not_initialized)
{

}

ShellData::ShellData(
    int idx,
    Ref<GaussianBasisSet> basis,
    Ref<GaussianBasisSet> dfbasis
) : index(idx),
    basis_(basis),
    dfbasis_(dfbasis),
    bfoff(this, &ShellData::get_bfoff),
    nbf(this, &ShellData::get_nbf),
    center(this, &ShellData::get_center),
    atom_bfoff(this, &ShellData::get_atom_bfoff),
    atom_shoff(this, &ShellData::get_atom_shoff),
    atom_nsh(this, &ShellData::get_atom_nsh),
    atom_nbf(this, &ShellData::get_atom_nbf),
    bfoff_in_atom(this, &ShellData::get_bfoff_in_atom),
    shoff_in_atom(this, &ShellData::get_shoff_in_atom),
    atom_last_function(this, &ShellData::get_atom_last_function),
    atom_last_shell(this, &ShellData::get_atom_last_shell),
    last_function(this, &ShellData::get_last_function),
    atom_dfshoff(this, &ShellData::get_atom_dfshoff),
    atom_dfbfoff(this, &ShellData::get_atom_dfbfoff),
    atom_dfnbf(this, &ShellData::get_atom_dfnbf),
    atom_dfnsh(this, &ShellData::get_atom_dfnsh),
    atom_df_last_function(this, &ShellData::get_atom_df_last_function),
    atom_df_last_shell(this, &ShellData::get_atom_df_last_shell)
{
  assert(basis_->nshell() > 0);
}

ShellData::ShellData(
    const ShellData& other
) : index(other.index),
    basis_(other.basis_),
    dfbasis_(other.dfbasis_),
    bfoff(this, &ShellData::get_bfoff),
    nbf(this, &ShellData::get_nbf),
    center(this, &ShellData::get_center),
    atom_bfoff(this, &ShellData::get_atom_bfoff),
    atom_shoff(this, &ShellData::get_atom_shoff),
    atom_nsh(this, &ShellData::get_atom_nsh),
    atom_nbf(this, &ShellData::get_atom_nbf),
    bfoff_in_atom(this, &ShellData::get_bfoff_in_atom),
    shoff_in_atom(this, &ShellData::get_shoff_in_atom),
    atom_last_function(this, &ShellData::get_atom_last_function),
    atom_last_shell(this, &ShellData::get_atom_last_shell),
    last_function(this, &ShellData::get_last_function),
    atom_dfshoff(this, &ShellData::get_atom_dfshoff),
    atom_dfbfoff(this, &ShellData::get_atom_dfbfoff),
    atom_dfnbf(this, &ShellData::get_atom_dfnbf),
    atom_dfnsh(this, &ShellData::get_atom_dfnsh),
    atom_df_last_function(this, &ShellData::get_atom_df_last_function),
    atom_df_last_shell(this, &ShellData::get_atom_df_last_shell)
{
  assert(basis_->nshell() > 0);
}

ShellData
sc::shell_data(
    Ref<GaussianBasisSet> basis,
    int ish, Ref<GaussianBasisSet> dfbasis
)
{
  return ShellData(ish, basis, dfbasis);
}

ShellData&
ShellData::operator=(const ShellData& other)
{
  index = other.index;
  basis_ = other.basis_;
  dfbasis_ = other.dfbasis_;
  bfoff.set_target_getter_setter(this, &ShellData::get_bfoff);
  nbf.set_target_getter_setter(this, &ShellData::get_nbf);
  center.set_target_getter_setter(this, &ShellData::get_center);
  atom_bfoff.set_target_getter_setter(this, &ShellData::get_atom_bfoff);
  atom_shoff.set_target_getter_setter(this, &ShellData::get_atom_shoff);
  atom_nsh.set_target_getter_setter(this, &ShellData::get_atom_nsh);
  atom_nbf.set_target_getter_setter(this, &ShellData::get_atom_nbf);
  bfoff_in_atom.set_target_getter_setter(this, &ShellData::get_bfoff_in_atom);
  shoff_in_atom.set_target_getter_setter(this, &ShellData::get_shoff_in_atom);
  atom_last_function.set_target_getter_setter(this, &ShellData::get_atom_last_function);
  atom_last_shell.set_target_getter_setter(this, &ShellData::get_atom_last_shell);
  last_function.set_target_getter_setter(this, &ShellData::get_last_function);
  atom_dfshoff.set_target_getter_setter(this, &ShellData::get_atom_dfshoff);
  atom_dfbfoff.set_target_getter_setter(this, &ShellData::get_atom_dfbfoff);
  atom_dfnbf.set_target_getter_setter(this, &ShellData::get_atom_dfnbf);
  atom_dfnsh.set_target_getter_setter(this, &ShellData::get_atom_dfnsh);
  atom_df_last_function.set_target_getter_setter(this, &ShellData::get_atom_df_last_function);
  atom_df_last_shell.set_target_getter_setter(this, &ShellData::get_atom_df_last_shell);
  return *this;
}

//============================================================================//
// BasisFunctionData

BasisFunctionData::BasisFunctionData(
    int idx,
    const Ref<GaussianBasisSet>& basis,
    const Ref<GaussianBasisSet>& dfbasis
) : shell(&(basis->shell(basis->function_to_shell(idx)))),
    index(idx)
{
  shell_index = basis->function_to_shell(index);
  center = basis->shell_to_center(shell_index);
  shell_bfoff = basis->shell_to_function(shell_index);
  shell_nbf = shell->nfunction();
  atom_shoff = basis->shell_on_center(center, 0);
  atom_bfoff = basis->shell_to_function(atom_shoff);
  atom_nsh = basis->nshell_on_center(center);
  atom_nbf = basis->nbasis_on_center(center);
  bfoff_in_atom = index - atom_bfoff;
  bfoff_in_shell = index - shell_bfoff;
  shoff_in_atom = shell_index - atom_shoff;
  atom_last_function = atom_bfoff + atom_nbf - 1;
  atom_last_shell = atom_shoff + atom_nsh - 1;

  if(dfbasis.nonnull()){
    atom_dfshoff = dfbasis->shell_on_center(center, 0);
    atom_dfbfoff = dfbasis->shell_to_function(atom_dfshoff);
    atom_dfnsh = dfbasis->nshell_on_center(center);
    atom_dfnbf = dfbasis->nbasis_on_center(center);
    atom_df_last_function = atom_dfbfoff + atom_dfnbf - 1;
    atom_df_last_shell = atom_dfshoff + atom_dfnsh - 1;
  }
}

const BasisFunctionData
sc::function_data(
    const Ref<GaussianBasisSet>& basis,
    int ish, const Ref<GaussianBasisSet>& dfbasis
)
{
  return BasisFunctionData(ish, basis, dfbasis);
}

//============================================================================//
// shell_iterator_wrapper

shell_iter_wrapper::shell_iter_wrapper(
    const Ref<GaussianBasisSet>& basis,
    const Ref<GaussianBasisSet>& dfbasis,
    int first_shell,
    int last_shell
) : basis_(basis),
    dfbasis_(dfbasis),
    first_shell_(first_shell),
    last_shell_(last_shell)
{

}

ShellIter
shell_iter_wrapper::begin() const
{
  if(dfbasis_.nonnull()){
    ShellIter rv(basis_, dfbasis_, first_shell_);
    return rv;
  }
  else{
    ShellIter rv(basis_, first_shell_);
    return rv;
  }
}

ShellIter
shell_iter_wrapper::end() const
{
  if(dfbasis_.nonnull()){
    ShellIter rv(basis_, dfbasis_, last_shell_ == -1 ? basis_->nshell() : last_shell_ + 1);
    return rv;
  }
  else{
    ShellIter rv(basis_, last_shell_ == -1 ? basis_->nshell() : last_shell_ + 1);
    return rv;
  }
}

//============================================================================//
// function_iterator_wrapper

function_iter_wrapper::function_iter_wrapper(
    const Ref<GaussianBasisSet>& basis,
    const Ref<GaussianBasisSet>& dfbasis,
    int first_function,
    int last_function
) : basis_(basis),
    dfbasis_(dfbasis),
    first_function_(first_function),
    last_function_(last_function)
{

}

BasisFunctionIter
function_iter_wrapper::begin() const
{
  if(dfbasis_.nonnull()){
    BasisFunctionIter rv(basis_, dfbasis_, first_function_);
    return rv;
  }
  else{
    BasisFunctionIter rv(basis_, first_function_);
    return rv;
  }
}

BasisFunctionIter
function_iter_wrapper::end() const
{
  if(dfbasis_.nonnull()){
    BasisFunctionIter rv(basis_, dfbasis_, last_function_ == -1 ? basis_->nbasis() : last_function_ + 1);
    return rv;
  }
  else{
    BasisFunctionIter rv(basis_, last_function_ == -1 ? basis_->nbasis() : last_function_ + 1);
    return rv;
  }
}

//============================================================================//
// shell_iterator

const shell_iter_wrapper
sc::shell_iterator(
    const Ref<GaussianBasisSet>& basis,
    const Ref<GaussianBasisSet>& dfbasis
)
{
  return shell_iter_wrapper(basis, dfbasis, 0, -1);
}

const shell_iter_wrapper
sc::shell_iterator(
    const Ref<GaussianBasisSet>& basis,
    int last_shell
)
{
  return shell_iter_wrapper(basis, 0, 0, last_shell);
}

const shell_iter_wrapper
sc::shell_iterator(
    const Ref<GaussianBasisSet>& basis,
    int first_shell,
    int last_shell
)
{
  return shell_iter_wrapper(basis, 0, first_shell, last_shell);
}

const shell_iter_wrapper
sc::shell_iterator(
    const Ref<GaussianBasisSet>& basis,
    const Ref<GaussianBasisSet>& dfbasis,
    int last_shell
)
{
  return shell_iter_wrapper(basis, dfbasis, 0, last_shell);
}

const shell_iter_wrapper
sc::shell_iterator(
    const Ref<GaussianBasisSet>& basis,
    const Ref<GaussianBasisSet>& dfbasis,
    int first_shell,
    int last_shell
)
{
  return shell_iter_wrapper(basis, dfbasis, first_shell, last_shell);
}

//============================================================================//
// function_iterator

const function_iter_wrapper
sc::function_iterator(
    const Ref<GaussianBasisSet>& basis,
    const Ref<GaussianBasisSet>& dfbasis
)
{
  return function_iter_wrapper(basis, dfbasis, 0, -1);
}

const function_iter_wrapper
sc::function_iterator(
    const Ref<GaussianBasisSet>& basis,
    int last_function,
    const Ref<GaussianBasisSet>& dfbasis
)
{
  return function_iter_wrapper(basis, dfbasis, 0, last_function);
}

const function_iter_wrapper
sc::function_iterator(
    const Ref<GaussianBasisSet>& basis,
    int first_function,
    int last_function,
    const Ref<GaussianBasisSet>& dfbasis
)
{
  return function_iter_wrapper(basis, dfbasis, first_function, last_function);
}

const function_iter_wrapper
sc::function_iterator(
    const Ref<GaussianBasisSet>& basis,
    const Ref<GaussianBasisSet>& dfbasis,
    int first_function,
    int last_function
)
{
  return function_iter_wrapper(basis, dfbasis, first_function, last_function);
}

const function_iter_wrapper
sc::function_iterator(
    const Ref<GaussianBasisSet>& basis,
    const Ref<GaussianBasisSet>& dfbasis,
    int last_function
)
{
  return function_iter_wrapper(basis, dfbasis, 0, last_function);
}

//============================================================================//

const shell_iter_wrapper
sc::iter_shells_on_center(
    const Ref<GaussianBasisSet>& basis,
    int center,
    const Ref<GaussianBasisSet>& dfbasis
)
{
  ShellData sh(basis->shell_on_center(center, 0), basis, dfbasis);
  return shell_iter_wrapper(
      basis, dfbasis, sh.atom_shoff, sh.atom_last_shell
  );
}

const function_iter_wrapper
sc::iter_functions_on_center(
    const Ref<GaussianBasisSet>& basis,
    int center,
    const Ref<GaussianBasisSet>& dfbasis
)
{
  ShellData sh(basis->shell_on_center(center, 0), basis, dfbasis);
  return function_iter_wrapper(
      basis, dfbasis, sh.atom_bfoff, sh.atom_last_function
  );
}

const function_iter_wrapper
sc::iter_functions_in_shell(
    const Ref<GaussianBasisSet>& basis,
    int shell,
    const Ref<GaussianBasisSet>& dfbasis
)
{
  ShellData sh(shell, basis, dfbasis);
  return function_iter_wrapper(
      basis, dfbasis, sh.bfoff, sh.last_function
  );
}

