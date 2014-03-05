//
// compute_j.cc
//
// Copyright (C) 2014 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: Feb 13, 2014
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
#include <util/misc/xmlwriter.h>

#include "cadfclhf.h"

using namespace sc;

RefSCMatrix
CADFCLHF::compute_J()
{
  /*=======================================================================================*/
  /* Setup                                                 		                        {{{1 */ #if 1 // begin fold
  //----------------------------------------//
  // Convenience variables
  Timer timer("compute J");
  const int me = scf_grp_->me();
  const int n_node = scf_grp_->n();
  const Ref<GaussianBasisSet>& obs = gbs_;
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
    for(int ithr = 0; ithr < nthread_; ++ithr) {
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
              auto& cpair = coefs_[{mu, nu}];
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
  if(xml_debug_ and scf_grp_->me() == 0) {
    write_as_xml("C_tilde", C_tilde);
    //end_xml_context("compute_J"); end_xml_context("compute_fock"); assert(scf_grp_->me() != 0);
  }
  timer.exit("compute C_tilde");

  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Form g_tilde                                          		                        {{{1 */ #if 1 // begin fold
  //----------------------------------------//
  // TODO Thread this
  timer.enter("compute g_tilde");
  auto X_g_Y = ints_to_eigen_threaded(
    ShellBlockData<>(dfbs_), ShellBlockData<>(dfbs_),
    eris_2c_,
    coulomb_oper_type_
  );
  Eigen::VectorXd g_tilde(dfnbf);
  g_tilde = (*X_g_Y) * C_tilde;
  if(xml_debug_) {
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
    for(int ithr = 0; ithr < nthread_; ++ithr) {
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
          for(auto&& Xblk : shell_block_range(dfbs_)){
            auto g_part = ints_to_eigen(
                ish, jsh, Xblk,
                eris_3c_[ithr],
                coulomb_oper_type_
            );
            for(auto&& mu : function_range(ish)) {
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
  scf_grp_->sum((double*)d_tilde.data(), dfnbf);
  if(xml_debug_) {
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
    for(int ithr = 0; ithr < nthread_; ++ithr) {
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
          for(auto&& mu : function_range(ish)){
            for(auto&& nu : function_range(jsh)){
              //----------------------------------------//
              auto cpair = coefs_[{mu, nu}];
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
  scf_grp_->sum((double*)J.data(), nbf*nbf);
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


