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
#include <util/misc/hash.h>

#include "cadfclhf.h"

using namespace sc;
#define DEBUG_J_INTERMEDIATES 0

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
  const int natom = gbs_->ncenter();
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
  RowMatrix J(nbf, nbf);
  J = RowMatrix::Zero(nbf, nbf);
  VectorMap j(J.data(), nbf*nbf);
  //----------------------------------------//
  // reset the iteration over local pairs
  local_pairs_spot_ = 0;

  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Form C_tilde and d_tilde                              		                        {{{1 */ #if 1 // begin fold
  //----------------------------------------//

  timer.enter("compute C_tilde");
  Eigen::VectorXd C_tilde(dfnbf);
  C_tilde = Eigen::VectorXd::Zero(dfnbf);
  RowMatrix C_t_ex;
  if(exact_diagonal_J_) {
    C_t_ex.resize(natom, dfnbf);
    C_t_ex = RowMatrix::Zero(natom, dfnbf);
  }

  //----------------------------------------//

  {

    boost::mutex C_tilde_mutex;
    boost::thread_group compute_threads;
    // Loop over number of threads
    for(int ithr = 0; ithr < nthread_; ++ithr) {
      // ...and create each thread that computes pairs
      compute_threads.create_thread([&,ithr](){
        /*-----------------------------------------------------*/
        /* C_tilde compute thread                         {{{2 */ #if 2 // begin fold

        decltype(C_tilde) Ct(dfnbf);
        Ct = decltype(C_tilde)::Zero(dfnbf);

        decltype(C_t_ex) Ct_ex_thr;
        if(exact_diagonal_J_) {
          Ct_ex_thr.resize(natom, dfnbf);
          Ct_ex_thr = decltype(C_t_ex)::Zero(natom, dfnbf);
        }

        ShellData ish, jsh;
        //----------------------------------------//
        while(get_shell_pair(ish, jsh, sig_pairs_J_ ? SignificantPairs : AllPairs)){
          //----------------------------------------//
          // Permutation prefactor
          double pf = (ish == jsh) ? 2.0 : 4.0;
          //----------------------------------------//
          for(auto&& mu : function_range(ish)) {
            for(auto&& nu : function_range(jsh)) {

              auto& cpair = coefs_[{mu, nu}];
              auto& Ca = *(cpair.first);
              auto& Cb = *(cpair.second);
              Ct.segment(ish.atom_dfbfoff, ish.atom_dfnbf) += pf * D(mu, nu) * Ca;
              if(ish.center != jsh.center) {
                Ct.segment(jsh.atom_dfbfoff, jsh.atom_dfnbf) += pf * D(mu, nu) * Cb;
              }

              // The exact diagonal part
              if(exact_diagonal_J_ and nu <= mu) {
                const double epf = mu == nu ? 1.0 : 2.0;
                Ct_ex_thr.row(jsh.center).segment(ish.atom_dfbfoff, ish.atom_dfnbf) += epf * D(mu, nu) * Ca;
                if(ish.center != jsh.center) {
                  Ct_ex_thr.row(ish.center).segment(jsh.atom_dfbfoff, jsh.atom_dfnbf) += epf * D(mu, nu) * Cb;
                }
              }

            } // end loop over mu
          } // end loop over nu
        } // end while get shell pair
        //----------------------------------------//
        // add our contribution to the node level C_tilde
        boost::lock_guard<boost::mutex> Ctlock(C_tilde_mutex);
        C_tilde += Ct;
        if(exact_diagonal_J_) {
          C_t_ex += Ct_ex_thr;
        }
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

  if(exact_diagonal_J_) {
    scf_grp_->sum(C_t_ex.data(), dfnbf * natom);
  }

  if(xml_debug_ and scf_grp_->me() == 0) {
    write_as_xml("C_tilde", C_tilde);
    if(exact_diagonal_J_) write_as_xml("C_t_ex", C_t_ex);
    //end_xml_context("compute_J"); end_xml_context("compute_fock"); assert(scf_grp_->me() != 0);
  }
  timer.exit("compute C_tilde");

  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Form g_tilde                                          		                        {{{1 */ #if 1 // begin fold
  //----------------------------------------//
  // TODO Thread this
  timer.enter("compute g_tilde");

  Eigen::VectorXd g_tilde(dfnbf);
  const auto& g2 = *g2_full_ptr_;

  g_tilde = g2 * C_tilde;

  std::unordered_map<std::pair<int, int>, RowMatrix, sc::hash<std::pair<int, int>>> W;
  if(exact_diagonal_J_) {
    timer.enter("exact diagonal: compute local W");
    local_pairs_spot_ = 0;
    ShellData ish, jsh;
    while(get_shell_pair(ish, jsh, sig_pairs_J_ ? SignificantPairs : AllPairs)){
      W[{ish, jsh}].resize(ish.nbf*jsh.nbf, dfnbf);
      W[{ish, jsh}] = RowMatrix::Zero(ish.nbf*jsh.nbf, dfnbf);
      W[{jsh, ish}].resize(ish.nbf*jsh.nbf, dfnbf);
      W[{jsh, ish}] = RowMatrix::Zero(ish.nbf*jsh.nbf, dfnbf);
    }

    local_pairs_spot_ = 0;
    // TODO This part needs to be optimized significantly!  At the very least, it needs to be threaded
    //do_threaded(nthread_, [&](int ithr){
      while(get_shell_pair(ish, jsh, sig_pairs_J_ ? SignificantPairs : AllPairs)){
        auto& Wij = W[{ish, jsh}];
        auto& Wji = W[{jsh, ish}];
        for(auto&& mu : function_range(ish)) {
          for(auto&& nu : function_range(jsh)) {

            auto& cpair = coefs_[{mu, nu}];
            auto& Ca = *(cpair.first);
            auto& Cb = *(cpair.second);
            Wij.row(mu.bfoff_in_shell*jsh.nbf + nu.bfoff_in_shell) +=
                Ca.transpose() * g2.middleRows(ish.atom_dfbfoff, ish.atom_dfnbf);
            if(ish != jsh) {
              if(ish.center != jsh.center) {
                Wji.row(nu.bfoff_in_shell*ish.nbf + mu.bfoff_in_shell) +=
                    Cb.transpose() * g2.middleRows(jsh.atom_dfbfoff, jsh.atom_dfnbf);
              }
              else{
                Wji.row(nu.bfoff_in_shell*ish.nbf + mu.bfoff_in_shell) +=
                    Ca.transpose() * g2.middleRows(jsh.atom_dfbfoff, jsh.atom_dfnbf);
              }
            }

          }
        }
      }

    //});
    timer.exit();
  }

  if(xml_debug_) {
    write_as_xml("g_tilde", g_tilde);
    write_as_xml("g2", g2);
    if(exact_diagonal_J_) {
      local_pairs_spot_ = 0;
      ShellData ish, jsh;
      while(get_shell_pair(ish, jsh, sig_pairs_J_ ? SignificantPairs : AllPairs)){
        auto& Wij = W[{ish, jsh}];
        auto& Wji = W[{jsh, ish}];
        for(auto&& mu : function_range(ish)) {
          for(auto&& nu : function_range(jsh)) {
            write_as_xml("W_j", Wij.row(mu.bfoff_in_shell*jsh.nbf + nu.bfoff_in_shell),
                std::map<std::string, int>{
                  {"ao_index1", mu},
                  {"ao_index2", nu}
            });
            if(ish != jsh) {
              write_as_xml("W_j", Wji.row(nu.bfoff_in_shell*ish.nbf + mu.bfoff_in_shell),
                  std::map<std::string, int>{
                    {"ao_index1", nu},
                    {"ao_index2", mu}
              });
            }
          }
        }

      }
    }
  }
  timer.exit("compute g_tilde");
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Form dtilde and add in Ctilde contribution to J       		                        {{{1 */ #if 1 // begin fold
  timer.enter("compute d_tilde");
  Eigen::VectorXd d_tilde(dfnbf);
  d_tilde = Eigen::VectorXd::Zero(dfnbf);

  RowMatrix d_t_ex;
  if(exact_diagonal_J_) {
    d_t_ex.resize(natom, dfnbf);
    d_t_ex = RowMatrix::Zero(natom, dfnbf);
  }
  //----------------------------------------//
  {
    boost::mutex tmp_mutex;
    boost::thread_group compute_threads;
    MultiThreadTimer mt_timer("threaded part", nthread_);
    //----------------------------------------//
    // reset the iteration over local pairs
    local_pairs_spot_ = 0;
    // Loop over number of threads
    for(int ithr = 0; ithr < nthread_; ++ithr) {
      // ...and create each thread that computes pairs
      compute_threads.create_thread([&,ithr](){
        auto ints_timer = mt_timer.get_subtimer("compute_ints", ithr);
        auto contract_timer = mt_timer.get_subtimer("contract", ithr);
        auto ex_timer = mt_timer.get_subtimer("exact diagonal", ithr);
        /*-----------------------------------------------------*/
        /* d_tilde compute and C_tilde contract thread    {{{2 */ #if 2 // begin fold
        Eigen::VectorXd dt(dfnbf);
        dt = Eigen::VectorXd::Zero(dfnbf);
        Eigen::MatrixXd jpart(nbf, nbf);
        jpart = Eigen::MatrixXd::Zero(nbf, nbf);
        RowMatrix dt_ex_thr;
        if(exact_diagonal_J_) {
          dt_ex_thr.resize(natom, dfnbf);
          dt_ex_thr = RowMatrix::Zero(natom, dfnbf);
        }
        //----------------------------------------//
        ShellData ish, jsh;
        while(get_shell_pair(ish, jsh, sig_pairs_J_ ? SignificantPairs : AllPairs)){
          //----------------------------------------//
          double perm_fact = (ish == jsh) ? 2.0 : 4.0;
          //----------------------------------------//
          // Note:  SameCenter shell_block requirement is the default
          for(auto&& Xblk : shell_block_range(dfbs_)){

            TimerHolder subtimer(ints_timer);
            auto g_part = ints_to_eigen(
                ish, jsh, Xblk,
                eris_3c_[ithr],
                coulomb_oper_type_
            );
            const auto& g3 = *g_part;

            subtimer.change(contract_timer);
            for(auto&& mu : function_range(ish)) {
              dt.segment(Xblk.bfoff, Xblk.nbf).transpose() += perm_fact * D.row(mu).segment(jsh.bfoff, jsh.nbf)
                  * g_part->middleRows(mu.bfoff_in_shell*jsh.nbf, jsh.nbf);
              //----------------------------------------//
              jpart.row(mu).segment(jsh.bfoff, jsh.nbf).transpose() +=
                  g_part->middleRows(mu.bfoff_in_shell*jsh.nbf, jsh.nbf)
                  * C_tilde.segment(Xblk.bfoff, Xblk.nbf);
            }

            if(exact_diagonal_J_) {

              subtimer.change(ex_timer);

              const auto& Wij = W[{ish, jsh}];
              const auto& Wji = W[{jsh, ish}];


              if(Xblk.center == ish.center) {

                for(auto&& mu : function_range(ish)) {
                  for(auto&& nu : function_range(jsh)) {
                    // TODO remove one loop
                    dt_ex_thr.row(jsh.center).segment(ish.atom_dfbfoff + Xblk.bfoff_in_atom, Xblk.nbf) += 0.5 * perm_fact *
                        D(mu, nu) * g3.row(mu.bfoff_in_shell*jsh.nbf + nu.bfoff_in_shell);
                  }
                }

                for(auto&& mu : function_range(ish)) {
                  // TODO Remove loops where possible (may have to reconsider how Wij is stored)
                  const double delta_factor = ish.center != jsh.center ? 2.0 : 1.0;
                  jpart.row(mu).segment(jsh.bfoff, jsh.nbf).transpose() -= delta_factor
                      * g3.middleRows(mu.bfoff_in_shell*jsh.nbf, jsh.nbf)
                      * C_t_ex.row(jsh.center).segment(Xblk.bfoff, Xblk.nbf).transpose();

                  for(auto&& nu : function_range(jsh)) {
                    jpart(mu, nu) += delta_factor
                        * Wij.row(mu.bfoff_in_shell*jsh.nbf + nu.bfoff_in_shell).segment(Xblk.bfoff, Xblk.nbf)
                        * C_t_ex.row(jsh.center).segment(Xblk.bfoff, Xblk.nbf).transpose();
                    if(ish.center != jsh.center) {
                      jpart(mu, nu) += delta_factor
                          * Wji.row(nu.bfoff_in_shell*ish.nbf + mu.bfoff_in_shell).segment(Xblk.bfoff, Xblk.nbf)
                          * C_t_ex.row(jsh.center).segment(Xblk.bfoff, Xblk.nbf).transpose();
                    }
                  }

                }

              }
              else if(Xblk.center == jsh.center) { // and thus ish.center != jsh.center

                for(auto&& mu : function_range(ish)) {
                  for(auto&& nu : function_range(jsh)) {
                    // TODO remove one loop
                    dt_ex_thr.row(ish.center).segment(jsh.atom_dfbfoff + Xblk.bfoff_in_atom, Xblk.nbf) += 2.0 *
                        D(mu, nu) * g3.row(mu.bfoff_in_shell*jsh.nbf + nu.bfoff_in_shell);
                  }
                }

                for(auto&& mu : function_range(ish)) {
                  // TODO Remove loops where possible
                  jpart.row(mu).segment(jsh.bfoff, jsh.nbf).transpose() -= 2.0 *
                      g3.middleRows(mu.bfoff_in_shell*jsh.nbf, jsh.nbf)
                      * C_t_ex.row(ish.center).segment(Xblk.bfoff, Xblk.nbf).transpose();
                  for(auto&& nu : function_range(jsh)) {
                    jpart(mu, nu) += 2.0 * Wij.row(mu.bfoff_in_shell*jsh.nbf + nu.bfoff_in_shell).segment(Xblk.bfoff, Xblk.nbf)
                        * C_t_ex.row(ish.center).segment(Xblk.bfoff, Xblk.nbf).transpose();
                    jpart(mu, nu) += 2.0 * Wji.row(nu.bfoff_in_shell*ish.nbf + mu.bfoff_in_shell).segment(Xblk.bfoff, Xblk.nbf)
                        * C_t_ex.row(ish.center).segment(Xblk.bfoff, Xblk.nbf).transpose();
                  }
                }

              }
            }

          } // end loop over kshbf
          // only constructing the lower triangle of the J matrix, so zero the strictly upper part
          jpart.triangularView<Eigen::StrictlyUpper>().setZero();
        } // end while get_shell_pair
        //----------------------------------------//
        // add our contribution to the node level d_tilde
        boost::lock_guard<boost::mutex> tmp_lock(tmp_mutex);
        d_tilde += dt;
        if(exact_diagonal_J_) {
          d_t_ex += dt_ex_thr;
        }
        J += jpart;
        /*******************************************************/ #endif //2}}}
        /*-----------------------------------------------------*/
      }); // end create_thread
    } // end enumeration of threads
    compute_threads.join_all();
    mt_timer.exit();
    timer.insert(mt_timer);
    if(print_iteration_timings_) mt_timer.print(ExEnv::out0(), 12, 45);
  } // compute_threads is destroyed here
  //----------------------------------------//
  // Global sum d_tilde
  scf_grp_->sum((double*)d_tilde.data(), dfnbf);
  if(exact_diagonal_J_) {
    scf_grp_->sum((double*)d_t_ex.data(), natom * dfnbf);
  }
  if(xml_debug_) {
    write_as_xml("d_tilde", d_tilde);
    if(exact_diagonal_J_) {
      write_as_xml("d_t_ex", d_t_ex);
    }
  }
  timer.exit("compute d_tilde");
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Add in first and third term contributions to J       		                        {{{1 */ #if 1 // begin fold
  {
    timer.enter("J contributions");
    boost::mutex tmp_mutex;
    boost::thread_group compute_threads;
    MultiThreadTimer mt_timer("threaded part", nthread_);
    //----------------------------------------//
    // reset the iteration over local pairs
    local_pairs_spot_ = 0;
    // Loop over number of threads
    for(int ithr = 0; ithr < nthread_; ++ithr) {
      // ...and create each thread that computes pairs
      compute_threads.create_thread([&,ithr](){
        auto contract_timer = mt_timer.get_subtimer("contract", ithr);
        auto ex_timer = mt_timer.get_subtimer("exact diagonal", ithr);
        Eigen::MatrixXd jpart(nbf, nbf);
        jpart = Eigen::MatrixXd::Zero(nbf, nbf);
        //----------------------------------------//
        ShellData ish, jsh;
        while(get_shell_pair(ish, jsh, sig_pairs_J_ ? SignificantPairs : AllPairs)){
          /*-----------------------------------------------------*/
          /* first and third terms compute thread           {{{2 */ #if 2 // begin fold
          //----------------------------------------//
          TimerHolder subtimer(contract_timer);
          for(auto&& mu : function_range(ish)){
            for(auto&& nu : function_range(jsh)){
              //----------------------------------------//
              auto cpair = coefs_[{mu, nu}];
              auto& Ca = *(cpair.first);
              auto& Cb = *(cpair.second);
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
              if(exact_diagonal_J_) {
                jpart(mu, nu) -= 2.0 * d_t_ex.row(jsh.center).segment(ish.atom_dfbfoff, ish.atom_dfnbf) * Ca;
                if(ish.center != jsh.center) {
                  jpart(mu, nu) -= 2.0 * d_t_ex.row(ish.center).segment(jsh.atom_dfbfoff, jsh.atom_dfnbf) * Cb;
                }
              }
            } // end loop over nu
          } // end loop over mu
          /*******************************************************/ #endif //2}}}
          /*-----------------------------------------------------*/
          /* Add back in the exact diagonal integrals       {{{2 */ #if 2 // begin fold
          //----------------------------------------//
          if(exact_diagonal_J_) {
            subtimer.change(ex_timer);
            double epf = (ish.center == jsh.center) ? 2.0 : 4.0;
            for(auto&& ksh : iter_shells_on_center(obs, ish.center)) {
              for(auto&& lsh : iter_shells_on_center(obs, jsh.center)) {
                if(not is_sig_pair(ksh, lsh)) continue;

                auto g4_ptr = ints_to_eigen(ish, jsh, ksh, lsh, tbis_[ithr], coulomb_oper_type_);
                const auto& g4 = *g4_ptr;
                for(auto&& mu : function_range(ish)) {
                  for(auto&& rho : function_range(ksh)) {
                    // TODO More Vectorization
                    jpart.row(mu).segment(jsh.bfoff, jsh.nbf).transpose() += epf *
                        g4.middleRows(mu.bfoff_in_shell*jsh.nbf, jsh.nbf).middleCols(rho.bfoff_in_shell*lsh.nbf, lsh.nbf)
                        * d.segment(rho*nbf + lsh.bfoff, lsh.nbf);
                  }
                }

              }
            }
          }

          /*******************************************************/ #endif //2}}}
          /*-----------------------------------------------------*/

        } // end while get_shell_pair
        // Sum the thread's contributions to the node-level J
        boost::lock_guard<boost::mutex> tmp_lock(tmp_mutex);
        J += jpart;
      }); // end create_thread
    } // end enumeration of threads
    compute_threads.join_all();
    mt_timer.exit();
    timer.insert(mt_timer);
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
  //ExEnv::out0() << indent << "J checksum: " << scprintf("%20.15f", J.sum()) << std::endl;
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


