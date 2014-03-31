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
#define DEBUG_J_INTERMEDIATES 1

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
  //----------------------------------------//
  #if DEBUG_J_INTERMEDIATES
  DEBUG_DELETE_THIS
  DEBUG_DELETE_THIS if(xml_debug_) {
  DEBUG_DELETE_THIS   for(auto&& ish : shell_range(gbs_)) {
  DEBUG_DELETE_THIS     for(auto&& jsh : shell_range(gbs_)) {
  DEBUG_DELETE_THIS       std::vector<RowMatrix> g4_out(ish.nbf*jsh.nbf);
  DEBUG_DELETE_THIS       for(auto&& mu : function_range(ish))
  DEBUG_DELETE_THIS         for(auto&& nu : function_range(jsh))
  DEBUG_DELETE_THIS           g4_out[mu.off*jsh.nbf + nu.off].resize(nbf, nbf);
  DEBUG_DELETE_THIS       for(auto&& ksh : shell_range(gbs_)) {
  DEBUG_DELETE_THIS         for(auto&& lsh : shell_range(gbs_)) {
  DEBUG_DELETE_THIS           auto g4_ptr = ints_to_eigen(ish, jsh, ksh, lsh, tbis_[0], coulomb_oper_type_);
  DEBUG_DELETE_THIS           const auto& g4 = *g4_ptr;

  DEBUG_DELETE_THIS           for(auto&& mu : function_range(ish)) {
  DEBUG_DELETE_THIS             for(auto&& nu : function_range(jsh)) {
  DEBUG_DELETE_THIS               for(auto&& rho : function_range(ksh)) {
  DEBUG_DELETE_THIS                 for(auto&& sigma : function_range(lsh)) {
  DEBUG_DELETE_THIS                   g4_out[mu.off*jsh.nbf + nu.off](rho, sigma) = g4(mu.off*jsh.nbf + nu.off, rho.off*lsh.nbf + sigma.off);
  DEBUG_DELETE_THIS                 }
  DEBUG_DELETE_THIS               }
  DEBUG_DELETE_THIS             }
  DEBUG_DELETE_THIS           }

  DEBUG_DELETE_THIS         }
  DEBUG_DELETE_THIS       }

  DEBUG_DELETE_THIS       for(auto&& mu : function_range(ish)) {
  DEBUG_DELETE_THIS         for(auto&& nu : function_range(jsh)) {
  DEBUG_DELETE_THIS           write_as_xml("g4", g4_out[mu.off*jsh.nbf + nu.off], attrs<int>{
  DEBUG_DELETE_THIS             {"ao_index1", mu}, {"ao_index2", nu}
  DEBUG_DELETE_THIS           });
  DEBUG_DELETE_THIS         }
  DEBUG_DELETE_THIS       }

  DEBUG_DELETE_THIS     }
  DEBUG_DELETE_THIS   }
  DEBUG_DELETE_THIS }
  DEBUG_DELETE_THIS
  #endif
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
        while(get_shell_pair(ish, jsh, SignificantPairs)){
          //----------------------------------------//
          // Permutation prefactor
          double pf = (ish == jsh) ? 2.0 : 4.0;
          //----------------------------------------//
          //for(int mu = ish.bfoff; mu <= ish.last_function; ++mu){
          //  for(int nu = jsh.bfoff; nu <= jsh.last_function; ++nu){
          for(auto&& mu : function_range(ish)) {
            for(auto&& nu : function_range(jsh)) {
              //if(nu > mu) continue;
            //for(int nu = jsh.bfoff; nu <= mu?; ++nu){

              auto& cpair = coefs_[{mu, nu}];
              VectorMap& Ca = *(cpair.first);
              VectorMap& Cb = *(cpair.second);
              Ct.segment(ish.atom_dfbfoff, ish.atom_dfnbf) += pf * D(mu, nu) * Ca;
              if(ish.center != jsh.center) {
                Ct.segment(jsh.atom_dfbfoff, jsh.atom_dfnbf) += pf * D(mu, nu) * Cb;
              }

              // The exact diagonal part
              if(exact_diagonal_J_ and nu <= mu) {
                const double epf = mu == nu ? 1.0 : 2.0;
                //for(auto&& Xsh : iter_shells_on_centers(dfbs_, ish.center, jsh.center)) {
                //  for(auto&& X : function_range(Xsh)) {
                //    if(Xsh.center == ish.center) {
                //      Ct_ex_thr(jsh.center, X) += Ca(X.bfoff_in_atom) * D(mu, nu) * epf;
                //    }
                //    else {
                //      Ct_ex_thr(ish.center, X) += Cb(X.bfoff_in_atom) * D(mu, nu) * epf;
                //    }
                //  }
                //}
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
  auto X_g_Y = ints_to_eigen_threaded(
    ShellBlockData<>(dfbs_), ShellBlockData<>(dfbs_),
    eris_2c_,
    coulomb_oper_type_
  );
  Eigen::VectorXd g_tilde(dfnbf);
  const auto& g2 = *X_g_Y;

  g_tilde = g2 * C_tilde;

  std::unordered_map<std::pair<int, int>, RowMatrix, sc::hash<std::pair<int, int>>> W;
  if(exact_diagonal_J_) {
    timer.enter("exact diagonal: compute local W");
    local_pairs_spot_ = 0;
    ShellData ish, jsh;
    while(get_shell_pair(ish, jsh, SignificantPairs)){
      W[{ish, jsh}].resize(ish.nbf*jsh.nbf, dfnbf);
      W[{ish, jsh}] = RowMatrix::Zero(ish.nbf*jsh.nbf, dfnbf);
      W[{jsh, ish}].resize(ish.nbf*jsh.nbf, dfnbf);
      W[{jsh, ish}] = RowMatrix::Zero(ish.nbf*jsh.nbf, dfnbf);
    }

    local_pairs_spot_ = 0;
    //do_threaded(nthread_, [&](int ithr){
      while(get_shell_pair(ish, jsh, SignificantPairs)){
        auto& Wij = W[{ish, jsh}];
        auto& Wji = W[{jsh, ish}];
        for(auto&& mu : function_range(ish)) {
          for(auto&& nu : function_range(jsh)) {

            auto& cpair = coefs_[{mu, nu}];
            VectorMap& Ca = *(cpair.first);
            VectorMap& Cb = *(cpair.second);
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
      while(get_shell_pair(ish, jsh, SignificantPairs)){
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
  #if DEBUG_J_INTERMEDIATES
  Eigen::MatrixXd Jex2(nbf, nbf);
  Jex2 = Eigen::MatrixXd::Zero(nbf, nbf);
  Eigen::MatrixXd Jex2_1(nbf, nbf);
  Jex2_1 = Eigen::MatrixXd::Zero(nbf, nbf);
  #endif

  RowMatrix d_t_ex;
  if(exact_diagonal_J_) {
    d_t_ex.resize(natom, dfnbf);
    d_t_ex = RowMatrix::Zero(natom, dfnbf);
  }
  #if DEBUG_J_INTERMEDIATES
  DEBUG_DELETE_THIS std::vector<Eigen::MatrixXd> g3_out(nbf);
  DEBUG_DELETE_THIS for(auto&& mu : function_range(gbs_)) { g3_out[mu].resize(nbf, dfnbf); }
  #endif
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
        Eigen::MatrixXd dt_ex_thr(natom, dfnbf);
        dt_ex_thr = Eigen::MatrixXd::Zero(natom, dfnbf);
        //----------------------------------------//
        ShellData ish, jsh;
        while(get_shell_pair(ish, jsh, SignificantPairs)){
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

            #if DEBUG_J_INTERMEDIATES
            DEBUG_DELETE_THIS {
            DEBUG_DELETE_THIS  for(auto&& mu : function_range(ish)) {
            DEBUG_DELETE_THIS    for(auto&& nu : function_range(jsh)) {
            DEBUG_DELETE_THIS      int block_offset = 0;
            DEBUG_DELETE_THIS      for(auto&& X : function_range(Xblk)) {
            DEBUG_DELETE_THIS        g3_out[mu](nu, X) = g3(mu.bfoff_in_shell*jsh.nbf + nu.bfoff_in_shell, block_offset);
            DEBUG_DELETE_THIS        g3_out[nu](mu, X) = g3(mu.bfoff_in_shell*jsh.nbf + nu.bfoff_in_shell, block_offset);
            DEBUG_DELETE_THIS        ++block_offset;
            DEBUG_DELETE_THIS      }
            DEBUG_DELETE_THIS    }
            DEBUG_DELETE_THIS  }
            DEBUG_DELETE_THIS }
            #endif

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

                //Eigen::MatrixXd part(ish.nbf*jsh.nbf, Xblk.nbf);
                //part = g3 - Wij.middleCols(Xblk.bfoff_in_atom, Xblk.nbf);
                //if(ish.center != jsh.center) {
                //  for(auto&& mu : function_range(ish)) {
                //    for(auto&& rho : function_range(jsh)) {
                //      part.row(mu.bfoff_in_shell*jsh.nbf + rho.bfoff_in_shell) -=
                //        Wji.row(rho.bfoff_in_shell*ish.nbf + mu.bfoff_in_shell).middleCols(Xblk.bfoff_in_atom, Xblk.nbf);
                //    }
                //  }
                //}

                for(auto&& mu : function_range(ish)) {
                  // TODO Remove loops where possible (may have to reconsider how Wij is stored)
                  const double delta_factor = ish.center != jsh.center ? 2.0 : 1.0;
                  jpart.row(mu).segment(jsh.bfoff, jsh.nbf).transpose() -= delta_factor
                      * g3.middleRows(mu.bfoff_in_shell*jsh.nbf, jsh.nbf)
                      * C_t_ex.row(jsh.center).segment(Xblk.bfoff, Xblk.nbf).transpose();
                  #if DEBUG_J_INTERMEDIATES
                  Jex2.row(mu).segment(jsh.bfoff, jsh.nbf).transpose() -= delta_factor
                      * g3.middleRows(mu.bfoff_in_shell*jsh.nbf, jsh.nbf)
                      * C_t_ex.row(jsh.center).segment(Xblk.bfoff, Xblk.nbf).transpose();
                  Jex2_1.row(mu).segment(jsh.bfoff, jsh.nbf).transpose() -= delta_factor
                      * g3.middleRows(mu.bfoff_in_shell*jsh.nbf, jsh.nbf)
                      * C_t_ex.row(jsh.center).segment(Xblk.bfoff, Xblk.nbf).transpose();
                  #endif
                  for(auto&& nu : function_range(jsh)) {
                    jpart(mu, nu) += delta_factor
                        * Wij.row(mu.bfoff_in_shell*jsh.nbf + nu.bfoff_in_shell).segment(Xblk.bfoff, Xblk.nbf)
                        * C_t_ex.row(jsh.center).segment(Xblk.bfoff, Xblk.nbf).transpose();
                    #if DEBUG_J_INTERMEDIATES
                    Jex2(mu, nu) += delta_factor
                        * Wij.row(mu.bfoff_in_shell*jsh.nbf + nu.bfoff_in_shell).segment(Xblk.bfoff, Xblk.nbf)
                        * C_t_ex.row(jsh.center).segment(Xblk.bfoff, Xblk.nbf).transpose();
                    #endif
                    if(ish.center != jsh.center) {
                      jpart(mu, nu) += delta_factor
                          * Wji.row(nu.bfoff_in_shell*ish.nbf + mu.bfoff_in_shell).segment(Xblk.bfoff, Xblk.nbf)
                          * C_t_ex.row(jsh.center).segment(Xblk.bfoff, Xblk.nbf).transpose();
                      #if DEBUG_J_INTERMEDIATES
                      Jex2(mu, nu) += delta_factor
                          * Wji.row(nu.bfoff_in_shell*ish.nbf + mu.bfoff_in_shell).segment(Xblk.bfoff, Xblk.nbf)
                          * C_t_ex.row(jsh.center).segment(Xblk.bfoff, Xblk.nbf).transpose();
                      #endif
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

                //Eigen::MatrixXd part(ish.nbf*jsh.nbf, Xblk.nbf);
                //part = g3 - Wij.middleCols(Xblk.bfoff_in_atom, Xblk.nbf);
                //for(auto&& mu : function_range(ish)) {
                //  for(auto&& rho : function_range(jsh)) {
                //    part.row(mu.bfoff_in_shell*jsh.nbf + rho.bfoff_in_shell) -=
                //        Wji.row(rho.bfoff_in_shell*ish.nbf + mu.bfoff_in_shell).middleCols(Xblk.bfoff_in_atom, Xblk.nbf);
                //  }
                //}

                for(auto&& mu : function_range(ish)) {
                  // TODO Remove loops where possible
                  jpart.row(mu).segment(jsh.bfoff, jsh.nbf).transpose() -= 2.0 *
                      g3.middleRows(mu.bfoff_in_shell*jsh.nbf, jsh.nbf)
                      * C_t_ex.row(ish.center).segment(Xblk.bfoff, Xblk.nbf).transpose();
                  #if DEBUG_J_INTERMEDIATES
                  Jex2.row(mu).segment(jsh.bfoff, jsh.nbf).transpose() -= 2.0 *
                      g3.middleRows(mu.bfoff_in_shell*jsh.nbf, jsh.nbf)
                      * C_t_ex.row(ish.center).segment(Xblk.bfoff, Xblk.nbf).transpose();
                  Jex2_1.row(mu).segment(jsh.bfoff, jsh.nbf).transpose() -= 2.0 *
                      g3.middleRows(mu.bfoff_in_shell*jsh.nbf, jsh.nbf)
                      * C_t_ex.row(ish.center).segment(Xblk.bfoff, Xblk.nbf).transpose();
                  #endif
                  for(auto&& nu : function_range(jsh)) {
                    jpart(mu, nu) += 2.0 * Wij.row(mu.bfoff_in_shell*jsh.nbf + nu.bfoff_in_shell).segment(Xblk.bfoff, Xblk.nbf)
                        * C_t_ex.row(ish.center).segment(Xblk.bfoff, Xblk.nbf).transpose();
                    jpart(mu, nu) += 2.0 * Wji.row(nu.bfoff_in_shell*ish.nbf + mu.bfoff_in_shell).segment(Xblk.bfoff, Xblk.nbf)
                        * C_t_ex.row(ish.center).segment(Xblk.bfoff, Xblk.nbf).transpose();
                    #if DEBUG_J_INTERMEDIATES
                    Jex2(mu, nu) += 2.0 * Wij.row(mu.bfoff_in_shell*jsh.nbf + nu.bfoff_in_shell).segment(Xblk.bfoff, Xblk.nbf)
                        * C_t_ex.row(ish.center).segment(Xblk.bfoff, Xblk.nbf).transpose();
                    Jex2(mu, nu) += 2.0 * Wji.row(nu.bfoff_in_shell*ish.nbf + mu.bfoff_in_shell).segment(Xblk.bfoff, Xblk.nbf)
                        * C_t_ex.row(ish.center).segment(Xblk.bfoff, Xblk.nbf).transpose();
                    #endif
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
        d_t_ex += dt_ex_thr;
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
  if(xml_debug_) {
    write_as_xml("d_tilde", d_tilde);
    if(exact_diagonal_J_) {
      write_as_xml("d_t_ex", d_t_ex);
      #if DEBUG_J_INTERMEDIATES
      for(auto&& mu : function_range(gbs_)) {
        write_as_xml("g3", g3_out[mu], std::map<std::string, int>{ { "ao_index1", mu.index } });
      }
      write_as_xml("Jex2", Jex2.triangularView<Eigen::Lower>());
      write_as_xml("Jex2_1", Jex2_1.triangularView<Eigen::Lower>());
      #endif
    }
  }
  timer.exit("compute d_tilde");
  #if DEBUG_J_INTERMEDIATES
  Eigen::MatrixXd Jex(nbf, nbf);
  Jex = Eigen::MatrixXd::Zero(nbf, nbf);
  Eigen::MatrixXd Jex3(nbf, nbf);
  Jex3 = Eigen::MatrixXd::Zero(nbf, nbf);
  #endif
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
        while(get_shell_pair(ish, jsh, SignificantPairs)){
          /*-----------------------------------------------------*/
          /* first and third terms compute thread           {{{2 */ #if 2 // begin fold
          //----------------------------------------//
          TimerHolder subtimer(contract_timer);
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
              if(exact_diagonal_J_) {
                jpart(mu, nu) -= 2.0 * d_t_ex.row(jsh.center).segment(ish.atom_dfbfoff, ish.atom_dfnbf) * Ca;
                #if DEBUG_J_INTERMEDIATES
                Jex3(mu, nu) -= 2.0 * d_t_ex.row(jsh.center).segment(ish.atom_dfbfoff, ish.atom_dfnbf) * Ca;
                #endif
                if(ish.center != jsh.center) {
                  jpart(mu, nu) -= 2.0 * d_t_ex.row(ish.center).segment(jsh.atom_dfbfoff, jsh.atom_dfnbf) * Cb;
                  #if DEBUG_J_INTERMEDIATES
                  Jex3(mu, nu) -= 2.0 * d_t_ex.row(ish.center).segment(jsh.atom_dfbfoff, jsh.atom_dfnbf) * Cb;
                  #endif
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
                auto g4_ptr = ints_to_eigen(ish, jsh, ksh, lsh, tbis_[ithr], coulomb_oper_type_);
                const auto& g4 = *g4_ptr;
                for(auto&& mu : function_range(ish)) {
                  for(auto&& rho : function_range(ksh)) {
                    jpart.row(mu).segment(jsh.bfoff, jsh.nbf).transpose() += epf *
                        g4.middleRows(mu.bfoff_in_shell*jsh.nbf, jsh.nbf).middleCols(rho.bfoff_in_shell*lsh.nbf, lsh.nbf)
                        * d.segment(rho*nbf + lsh.bfoff, lsh.nbf);

                    #if DEBUG_J_INTERMEDIATES
                    Jex.row(mu).segment(jsh.bfoff, jsh.nbf).transpose() += epf *
                        g4.middleRows(mu.bfoff_in_shell*jsh.nbf, jsh.nbf).middleCols(rho.bfoff_in_shell*lsh.nbf, lsh.nbf)
                        * d.segment(rho*nbf + lsh.bfoff, lsh.nbf);
                    // TODO Vectorize
                    //for(auto&& nu : function_range(jsh)) {
                    //  if(nu > mu) continue;
                    //  for(auto&& sigma : function_range(lsh)) {
                    //    jpart(mu, nu) += epf *
                    //        g4(mu.bfoff_in_shell*jsh.nbf + nu.bfoff_in_shell, rho.bfoff_in_shell*lsh.nbf + sigma.bfoff_in_shell)
                    //        * D(rho, sigma);
                    //    Jex(mu, nu) += epf *
                    //        g4(mu.bfoff_in_shell*jsh.nbf + nu.bfoff_in_shell, rho.bfoff_in_shell*lsh.nbf + sigma.bfoff_in_shell)
                    //        * D(rho, sigma);
                    //  }
                    //}
                    #endif

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
  if(xml_debug_) {
    #if DEBUG_J_INTERMEDIATES
    if(exact_diagonal_J_) {
      write_as_xml("Jex", Jex.triangularView<Eigen::Lower>());
      write_as_xml("Jex3", Jex3.triangularView<Eigen::Lower>());
    }
    #endif
  }
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


