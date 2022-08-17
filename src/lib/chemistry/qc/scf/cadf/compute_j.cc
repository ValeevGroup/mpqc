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
typedef unsigned long int uli;

RefSCMatrix
CADFCLHF::compute_J()
{
  auto result = (cadf_K_only_) ? compute_J_df() : compute_J_cadf();
  return result;
}

RefSCMatrix
CADFCLHF::compute_J_df()
{
  Timer timer("compute J(DF)"); timer.enter("build");
  // transform the density difference to the AO basis
  RefSymmSCMatrix dd = cl_dens_diff_;
  Ref<PetiteList> pl = integral()->petite_list();
  cl_dens_diff_ = pl->to_AO_basis(dd);

  double gmat_accuracy = 1e-15;
  if (min_orthog_res() < 1.0) { gmat_accuracy *= min_orthog_res(); }

  Ref<OrbitalSpaceRegistry> oreg = world_->moints_runtime()->factory()->orbital_registry();
  Ref<AOSpaceRegistry> aoreg = world_->moints_runtime()->factory()->ao_registry();
  std::string aospace_id = new_unique_key(oreg);
  if (aoreg->key_exists(basis()) == false) {
    Ref<OrbitalSpace> aospace = new AtomicOrbitalSpace(aospace_id, "CADFCLHF AO basis set", basis(), integral());
    aoreg->add(basis(), aospace);
    MPQC_ASSERT(oreg->key_exists(aospace_id) == false); // should be ensured by using new_unique_key
    oreg->add(make_keyspace_pair(aospace));
  }
  // feed the spin densities to the builder, cl_dens_diff_ includes total density right now, so halve it
  RefSymmSCMatrix Pa = cl_dens_diff_.copy(); Pa.scale(0.5);
  RefSymmSCMatrix Pb = Pa;

  Ref<FockBuildRuntime> fb_rtime = world_->fockbuild_runtime();
  fb_rtime->set_densities(Pa, Pb);

  Ref<OrbitalSpace> aospace = aoreg->value(basis());
  RefSCMatrix G;

  const std::string jkey = ParsedOneBodyIntKey::key(aospace->id(),aospace->id(),std::string("J"));
  RefSCMatrix J = fb_rtime->get(jkey);
  timer.exit("build");
  return J.copy();
}

RefSCMatrix
CADFCLHF::compute_J_cadf()
{
  /*=======================================================================================*/
  /* Setup                                                 		                        {{{1 */ #if 1 // begin fold
  //----------------------------------------//
  // Convenience variables

  Timer timer("compute J");
  const int me = scf_grp_->me();
  const int n_node = scf_grp_->n();
  const int natom = gbs_->ncenter();
  const uli nbf = gbs_->nbasis();
  const uli dfnbf = dfbs_->nbasis();
  //----------------------------------------//
  // Get the density in an Eigen::Map form
  //double *D_ptr = allocate<double>(nbf*nbf);
  double* __restrict__ D_ptr = new double[nbf*nbf];
  D_.convert(D_ptr);
  // Matrix and vector wrappers, for convenience
  Eigen::Map<Eigen::VectorXd> d(D_ptr, nbf*nbf);
  Eigen::Map<ColMatrix> D(D_ptr, nbf, nbf);
  //----------------------------------------//
  double* __restrict__ J_data = new double[nbf*nbf];
  Eigen::Map<RowMatrix> J(J_data, nbf, nbf);
  J.setZero();
  //----------------------------------------//
  // reset the iteration over local pairs
  local_pairs_spot_ = 0;

  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Form C_tilde and d_tilde                              		                        {{{1 */ #if 1 // begin fold
  //----------------------------------------//

  timer.enter("compute C_tilde");
  Eigen::VectorXd C_tilde(dfnbf);
  C_tilde.setZero();
  RowMatrix C_t_ex;
  if(exact_diagonal_J_) {
    C_t_ex.resize(natom, dfnbf);
    C_t_ex.setZero();
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
        Ct.setZero(dfnbf);

        decltype(C_t_ex) Ct_ex_thr;
        if(exact_diagonal_J_) {
          Ct_ex_thr.resize(natom, dfnbf);
          Ct_ex_thr.setZero(natom, dfnbf);
        }

        ShellData ish, jsh;
        //----------------------------------------//
        while(get_shell_pair(ish, jsh, SignificantPairs)){
          //----------------------------------------//
          // Permutation prefactor
          double pf = (ish == jsh) ? 2.0 : 4.0;
          //----------------------------------------//
          for(auto&& mu : function_range(ish)) {
            for(auto&& nu : function_range(jsh)) {

              auto& cpair = coefs_[{mu, nu}];
              auto& Ca = *(cpair.first);
              auto& Cb = *(cpair.second);
              Ct.segment(ish.atom_dfbfoff, ish.atom_dfnbf) += pf * D.coeff(mu, nu) * Ca;
              if(ish.center != jsh.center) {
                Ct.segment(jsh.atom_dfbfoff, jsh.atom_dfnbf) += pf * D.coeff(mu, nu) * Cb;
              }

              // The exact diagonal part
              if(exact_diagonal_J_ and nu <= mu) {
                const double epf = mu == nu ? 1.0 : 2.0;
                Ct_ex_thr.row(jsh.center).segment(ish.atom_dfbfoff, ish.atom_dfnbf) += epf * D.coeff(mu, nu) * Ca;
                if(ish.center != jsh.center) {
                  Ct_ex_thr.row(ish.center).segment(jsh.atom_dfbfoff, jsh.atom_dfnbf) += epf * D.coeff(mu, nu) * Cb;
                }
              }

            } // end loop over mu
          } // end loop over nu
        } // end while get shell pair
        //----------------------------------------//
        // add our contribution to the node level C_tilde
        boost::lock_guard<boost::mutex> Ctlock(C_tilde_mutex);
        C_tilde.noalias() += Ct;
        if(exact_diagonal_J_) {
          C_t_ex.noalias() += Ct_ex_thr;
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
  timer.enter("compute g_tilde");

  Eigen::VectorXd g_tilde(dfnbf);
  const auto& g2 = *g2_full_ptr_;

  g_tilde = g2 * C_tilde;

  std::unordered_map<std::pair<int, int>, RowMatrix, sc::hash<std::pair<int, int>>> W;
  if(exact_diagonal_J_) {
    timer.enter("exact diagonal: compute local W");
    local_pairs_spot_ = 0;
    ShellData ish, jsh;
    while(get_shell_pair(ish, jsh, SignificantPairs)){
      W[{ish, jsh}].resize(ish.nbf*jsh.nbf, dfnbf);
      W[{ish, jsh}].setZero();
      W[{jsh, ish}].resize(ish.nbf*jsh.nbf, dfnbf);
      W[{jsh, ish}].setZero();
    }

    local_pairs_spot_ = 0;
    // TODO This part needs to be optimized significantly!  At the very least, it needs to be threaded
    //do_threaded(nthread_, [&](int ithr){
      while(get_shell_pair(ish, jsh, SignificantPairs)){
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
        /*-----------------------------------------------------*/
        /* d_tilde compute and C_tilde contract thread    {{{2 */ #if 2 // begin fold
        mt_timer.enter("misc", ithr);
        Eigen::VectorXd dt(dfnbf);
        dt.setZero();

        RowMatrix dt_ex_thr;
        if(exact_diagonal_J_) {
          dt_ex_thr.resize(natom, dfnbf);
          dt_ex_thr.setZero();
        }

        double* __restrict__ j_intbuff = new double[max_fxn_obs_j_ish_*max_fxn_obs_j_jsh_*
          // The extra max_fxn_dfbs_ is just to be safe, since the block size is only a "target"
          std::min(dfbs_->nbasis(), (unsigned int)(DEFAULT_TARGET_BLOCK_SIZE + max_fxn_dfbs_))
        ];

        //----------------------------------------//
        mt_timer.change("loop significant ij", ithr);
        auto ints_timer = mt_timer.get_subtimer("compute_ints", ithr);
        auto contract_timer = mt_timer.get_subtimer("contract", ithr);
        auto ex_timer = mt_timer.get_subtimer("exact diagonal", ithr);
        ShellData ish, jsh;

        while(get_shell_pair(ish, jsh, SignificantPairs)){
          //----------------------------------------//
          double perm_fact = (ish == jsh) ? 2.0 : 4.0;
          //----------------------------------------//
          // Note:  SameCenter shell_block requirement is the default
          for(auto&& Xblk : shell_block_range(dfbs_, gbs_, 0, NoLastIndex, exact_diagonal_J_ ? SameCenter : NoRestrictions)){

            TimerHolder subtimer(ints_timer);
            auto g3 = ints_to_eigen_map(
                ish, jsh, Xblk,
                eris_3c_[ithr], coulomb_oper_type_,
                j_intbuff
            );

            subtimer.change(contract_timer);
            for(auto&& mu : function_range(ish)) {
              dt.segment(Xblk.bfoff, Xblk.nbf).transpose() += perm_fact * D.row(mu).segment(jsh.bfoff, jsh.nbf)
                  * g3.middleRows(mu.bfoff_in_shell*jsh.nbf, jsh.nbf);
              //----------------------------------------//
              J.block(mu, jsh.bfoff, 1, jsh.nbf).transpose() +=
                  g3.middleRows(mu.bfoff_in_shell*jsh.nbf, jsh.nbf)
                  * C_tilde.segment(Xblk.bfoff, Xblk.nbf);
            }

            /*-----------------------------------------------------*/
            /* Exact diagonal part                            {{{3 */ #if 3 // begin fold
            if(exact_diagonal_J_) {
              subtimer.change(ex_timer);

              const auto& Wij = W[{ish, jsh}];
              const auto& Wji = W[{jsh, ish}];

              if(Xblk.center == ish.center) {

                for(auto&& mu : function_range(ish)) {
                  for(auto&& nu : function_range(jsh)) {
                    // TODO remove one loop
                    dt_ex_thr.row(jsh.center).segment(ish.atom_dfbfoff + Xblk.bfoff_in_atom, Xblk.nbf) += 0.5 * perm_fact *
                        D.coeff(mu, nu) * g3.row(mu.bfoff_in_shell*jsh.nbf + nu.bfoff_in_shell);
                  }
                }

                for(auto&& mu : function_range(ish)) {
                  // TODO Remove loops where possible (may have to reconsider how Wij is stored)
                  const double delta_factor = ish.center != jsh.center ? 2.0 : 1.0;
                  J.block(mu, jsh.bfoff, 1, jsh.nbf).transpose() -= delta_factor
                      * g3.middleRows(mu.bfoff_in_shell*jsh.nbf, jsh.nbf)
                      * C_t_ex.row(jsh.center).segment(Xblk.bfoff, Xblk.nbf).transpose();

                  for(auto&& nu : function_range(jsh)) {
                    J.coeffRef(mu, nu) += delta_factor
                        * Wij.row(mu.bfoff_in_shell*jsh.nbf + nu.bfoff_in_shell).segment(Xblk.bfoff, Xblk.nbf)
                        * C_t_ex.row(jsh.center).segment(Xblk.bfoff, Xblk.nbf).transpose();
                    if(ish.center != jsh.center) {
                      J.coeffRef(mu, nu) += delta_factor
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
                        D.coeff(mu, nu) * g3.row(mu.bfoff_in_shell*jsh.nbf + nu.bfoff_in_shell);
                  }
                }

                for(auto&& mu : function_range(ish)) {
                  // TODO Remove loops where possible
                  J.block(mu, jsh.bfoff, 1, jsh.nbf).transpose() -= 2.0 *
                      g3.middleRows(mu.bfoff_in_shell*jsh.nbf, jsh.nbf)
                      * C_t_ex.row(ish.center).segment(Xblk.bfoff, Xblk.nbf).transpose();
                  for(auto&& nu : function_range(jsh)) {
                    J.coeffRef(mu, nu) += 2.0 * Wij.row(mu.bfoff_in_shell*jsh.nbf + nu.bfoff_in_shell).segment(Xblk.bfoff, Xblk.nbf)
                        * C_t_ex.row(ish.center).segment(Xblk.bfoff, Xblk.nbf).transpose();
                    J.coeffRef(mu, nu) += 2.0 * Wji.row(nu.bfoff_in_shell*ish.nbf + mu.bfoff_in_shell).segment(Xblk.bfoff, Xblk.nbf)
                        * C_t_ex.row(ish.center).segment(Xblk.bfoff, Xblk.nbf).transpose();
                  }
                }

              }
            } // end if exact_diagonal_J_
            /*******************************************************/ #endif //3}}}
            /*-----------------------------------------------------*/

          } // end loop over Xblk
        } // end while get_shell_pair
        mt_timer.change("misc", ithr);
        // only constructing the lower triangle of the J matrix, so zero the strictly upper part
        delete[] j_intbuff;
        //----------------------------------------//
        // add our contribution to the node level d_tilde
        mt_timer.change("sum thread contributions", ithr);
        {
          boost::lock_guard<boost::mutex> tmp_lock(tmp_mutex);
          d_tilde.noalias() += dt;
          if(exact_diagonal_J_) {
            d_t_ex.noalias() += dt_ex_thr;
          }
        }
        mt_timer.exit(ithr);
        /*******************************************************/ #endif //2}}}
        /*-----------------------------------------------------*/
      }); // end create_thread
    } // end enumeration of threads
    compute_threads.join_all();
    // Some diagonal shells might have basis functions that spill over into the
    //   upper triangle.  Zero these values to avoid double counting.
    J.triangularView<Eigen::StrictlyUpper>().setZero();
    mt_timer.exit();
    timer.insert(mt_timer);
    if(print_iteration_timings_) mt_timer.print(ExEnv::out0(), 12, 45);
  } // compute_threads is destroyed here
  //----------------------------------------//
  // Global sum d_tilde
  timer.enter("global sum d_tilde");
  scf_grp_->sum((double*)d_tilde.data(), dfnbf);
  if(exact_diagonal_J_) {
    scf_grp_->sum((double*)d_t_ex.data(), natom * dfnbf);
  }
  timer.exit();
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
              auto& Ca = *(cpair.first);
              auto& Cb = *(cpair.second);
              //----------------------------------------//
              // First term contribution
              J.coeffRef(mu, nu) += d_tilde.segment(ish.atom_dfbfoff, ish.atom_dfnbf).transpose() * Ca;
              if(ish.center != jsh.center){
                J.coeffRef(mu, nu) += d_tilde.segment(jsh.atom_dfbfoff, jsh.atom_dfnbf).transpose() * Cb;
              }
              //----------------------------------------//
              // Third term contribution
              J.coeffRef(mu, nu) -= g_tilde.segment(ish.atom_dfbfoff, ish.atom_dfnbf).transpose() * Ca;
              if(ish.center != jsh.center){
                J.coeffRef(mu, nu) -= g_tilde.segment(jsh.atom_dfbfoff, jsh.atom_dfnbf).transpose() * Cb;
              }
              //----------------------------------------//
              if(exact_diagonal_J_) {
                J.coeffRef(mu, nu) -= 2.0 * d_t_ex.row(jsh.center).segment(ish.atom_dfbfoff, ish.atom_dfnbf) * Ca;
                if(ish.center != jsh.center) {
                  J.coeffRef(mu, nu) -= 2.0 * d_t_ex.row(ish.center).segment(jsh.atom_dfbfoff, jsh.atom_dfnbf) * Cb;
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
            for(auto&& ksh : iter_shells_on_center(gbs_, ish.center)) {
              for(auto&& lsh : iter_shells_on_center(gbs_, jsh.center)) {
                if(not is_sig_pair(ksh, lsh)) continue;

                auto g4_ptr = ints_to_eigen(ish, jsh, ksh, lsh, tbis_[ithr], coulomb_oper_type_);
                const auto& g4 = *g4_ptr;
                for(auto&& mu : function_range(ish)) {
                  for(auto&& rho : function_range(ksh)) {
                    // TODO More Vectorization
                    J.block(mu, jsh.bfoff, 1, jsh.nbf).transpose() += epf *
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
  scf_grp_->sum(J_data, nbf*nbf);
  //----------------------------------------//
  // Fill in the upper triangle of J
  for(int mu = 0; mu < nbf; ++mu) {
    for(int nu = 0; nu < mu; ++nu) {
      J.coeffRef(nu, mu) = J.coeff(mu, nu);
    }
  }
  if(debug_coulomb_energy_) {
    double jenergy = (J.array() * D.array()).sum();
    ExEnv::out0() << indent << "Coulomb energy: " << std::setprecision(12) << jenergy << std::endl;
  }
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Transfer J to a RefSCMatrix                           		                        {{{1 */ #if 1 // begin fold
  Ref<Integral> localints = integral()->clone();
  localints->set_basis(gbs_);
  Ref<PetiteList> pl = localints->petite_list();
  RefSCDimension obsdim = pl->AO_basisdim();
  RefSCMatrix result(
      obsdim,
      obsdim,
      gbs_->so_matrixkit()
  );
  result.assign(J_data);
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Clean up                                             		                        {{{1 */ #if 1 // begin fold
  //----------------------------------------//
  //deallocate(D_ptr);
  delete[] J_data;
  delete[] D_ptr;
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  return result;
}


