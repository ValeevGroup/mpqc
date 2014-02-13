//
// coefficients.cc
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


#include <util/misc/xmlwriter.h>

#include "cadfclhf.h"

using namespace sc;

typedef std::pair<int, int> IntPair;

void
CADFCLHF::compute_coefficients()
{                                                                                          //latex `\label{sc:compute_coefficients}`
  /*=======================================================================================*/
  /* Setup                                                		                        {{{1 */ #if 1 // begin fold
  //---------------------------------------------------------------------------------------//
  // References for speed
  Timer timer("compute coefficients");
  const Ref<GaussianBasisSet>& obs = gbs_;
  // Constants for convenience
  const int nbf = obs->nbasis();
  const int dfnbf = dfbs_->nbasis();
  const int natom = obs->ncenter();
  /*-----------------------------------------------------*/
  /* Initialize coefficient memory                  {{{2 */ #if 2 //latex `\label{sc:coefmem}`
  // Now initialize the coefficient memory.
  // First, compute the amount of memory needed
  // Coefficients will be stored jbf <= ibf
  int ncoefs = 0;                                                                          //latex `\label{sc:coefcountbegin}`
  for(auto ibf : function_range(obs, dfbs_)){
    for(auto jbf : function_range(obs, dfbs_, 0, ibf)){
      ncoefs += ibf.atom_dfnbf;
      if(ibf.center != jbf.center){
        ncoefs += jbf.atom_dfnbf;
      }
    }
  }                                                                                        //latex `\label{sc:coefcountend}`
  coefficients_data_ = allocate<double>(ncoefs);                                           //latex `\label{sc:coefalloc}`
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
  /********************************************************/ #endif //2}}} //latex `\label{sc:coefmemend}`
  /*-----------------------------------------------------*/
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Compute the coefficients in threads                   		                        {{{1 */ #if 1 //latex `\label{sc:coefloop}`
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
        while(get_shell_pair(ish, jsh)){                                                   //latex `\label{sc:coefgetpair}`
          const int dfnbfAB = ish.center == jsh.center ? ish.atom_dfnbf : ish.atom_dfnbf + jsh.atom_dfnbf;
          //----------------------------------------//
          std::shared_ptr<Decomposition> decomp = get_decomposition(                       //latex `\label{sc:coefgetdecomp}`
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
          for(auto mu : function_range(ish)){                                             //latex `\label{sc:coefsolveloop}`
            // Since coefficients are only stored jbf <= ibf,
            //   we need to figure out how many jbf's to do
            for(auto nu : function_range(jsh)){
              if(ish == jsh and nu > mu) continue;
              const int ijbf = mu.bfoff_in_shell*jsh.nbf + nu.bfoff_in_shell;
              IntPair mu_nu(mu, nu);
              Eigen::VectorXd Ctmp(dfnbfAB);
              Ctmp = decomp->solve(ij_M_X.row(ijbf).transpose());                         //latex `\label{sc:coefsolve}`
              assert((coefs_.find(mu_nu) != coefs_.end())
                  || ((cout << "couldn't find coefficients for " << mu << ", " << nu << endl), false));
              *(coefs_[mu_nu].first) = Ctmp.head(ish.atom_dfnbf);
              if(ish.center != jsh.center){
                *(coefs_[mu_nu].second) = Ctmp.tail(jsh.atom_dfnbf);
              }
            } // end jbf loop
          } // end ibf loop                                                               //latex `\label{sc:coefsolveloopend}`
        } // end while get_shell_pair                                                     //latex `\label{sc:coefendgetpair}`
        /********************************************************/ #endif //2}}}
        /*-----------------------------------------------------*/
      });  // End threaded compute function and compute_threads.create_thread() call
    } // end loop over number of threads
    compute_threads.join_all();
  } // thread_group compute_threads is destroyed and goes out of scope, threads are destroyed (?)
  /*****************************************************************************************/ #endif //1}}} //latex `\label{sc:coefloopend}`
  /*=======================================================================================*/
  /* Global sum coefficient memory                        		                        {{{1 */ #if 1 // begin fold
  //---------------------------------------------------------------------------------------//
  scf_grp_->sum(coefficients_data_, ncoefs);                                               //latex `\label{sc:coefsum}`
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Store the transpose coefficients                     		                        {{{1 */ #if 1 //latex `\label{sc:coeftrans}`
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
  /*****************************************************************************************/ #endif //1}}} //latex `\label{sc:coeftransend}`
  /*=======================================================================================*/
  /* Debugging output                                     		                        {{{1 */ #if 1 // begin fold
  if(xml_debug_) {
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
    //----------------------------------------//
    // Compute Cmaxes_

    Cmaxes_.resize(dfbs_->nshell());
    for(auto Xsh : shell_range(dfbs_)){

      Cmaxes_[Xsh].resize(obs->nshell());
      auto& Cmaxes_X = Cmaxes_[Xsh];

      for(auto lsh : shell_range(obs)) {

        int max_index;
        double max_val;

        if(lsh.center == Xsh.center) {
          // We might be able to get away with only nu on same center as X and sigma,
          //  but for now approach it more rigorously
          if(use_norms_nu_){
            max_val = C_trans_frob_[Xsh].row(lsh.shoff_in_atom).norm();
            max_index = -1;
          }
          else {
            max_val = C_trans_frob_[Xsh].row(lsh.shoff_in_atom).maxCoeff(&max_index);
          }
        }
        else { // different centers...
          if(use_norms_nu_){
            max_val = C_trans_frob_[Xsh].col(lsh).norm();
            max_index = -1;
          }
          else{
            max_val = C_trans_frob_[Xsh].col(lsh).maxCoeff(&max_index);
            max_index += lsh.atom_shoff;
          }
        }
        // Note:  the index attribute is unused and deprecated.
        // TODO change this to just a double rather than an index and value pair
        Cmaxes_X[lsh].index = max_index;
        Cmaxes_X[lsh].value = max_val;

      } // end loop over lsh

    } // end loop over Xsh

    //----------------------------------------//
  } // end if do_linK
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Clean up                                              		                        {{{1 */ #if 1 // begin fold
  //---------------------------------------------------------------------------------------//
  have_coefficients_ = true;                                                               //latex `\label{sc:coefflag}`
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
}

//////////////////////////////////////////////////////////////////////////////////
