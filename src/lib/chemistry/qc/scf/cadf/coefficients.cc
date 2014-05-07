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
#include <util/container/conc_cache.h>

#include "cadfclhf.h"
#include "assignments.h"
#include "treemat.h"
#include "ordered_shells.h"

using namespace sc;

typedef std::pair<int, int> IntPair;
typedef uint64_t uli;
typedef unsigned int uint;


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
  const cadf::Node& my_part = atom_pair_assignments_k_->my_assignments(scf_grp_->me());

  /*-----------------------------------------------------*/
  /* Initialize coefficient memory                  {{{2 */ #if 2 //latex `\label{sc:coefmem}`
  // Now initialize the coefficient memory.
  // First, compute the amount of memory needed

  timer.enter("01 - init coef memory");

  uli ncoefs = 0;                                                                          //latex `\label{sc:coefcountbegin}`

  std::vector<uli> offsets;
  std::vector<uli> block_sizes;
  uli ioffset = 0;

  if(distribute_coefficients_) {
    local_pairs_spot_ = 0;
    ShellData ish, jsh;
    while(get_shell_pair(ish, jsh, SignificantPairs)) {
      ncoefs += ish.nbf * jsh.nbf * (ish.atom_dfnbf + (ish.center == jsh.center ? 0 : jsh.atom_dfnbf));
    }

    ExEnv::out0() << indent << "Need " << data_size_to_string(ncoefs*sizeof(double))
                  << " for main coefficients storage." << std::endl;
    memory_used_ += ncoefs*sizeof(double);
    ExEnv::out0() << indent << "Total memory usage is at least "
                  << data_size_to_string(memory_used_) << std::endl;

    coefficients_data_ = new double[ncoefs];

    memset(coefficients_data_, 0, ncoefs * sizeof(double));

    double* spot = coefficients_data_;
    local_pairs_spot_ = 0;
    while(get_shell_pair(ish, jsh, SignificantPairs)) {
      for(auto&& mu : function_range(ish)) {
        for(auto&& nu : function_range(jsh)) {
          if(nu > mu) continue;
          CoefContainer coefs_a = make_shared<CoefView>(
              spot,
              ish.atom_dfnbf,
              Eigen::Stride<1, Eigen::Dynamic>(1, 1)
          );
          CoefContainer coefs_b = make_shared<CoefView>(
              spot + (ish.center == jsh.center ? 0 : ish.atom_dfnbf),
              jsh.atom_dfnbf,
              Eigen::Stride<1, Eigen::Dynamic>(1, 1)
          );

          coefs_.emplace(std::make_pair(int(mu), int(nu)), std::make_pair(coefs_a, coefs_b));
          if(mu != nu) {
            coefs_.emplace(std::make_pair(int(nu), int(mu)), std::make_pair(coefs_b, coefs_a));
          }
          spot += ish.atom_dfnbf;
          if(mu.center != nu.center) spot += jsh.atom_dfnbf;
        }
      }
    }

  }
  else {
    for(int iatom = 0; iatom < natom; ++iatom) {
      const uli atom_nbf = gbs_->nbasis_on_center(iatom);
      const uli atom_dfnbf = dfbs_->nbasis_on_center(iatom);
      const uli block_size = atom_nbf * nbf * atom_dfnbf;
      offsets.push_back(ioffset);
      ioffset += block_size;
      block_sizes.push_back(block_size);
    }
    ncoefs = ioffset;

    ExEnv::out0() << indent << "Need " << data_size_to_string(ncoefs*sizeof(double))
                  << " for main coefficients storage." << std::endl;
    memory_used_ += ncoefs*sizeof(double);
    ExEnv::out0() << indent << "Total memory usage is at least "
                  << data_size_to_string(memory_used_) << std::endl;

    coefficients_data_ = new double[ncoefs];                                           //latex `\label{sc:coefalloc}`

    memset(coefficients_data_, 0, ncoefs * sizeof(double));

    coefs_transpose_blocked_.reserve(natom);
    coefs_transpose_blocked_other_.reserve(natom);
    coefs_transpose_.reserve(dfnbf);

    for(int iatom = 0; iatom < natom; ++iatom){

      const int atom_nbf = gbs_->nbasis_on_center(iatom);
      coefs_transpose_blocked_.emplace_back(
          coefficients_data_ + offsets[iatom],
          dfbs_->nbasis_on_center(iatom),
          atom_nbf * nbf
      );
      coefs_transpose_blocked_other_.emplace_back(
          coefficients_data_ + offsets[iatom],
          dfbs_->nbasis_on_center(iatom) * atom_nbf,
          nbf
      );

      // The transpose coefficients
      auto& cblock = coefs_transpose_blocked_[iatom];
      for(auto&& X : iter_functions_on_center(dfbs_, iatom)) {
        double* offset = cblock.data() + X.bfoff_in_atom * nbf * atom_nbf;
        coefs_transpose_.emplace_back(
            offset,
            atom_nbf, nbf
        );

        memory_used_ += sizeof(Eigen::Map<RowMatrix>);
      }

    }

    for(auto mu : function_range(gbs_, dfbs_)){
      for(auto rho : function_range(gbs_, dfbs_)){
        double *data_spot_a = coefficients_data_ + offsets[mu.center] + mu.bfoff_in_atom*nbf + rho;
        double *data_spot_b = coefficients_data_ + offsets[rho.center] + rho.bfoff_in_atom*nbf + mu;
        IntPair ij(mu, rho);
        CoefContainer coefs_a = make_shared<CoefView>(
            data_spot_a, mu.atom_dfnbf, Eigen::Stride<1, Eigen::Dynamic>(1, mu.atom_nbf * nbf)
        );
        CoefContainer coefs_b = make_shared<CoefView>(
            data_spot_b, rho.atom_dfnbf, Eigen::Stride<1, Eigen::Dynamic>(1, rho.atom_nbf * nbf)
        );
        coefs_.emplace(ij, std::make_pair(coefs_a, coefs_b));
      }
    }
  }


  /********************************************************/ #endif //2}}} //latex `\label{sc:coefmemend}`
  /*-----------------------------------------------------*/

  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Compute the coefficients in threads                   		                        {{{1 */ #if 1 //latex `\label{sc:coefloop}`
  //---------------------------------------------------------------------------------------//

  ExEnv::out0() << indent << "Computing coefficients" << std::endl;
  // reset the iteration over local pairs
  local_pairs_spot_ = 0;
  timer.change("02 - compute coefficients");
  {
    boost::thread_group compute_threads;
    // Loop over number of threads
    for(int ithr = 0; ithr < nthread_; ++ithr) {
      // ...and create each thread that computes pairs
      compute_threads.create_thread([&,ithr](){

        ShellData ish, jsh;
        while(get_shell_pair(ish, jsh, SignificantPairs)){                                 //latex `\label{sc:coefgetpair}`

          std::vector<CoefView> coefsA;
          std::vector<CoefView> coefsB;

          for(auto&& mu : function_range(ish)) {
            for(auto&& nu : function_range(jsh)) {
              auto& cmn = coefs_[{mu, nu}];
              coefsA.push_back(*cmn.first);
              if(ish.center != jsh.center) {
                coefsB.push_back(*cmn.second);
              }
            }
          }
          get_coefs_ish_jsh(ish, jsh, ithr, coefsA, coefsB);

          if(ish != jsh and !distribute_coefficients_){
            coefsA.clear();
            coefsB.clear();
            for(auto&& nu : function_range(jsh)) {
              for(auto&& mu : function_range(ish)) {
                auto& cnm = coefs_[{nu, mu}];
                if(ish.center != jsh.center) {
                  coefsA.push_back(*cnm.second);
                }
                coefsB.push_back(*cnm.first);
              }
            }
            get_coefs_ish_jsh(jsh, ish, ithr, coefsB, coefsA);
          }

        } // end while get_shell_pair                                                     //latex `\label{sc:coefendgetpair}`

      });  // End threaded compute function and compute_threads.create_thread() call
    } // end loop over number of threads
    compute_threads.join_all();
  } // thread_group compute_threads is destroyed and goes out of scope

  /*****************************************************************************************/ #endif //1}}} //latex `\label{sc:coefloopend}`
  /*=======================================================================================*/
  /* Global sum coefficient memory                        		                        {{{1 */ #if 1 // begin fold
  //---------------------------------------------------------------------------------------//

  ExEnv::out0() << indent << "Distributing coefficients" << std::endl;
  timer.change("03 - global sum coefficient memory");
  if(not distribute_coefficients_) {
    if(ncoefs * sizeof(double) < std::numeric_limits<int>::max()) {
      scf_grp_->sum(coefficients_data_, ncoefs);                                               //latex `\label{sc:coefsum}`
    }
    else {
      int chunk_size = std::numeric_limits<int>::max() / sizeof(double);
      for(int ichunk = 0; ichunk < ncoefs / chunk_size; ++ichunk) {
        scf_grp_->sum(coefficients_data_ + ichunk*chunk_size, chunk_size);
      }
      if(ncoefs % chunk_size > 0) {
        scf_grp_->sum(coefficients_data_ + (ncoefs / chunk_size) * chunk_size, ncoefs % chunk_size);
      }
    }
  }

  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Compute the distributed coefficients                  		                        {{{1 */ #if 1 //latex `\label{sc:coeftrans}`
  //---------------------------------------------------------------------------------------//
  if(distribute_coefficients_) {
    timer.change("04 - distributed coefficients");
    uli n_Xcoefs = 0;
    uli ncoefs_dist = my_part.bin->obs_ncoefs + my_part.bin->dfbs_ncoefs;
    try {
      if(ncoefs_dist * sizeof(double) > std::numeric_limits<int>::max()) {
        dist_coefs_data_ = new double[my_part.bin->obs_ncoefs];
        dist_coefs_data_df_ = new double[my_part.bin->dfbs_ncoefs];
      }
      else {
        dist_coefs_data_ = new double[ncoefs_dist];
        dist_coefs_data_df_ = dist_coefs_data_ + my_part.bin->obs_ncoefs;
      }
      memory_used_ += ncoefs_dist * sizeof(double);
    }
    catch(std::bad_alloc& e) {
      ExEnv::outn() << "Failed to allocate " << data_size_to_string(ncoefs_dist*sizeof(double))
                    << " for distributed coefficients on node " << scf_grp_->me() << std::endl;
      ExEnv::outn() << "Before allocation, memory in use was at least " << data_size_to_string(memory_used_.load()) << std::endl;
      throw;
    }
    memset(dist_coefs_data_, 0, my_part.bin->obs_ncoefs*sizeof(double));
    memset(dist_coefs_data_df_, 0, my_part.bin->dfbs_ncoefs*sizeof(double));
    ExEnv::out0() << indent << "Allocated " << data_size_to_string(ncoefs_dist*sizeof(double))
                  << " (on node 0) for distributed coefficients." << std::endl;
    ExEnv::out0() << indent << "Total memory usage is now at least " << data_size_to_string(memory_used_.load()) << std::endl;

    // Initialize the Eigen::Maps of data parts
    for(auto&& obs_shell : my_part.bin->assigned_obs_shells) {
      ShellData ish(obs_shell->index, gbs_, dfbs_);
      coefs_mu_X.emplace(
          std::piecewise_construct,
          std::forward_as_tuple(ish),
          std::forward_as_tuple(
              dist_coefs_data_ + my_part.bin->obs_coef_offsets.at(ish),
              ish.nbf, nbf*ish.atom_dfnbf
          )
      );
      coefs_mu_X_other.emplace(
          std::piecewise_construct,
          std::forward_as_tuple(ish),
          std::forward_as_tuple(
              dist_coefs_data_ + my_part.bin->obs_coef_offsets.at(ish),
              ish.nbf*nbf, ish.atom_dfnbf
          )
      );
    }
    for(auto&& dfbs_atom : my_part.bin->assigned_dfbs_atoms) {
      ShellBlockData<> Xblk = ShellBlockData<>::atom_block(dfbs_atom->index, dfbs_, gbs_);
      coefs_X_nu.emplace(
          std::piecewise_construct,
          std::forward_as_tuple(Xblk.center),
          std::forward_as_tuple(
              dist_coefs_data_df_ + my_part.bin->dfbs_coef_offsets.at(Xblk.center),
              Xblk.atom_nbf, Xblk.atom_obsnbf*nbf
          )
      );
      coefs_X_nu_other.emplace(
          std::piecewise_construct,
          std::forward_as_tuple(Xblk.center),
          std::forward_as_tuple(
              dist_coefs_data_df_ + my_part.bin->dfbs_coef_offsets.at(Xblk.center),
              Xblk.atom_nbf*Xblk.atom_obsnbf, nbf
          )
      );
    }

    do_threaded(nthread_, [&](int ithr) {
      std::vector<Eigen::Map<Eigen::VectorXd>> empty;
      for(auto&& obs_shell : thread_over_range(
          my_part.compute_coef_items[false], ithr, nthread_
      )) {
        // Do the C_mu_X part first
        ShellData ish(obs_shell->index, gbs_, dfbs_);
        for(auto&& jsh : iter_significant_partners(ish)) {
          std::vector<Eigen::Map<Eigen::VectorXd>> coefs;
          for(auto&& mu : function_range(ish)) {
            for(auto&& rho : function_range(jsh)) {
              coefs.emplace_back(
                  // NOTE: POSSIBLE DATA CONCURRENCY ISSUES!
                  coefs_mu_X.at(ish).data()
                    + mu.off*nbf*ish.atom_dfnbf + rho*ish.atom_dfnbf,
                  ish.atom_dfnbf
              );
            }
          }
          get_coefs_ish_jsh(ish, jsh, ithr, coefs, empty);
        }
      }
    });

    // Node-row-wise sum of mu coefficients
    {
      Ref<MessageGrp> mu_grp = scf_grp_->split(my_part.bin->obs_row_id);
      const uli ncfs = my_part.bin->obs_ncoefs;
      if(ncfs * sizeof(double) < std::numeric_limits<int>::max()) {
        mu_grp->sum(dist_coefs_data_, ncfs);
      }
      else {
        int chunk_size = std::numeric_limits<int>::max() / sizeof(double);
        for(int ichunk = 0; ichunk < ncfs / chunk_size; ++ichunk) {
          mu_grp->sum(dist_coefs_data_ + ichunk*chunk_size, chunk_size);
        }
        if(ncfs % chunk_size > 0) {
          mu_grp->sum(dist_coefs_data_ + (ncfs / chunk_size) * chunk_size,
              ncfs % chunk_size
          );
        }
      }
    } // mu_grp is deleted

    sc::SCFormIO::init_mp(scf_grp_->me());

    std::vector<CoefView> empty_df;
    do_threaded(nthread_, [&](int ithr) {
      for(auto&& dfbs_atom : thread_over_range(
          my_part.compute_coef_items[true],
          ithr, nthread_
      )) {

        // Now do the C_X_mu part
        ShellBlockData<> Xblk = ShellBlockData<>::atom_block(dfbs_atom->index, gbs_, dfbs_);
        for(auto&& ish : iter_shells_on_center(gbs_, Xblk.center, dfbs_)) {
          for(auto&& jsh : iter_significant_partners(ish)) {
            std::vector<CoefView> coefs;
            for(auto&& mu : function_range(ish)) {
              for(auto&& rho : function_range(jsh)) {
                coefs.emplace_back(
                    // NOTE: POSSIBLE DATA CONCURRENCY ISSUES!
                    coefs_X_nu.at(Xblk.center).data() + mu.bfoff_in_atom*nbf + rho,
                    ish.atom_dfnbf,
                    Eigen::Stride<1, Eigen::Dynamic>(1, mu.atom_nbf * nbf)
                );
              }
            }
            get_coefs_ish_jsh(ish, jsh, ithr, coefs, empty_df);
          }
        }
      }
    });

    // Node-row-wise sum of X coefficients
    {
      Ref<MessageGrp> X_grp = scf_grp_->split(my_part.bin->dfbs_row_id);
      const uli ncfs = my_part.bin->dfbs_ncoefs;
      if(ncfs * sizeof(double) < std::numeric_limits<int>::max()) {
        X_grp->sum(dist_coefs_data_df_, ncfs);
      }
      else {
        int chunk_size = std::numeric_limits<int>::max() / sizeof(double);
        for(int ichunk = 0; ichunk < ncfs / chunk_size; ++ichunk) {
          X_grp->sum(dist_coefs_data_df_ + ichunk*chunk_size, chunk_size);
        }
        if(ncfs % chunk_size > 0) {
          X_grp->sum(dist_coefs_data_df_ + (ncfs / chunk_size) * chunk_size,
              ncfs % chunk_size
          );
        }
      }
    } // X_grp is deleted

    sc::SCFormIO::init_mp(scf_grp_->me());
  }

  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Store the transpose and blocked coefficients          		                        {{{1 */ #if 1 //latex `\label{sc:coeftrans}`
  //---------------------------------------------------------------------------------------//

  // TODO Get rid of this!!!
  timer.change("05 - store blocked coefs");
  if(store_coefs_transpose_ and not distribute_coefficients_) {
    coef_block_offsets_.resize(natom);
    for(int iatom = 0; iatom < natom; ++iatom) {

      const int atom_nbf = gbs_->nbasis_on_center(iatom);
      const int atom_dfnbf = dfbs_->nbasis_on_center(iatom);

      timer.enter("count columns");
      // First, figure out the number of columns in the matrix
      int ncol = 0;
      coef_block_offsets_[iatom].resize(natom);
      int current_atom = -1;
      for(auto&& jsh : shell_range(gbs_, dfbs_)) {
        if(current_atom != jsh.center) {
          current_atom = jsh.center;
          coef_block_offsets_[iatom][current_atom] = ncol;
        }
        ncol += jsh.nbf * jsh.atom_dfnbf;
        if(jsh.center != iatom) {
          ncol += jsh.nbf * atom_dfnbf;
        }
      }

      timer.change("compute coefficients");
      coefs_blocked_.emplace_back(atom_nbf, ncol);
      auto& cblock = coefs_blocked_.back();
      cblock = RowMatrix::Zero(atom_nbf, ncol);

      // Now transfer the coefficients
      for(auto&& mu : iter_functions_on_center(gbs_, iatom, dfbs_)) {
        int offset = 0;
        for(auto&& rho : function_range(gbs_, dfbs_)) {
          if(rho <= mu) {
            auto& cpair = coefs_[{mu, rho}];
            cblock.row(mu.bfoff_in_atom).segment(offset, mu.atom_dfnbf) = *(cpair.first);
            offset += mu.atom_dfnbf;
            if(mu.center != rho.center) {
              cblock.row(mu.bfoff_in_atom).segment(offset, rho.atom_dfnbf) = *(cpair.second);
              offset += rho.atom_dfnbf;
            }
          }
          else{
            auto& cpair = coefs_[{rho, mu}];
            if(mu.center == rho.center) {
              cblock.row(mu.bfoff_in_atom).segment(offset, mu.atom_dfnbf) = *(cpair.first);
              offset += mu.atom_dfnbf;
            }
            else {
              cblock.row(mu.bfoff_in_atom).segment(offset, mu.atom_dfnbf) = *(cpair.second);
              offset += mu.atom_dfnbf;
              cblock.row(mu.bfoff_in_atom).segment(offset, rho.atom_dfnbf) = *(cpair.first);
              offset += rho.atom_dfnbf;
            }
          }

        }
      }

      timer.exit();

    }
  }
  /*****************************************************************************************/ #endif //1}}} //latex `\label{sc:coeftransend}`
  /*=======================================================================================*/
  /* Debugging output                                     		                        {{{1 */ #if 1 // begin fold
  // debugging not implemented for sparse yet
  if(xml_debug_) {
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
        write_as_xml(
            "coefficient_vector",
            coefs,
            std::map<std::string, int>{
              { "ao_index1", nu },
              { "ao_index2", mu }
            }
        );
      }
    }
    //end_xml_context("df_coefficients");
  }
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Make the CADF-LinK lists                                                         {{{1 */ #if 1 // begin fold
  timer.change("06 - LinK coef lists");

  if(do_linK_) {

    //----------------------------------------//
    // Compute the Frobenius norm of C_transpose_ blocks
    timer.change("07 - C_trans_frob");

    // Since C_trans_frob is only used here, we don't need it to be a class-scope member now.
    std::vector<Eigen::MatrixXd> C_trans_frob(dfbs_->nshell());

    if(distribute_coefficients_) {
      for(auto Xatom : my_part.bin->assigned_dfbs_atoms) {
        for(auto&& Xsh : iter_shells_on_center(dfbs_, Xatom->index, gbs_)) {

          resize_and_zero_matrix(C_trans_frob[Xsh], Xsh.atom_obsnsh, obs->nshell());

          for(auto&& ish : iter_shells_on_center(obs, Xsh.center)) {
            for(auto&& jsh : shell_range(obs)) {
              for(auto&& mu : function_range(ish)) {
                C_trans_frob[Xsh](ish.shoff_in_atom, jsh) += coefs_X_nu.at(Xsh.center).block(
                    Xsh.bfoff_in_atom, mu.bfoff_in_atom*nbf + jsh.bfoff,
                    Xsh.nbf, jsh.nbf
                ).squaredNorm();
              }
            }
          }

          C_trans_frob[Xsh] = C_trans_frob[Xsh].array().sqrt();
        }
      }
    }
    else {
      for(auto&& Xsh : shell_range(dfbs_, obs)) {

        resize_and_zero_matrix(C_trans_frob[Xsh], Xsh.atom_obsnsh, obs->nshell());

        for(auto&& ish : iter_shells_on_center(obs, Xsh.center)) {
          for(auto&& jsh : shell_range(obs)) {
            for(auto&& mu : function_range(ish)) {
              C_trans_frob[Xsh](ish.shoff_in_atom, jsh) += coefs_transpose_blocked_[Xsh.center].block(
                  Xsh.bfoff_in_atom, mu.bfoff_in_atom*nbf + jsh.bfoff,
                  Xsh.nbf, jsh.nbf
              ).squaredNorm();
            }
          } // end loop over jsh
        } // end loop over ish

        C_trans_frob[Xsh] = C_trans_frob[Xsh].array().sqrt();

      } // end loop over Xsh
    }
    //----------------------------------------//

    // Compute Cmaxes_
    timer.change("08 - Cbar");
    resize_and_zero_matrix(C_bar_, gbs_->nshell(), dfbs_->nshell());
    resize_and_zero_matrix(C_bar_mine_, gbs_->nshell(), my_part.bin->dfnsh());
    if(distribute_coefficients_) {
      resize_and_zero_matrix(C_underbar_, gbs_->nshell(), dfbs_->nshell());
    }
    // TODO switch C_dfsame_ to offset indices
    if(use_norms_nu_) {
      uint Xsh_off = 0;
      for(auto&& Xsh_index : my_part.bin->assigned_dfbs_shells) {
        ShellData Xsh(Xsh_index, dfbs_, gbs_);
        {
          const auto& sqnorm1 = C_trans_frob[Xsh].colwise().squaredNorm();
          if(distribute_coefficients_) {
            C_bar_.col(Xsh) = sqnorm1;
            C_underbar_.col(Xsh) = sqnorm1;
          }
          C_bar_mine_.col(Xsh_off) = sqnorm1;
        }
        {
          const auto& sqnorm2 = C_trans_frob[Xsh].rowwise().squaredNorm();
          if(distribute_coefficients_) C_bar_.col(Xsh).segment(Xsh.atom_obsshoff, Xsh.atom_obsnsh) += sqnorm2;
          C_bar_mine_.col(Xsh_off).segment(Xsh.atom_obsshoff, Xsh.atom_obsnsh) += sqnorm2;
        }
        {
          const auto& sqnorm3 = C_trans_frob[Xsh].middleCols(
              Xsh.atom_obsshoff, Xsh.atom_obsnsh
          ).rowwise().squaredNorm();
          if(distribute_coefficients_) C_bar_.col(Xsh).segment(Xsh.atom_obsshoff, Xsh.atom_obsnsh) -= sqnorm3;
          C_bar_mine_.col(Xsh_off).segment(Xsh.atom_obsshoff, Xsh.atom_obsnsh) -= sqnorm3;
        }
        ++Xsh_off;
      }
      if(distribute_coefficients_) {
        C_bar_ = C_bar_.array().sqrt();
        C_underbar_ = C_underbar_.array().sqrt();
      }
      C_bar_mine_ = C_bar_mine_.array().sqrt();

      if(distribute_coefficients_) {
        // Node-row-wise sum of X parts of C_bar_
        {
          Ref<MessageGrp> mu_grp = scf_grp_->split(my_part.bin->obs_row_id);
          mu_grp->sum(C_bar_.data(), gbs_->nshell() * dfbs_->nshell());
          mu_grp->sum(C_underbar_.data(), dfbs_->nshell() * gbs_->nshell());
        } // mu_grp is deleted
        sc::SCFormIO::init_mp(scf_grp_->me());
      }
      else {
        for(ShellData&& Xsh : shell_range(dfbs_, gbs_)) {
          C_bar_.col(Xsh) = C_trans_frob[Xsh].colwise().squaredNorm();
          C_bar_.col(Xsh).segment(Xsh.atom_obsshoff, Xsh.atom_obsnsh) +=
              C_trans_frob[Xsh].rowwise().squaredNorm();
          C_bar_.col(Xsh).segment(Xsh.atom_obsshoff, Xsh.atom_obsnsh) -=
              C_trans_frob[Xsh].middleCols(Xsh.atom_obsshoff, Xsh.atom_obsnsh).rowwise().squaredNorm();
        }
        C_bar_ = C_bar_.array().sqrt();
      }
      if(distribute_coefficients_) {
        C_dfsame_ = std::make_shared<typename decltype(C_dfsame_)::element_type>(C_underbar_, gbs_);
        for(auto&& Xsh_index : my_part.bin->assigned_dfbs_shells) {
          auto& L_C_under_X = L_C_under[Xsh_index];
          L_C_under_X.acquire_and_sort(
              C_underbar_.data() + Xsh_index * gbs_->nshell(),
              gbs_->nshell(), 0.0, true
          );
          L_C_under_X.set_basis(gbs_, dfbs_);
        }
      }
    }
    else {
      if(distribute_coefficients_) {
        if(screen_B_) {
          throw FeatureNotImplemented("use_norms_nu = false with distributed coefficients and screen_B", __FILE__, __LINE__, class_desc());
        }
        throw FeatureNotImplemented("use_norms_nu = false with distributed coefficients", __FILE__, __LINE__, class_desc());

      }

      // TODO just do this in the above loop similarly
      for(auto Xsh : shell_range(dfbs_)){
        for(auto lsh : shell_range(obs)) {

          double max_val;

          if(lsh.center == Xsh.center) {
            // We might be able to get away with only nu on same center as X and sigma,
            //  but for now approach it more rigorously
            max_val = C_trans_frob[Xsh].row(lsh.shoff_in_atom).maxCoeff();
          }
          else { // different centers...
            max_val = C_trans_frob[Xsh].col(lsh).maxCoeff();
          }

          C_bar_(lsh, Xsh) = max_val;

        } // end loop over lsh
      } // end loop over Xsh

    }

    //----------------------------------------//
  } // end if do_linK
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Sparse coefficients                                   		                        {{{1 */ #if 1 // begin fold
  //---------------------------------------------------------------------------------------//

  //if(use_sparse_) {
  //  std::vector<Eigen::Triplet<double>> triplets;
  //  for(auto&& ckvpair : coefs_) {
  //    const BasisFunctionData mu(ckvpair.first.first, gbs_, dfbs_);
  //    const BasisFunctionData nu(ckvpair.first.second, gbs_, dfbs_);
  //    const auto& Ca = *ckvpair.second.first;
  //    const auto& Cb = *ckvpair.second.second;
  //    for(int Xoff = 0; Xoff < mu.atom_dfnbf; ++Xoff) {
  //      const double val = Ca[Xoff];
  //      triplets.emplace_back(mu.index, nu.index*dfnbf + mu.atom_dfnbf + Xoff, val);
  //    }
  //    coefs_sp_.resize(nbf, nbf*dfnbf);

  //  }
  //}



  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /*=======================================================================================*/
  /* Clean up                                              		                        {{{1 */ #if 1 // begin fold
  //---------------------------------------------------------------------------------------//
  decomps_->clear();
  have_coefficients_ = true;                                                               //latex `\label{sc:coefflag}`
  timer.exit();
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
}

//////////////////////////////////////////////////////////////////////////////////
