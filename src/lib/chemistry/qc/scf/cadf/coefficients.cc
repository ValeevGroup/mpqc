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
#include "assignments.h"

using namespace sc;

typedef std::pair<int, int> IntPair;
typedef unsigned long uli;


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

  timer.enter("compute (X|Y)");
  ExEnv::out0() << indent << "Computing two center integrals" << std::endl;
  g2_full_ptr_ = ints_to_eigen_threaded(
      ShellBlockData<>(dfbs_),
      ShellBlockData<>(dfbs_),
      eris_2c_, coulomb_oper_type_
  );
  memory_used_ += sizeof(TwoCenterIntContainer) + dfnbf*dfnbf*sizeof(double);
  const auto& g2 = *g2_full_ptr_;

  schwarz_df_.resize(dfbs_->nshell());
  for(auto&& Xsh : shell_range(dfbs_)) {
    schwarz_df_(Xsh) = g2.block(Xsh.bfoff, Xsh.bfoff, Xsh.nbf, Xsh.nbf).norm();
  }

  timer.exit();


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
  // TODO MessageGrp takes an int here, and if ncoefs is large, it needs to take a long
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
    uli n_Xcoefs = 0;
    const cadf::Node& my_part = atom_pair_assignments_k_->my_assignments(scf_grp_->me());
    uli ncoefs_dist = my_part.bin->obs_ncoefs + my_part.bin->dfbs_ncoefs;
    dist_coefs_data_ = new double[ncoefs_dist];
    memset(dist_coefs_data_, 0, ncoefs_dist*sizeof(double));

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
    }
    for(auto&& dfbs_atom : my_part.bin->assigned_dfbs_atoms) {
      ShellBlockData<> Xblk = ShellBlockData<>::atom_block(dfbs_atom->index, dfbs_, gbs_);
      coefs_X_nu.emplace(
          std::piecewise_construct,
          std::forward_as_tuple(Xblk.center),
          std::forward_as_tuple(
              dist_coefs_data_ + my_part.bin->dfbs_coef_offsets.at(Xblk.center),
              Xblk.atom_nbf, Xblk.atom_obsnbf*nbf
          )
      );
      coefs_X_nu_other.emplace(
          std::piecewise_construct,
          std::forward_as_tuple(Xblk.center),
          std::forward_as_tuple(
              dist_coefs_data_ + my_part.bin->dfbs_coef_offsets.at(Xblk.center),
              Xblk.atom_nbf*Xblk.atom_obsnbf, nbf
          )
      );
    }

    // TODO threads
    std::vector<Eigen::Map<Eigen::VectorXd>> empty;
    //assert(scf_grp_->n() > 1 || my_part.bin->compute_coef_items[false].size() > 0);
    //assert(scf_grp_->n() > 1 || my_part.compute_coef_items[false].size() > 0);
    for(auto&& obs_shell : my_part.compute_coef_items[false]) {

      //DUMP(obs_shell->index);
      // Do the C_mu_X part first
      ShellData ish(obs_shell->index, gbs_, dfbs_);
      for(auto&& jsh : iter_significant_partners(ish)) {
        std::vector<Eigen::Map<Eigen::VectorXd>> coefs;
        for(auto&& mu : function_range(ish)) {
          for(auto&& rho : function_range(jsh)) {
            coefs.emplace_back(
                coefs_mu_X.at(ish).data()
                  + mu.off*nbf*ish.atom_dfnbf + rho*ish.atom_dfnbf,
                ish.atom_dfnbf
            );
          }
        }
        //DUMP2(ish.center, jsh.center)
        //if(ish.center >= jsh.center)
          get_coefs_ish_jsh(ish, jsh, 0, coefs, empty);
        //else
        //  get_coefs_ish_jsh(jsh, ish, 0, empty, coefs);

        //DEBUG_DELETE_THIS
        //bool failed = false;
        //for(auto&& mu : function_range(ish)) {
        //  for(auto&& rho : function_range(jsh)) {
        //    for(auto&& X : iter_functions_on_center(dfbs_, ish.center)) {
        //      if(fabs(coefs_mu_X.at(ish)(mu.off, rho*ish.atom_dfnbf + X.bfoff_in_atom)
        //          - coefs_transpose_[X](mu.bfoff_in_atom, rho)) > 1e-14
        //      ) {
        //        ExEnv::out0() << "C_{" << mu.index << ", " << rho.index << "}^{" << X.index << "} is not correct. ("
        //            << coefs_mu_X.at(ish)(mu.off, rho*ish.atom_dfnbf + X.bfoff_in_atom)
        //            << " != "
        //            << coefs_transpose_[X](mu.bfoff_in_atom, rho) << ")" << std::endl;
        //        failed = true;
        //      }
        //    }
        //  }
        //}
        //assert(!failed);
        //DEBUG_DELETE_THIS
      }
    }

    // Node-row-wise sum of mu coefficients
    {
      Ref<MessageGrp> mu_grp = scf_grp_->split(my_part.bin->obs_row_id);
      mu_grp->sum(dist_coefs_data_, my_part.bin->obs_ncoefs);
    } // mu_grp is deleted

    sc::SCFormIO::init_mp(scf_grp_->me());
    //assert(scf_grp_->n() > 1 || my_part.bin->assigned_obs_shells.size() == gbs_->nshell());
    //assert(scf_grp_->n() > 1 || my_part.bin->assigned_dfbs_atoms.size() == gbs_->ncenter());
    //assert(scf_grp_->n() > 1 || my_part.pairs.size() == gbs_->ncenter() * gbs_->nshell());

    std::vector<CoefView> empty_df;
    //assert(scf_grp_->n() > 1 || my_part.bin->compute_coef_items[true].size() > 0);
    //assert(scf_grp_->n() > 1 || my_part.compute_coef_items[true].size() == gbs_->ncenter());
    for(auto&& dfbs_atom : my_part.compute_coef_items[true]) {

      // Now do the C_X_mu part
      ShellBlockData<> Xblk = ShellBlockData<>::atom_block(dfbs_atom->index, gbs_, dfbs_);
      for(auto&& ish : iter_shells_on_center(gbs_, Xblk.center, dfbs_)) {
        for(auto&& jsh : iter_significant_partners(ish)) {
          std::vector<CoefView> coefs;
          for(auto&& mu : function_range(ish)) {
            for(auto&& rho : function_range(jsh)) {
              coefs.emplace_back(
                  coefs_X_nu.at(Xblk.center).data() + mu.bfoff_in_atom*nbf + rho,
                  ish.atom_dfnbf,
                  Eigen::Stride<1, Eigen::Dynamic>(1, mu.atom_nbf * nbf)
              );
            }
          }
          get_coefs_ish_jsh(ish, jsh, 0, coefs, empty_df);
          //DEBUG_DELETE_THIS
          //bool failed = false;
          //for(auto&& mu : function_range(ish)) {
          //  for(auto&& rho : function_range(jsh)) {
          //    for(auto&& X : iter_functions_on_center(dfbs_, ish.center)) {
          //      if(fabs(coefs_X_nu.at(Xblk.center)(X.bfoff_in_atom, mu.bfoff_in_atom*nbf + rho)
          //          - coefs_transpose_[X](mu.bfoff_in_atom, rho)) > 1e-14
          //      ) {
          //        ExEnv::out0() << "C_{" << mu.index << ", " << rho.index << "}^{" << X.index << "} is not correct. ("
          //            << coefs_X_nu.at(Xblk.center)(X.bfoff_in_atom, mu.bfoff_in_atom*nbf + rho)
          //            << " != "
          //            << coefs_transpose_[X](mu.bfoff_in_atom, rho) << ")" << std::endl;
          //        failed = true;
          //      }
          //    }
          //  }
          //}
          //assert(!failed);
          //DEBUG_DELETE_THIS
        }
      }
    }

    // Node-row-wise sum of X coefficients
    {
      Ref<MessageGrp> X_grp = scf_grp_->split(my_part.bin->dfbs_row_id);
      X_grp->sum(dist_coefs_data_ + my_part.bin->obs_ncoefs, my_part.bin->dfbs_ncoefs);
    } // X_grp is deleted

    sc::SCFormIO::init_mp(scf_grp_->me());
  }

  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Store the transpose and blocked coefficients          		                        {{{1 */ #if 1 //latex `\label{sc:coeftrans}`
  //---------------------------------------------------------------------------------------//

  // TODO Get rid of this!!!
  timer.change("05 - store blocked coefs");
  if(store_coefs_transpose_) {
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
    //begin_xml_context(
    //    "df_coefficients",
    //    "compute_.xml"
    //);
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

  //schwarz_df_.resize(dfbs_->nshell());
  //for(auto Xsh : shell_range(dfbs_)){
  //  auto g_XX_ptr = ints_to_eigen(
  //      Xsh, Xsh, eris_2c_[0], coulomb_oper_type_
  //  );
  //  auto& g_XX = *g_XX_ptr;
  //  double frob_norm = 0.0;
  //  for(auto X : function_range(Xsh)) {
  //    frob_norm += g_XX(X.bfoff_in_shell, X.bfoff_in_shell);
  //  }
  //  schwarz_df_[Xsh] = sqrt(frob_norm);
  //}

  if(do_linK_) {

    /*
    do_threaded(nthread_, [&](int ithr) {
      for(auto&& lsh : thread_over_range(shell_range(obs), ithr, nthread_)) {
        for(auto&& ksh : L_schwarz[lsh]) {
          double Ct = 0.0;
          //----------------------------------------//
          for(auto&& sigma : function_range(lsh)) {
            for(auto&& nu : function_range(ksh)) {
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

              for(auto&& Xsh : iter_shells_on_center(dfbs_, first.center)){
                auto& C = *cpair.first;
                for(auto&& X : function_range(Xsh)) {
                  const double CX = C[X.bfoff_in_atom];
                  Ct += CX * CX * g2(X, X);
                }
              }

              if(first.center != second.center) {
                for(auto&& Xsh : iter_shells_on_center(dfbs_, second.center)){
                  auto& C = *cpair.second;
                  for(auto&& X : function_range(Xsh)) {
                    const double CX = C[X.bfoff_in_atom];
                    Ct += CX * CX * g2(X, X);
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
    });

    do_threaded(nthread_, [&](int ithr){
      auto L_coefs_iter = L_coefs.begin();
      const auto& L_coefs_end = L_coefs.end();
      L_coefs_iter.advance(ithr);
      while(L_coefs_iter != L_coefs_end) {
        L_coefs_iter->second.sort();
        L_coefs_iter.advance(nthread_);
      }
    });
    */

    //----------------------------------------//
    // Compute the Frobenius norm of C_transpose_ blocks
    timer.change("07 - C_trans_frob");
    if(distribute_coefficients_) {
      C_trans_frob_.resize(dfbs_->nshell());
      const cadf::Node& my_part = atom_pair_assignments_k_->my_assignments(scf_grp_->me());
      for(auto Xatom : my_part.bin->assigned_dfbs_atoms) {
        //ShellBlockData<> Xblk = ShellBlockData<>::atom_block(Xatom->index, dfbs_, gbs_);
        for(auto&& Xsh : iter_shells_on_center(dfbs_, Xatom->index, gbs_)) {

          resize_and_zero_matrix(C_trans_frob_[Xsh], Xsh.atom_obsnsh, obs->nshell());

          for(auto&& ish : iter_shells_on_center(obs, Xsh.center)) {
            for(auto&& jsh : shell_range(obs)) {
              C_trans_frob_[Xsh](ish.shoff_in_atom, jsh) += coefs_X_nu.at(Xsh.center).block(
                  Xsh.bfoff_in_atom, ish.bfoff_in_atom*nbf + jsh.bfoff,
                  Xsh.nbf, jsh.nbf
              ).squaredNorm();
            }
          }

          C_trans_frob_[Xsh] = C_trans_frob_[Xsh].array().sqrt();

        }
      }
    }
    else {
      C_trans_frob_.resize(dfbs_->nshell());
      for(auto Xsh : shell_range(dfbs_, obs)) {

        resize_and_zero_matrix(C_trans_frob_[Xsh], Xsh.atom_obsnsh, obs->nshell());

        for(auto&& ish : iter_shells_on_center(obs, Xsh.center)) {
          for(auto&& jsh : shell_range(obs)) {
            for(auto&& mu : function_range(ish)) {
              C_trans_frob_[Xsh](ish.shoff_in_atom, jsh) += coefs_transpose_blocked_[Xsh.center].block(
                  Xsh.bfoff_in_atom, mu.bfoff_in_atom*nbf + jsh.bfoff,
                  Xsh.nbf, jsh.nbf
              ).squaredNorm();
            }
          } // end loop over jsh
        } // end loop over ish

        C_trans_frob_[Xsh] = C_trans_frob_[Xsh].array().sqrt();

      } // end loop over Xsh
    }
    //----------------------------------------//

    // Compute Cmaxes_
    timer.change("08 - Cbar");
    resize_and_zero_matrix(C_bar_, gbs_->nshell(), dfbs_->nshell());
    if(use_norms_nu_) {
      if(distribute_coefficients_) {
        const cadf::Node& my_part = atom_pair_assignments_k_->my_assignments(scf_grp_->me());
        for(auto Xatom : my_part.bin->assigned_dfbs_atoms) {
          for(auto&& Xsh : iter_shells_on_center(dfbs_, Xatom->index, gbs_)) {
            C_bar_.col(Xsh) = C_trans_frob_[Xsh].colwise().squaredNorm();
            C_bar_.col(Xsh).segment(Xsh.atom_obsshoff, Xsh.atom_obsnsh) +=
                C_trans_frob_[Xsh].rowwise().squaredNorm();
            C_bar_.col(Xsh).segment(Xsh.atom_obsshoff, Xsh.atom_obsnsh) -=
                C_trans_frob_[Xsh].middleCols(Xsh.atom_obsshoff, Xsh.atom_obsnsh).rowwise().squaredNorm();
          }
          C_bar_ = C_bar_.array().sqrt();
        }

        // Node-row-wise sum of X parts of C_bar_
        {
          Ref<MessageGrp> X_grp = scf_grp_->split(my_part.bin->dfbs_row_id);
          X_grp->sum(C_bar_.data(), gbs_->nshell() * dfbs_->nshell());
        } // X_grp is deleted
        sc::SCFormIO::init_mp(scf_grp_->me());

      }
      else {
        for(auto&& Xsh : shell_range(dfbs_, gbs_)) {
          C_bar_.col(Xsh) = C_trans_frob_[Xsh].colwise().squaredNorm();
          C_bar_.col(Xsh).segment(Xsh.atom_obsshoff, Xsh.atom_obsnsh) +=
              C_trans_frob_[Xsh].rowwise().squaredNorm();
          C_bar_.col(Xsh).segment(Xsh.atom_obsshoff, Xsh.atom_obsnsh) -=
              C_trans_frob_[Xsh].middleCols(Xsh.atom_obsshoff, Xsh.atom_obsnsh).rowwise().squaredNorm();
        }
        C_bar_ = C_bar_.array().sqrt();
      }
    }
    else {
      if(distribute_coefficients_) {
        throw FeatureNotImplemented("use_norms_nu = false with distributed coefficients", __FILE__, __LINE__, class_desc());
      }

      // TODO just do this in the above loop similarly
      for(auto Xsh : shell_range(dfbs_)){
        for(auto lsh : shell_range(obs)) {

          double max_val;

          if(lsh.center == Xsh.center) {
            // We might be able to get away with only nu on same center as X and sigma,
            //  but for now approach it more rigorously
            max_val = C_trans_frob_[Xsh].row(lsh.shoff_in_atom).maxCoeff();
          }
          else { // different centers...
            max_val = C_trans_frob_[Xsh].col(lsh).maxCoeff();
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
  have_coefficients_ = true;                                                               //latex `\label{sc:coefflag}`
  timer.exit();
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
}

//////////////////////////////////////////////////////////////////////////////////
