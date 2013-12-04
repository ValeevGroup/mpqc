//
// cadfclhf.h
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


#ifndef _chemistry_qc_scf_cadfclhf_h
#define _chemistry_qc_scf_cadfclhf_h

#include <boost/thread/thread.hpp>
#include <boost/range/irange.hpp>
#include <chemistry/qc/scf/clhf.h>
#include <util/container/conc_cache.h>
#include <atomic>
#include <future>
#include <Eigen/Dense>

typedef typename boost::integer_range<int> int_range;


namespace sc {

// Forward declarations
struct ShellData;

/**
 * A specialization of CLHF that uses concentric atomic
 *   density fitting to build fock matrices
 */
class CADFCLHF: public CLHF {

  protected:
    void ao_fock(double accuracy);
    void reset_density();

  public:

    // This will later be changed to simple double* to allow for direct BLAS calls
    typedef std::shared_ptr<Eigen::Map<Eigen::VectorXd>> CoefContainer;

    typedef std::map<std::pair<int, int>, std::pair<CoefContainer, CoefContainer>> CoefMap;
    typedef Eigen::HouseholderQR<Eigen::MatrixXd> Decomposition;
    // Shared mutexes take up too much memory to have one per shell
    //   In really large cases, we may have to switch to row-locked caches,
    //   which only create k*N mutexes but are much less efficient.
    typedef ElementLockedExclusiveReadCache<
        std::shared_ptr<Eigen::MatrixXd>,
        int, int, TwoBodyOper::type
    > TwoCenterIntCache;
    typedef ElementLockedSharedReadCache<
        std::shared_ptr<Eigen::MatrixXd>,
        int, int, int, TwoBodyOper::type
    > ThreeCenterIntCache;
    // In the case of the Natom x Natom decomposition cache, it's probably
    //   safe to use shared mutexes.
    typedef ElementLockedExclusiveReadCache<std::shared_ptr<Decomposition>, int, int> DecompositionCache;

    CADFCLHF(StateIn&);

    /** Accepts all keywords of CLHF class + the following keywords:
        <table border="1">

          <tr><td>%Keyword<td>Type<td>Default<td>Description

          </table>

     */
    CADFCLHF(const Ref<KeyVal>&);
    ~CADFCLHF();

    void save_data_state(StateOut&);

  private:


    enum {
      NoMorePairs = -1
    };

    bool is_master() {
      if(dynamic_){
        // could be improved to have multiple masters
        return scf_grp_->me() == 0;
      }
      else{
        return false;
      }
    }

    RefSCMatrix D_;

    // For now, just do static load balancing
    bool dynamic_ = false;

    TwoBodyOper::type metric_oper_type_;

    // Convenience variable for better code readibility
    static const TwoBodyOper::type coulomb_oper_type_ = TwoBodyOper::eri;

    bool get_shell_pair(ShellData& mu, ShellData& nu);

    RefSCMatrix compute_J();

    RefSCMatrix compute_K();

    void compute_coefficients();

    void initialize();

    void init_threads();

    /// returns shell ints in inbf x jnbf Eigen Matrix pointer
    std::shared_ptr<Eigen::MatrixXd> ints_to_eigen(
        int ish, int jsh,
        Ref<TwoBodyTwoCenterInt> ints,
        TwoBodyOper::type ints_type
    );

    /// Range version of the above
    std::shared_ptr<Eigen::MatrixXd> ints_to_eigen(
        int_range&& ishs,
        int_range&& jshs,
        Ref<TwoBodyTwoCenterInt> ints,
        TwoBodyOper::type ints_type
    );

    /// returns ints for shell in (inbf, jnbf) x kdfnbf matrix in chemists' notation
    std::shared_ptr<Eigen::MatrixXd> ints_to_eigen(
        int ish, int jsh, int ksh,
        Ref<TwoBodyThreeCenterInt> ints,
        TwoBodyOper::type ints_type
    );

    std::shared_ptr<Decomposition> get_decomposition(int ish, int jsh, Ref<TwoBodyTwoCenterInt> ints);

    // Non-blocked version of cl_gmat_
    RefSymmSCMatrix gmat_;

    // Whether or not the MessageGrp is an instance of MPIMessageGrp
    bool using_mpi_;

    // Whether or not we've computed the fitting coefficients yet (only needs to be done on first iteration)
    bool have_coefficients_;

    // The density fitting auxiliary basis
    Ref<GaussianBasisSet> dfbs_ = 0;

    // Where are the shell pairs being evaluated?
    std::map<std::pair<int, int>, int> pair_assignments_;

    // What shell pairs are being evaluated on the current node?
    std::vector<std::pair<int, int> > local_pairs_;

    // Where are we in the iteration over the local_pairs_?
    std::atomic<int> local_pairs_spot_;

    // Is the pair fetch process in progress?  Only used when dynamic_ is true
    std::atomic<bool> is_fetching_pairs_;

    // TwoBodyThreeCenterInt integral objects for each thread
    std::vector<Ref<TwoBodyThreeCenterInt>> eris_3c_;
    std::vector<Ref<TwoBodyThreeCenterInt>> metric_ints_3c_;

    // TwoBodyTwoCenterInt integral objects for each thread
    std::vector<Ref<TwoBodyTwoCenterInt>> eris_2c_;
    std::vector<Ref<TwoBodyTwoCenterInt>> metric_ints_2c_;

    // Coefficients storage.  Not accessed directly
    double* coefficients_data_ = 0;

    CoefMap coefs_;

    DecompositionCache decomps_;
    TwoCenterIntCache ints2_;
    ThreeCenterIntCache ints3_;

    static ClassDesc cd_;

};

} // end namespace sc

#endif /* _chemistry_qc_scf_cadfclhf_h */
