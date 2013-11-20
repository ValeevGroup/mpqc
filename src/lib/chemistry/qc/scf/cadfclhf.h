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

#include <chemistry/qc/scf/clhf.h>
#include <util/container/conc_cache.h>
#include <atomic>
#include <future>
#include <Eigen/Dense>


namespace sc {

// Forward declaration
class CADFCLHF;

namespace detail {

  struct compute_task {
    public:
      compute_task(Ref<CADFCLHF> wfn, int ithr) : wfn_(wfn), ithr_(ithr) { }
    protected:
      Ref<CADFCLHF> wfn_;
      int ithr_;
  };

  struct compute_J_task : compute_task {
    void operator()();
  };

  struct compute_K_task : compute_task {
    void operator()();
  };

  struct compute_coefs_task : compute_task {
    void operator()();
  };
}

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
    typedef Eigen::Map<Eigen::VectorXd> CoefContainer;

    typedef std::map<std::pair<int, int>, std::pair<CoefContainer, CoefContainer>> CoefMap;
    typedef std::shared_ptr<Eigen::MatrixXd> IntContainer2;
    typedef std::vector<std::vector<Eigen::VectorXd>> IntContainer3;
    typedef Eigen::HouseholderQR<Eigen::MatrixXd> Decomposition;
    typedef ConcurrentCache<std::shared_ptr<Decomposition>, int, int> DecompositionCache;
    typedef ConcurrentCache<std::shared_ptr<Eigen::MatrixXd>, int, int> TwoCenterIntCache;

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

    // For now, just do static load balancing
    bool dynamic_ = false;

    bool get_shell_pair(int& mu, int& nu);

    RefSCMatrix compute_J();

    RefSCMatrix compute_K();

    void compute_coefficients();

    void initialize();

    void init_threads();

    /// returns shell ints in inbf x jnbf Eigen Matrix pointer
    std::shared_ptr<Eigen::MatrixXd> ints_to_eigen(int ish, int jsh, Ref<TwoBodyTwoCenterInt> ints);

    /// returns ints for shell in (inbf, jnbf) x kdfnbf matrix
    std::shared_ptr<Eigen::MatrixXd> ints_to_eigen(int ish, int jsh, int ksh, Ref<TwoBodyThreeCenterInt> ints);

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

    static ClassDesc cd_;

    friend struct detail::compute_J_task;
    friend struct detail::compute_K_task;
    friend struct detail::compute_coefs_task;

};

} // end namespace sc

#endif /* _chemistry_qc_scf_cadfclhf_h */
