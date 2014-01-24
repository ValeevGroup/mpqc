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
#include <boost/functional/hash.hpp>
#include <chemistry/qc/scf/clhf.h>
#include <util/container/conc_cache.h>
#include <atomic>
#include <future>
#include <functional>
#include <type_traits>
#include <Eigen/Dense>
#include <util/misc/property.h>
#include <chemistry/qc/scf/cadf_iters.h>
#include <boost/static_assert.hpp>

typedef typename boost::integer_range<int> int_range;


namespace sc {

template <typename... Args>
void
do_threaded(int nthread, Args... args, const std::function<void(int, Args...)>& f){
  boost::thread_group compute_threads;
  // Loop over number of threads
  for(int ithr = 0; ithr < nthread; ++ithr) {
    // create each thread that runs f
    compute_threads.create_thread([&, ithr](){
      // run the work
      f(ithr, args...);
    });
  }
  // join the created threads
  compute_threads.join_all();
}

//============================================================================//
//============================================================================//
//============================================================================//

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
    typedef ConcurrentCache<
        std::shared_ptr<Eigen::MatrixXd>,
        int, int, TwoBodyOper::type
    > TwoCenterIntCache;
    typedef ConcurrentCache<
        std::shared_ptr<Eigen::MatrixXd>,
        int, int, int, TwoBodyOper::type
    > ThreeCenterIntCache;
    typedef ConcurrentCache<
        double,
        int, int, int, int, TwoBodyOper::type
    > FourCenterMaxIntCache;
    typedef ConcurrentCache<std::shared_ptr<Decomposition>, int, int> DecompositionCache;

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

    typedef enum {
      AllPairs = 0,
      SignificantPairs = 1,
      ExchangeOuterLoopPairs = 2
    } PairSet;

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

    bool get_shell_pair(ShellData& mu, ShellData& nu, PairSet pset = AllPairs);

    void loop_shell_pairs_threaded(PairSet pset,
        const std::function<void(int, const ShellData&, const ShellData&)>& f
    );

    void loop_shell_pairs_threaded(
        const std::function<void(int, const ShellData&, const ShellData&)>& f
    );


    RefSCMatrix compute_J();

    RefSCMatrix compute_K();

    void compute_coefficients();

    void initialize();

    void init_significant_pairs();

    void init_threads();

    /// returns shell ints in inbf x jnbf Eigen Matrix pointer
    std::shared_ptr<Eigen::MatrixXd> ints_to_eigen(
        int ish, int jsh,
        Ref<TwoBodyTwoCenterInt>& ints,
        TwoBodyOper::type ints_type
    );

    template <typename ShellRange>
    std::shared_ptr<Eigen::MatrixXd> ints_to_eigen(
        const ShellBlockData<ShellRange>& ish,
        const ShellBlockData<ShellRange>& jsh,
        Ref<TwoBodyTwoCenterInt>& ints,
        TwoBodyOper::type ints_type
    );

    template <typename ShellRange>
    std::shared_ptr<Eigen::MatrixXd> ints_to_eigen(
        const ShellBlockData<ShellRange>& ish, const ShellData& jsh,
        Ref<TwoBodyTwoCenterInt>& ints,
        TwoBodyOper::type ints_type
    );

    /// returns ints for shell in (inbf, jnbf) x kdfnbf matrix in chemists' notation
    std::shared_ptr<Eigen::MatrixXd> ints_to_eigen(
        int ish, int jsh, int ksh,
        Ref<TwoBodyThreeCenterInt>& ints,
        TwoBodyOper::type ints_type
    );

    /// returns ints for shell in (inbf, jnbf) x kblk.nbf matrix in chemists' notation
    template <typename ShellRange>
    std::shared_ptr<Eigen::MatrixXd> ints_to_eigen(
        const ShellData& ish, const ShellData& jsh,
        const ShellBlockData<ShellRange>& kblk,
        Ref<TwoBodyThreeCenterInt>& ints,
        TwoBodyOper::type ints_type
    );

    /// returns ints for shell in (inbf, jblk.nbf) x kblk.nbf matrix in chemists' notation
    template <typename ShellRange1, typename ShellRange2>
    std::shared_ptr<Eigen::MatrixXd> ints_to_eigen(
        const ShellBlockData<ShellRange1>& ish,
        const ShellData& jsh,
        const ShellBlockData<ShellRange2>& kblk,
        Ref<TwoBodyThreeCenterInt>& ints,
        TwoBodyOper::type ints_type
    );

    std::shared_ptr<Decomposition> get_decomposition(int ish, int jsh, Ref<TwoBodyTwoCenterInt> ints);

    // Non-blocked version of cl_gmat_
    RefSymmSCMatrix gmat_;

    /// The threshold for including pairs in the pair list using half of the schwarz bound
    double pair_screening_thresh_;
    double density_screening_thresh_;
    double full_screening_thresh_;
    double coef_screening_thresh_;

    bool print_screening_stats_;

    // Whether or not the MessageGrp is an instance of MPIMessageGrp
    bool using_mpi_;

    // Whether or not we've computed the fitting coefficients yet (only needs to be done on first iteration)
    bool have_coefficients_;

    // The density fitting auxiliary basis
    Ref<GaussianBasisSet> dfbs_ = 0;

    // Where are the shell pairs being evaluated?
    std::map<PairSet, std::map<std::pair<int, int>, int>> pair_assignments_;

    // Pair assignments for K
    std::map<
      PairSet,
      std::map<
        std::pair<int, ShellBlockSkeleton<>>,
        int
      >
    > pair_assignments_k_;

    // What pairs are being evaluated on the current node?
    std::map<PairSet, std::vector<std::pair<int, int>>> local_pairs_;

    // What pairs are being evaluated on the current node?
    std::map<
      PairSet,
      std::vector<
        std::pair<int, ShellBlockSkeleton<>>
      >
    > local_pairs_k_;

    // List of the permutationally unique pairs with half-schwarz bounds larger than pair_thresh_
    std::vector<std::pair<int, int>> sig_pairs_;

    std::vector<std::set<int>> sig_partners_;
    std::vector<std::set<ShellBlockSkeleton<>>> sig_blocks_;

    // The same as sig_pairs_, but organized differently
    std::vector<std::vector<int>> shell_to_sig_shells_;

    std::vector<double> max_schwarz_;
    std::vector<double> schwarz_df_;

    // Where are we in the iteration over the local_pairs_?
    std::atomic<int> local_pairs_spot_;

    Eigen::MatrixXd schwarz_frob_;
    std::vector<Eigen::MatrixXd> C_trans_frob_;

    // TwoBodyThreeCenterInt integral objects for each thread
    std::vector<Ref<TwoBodyThreeCenterInt>> eris_3c_;
    std::vector<Ref<TwoBodyThreeCenterInt>> metric_ints_3c_;

    // TwoBodyTwoCenterInt integral objects for each thread
    std::vector<Ref<TwoBodyTwoCenterInt>> eris_2c_;
    std::vector<Ref<TwoBodyTwoCenterInt>> metric_ints_2c_;

    // Integrals computed locally this iteration
    std::atomic<int> ints_computed_locally_;
    int ints_computed_;

    // Coefficients storage.  Not accessed directly
    double* coefficients_data_ = 0;

    CoefMap coefs_;

    std::vector<Eigen::MatrixXd> coefs_transpose_;

    DecompositionCache decomps_;
    TwoCenterIntCache ints2_;
    ThreeCenterIntCache ints3_;
    FourCenterMaxIntCache ints4maxes_;

    bool threads_initialized_ = false;

    static ClassDesc cd_;

    typedef typename decltype(sig_partners_)::value_type::iterator sig_partners_iter_t;
    typedef range_of<ShellData, sig_partners_iter_t> sig_partners_range_t;

    sig_partners_range_t
    iter_significant_partners(const ShellData& ish){
      const auto& sig_parts = sig_partners_[ish];
      return boost::make_iterator_range(
          basis_element_iterator<ShellData, sig_partners_iter_t>(
              ish.basis, ish.dfbasis, sig_parts.begin()
          ),
          basis_element_iterator<ShellData, sig_partners_iter_t>(
              ish.basis, ish.dfbasis, sig_parts.end()
          )
      );
    }

    bool do_linK_ = false;

    // CADF-LinK lists
    template <typename T> struct hash_;
    template<typename A, typename B>
    struct hash_<std::pair<A, B>>{
      std::size_t operator()(const std::pair<A, B>& val) const {
        std::size_t seed = 0;
        boost::hash_combine(seed, val.first);
        boost::hash_combine(seed, val.second);
        return seed;
      }
    };
    typedef madness::ConcurrentHashMap<int, OrderedShellList> IndexListMap;
    typedef madness::ConcurrentHashMap<
        std::pair<int, int>,
        OrderedShellList,
        hash_<std::pair<int, int>>
    > IndexListMap2;
    //typedef ConcurrentMap<OrderedShellList, int, int, int> LinKListCache3;

    IndexListMap L_schwarz;
    IndexListMap L_coefs;
    IndexListMap L_D;
    IndexListMap2 L_3;
    /*
    LinKListCache2 L_4;
    std::set<std::pair<int, int>> L_4_keys;
    LinKListCache2 L_K;
    std::set<std::pair<int, int>> L_K_keys;
    LinKListCache2 L_5;
    */

};

template <typename ShellRange>
std::shared_ptr<Eigen::MatrixXd>
CADFCLHF::ints_to_eigen(
    const ShellBlockData<ShellRange>& iblk,
    const ShellBlockData<ShellRange>& jblk,
    Ref<TwoBodyTwoCenterInt>& ints,
    TwoBodyOper::type int_type
){
  auto rv = std::make_shared<Eigen::MatrixXd>(iblk.nbf, jblk.nbf);
  for(auto ish : shell_range(iblk)) {
    for(auto jsh : shell_range(jblk)) {
      const auto& ints_ptr = ints_to_eigen(ish, jsh, ints, int_type);
      rv->block(ish.bfoff - iblk.bfoff, jsh.bfoff - jblk.bfoff, ish.nbf, jsh.nbf) = *ints_ptr;
    }
  }
  return rv;
}

template <typename ShellRange>
std::shared_ptr<Eigen::MatrixXd>
CADFCLHF::ints_to_eigen(
    const ShellBlockData<ShellRange>& iblk,
    const ShellData& jsh,
    Ref<TwoBodyTwoCenterInt>& ints,
    TwoBodyOper::type int_type
){
  auto rv = std::make_shared<Eigen::MatrixXd>((const int)iblk.nbf, (const int)jsh.nbf);
  for(auto ish : shell_range(iblk)) {
    const auto& ints_ptr = ints_to_eigen(ish, jsh, ints, int_type);
    rv->block(ish.bfoff - iblk.bfoff, 0, ish.nbf, jsh.nbf) = *ints_ptr;
  }
  return rv;
}

template <typename ShellRange>
std::shared_ptr<Eigen::MatrixXd>
CADFCLHF::ints_to_eigen(
    const ShellData& ish, const ShellData& jsh,
    const ShellBlockData<ShellRange>& Xblk,
    Ref<TwoBodyThreeCenterInt>& ints,
    TwoBodyOper::type int_type
){
  auto rv = std::make_shared<Eigen::MatrixXd>(
      ish.nbf * jsh.nbf,
      Xblk.nbf
  );
  for(auto Xsh : shell_range(Xblk)) {
    out_assert(ints->basis3()->nbasis(), ==, Xsh.basis->nbasis());
    const auto& ints_ptr = ints_to_eigen(ish, jsh, Xsh, ints, int_type);
    rv->middleCols(Xsh.bfoff - Xblk.bfoff, Xsh.nbf) = *ints_ptr;
  }
  return rv;
}

template <typename ShellRange1, typename ShellRange2>
std::shared_ptr<Eigen::MatrixXd>
CADFCLHF::ints_to_eigen(
    const ShellBlockData<ShellRange1>& iblk,
    const ShellData& jsh,
    const ShellBlockData<ShellRange2>& Xblk,
    Ref<TwoBodyThreeCenterInt>& ints,
    TwoBodyOper::type int_type
){
  auto rv = std::make_shared<Eigen::MatrixXd>(
      iblk.nbf * jsh.nbf,
      Xblk.nbf
  );
  for(auto ish: shell_range(iblk)) {
    for(auto Xsh : shell_range(Xblk)) {
      const auto& ints_ptr = ints_to_eigen(ish, jsh, Xsh, ints, int_type);
      rv->block(
          (ish.bfoff - iblk.bfoff) * jsh.nbf, Xsh.bfoff - Xblk.bfoff,
          ish.nbf*jsh.nbf, Xsh.nbf
      ) = *ints_ptr;
    }
  }
  return rv;
}

} // end namespace sc

#endif /* _chemistry_qc_scf_cadfclhf_h */
