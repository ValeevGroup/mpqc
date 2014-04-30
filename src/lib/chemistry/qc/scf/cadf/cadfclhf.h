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

#define USE_SPARSE 0

// Standard library includes
#include <atomic>
#include <future>
#include <functional>
#include <type_traits>
#include <unordered_set>
#include <unordered_map>

// Eigen includes
#include <Eigen/Dense>
#include <Eigen/Sparse>

// MPQC includes
#include <chemistry/qc/scf/clhf.h>
#include <util/container/conc_cache_fwd.h>
#include <util/container/thread_wrap.h>
#include <util/misc/hash.h>

#include "iters.h"

#ifdef ECLIPSE_PARSER_ONLY
#include <util/misc/sharedptr.h>
#endif

// Allows easy switching between std::shared_ptr and e.g. boost::shared_ptr
template<typename... Types>
using shared_ptr = std::shared_ptr<Types...>;
using std::make_shared;

typedef typename boost::integer_range<int> int_range;
typedef unsigned long long ull;

#define USE_INTEGRAL_CACHE 0

namespace sc {


// Forward Declarations

class XMLWriter;
namespace cadf { class AssignmentGrid; }

//============================================================================//

inline void
do_threaded(int nthread, const std::function<void(int)>& f){
  boost::thread_group compute_threads;
  // Loop over number of threads
  for(int ithr = 0; ithr < nthread; ++ithr) {
    // create each thread that runs f
    compute_threads.create_thread([&, ithr](){
      // run the work
      f(ithr);
    });
  }
  // join the created threads
  compute_threads.join_all();
}

//============================================================================//
//============================================================================//
//============================================================================//

/// @addtogroup ChemistryElectronicStructureOneBodyHFCADF
/// @{

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
    typedef Eigen::Map<Eigen::VectorXd, Eigen::Unaligned, Eigen::Stride<1, Eigen::Dynamic>> CoefView;
    typedef shared_ptr<CoefView> CoefContainer;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> TwoCenterIntContainer;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ThreeCenterIntContainer;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> FourCenterIntContainer;
    typedef shared_ptr<TwoCenterIntContainer> TwoCenterIntContainerPtr;
    typedef shared_ptr<ThreeCenterIntContainer> ThreeCenterIntContainerPtr;
    typedef shared_ptr<FourCenterIntContainer> FourCenterIntContainerPtr;

    typedef std::unordered_map<
        std::pair<int, int>,
        std::pair<CoefContainer, CoefContainer>,
        sc::hash<std::pair<int, int>>
    > CoefMap;
    typedef Eigen::HouseholderQR<Eigen::MatrixXd> Decomposition;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMatrix;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> ColMatrix;
    typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SparseRowMatrix;
    typedef Eigen::Map<RowMatrix, Eigen::Unaligned, Eigen::OuterStride<>> StridedRowMap;

#if USE_INTEGRAL_CACHE

    typedef ConcurrentCacheWithSymmetry<
        TwoCenterIntContainerPtr,
        // 3 total keys, allow transpose of first two (0, 1)
        SingleTranspositionKeySymmetry<3, 0, 1>,
        int, int, TwoBodyOper::type
    > TwoCenterIntCache;

    typedef ConcurrentCacheWithSymmetry<
        ThreeCenterIntContainerPtr,
        // 4 total keys, allow transpose of first two
        SingleTranspositionKeySymmetry<4, 0, 1>,
        int, int, int, TwoBodyOper::type
    > ThreeCenterIntCache;
#endif

    // TODO Get rid of this
    typedef ConcurrentCache<
        double,
        int, int, int, int, TwoBodyOper::type
    > FourCenterMaxIntCache;
    typedef ConcurrentCache<shared_ptr<Decomposition>, int, int> DecompositionCache;

    CADFCLHF(StateIn&);

    /** Accepts all keywords of CLHF class + the following keywords:
        <table border="1">

          <tr><td>%Keyword<td>Type<td>Default<td>Description

          </table>

     */
    CADFCLHF(const Ref<KeyVal>&);

    ~CADFCLHF();

    void save_data_state(StateOut&);

    class ScreeningStatistics {

      public:

        typedef std::atomic_uint_fast64_t accumulate_t;
        typedef decltype(accumulate_t().load()) count_t;

        struct Iteration {

          public:

            Iteration() = default;

            Iteration(Iteration&& other)
              : int_screening_values(other.int_screening_values),
                int_actual_values(other.int_actual_values),
                int_distance_factors(other.int_distance_factors),
                int_distances(other.int_distances),
                int_indices(other.int_indices),
                int_ams(other.int_ams)
            {
              K_3c_needed.store(other.K_3c_needed.load());
              K_3c_needed_fxn.store(other.K_3c_needed_fxn.load());
              K_3c_dist_screened.store(other.K_3c_dist_screened.load());
              K_3c_dist_screened_fxn.store(other.K_3c_dist_screened_fxn.load());
              K_3c_underestimated.store(other.K_3c_underestimated.load());
              K_3c_underestimated_fxn.store(other.K_3c_underestimated_fxn.load());
              K_3c_perfect.store(other.K_3c_perfect.load());
              K_3c_perfect_fxn.store(other.K_3c_perfect_fxn.load());

              parent = other.parent;
              is_first = other.is_first;
            }

            accumulate_t K_3c_needed = { 0 };
            accumulate_t K_3c_needed_fxn = { 0 };
            accumulate_t K_3c_dist_screened = { 0 };
            accumulate_t K_3c_dist_screened_fxn = { 0 };
            accumulate_t K_3c_underestimated = { 0 };
            accumulate_t K_3c_underestimated_fxn = { 0 };
            accumulate_t K_3c_perfect = { 0 };
            accumulate_t K_3c_perfect_fxn = { 0 };

            ScreeningStatistics* parent;

            ThreadReplicated<std::vector<double>> int_screening_values;
            ThreadReplicated<std::vector<double>> int_actual_values;
            ThreadReplicated<std::vector<double>> int_distance_factors;
            ThreadReplicated<std::vector<double>> int_distances;
            ThreadReplicated<std::vector<std::tuple<int, int, int>>> int_indices;
            ThreadReplicated<std::vector<std::tuple<int, int, int>>> int_ams;

            void set_nthread(int nthr) {
              int_screening_values.set_nthread(nthr);
              int_actual_values.set_nthread(nthr);
              int_distance_factors.set_nthread(nthr);
              int_distances.set_nthread(nthr);
              int_indices.set_nthread(nthr);
              int_ams.set_nthread(nthr);
            }

            void global_sum(const Ref<MessageGrp>& msg) {
              auto accum_all = [&msg](accumulate_t& val) {
                long tmp = val.load();
                msg->sum(&tmp, 1);
                val.store(tmp);
              };

              accum_all(K_3c_needed);
              accum_all(K_3c_needed_fxn);
              accum_all(K_3c_dist_screened);
              accum_all(K_3c_dist_screened_fxn);
              accum_all(K_3c_underestimated);
              accum_all(K_3c_underestimated_fxn);
              accum_all(K_3c_perfect);
              accum_all(K_3c_perfect_fxn);
            }

            bool is_first = false;

        };

        mutable accumulate_t sig_pairs = { 0 };
        mutable accumulate_t sig_pairs_fxn = { 0 };
        int print_level = 0;
        mutable bool xml_stats_saved = false;

        void global_sum(const Ref<MessageGrp>& msg) const {
          auto accum_all = [&msg](accumulate_t& val) {
            long tmp = val.load();
            msg->sum(&tmp, 1);
            val.store(tmp);
          };
          accum_all(sig_pairs);
          accum_all(sig_pairs_fxn);
        }

        std::vector<Iteration> iterations;

        ScreeningStatistics()
          : iterations(),
            print_level(0)
        { }

        Iteration& next_iteration() {
          iterations.emplace_back();
          iterations.back().is_first = iterations.size() == 1;
          iterations.back().parent = this;
          return iterations.back();
        }

        void print_summary(std::ostream& out,
            const Ref<GaussianBasisSet>& basis,
            const Ref<GaussianBasisSet>& dfbs,
            int print_level = -1
        ) const;
    };

  private:

    typedef enum {
      AllPairs = 0,
      SignificantPairs = 1,
      ExchangeOuterLoopPairs = 2
    } PairSet;

    enum {
      NoMorePairs = -1
    };

    /**
     * Flags and other members that basically correspond directly to KeyVal options
     */
    //@{
    /// Number of threads to run on (defaults to threadgrp_->nthread();
    int nthread_;
    /// The threshold for including pairs in the pair list using half of the schwarz bound
    double pair_screening_thresh_;
    /// currently unused
    double density_screening_thresh_;
    /// Threshold for including a given maximum contribution to the K matrix
    double full_screening_thresh_;
    /// currently unused
    double coef_screening_thresh_;
    /** Different thresh to use for distance-including screening.  Defaults to full_screening_thresh_.
     *  Integrals will be included only if the non-distance-including estimate is greter than
     *  full_screening_thresh_ and the distance-including estimate is greater than distance_screening_thresh_.
     */
    double distance_screening_thresh_;
    /// Whether or not to do LinK style screening
    bool do_linK_;
    /// Advanced LinK options
    bool linK_block_rho_;
    /// Currently does nothing
    bool linK_sorted_B_contraction_;
    /// Use 1/r^(lX+1) factor in screening
    bool linK_use_distance_;
    /// Screening statistics print level.  Higher levels may result in a slight slowdown
    int print_screening_stats_;
    /// Exponent to raise the the distance denominator exponent to (i.e. (lX+1)^damping_factor)
    double distance_damping_factor_;
    /// Print timings after each iteration.  Useful for benchmarking large systems without doing full calculations
    bool print_iteration_timings_ = false;
    /// Use norms for nu dimension when screening
    bool use_norms_nu_ = true;
    /// Use norms for sigma dimension when screening.  It's actually a sum over this dimension rather than a norm.
    bool use_norms_sigma_ = true;
    /// Print a lot of stuff in an XML debug file.  Not really for end users...
    bool xml_debug_ = false;
    /// Dump screening data to XML
    bool xml_screening_data_ = false;
    /// Use extents to damp R distances for screening purposes
    bool use_extents_ = false;
    /// Use the maximum extent for contracted pairs.  Defaults to using coefficient weighted average
    bool use_max_extents_ = true;
    /// Should we subtract the extents from the denominator, or just use an if statement to avoid non-well-separated screening?
    bool subtract_extents_ = true;
    /// The CFMM well-separatedness thresh
    double well_separated_thresh_ = 1e-8;
    /// The old way of distributing LinK list work
    bool all_to_all_L_3_ = false;
    /// Use the exact semidiagonal integrals in J
    bool exact_diagonal_J_ = false;
    /// Use the exact semidiagonal integrals in K
    bool exact_diagonal_K_ = false;
    /// Use sparse matrix structures wherever available and reasonable
    bool use_sparse_ = false;
    /// Debugging
    bool sig_pairs_J_ = true;
    /// B use buffer
    bool B_use_buffer_ = false;
    /// B use buffer
    size_t B_buffer_size_;
    /// Scale the LinK screening threshold for differential density iterations by the ratio of the Frobenius norm of the density relative to the previous iteration
    bool scale_screening_thresh_ = true;
    /// The minimum value to scale the screening threshold down to
    double full_screening_thresh_min_ = 1e-16;
    /// Full screening exponent for non-reset iterations
    double full_screening_expon_ = 1.0;
    /// Should we screen the B contraction?
    bool screen_B_ = false;
    /// Should we distribute the coefficients among nodes?  (Note: Some replication is still necessary)
    bool distribute_coefficients_ = true;
    /// Should we use distance factor when we screen the B contraction?
    double screen_B_use_distance_ = false;
    /// Screening thresh for B contraction; defaults to full_screening_thresh_
    double B_screening_thresh_;
    /// Store the coefficients in both storage orders.  Takes more memory but also faster
    bool store_coefs_transpose_ = false;
    /// Thread the schwarz computation.  Requires more memory for 4c int evaluators and doesn't increase speed much.
    bool thread_4c_ints_ = false;
    //@}

    ScreeningStatistics stats_;
    ScreeningStatistics::Iteration* iter_stats_;

    double prev_density_frob_;
    double prev_epsilon_;
    double prev_epsilon_dist_;

    TwoCenterIntContainerPtr g2_full_ptr_;

    int max_fxn_obs_ = 0;
    int max_fxn_obs_todo_ = 0;
    int max_fxn_atom_obs_ = 0;
    int max_fxn_dfbs_ = 0;
    int max_fxn_dfbs_todo_ = 0;
    int max_fxn_atom_dfbs_ = 0;
    int max_obs_atom_fxn_on_dfbs_center_todo_ = 0;
    int max_fxn_atom_dfbs_todo_ = 0;
    //int max_fxn_obs_assigned_ = 0;
    //int max_fxn_atom_dfbs_assigned_ = 0;

    RefSCMatrix D_;

    // For now, just do static load balancing
    bool dynamic_ = false;

    TwoBodyOper::type metric_oper_type_;

    // A (very) rough estimate of the minimum amount of memory that is in use by CADF at the present time
    std::atomic<size_t> memory_used_;

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

    double get_distance_factor(
        const ShellData& ish,
        const ShellData& jsh,
        const ShellData& Xsh
    ) const;

    double get_R(
        const ShellData& ish,
        const ShellData& jsh,
        const ShellData& Xsh
    ) const;

    void compute_coefficients();

    template<template<typename...> class container, typename map_type>
    void get_coefs_ish_jsh(
        const ShellData& ish,
        const ShellData& jsh,
        int ithr,
        container<map_type>& coefsA,
        container<map_type>& coefsB
    );

    void initialize();

    void init_significant_pairs();

    void init_threads();

    //void init_vector();

    void done_threads();

    void print(std::ostream& o) const;

    /// returns shell ints in inbf x jnbf Eigen Matrix pointer
    TwoCenterIntContainerPtr ints_to_eigen(
        int ish, int jsh,
        Ref<TwoBodyTwoCenterInt>& ints,
        TwoBodyOper::type ints_type
    );

    //void ints_to_buffer(
    //    int ish, int jsh,
    //    Ref<TwoBodyTwoCenterInt>& ints,
    //    TwoBodyOper::type ints_type,
    //    double* buffer
    //);


    template <typename ShellRange>
    TwoCenterIntContainerPtr ints_to_eigen(
        const ShellBlockData<ShellRange>& ish,
        const ShellBlockData<ShellRange>& jsh,
        Ref<TwoBodyTwoCenterInt>& ints,
        TwoBodyOper::type ints_type
    );

    template <typename ShellRange>
    TwoCenterIntContainerPtr ints_to_eigen_threaded(
        const ShellBlockData<ShellRange>& ish,
        const ShellBlockData<ShellRange>& jsh,
        std::vector<Ref<TwoBodyTwoCenterInt>>& ints_for_thread,
        TwoBodyOper::type ints_type
    );

    template <typename ShellRange>
    Eigen::Map<TwoCenterIntContainer>
    ints_to_eigen_map_threaded(
        const ShellBlockData<ShellRange>& iblk,
        const ShellBlockData<ShellRange>& jblk,
        std::vector<Ref<TwoBodyTwoCenterInt>>& ints_for_thread,
        TwoBodyOper::type int_type,
        double* buffer
    );

    template <typename ShellRange>
    ThreeCenterIntContainerPtr ints_to_eigen(
        const ShellBlockData<ShellRange>& ish, const ShellData& jsh,
        Ref<TwoBodyTwoCenterInt>& ints,
        TwoBodyOper::type ints_type
    );

    /// returns ints for shell in (inbf, jnbf) x kdfnbf matrix in chemists' notation
    ThreeCenterIntContainerPtr ints_to_eigen(
        int ish, int jsh, int ksh,
        Ref<TwoBodyThreeCenterInt>& ints,
        TwoBodyOper::type ints_type
    );

    /// returns ints for shell in (inbf, jnbf) x kblk.nbf matrix in chemists' notation
    template <typename ShellRange>
    ThreeCenterIntContainerPtr ints_to_eigen(
        const ShellData& ish, const ShellData& jsh,
        const ShellBlockData<ShellRange>& kblk,
        Ref<TwoBodyThreeCenterInt>& ints,
        TwoBodyOper::type ints_type
    );

    /// returns ints for shell in (inbf, jblk.nbf) x kblk.nbf matrix in chemists' notation
    template <typename ShellRange1, typename ShellRange2>
    ThreeCenterIntContainerPtr ints_to_eigen(
        const ShellBlockData<ShellRange1>& ish,
        const ShellData& jsh,
        const ShellBlockData<ShellRange2>& kblk,
        Ref<TwoBodyThreeCenterInt>& ints,
        TwoBodyOper::type ints_type
    );

    /// returns ints for shell in (iblk.nbf, jsh.nbf) x kblk.nbf matrix in chemists' notation
    template <typename ShellRange>
    ThreeCenterIntContainerPtr ints_to_eigen(
        const ShellBlockData<ShellRange>& ish,
        const ShellData& jsh,
        const ShellData& ksh,
        Ref<TwoBodyThreeCenterInt>& ints,
        TwoBodyOper::type ints_type
    );

    template <typename ShellRange>
    Eigen::Map<ThreeCenterIntContainer> ints_to_eigen_map(
        const ShellBlockData<ShellRange>& ish,
        const ShellData& jsh,
        const ShellData& ksh,
        Ref<TwoBodyThreeCenterInt>& ints,
        TwoBodyOper::type ints_type,
        double* buffer
    );

    template <typename ShellRange1, typename ShellRange2>
    Eigen::Map<ThreeCenterIntContainer> ints_to_eigen_map(
        const ShellBlockData<ShellRange1>& ish,
        const ShellData& jsh,
        const ShellBlockData<ShellRange2>& Xblk,
        Ref<TwoBodyThreeCenterInt>& ints,
        TwoBodyOper::type ints_type,
        double* buffer
    );

    void ints_to_buffer(
        int ish, int jsh, int ksh,
        int nbfi, int nbfj, int nbfk,
        Ref<TwoBodyThreeCenterInt>& ints,
        TwoBodyOper::type ints_type,
        double* buffer, int stride=-1
    );

    /// returns ints for shell in (iblk.nbf, jsh.nbf) x kblk.nbf matrix in chemists' notation
    template <typename ShellRange>
    ThreeCenterIntContainerPtr ints_to_eigen(
        const ShellData& ish,
        const ShellBlockData<ShellRange>& jblk,
        const ShellData& ksh,
        Ref<TwoBodyThreeCenterInt>& ints,
        TwoBodyOper::type ints_type
    );

    FourCenterIntContainerPtr ints_to_eigen(
        const ShellData& ish, const ShellData& jsh,
        const ShellData& ksh, const ShellData& lsh,
        Ref<TwoBodyInt>& ints, TwoBodyOper::type ints_type
    );

    shared_ptr<Decomposition> get_decomposition(int ish, int jsh, Ref<TwoBodyTwoCenterInt> ints);

    // Non-blocked version of cl_gmat_
    RefSymmSCMatrix gmat_;

    /**
     * Whether or not the density was reset this iteration.  Numerical stability
     *  might be improved by lowering the cutoffs on iterations where the
     *  density is differential rather than full.
     */
    bool density_reset_ = true;

    /**
     * Whether or not the MessageGrp is an instance of MPIMessageGrp; might
     *  be used for dynamic load balancing in the future.
     */
    bool using_mpi_;

    /**
     * Whether or not we've computed the fitting coefficients yet
     *  (since that only needs to be done on the first iteration)
     */
    bool have_coefficients_;

    // The density fitting auxiliary basis
    Ref<GaussianBasisSet> dfbs_ = 0;

    // Where are the shell pairs being evaluated?
    std::map<PairSet, std::map<std::pair<int, int>, int>> pair_assignments_;

    // Pair assignments for K
    std::map<std::pair<int, ShellBlockSkeleton<>>, int> pair_assignments_k_;

    shared_ptr<cadf::AssignmentGrid> atom_pair_assignments_k_;

    // What pairs are being evaluated on the current node?
    std::vector<std::pair<int, int>> local_pairs_all_;
    std::vector<std::pair<int, int>> local_pairs_sig_;

    // What pairs are being evaluated on the current node?
    std::vector<std::pair<int, ShellBlockSkeleton<>>> local_pairs_k_;
    std::unordered_set<
      std::pair<int, int>, sc::hash<std::pair<int, int>>
    > local_pairs_linK_;
    std::unordered_map<int, std::vector<int>> linK_local_map_;

    // List of the permutationally unique pairs with half-schwarz bounds larger than pair_thresh_
    std::vector<std::pair<int, int>> sig_pairs_;

    std::vector<std::set<int>> sig_partners_;

    // The same as sig_pairs_, but organized differently
    std::vector<std::vector<int>> shell_to_sig_shells_;

    std::vector<double> max_schwarz_;
    Eigen::VectorXd schwarz_df_;

    // Where are we in the iteration over the local_pairs_?
    std::atomic<int> local_pairs_spot_;

    Eigen::MatrixXd schwarz_frob_;
    std::vector<Eigen::MatrixXd> C_trans_frob_;

    Eigen::MatrixXd C_bar_;

    // TwoBodyThreeCenterInt integral objects for each thread
    std::vector<Ref<TwoBodyThreeCenterInt>> eris_3c_;
    std::vector<Ref<TwoBodyThreeCenterInt>> metric_ints_3c_;

    // TwoBodyTwoCenterInt integral objects for each thread
    std::vector<Ref<TwoBodyTwoCenterInt>> eris_2c_;
    std::vector<Ref<TwoBodyTwoCenterInt>> metric_ints_2c_;

    // Integrals computed locally this iteration
    std::atomic<long> ints_computed_locally_;
    long ints_computed_;

    /// List of atom centers as Eigen::Vector3d objects, for convenience
    std::vector<Eigen::Vector3d> centers_;

    /**
     * List of orbital basis set shell pair centers, computed using
     *  the Gaussian product theorem centers of primitive pairs weighted
     *  by the absolute value of the products of primitive coefficients
     */
    std::map<std::pair<int, int>, Eigen::Vector3d> pair_centers_;
    std::map<std::pair<int, int>, double> pair_extents_;
    std::vector<double> df_extents_;

    /// Coefficients storage.  Not accessed directly
    double* __restrict__ coefficients_data_ = 0;
    double* __restrict__ dist_coefs_data_ = 0;
    double* __restrict__ dist_coefs_data_df_ = 0;

    std::unordered_map<int, Eigen::Map<RowMatrix>> coefs_X_nu;
    std::unordered_map<int, Eigen::Map<RowMatrix>> coefs_X_nu_other;
    std::unordered_map<int, Eigen::Map<RowMatrix>> coefs_mu_X;
    std::unordered_map<int, Eigen::Map<RowMatrix>> coefs_mu_X_other;

    CoefMap coefs_;

    // for each atom, matrix of mu_a in atom, rho_b, X(ab)
    std::vector<RowMatrix> coefs_blocked_;
    std::vector<std::vector<int>> coef_block_offsets_;

    //std::vector<std::vector<ShellIndexWithValue>> Cmaxes_;

    std::vector<Eigen::Map<RowMatrix>> coefs_transpose_blocked_;
    std::vector<Eigen::Map<RowMatrix>> coefs_transpose_blocked_other_;
    //std::vector<std::vector<StridedRowMap>> coefs_t_shell_blocked_;
#if USE_SPARSE
    std::vector<SparseRowMatrix> coefs_transpose_;
#else
    //std::vector<RowMatrix> coefs_transpose_;
    std::vector<Eigen::Map<RowMatrix>> coefs_transpose_;
#endif

#if USE_INTEGRAL_CACHE
    shared_ptr<TwoCenterIntCache> ints2_;
    shared_ptr<ThreeCenterIntCache> ints3_;
#endif

    shared_ptr<DecompositionCache> decomps_;
    shared_ptr<FourCenterMaxIntCache> ints4maxes_;

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

    range_of_shell_blocks<sig_partners_iter_t>
    iter_significant_partners_blocked(
        const ShellData& ish,
        int requirements = Contiguous,
        int target_size = DEFAULT_TARGET_BLOCK_SIZE
    ){
      const auto& sig_parts = sig_partners_[ish];
      return boost::make_iterator_range(
          shell_block_iterator<sig_partners_iter_t>(
              sig_parts.begin(), sig_parts.end(),
              ish.basis, ish.dfbasis, requirements, target_size
          ),
          shell_block_iterator<sig_partners_iter_t>(
              sig_parts.end(), sig_parts.end(),
              ish.basis, ish.dfbasis, requirements, target_size
          )
      );
    }

    bool is_sig_pair(int ish, int jsh) const {
      const auto& parts = sig_partners_[ish];
      return parts.find(jsh) != parts.end();
    }

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

    IndexListMap L_schwarz;
    IndexListMap L_DC;
    IndexListMap2 L_3;
    IndexListMap2 L_B;


    /**
     * Used for Eigen::SparseMatrix alternative implementation
     */
    //@{

    typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SparseMatrix;
    typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SparseColMatrix;

    SparseMatrix coefs_sp_;


    //@}

}; // CADFCLHF

/// @}
// end of addtogroup ChemistryElectronicStructureOneBodyHFCADF


boost::property_tree::ptree&
write_xml(
    const CADFCLHF::ScreeningStatistics& obj,
    ptree& parent,
    const XMLWriter& writer
);

boost::property_tree::ptree&
write_xml(
    const CADFCLHF::ScreeningStatistics::Iteration& obj,
    ptree& parent,
    const XMLWriter& writer
);

// Borrowed from ConsumeableResources
inline std::string data_size_to_string(size_t t) {
  const int prec = 3; // print this many digits

  // determine m such that 1024^m <= t <= 1024^(m+1)
  char m = 0;
  double thousand_m = 1;
  double thousand_mp1 = 1024;
  while (t >= thousand_mp1) {
    ++m;
    thousand_m = thousand_mp1;
    thousand_mp1 *= 1024;
  }

  // determine units
  std::string unit;
  switch (m) {
    case 0: unit = "B"; break;
    case 1: unit = "kB"; break;
    case 2: unit = "MB"; break;
    case 3: unit = "GB"; break;
    case 4: unit = "TB"; break;
    case 5: unit = "PB"; break;
    case 6: unit = "EB"; break;
    case 7: unit = "ZB"; break;
    case 8: unit = "YB"; break;
    default: MPQC_ASSERT(false); break;
  }

  // compute normalized mantissa
  std::ostringstream oss;
  oss.precision(prec);
  oss << (double)t/(double)thousand_m << unit;
  return oss.str();
}

} // end namespace sc

#include "get_ints_impl.h"
#include "coefs_impl.h"

#endif /* _chemistry_qc_scf_cadfclhf_h */
