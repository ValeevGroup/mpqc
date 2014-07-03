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
#include "ordered_shells.h"
#include "treemat_fwd.h"

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
class ApproximatePairWriter;
namespace cadf {
  class AssignmentGrid;
  namespace assignments {
    struct Assignments;
  }
}

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
    double new_density();
    void init_vector();

  public:

    // typedefs

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
    typedef Eigen::Map<RowMatrix, Eigen::Unaligned, Eigen::OuterStride<>> StridedRowMap;
    template<typename T1, typename T2>
    using fast_map = std::unordered_map<T1, T2, sc::hash<T1>>;

    // TODO Get rid of this
    typedef ConcurrentCache<
        double,
        int, int, int, int, TwoBodyOper::type
    > FourCenterMaxIntCache;
    typedef ConcurrentCache<shared_ptr<Decomposition>, int, int> DecompositionCache;

    CADFCLHF(StateIn&);

    /// TODO Document KeyVal options
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
    } PairSet;

    /**
     * Flags and other members that basically correspond directly to KeyVal options
     */
    //@{
    /// Number of threads to run on (defaults to threadgrp_->nthread());
    int nthread_;
    /// The threshold for including pairs in the pair list using half of the schwarz bound
    double pair_screening_thresh_;
    /// Threshold for including a given maximum contribution to the K matrix
    double full_screening_thresh_;
    /// currently unused
    double coef_screening_thresh_;
    /** Different thresh to use for distance-including screening.  Defaults to full_screening_thresh_.
     *  Integrals will be included only if the non-distance-including estimate is greater than
     *  full_screening_thresh_ and the distance-including estimate is greater than distance_screening_thresh_.
     */
    double distance_screening_thresh_;
    /// Whether or not to do LinK style screening
    bool do_linK_;
    /** Reorder the density matrix to block all rho in a given $L^{(3})_{\mu X}$ together
     *  as one contraction.  Defaults to true if screen_B is true, false otherwise.
     */
    bool linK_block_rho_;
    /// Use 1/r^(lX+1) factor in screening
    bool linK_use_distance_ = true;
    /// Screening statistics print level.  Higher levels may result in a slight slowdown
    int print_screening_stats_;
    /** Exponent to raise the the distance denominator exponent to (i.e. (l_X+1)^damping_factor)
     *  Use of values other than 1.0 for this option is not recommended.
     */
    double distance_damping_factor_ = 1.0;
    /** Print timings after each iteration.  Useful for benchmarking large systems without doing
     *  full calculations.  Note that the times printed currently exclude some sections of the
     *  code, and thus it is better to use the timings printed at the end instead.
     */
    bool print_iteration_timings_ = false;
    /** Use norms for nu dimension when screening.  This should pretty much always be true
     *  (though it doesn't usually make much difference).  use_norms_nu = false is not
     *  implemented with distribute_coefficients = true
     */
    bool use_norms_nu_ = true;
    /** Use norms to get L_DC.  If false, an N_shell^3 contraction will
     *  be done, which only starts to dominate for VERY large systems.
     *  Defaults to true, but note that unless sigma_norms_chunk_by_atoms is
     *  explicitly set to false (which gives *substantially* worse screening),
     *  an N_shell^2*N_atom contraction will still be done.
     */
    bool use_norms_sigma_ = true;
    /** When using norms for L_DC, we can also do a N_shell^2*N_atom contraction
     *  that is slightly cheaper but screens better than the full column norms
     *  (only relevent if use_norms_sigma = true).  Defaults to true.  Mostly
     *  useful to make the scaling of the L_DC computation look better, even though
     *  it isn't a big part of the computation.
     */
    bool sigma_norms_chunk_by_atoms_ = true;
    /** Print a lot of stuff in an XML debug file.  Also, if compiled with assertions
     *  enabled, exit after the first iteration (by asserting false).  Not for end users.
     */
    bool xml_debug_ = false;
    /** Dump screening data to XML.  This was mostly for making plots of screening
     *  efficiency and accuracy.  Unless you're doing that, this should be false.
     */
    bool xml_screening_data_ = false;
    /** Use extents to damp R distances for screening purposes.  Defaults to true.
     *  Setting this to false can lead to slightly stronger screening, but also leads to
     *  less rigorous screening.
     */
    bool use_extents_ = true;
    /** Use the maximum extent for contracted pairs.  If false,  uses a coefficient-weighted
     *  average.  This has very little effect on the final result (though it may be a little
     *  more significant if the basis has lots of contracted functions, such as the ANO basis
     *  sets).  The safer option of use_max_extents = true is prefered and is thus the default.
     */
    bool use_max_extents_ = true;
    /** Use subtracted the extents from the denominator.  If false, an if statement is used to
     *  avoid non-well-separated screening.  Subtracting extents leads to less "tight" screening,
     *  but it is much safer, particularly for higher angular momentum.  Defaults to true.
     */
    bool subtract_extents_ = true;
    /** The CFMM well-separatedness thresh.  Defaults to 1e-1, which is the recommended value by
     *  Ochsenfeld et. al.
     */
    double well_separated_thresh_ = 1e-1;
    /**
     * Use the exact semidiagonal integrals in the coulomb matrix computation.
     * distribute_coefficients = true is not implemented for exact_diagonal_J = true.
     * Defaults to false.
     */
    bool exact_diagonal_J_ = false;
    /** Use the exact semidiagonal integrals in the exchange matrix computation.
     *  distribute_coefficients = true and/or screen_B = true is not implemented
     *  for exact_diagonal_K = true.  Defaults to false.  Code for this option is
     *  not heavily optimized.
     */
    bool exact_diagonal_K_ = false;
    /** Use a buffer to minimize the number of small contractions when computing the B intermediate
     *  in the exchange matrix.  Setting this to true is not compatible with screen_B = true and/or
     *  distribute_coefficients = true and/or linK_block_rho = true, and thus this option is deprecated.
     *  Almost always slows things down anyway.
     */
    bool B_use_buffer_ = false;
    /// The size of the B buffer.  Only relevant if B_use_buffer = true.
    size_t B_buffer_size_;
    // TODO individual options to scale other screening thresholds
    /** Scale the LinK screening threshold(s) for differential density iterations by the ratio of the
     *  Frobenius norm of the density relative to the previous iteration.  Defaults to true.
     *  This is much safer and should almost always be enabled, unless the screening thresholds are
     *  quite small to begin with.  Defaults to true.
     */
    bool scale_screening_thresh_ = true;
    // TODO individual screening threshold minima
    /** The minimum value to scale the screening threshold down to.
     *  It is usually a good idea to set this to something not too small so that the cost of
     *  later differential iterations doesn't get out of hand.  Defaults to 10^-3 * full_screening_thresh
     */
    double full_screening_thresh_min_;
    /** Full screening exponent for differential density iterations, i.e. the full_screening_thresh for
     *  differential density iterations will be full_screening_thresh^full_screening_expon.  This
     *  was the old way to deal with differential density iterations.  The new way is to use
     *  scale_screening_thresh = true, which is recommended.  Using values other than 1.0 here
     *  is not recommended.
     */
    double full_screening_expon_ = 1.0;
    /** Screen the B intermediate formation when building the exchange matrix.  This is recommended
     *  for large and most medium-sized molecules.  Ignored if do_linK = false.
     */
    bool screen_B_ = false;
    /** Distribute the coefficients among nodes if true.
     *  (Note: Some replication is still necessary.  The decrease in coefficient memory
     *  is currently only proportional to sqrt(p)/2).  Note that as currently implemented,
     *  the number of mpi tasks *must* be a perfect square for this to work.
     *  Not implemented for do_linK = false.
     */
    bool distribute_coefficients_ = true;
    /// Should we use distance factor when we screen the B contraction?
    double screen_B_use_distance_ = false;
    /** Screening thresh for B intermediate construction in the exchange
     *  matrix build.  Defaults to full_screening_thresh.
     */
    double B_screening_thresh_;
    /** Store the coefficients in both storage orders.  Takes more memory but also a little bit faster.
     *  Irrelevant if distribute_coefficients = true.  Setting this option is deprecated and may
     *  be removed in the near future.
     */
    bool store_coefs_transpose_ = false;
    /** Thread the schwarz matrix computation.  Requires more memory for 4c int evaluators
     *  and doesn't increase speed much.  Memory increase is pretty small also.  Defaults to
     *  false.
     */
    bool thread_4c_ints_ = false;
    /** Screening threshold for d_overtilde intermediate in the exchange matrix build.
     *  Only relevent if screen_B = true and distribute_coefficients = true.  Defaults
     *  to B_screening_thresh.
     */
    double d_over_screening_thresh_;
    /** Screening threshold for d_undertilde intermediate in the exchange matrix build.
     *  Only relevent if screen_B = true and distribute_coefficients = true.  Defaults
     *  to B_screening_thresh * 1e-1.
     */
    double d_under_screening_thresh_;
    /** Whether or not Frobenius norms should be used in computing the B screening
     *  lists.  It is actually faster, more accurate, and more efficient to set this
     *  to false, though setting it to true could have a similar effect to raising
     *  the B screening threshold.
     */
    bool use_norms_B_ = false;
    /** Whether or not to shuffle the L_3 and L_3_star key lists.  Shuffling may improve thread
     *  load balancing.
     */
    bool shuffle_L_3_keys_ = true;
    /** Whether or not to shuffle assignments list for the J computation.  Shuffling can improve thread
     *  and node load balancing.
     */
    bool shuffle_J_assignments_ = true;
    /** Match up orbitals with previous iteration to make sure we aren't including "attractive electron"
     *  orbitals.
     */
    bool match_orbitals_ = true;
    /** Maximum norm of difference for which orbitals other than the minimum difference on can be considered
     *  the same for the purposes of orbital matching
     */
    double match_orbitals_max_diff_ = 1e-4;
    /** Maximum norm of difference for which orbitals other than the minimum difference on can be considered
     *  the same for the purposes of orbital matching
     */
    double match_orbitals_max_homo_offset_ = 0.05;
    /** Whether or not to do an SVD of the possibly occupied orbitals to determine the best means of
     *  projecting out the bad orbitals
     */
    bool match_orbitals_use_svd_ = true;
    /** Compute and print the Coulomb energy for debugging purposes
     */
    bool debug_coulomb_energy_ = false;
    /** Compute and print the exchange energy for debugging purposes
     */
    bool debug_exchange_energy_ = false;
    /** Use new exchange algorithm that recomputes some coefficients on the
     *  fly.  Ignores a lot of other options.  See source code for details.
     */
    bool new_exchange_algorithm_ = true;
    /** Minimum number of DFBS atoms to be assigned to a given node.  Smaller
     *  numbers will hurt load balancing, but may be necessary to make the
     *  coefficients fit into memory.  Defaults to 3.
     */
    int min_atoms_per_node_ = 3;
    //@}

    std::shared_ptr<ScreeningStatistics> stats_;
    ScreeningStatistics::Iteration* iter_stats_;

    double prev_density_frob_;
    double prev_epsilon_;
    double prev_epsilon_dist_;
    double prev_epsilon_B_;
    double prev_epsilon_d_over_;
    double prev_epsilon_d_under_;

    RefDiagSCMatrix prev_evals_;
    RefSCMatrix prev_evecs_;
    std::vector<int> prev_occ_;

    int max_fxn_obs_ = 0;
    int max_fxn_obs_j_ish_ = 0;
    int max_fxn_obs_j_jsh_ = 0;
    int max_fxn_obs_todo_ = 0;
    int max_fxn_atom_obs_ = 0;
    int max_fxn_dfbs_ = 0;
    int max_fxn_dfbs_todo_ = 0;
    int max_fxn_atom_dfbs_ = 0;
    int max_obs_atom_fxn_on_dfbs_center_todo_ = 0;
    int max_fxn_atom_dfbs_todo_ = 0;
    //int max_fxn_obs_assigned_ = 0;
    //int max_fxn_atom_dfbs_assigned_ = 0;

    TwoCenterIntContainerPtr g2_full_ptr_;

    RefDiagSCMatrix core_evals_;
    RefSCMatrix core_evecs_;

    std::vector<int> orb_to_hcore_orb_;

    RefSCMatrix D_;

    // A (very) rough estimate of the minimum amount of memory that is in use by CADF at the present time
    std::atomic<size_t> memory_used_;

    class TrackedMemoryHolder {
        std::atomic<size_t>& used_ref_;
        size_t size_;
      public:
        TrackedMemoryHolder(
            std::atomic<size_t>& accum,
            size_t size
        ) : used_ref_(accum), size_(size)
        {
          used_ref_ += size_;
        }
        ~TrackedMemoryHolder() { used_ref_ -= size_; }
    };
    void consume_memory(size_t size) {
      memory_used_ += size;
    }
    void release_memory(size_t size) {
      memory_used_ -= size;
    }
    TrackedMemoryHolder hold_memory(size_t size) {
      return TrackedMemoryHolder(memory_used_, size);
    }

    // Convenience variable for better code readibility
    static const TwoBodyOper::type coulomb_oper_type_ = TwoBodyOper::eri;
    // Currently only metric_oper_type_ == coulomb_oper_type_ is supported
    TwoBodyOper::type metric_oper_type_;

    bool get_shell_pair(ShellData& mu, ShellData& nu, PairSet pset = AllPairs);

    Eigen::MatrixXd J_;

    RefSCMatrix compute_J();

    RefSCMatrix compute_K();

    RefSCMatrix new_compute_K();

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

    void done_threads();

    void print(std::ostream& o) const;

    /// returns shell ints in inbf x jnbf Eigen Matrix pointer
    TwoCenterIntContainerPtr ints_to_eigen(
        int ish, int jsh,
        Ref<TwoBodyTwoCenterInt>& ints,
        TwoBodyOper::type ints_type,
        const GaussianBasisSet* bs1=0,
        const GaussianBasisSet* bs2=0
    );

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
        TwoBodyOper::type ints_type,
        GaussianBasisSet* bs1 = 0,
        GaussianBasisSet* bs2 = 0,
        GaussianBasisSet* bs3 = 0
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
        const ShellData& ish,
        const ShellData& jsh,
        const ShellBlockData<ShellRange>& Xsh,
        Ref<TwoBodyThreeCenterInt>& ints,
        TwoBodyOper::type ints_type,
        double* __restrict__ buffer
    );

    template <typename MapType>
    void ints_to_eigen_map(
        const ShellData& ish,
        const ShellData& jsh,
        const ShellData& Xsh,
        Ref<TwoBodyThreeCenterInt>& ints,
        TwoBodyOper::type int_type,
        MapType& out_map
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

    /** Whether or not the density was reset this iteration.  Numerical stability
     *  is improved by lowering the cutoffs on iterations where the
     *  density is differential rather than full.
     */
    bool density_reset_ = true;

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
    shared_ptr<cadf::assignments::Assignments> assignments_new_k_;

    // What pairs are being evaluated on the current node?
    std::vector<std::pair<int, int>> local_pairs_all_;
    std::vector<std::pair<int, int>> local_pairs_sig_;

    // What pairs are being evaluated on the current node?
    std::vector<std::pair<int, ShellBlockSkeleton<>>> local_pairs_k_;
    std::unordered_set<
      std::pair<int, int>, sc::hash<std::pair<int, int>>
    > local_pairs_linK_;
    std::unordered_map<int, std::vector<int>> linK_local_map_;
    std::unordered_map<int, std::vector<int>> linK_local_map_ish_Xsh_;

    // List of the permutationally unique pairs with half-schwarz bounds larger than pair_thresh_
    std::vector<std::pair<int, int>> sig_pairs_;

    // Significant partners for a given shell index
    std::vector<std::set<int>> sig_partners_;

    // Significant partners in pre-blocked form
    std::vector<ContiguousShellBlockList> sig_partner_blocks_;

    /// (X|X)^{1/2} Schwarz integrals for the DF basis
    Eigen::VectorXd schwarz_df_;
    /// (X|X)^{1/2} Schwarz integrals for
    Eigen::VectorXd schwarz_df_mine_;

    // Where are we in the iteration over the local_pairs_?
    std::atomic<int> local_pairs_spot_;

    /// Shell-wise Frobenius norms of orbital basis pairs
    Eigen::MatrixXd schwarz_frob_;

    /// Various ways to take Frobenius norms of the coefficient tensor
    Eigen::MatrixXd C_bar_;
    Eigen::MatrixXd C_bar_mine_;
    Eigen::MatrixXd C_underbar_;
    shared_ptr<cadf::TreeMatrix<>> C_dfsame_;

    // TwoBodyThreeCenterInt integral objects for each thread
    std::vector<Ref<TwoBodyThreeCenterInt>> eris_3c_;
    std::vector<Ref<TwoBodyThreeCenterInt>> metric_ints_3c_;

    // TwoBodyTwoCenterInt integral objects for each thread
    std::vector<Ref<TwoBodyTwoCenterInt>> eris_2c_;
    std::vector<Ref<TwoBodyTwoCenterInt>> metric_ints_2c_;

    // Count of integrals computed locally this iteration
    std::atomic<long> ints_computed_locally_;

    /// List of atom centers as Eigen::Vector3d objects, for convenience
    std::vector<Eigen::Vector3d> centers_;

    /**
     * List of orbital basis set shell pair centers, computed using
     *  the Gaussian product theorem centers of primitive pairs weighted
     *  by the absolute value of the products of primitive coefficients
     */
    fast_map<std::pair<int, int>, Eigen::Vector3d> pair_centers_;

    // Extents of the various pairs
    fast_map<std::pair<int, int>, double> pair_extents_;
    std::vector<double> df_extents_;

    /// Coefficients storage.  Not accessed directly
    double* __restrict__ coefficients_data_ = 0;
    double* __restrict__ dist_coefs_data_ = 0;
    double* __restrict__ dist_coefs_data_df_ = 0;

    fast_map<int, Eigen::Map<RowMatrix>> coefs_X_nu;
    fast_map<int, Eigen::Map<RowMatrix>> coefs_X_nu_other;
    fast_map<int, Eigen::Map<RowMatrix>> coefs_mu_X;
    fast_map<int, Eigen::Map<RowMatrix>> coefs_mu_X_other;

    CoefMap coefs_;

    // for each atom, matrix of mu_a in atom, rho_b, X(ab)
    std::vector<RowMatrix> coefs_blocked_;
    std::vector<std::vector<int>> coef_block_offsets_;

    std::vector<Eigen::Map<RowMatrix>> coefs_transpose_blocked_;
    std::vector<Eigen::Map<RowMatrix>> coefs_transpose_blocked_other_;

    std::vector<Eigen::Map<RowMatrix>> coefs_transpose_;

    shared_ptr<DecompositionCache> decomps_;

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
    typedef madness::ConcurrentHashMap<int, OrderedShellList> IndexListMap;
    typedef madness::ConcurrentHashMap<
        std::pair<int, int>,
        OrderedShellList,
        sc::hash<std::pair<int, int>>
    > IndexListMap2;
    typedef madness::ConcurrentHashMap<
        std::pair<int, int>,
        std::vector<std::pair<uint64_t, uint64_t>>,
        sc::hash<std::pair<int, int>>
    > LinKRangeMap2;

    IndexListMap L_schwarz;
    IndexListMap L_DC;
    IndexListMap L_C_under;
    IndexListMap2 L_3;
    IndexListMap2 L_3_star;
    IndexListMap2 L_B;
    IndexListMap2 L_d_over;
    LinKRangeMap2 L_d_under_ranges;

    friend class ApproximatePairWriter;

    Ref<ApproximatePairWriter> pair_writer_;

}; // CADFCLHF

/// @}
// end of addtogroup ChemistryElectronicStructureOneBodyHFCADF

void
get_split_range_part(
    const Ref<MessageGrp>& msg,
    int full_begin, int full_end,
    int& out_begin, int& out_size
);

void
get_split_range_part(
    int me, int n,
    int full_begin, int full_end,
    int& out_begin, int& out_size
);

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
