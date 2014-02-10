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
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> TwoCenterIntContainer;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ThreeCenterIntContainer;
    typedef std::shared_ptr<TwoCenterIntContainer> TwoCenterIntContainerPtr;
    typedef std::shared_ptr<ThreeCenterIntContainer> ThreeCenterIntContainerPtr;

    typedef std::map<std::pair<int, int>, std::pair<CoefContainer, CoefContainer>> CoefMap;
    typedef Eigen::HouseholderQR<Eigen::MatrixXd> Decomposition;
    typedef ConcurrentCache<
        TwoCenterIntContainerPtr,
        int, int, TwoBodyOper::type
    > TwoCenterIntCache;
    typedef ConcurrentCache<
        ThreeCenterIntContainerPtr,
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

    class ScreeningStatistics {


      public:

        typedef std::atomic_uint_fast64_t accumulate_t;

        struct Iteration {

          public:
            // TODO count screening due to Schwarz pair screening vs. LinK screening

            Iteration() = default;

            Iteration(Iteration&& other)
            {
              K_3c_needed.store(other.K_3c_needed.load());
              K_3c_needed_fxn.store(other.K_3c_needed_fxn.load());
              K_3c_dist_screened.store(other.K_3c_dist_screened.load());
              K_3c_dist_screened_fxn.store(other.K_3c_dist_screened_fxn.load());
              K_3c_underestimated.store(other.K_3c_underestimated.load());
              K_3c_underestimated_fxn.store(other.K_3c_underestimated_fxn.load());
              K_3c_perfect.store(other.K_3c_perfect.load());
              K_3c_perfect_fxn.store(other.K_3c_perfect_fxn.load());
            }

            accumulate_t K_3c_needed = { 0 };
            accumulate_t K_3c_needed_fxn = { 0 };
            accumulate_t K_3c_dist_screened = { 0 };
            accumulate_t K_3c_dist_screened_fxn = { 0 };
            accumulate_t K_3c_underestimated = { 0 };
            accumulate_t K_3c_underestimated_fxn = { 0 };
            accumulate_t K_3c_perfect = { 0 };
            accumulate_t K_3c_perfect_fxn = { 0 };

            //accumulate_t ints_2c_needed = { 0 };
            //accumulate_t ints_2c_underestimated = { 0 };

        };

        std::vector<Iteration> iterations;

        ScreeningStatistics() : iterations() { }

        Iteration& next_iteration() {
          iterations.emplace_back();
          return iterations.back();
        }

        void print_summary(std::ostream& out,
            const Ref<GaussianBasisSet>& basis,
            const Ref<GaussianBasisSet>& dfbs,
            int print_level
        ) const
        {
          using std::endl;
          using std::setw;
          const auto& old_loc = out.getloc(); out.imbue(std::locale(""));
          out << indent << "CADFCLHF Screening Statistics" << endl;
          out << indent << "-----------------------------" << endl;
          const int total_3c = basis->nshell() * basis->nshell() * dfbs->nshell();
          const int total_3c_fxn = basis->nbasis() * basis->nbasis() * dfbs->nbasis();
          out << indent << "Total shell triplets: " << total_3c << endl
              << indent << "Total function triplets: " << total_3c_fxn << endl;
          int iteration = 1;
          out << incindent;
          out << indent << setw(38) << " "
              << setw(22) << std::internal <<  "Shell-wise"
              << setw(22) << std::internal <<  "Function-wise"
              << endl;
          out << decindent;

          for(auto& iter : iterations) {
            if(print_level > 1) {
              out << indent << "Iteration " << iteration++ << ":" << endl;
              auto pr_iter_stat = [&](const std::string& title, int sh_num, int fxn_num) {
                out << indent << setw(38) << std::left << title
                    << setw(14) << std::right << sh_num
                    << setw(7) << scprintf(" (%3.1f", 100.0 * double(sh_num)/double(total_3c)) << "%)"
                    << setw(14) << std::right << fxn_num
                    << setw(7) << scprintf(" (%3.1f", 100.0 * double(fxn_num)/double(total_3c_fxn)) << "%)"
                    << std::endl;
              };
              out << incindent;
              pr_iter_stat("K: 3c ints needed", iter.K_3c_needed, iter.K_3c_needed_fxn);
              pr_iter_stat("K: 3c ints screened by distance", iter.K_3c_dist_screened, iter.K_3c_dist_screened_fxn);
              if(print_level > 2) {
                pr_iter_stat("K: 3c ints underestimated", iter.K_3c_underestimated, iter.K_3c_underestimated_fxn);
                pr_iter_stat("K: 3c ints needed, \"perfect screening\"", iter.K_3c_perfect, iter.K_3c_perfect_fxn);
                pr_iter_stat("K: \"extra\" ints computed",
                    iter.K_3c_needed - iter.K_3c_perfect,
                    iter.K_3c_needed_fxn - iter.K_3c_perfect_fxn
                );
              }
              out << decindent;
            }
          }
          out.imbue(old_loc);
        }
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
    /// Full screening exponent for non-reset iterations
    double full_screening_expon_;
    /// Use 1/r^(lX+1) factor in screening
    bool linK_use_distance_;
    /// Screening statistics print level.  Higher levels may result in a slight slowdown
    int print_screening_stats_;
    /// Exponent to raise the the distance denominator exponent to (i.e. (lX+1)^damping_factor)
    double distance_damping_factor_;
    //@}

    ScreeningStatistics stats_;
    ScreeningStatistics::Iteration* iter_stats_;

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

    void done_threads();

    void print(std::ostream& o) const;

    /// returns shell ints in inbf x jnbf Eigen Matrix pointer
    TwoCenterIntContainerPtr ints_to_eigen(
        int ish, int jsh,
        Ref<TwoBodyTwoCenterInt>& ints,
        TwoBodyOper::type ints_type
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

    /// returns ints for shell in (iblk.nbf, jsh.nbf) x kblk.nbf matrix in chemists' notation
    template <typename ShellRange>
    ThreeCenterIntContainerPtr ints_to_eigen(
        const ShellData& ish,
        const ShellBlockData<ShellRange>& jblk,
        const ShellData& ksh,
        Ref<TwoBodyThreeCenterInt>& ints,
        TwoBodyOper::type ints_type
    );

    std::shared_ptr<Decomposition> get_decomposition(int ish, int jsh, Ref<TwoBodyTwoCenterInt> ints);

    // Non-blocked version of cl_gmat_
    RefSymmSCMatrix gmat_;

    bool density_reset_ = true;

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

    std::vector<Eigen::Vector3d> centers_;
    std::map<std::pair<int, int>, Eigen::Vector3d> pair_centers_;

    // Coefficients storage.  Not accessed directly
    double* coefficients_data_ = 0;

    CoefMap coefs_;

    std::vector<std::vector<ShellIndexWithValue>> Cmaxes_;

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
    IndexListMap L_coefs;
    IndexListMap L_D;
    IndexListMap L_DC;
    IndexListMap2 L_3;

};

template <typename ShellRange>
CADFCLHF::TwoCenterIntContainerPtr
CADFCLHF::ints_to_eigen(
    const ShellBlockData<ShellRange>& iblk,
    const ShellBlockData<ShellRange>& jblk,
    Ref<TwoBodyTwoCenterInt>& ints,
    TwoBodyOper::type int_type
){
  auto rv = std::make_shared<TwoCenterIntContainer>(iblk.nbf, jblk.nbf);
  int block_offset_i = 0;
  for(auto ish : shell_range(iblk)) {
    int block_offset_j = 0;
    for(auto jsh : shell_range(jblk)) {
      const auto& ints_ptr = ints_to_eigen(ish, jsh, ints, int_type);
      rv->block(block_offset_i, block_offset_j, ish.nbf, jsh.nbf) = *ints_ptr;
      block_offset_j += jsh.nbf;
    }
    block_offset_i += ish.nbf;
  }
  return rv;
}

template <typename ShellRange>
CADFCLHF::TwoCenterIntContainerPtr
CADFCLHF::ints_to_eigen_threaded(
    const ShellBlockData<ShellRange>& iblk,
    const ShellBlockData<ShellRange>& jblk,
    std::vector<Ref<TwoBodyTwoCenterInt>>& ints_for_thread,
    TwoBodyOper::type int_type
){
  auto rv = std::make_shared<TwoCenterIntContainer>(iblk.nbf, jblk.nbf);
  // non-contiguous form not implemented yet
  assert(iblk.is_contiguous());
  assert(jblk.is_contiguous());
  const int nthread = threadgrp_->nthread();
  do_threaded(nthread, [&](int ithr){
    ShellData ish, jsh;
    for(auto pair : threaded_shell_block_pair_range(iblk, jblk, ithr, nthread)){
      boost::tie(ish, jsh) = pair;
      const auto& ints_ptr = ints_to_eigen(ish, jsh, ints_for_thread[ithr], int_type);
      rv->block(
          ish.bfoff - iblk.bfoff, jsh.bfoff - jblk.bfoff,
          ish.nbf, jsh.nbf
      ) = *ints_ptr;
    }
  });
  return rv;
}

template <typename ShellRange>
CADFCLHF::ThreeCenterIntContainerPtr
CADFCLHF::ints_to_eigen(
    const ShellBlockData<ShellRange>& iblk,
    const ShellData& jsh,
    Ref<TwoBodyTwoCenterInt>& ints,
    TwoBodyOper::type int_type
){
  auto rv = std::make_shared<ThreeCenterIntContainer>((const int)iblk.nbf, (const int)jsh.nbf);
  for(auto ish : shell_range(iblk)) {
    const auto& ints_ptr = ints_to_eigen(ish, jsh, ints, int_type);
    rv->block(ish.bfoff - iblk.bfoff, 0, ish.nbf, jsh.nbf) = *ints_ptr;
  }
  return rv;
}

template <typename ShellRange>
CADFCLHF::ThreeCenterIntContainerPtr
CADFCLHF::ints_to_eigen(
    const ShellData& ish, const ShellData& jsh,
    const ShellBlockData<ShellRange>& Xblk,
    Ref<TwoBodyThreeCenterInt>& ints,
    TwoBodyOper::type int_type
){
  auto rv = std::make_shared<ThreeCenterIntContainer>(
      ish.nbf * jsh.nbf,
      Xblk.nbf
  );
  for(auto Xsh : shell_range(Xblk)) {
    const auto& ints_ptr = ints_to_eigen(ish, jsh, Xsh, ints, int_type);
    rv->middleCols(Xsh.bfoff - Xblk.bfoff, Xsh.nbf) = *ints_ptr;
  }
  return rv;
}

template <typename ShellRange1, typename ShellRange2>
CADFCLHF::ThreeCenterIntContainerPtr
CADFCLHF::ints_to_eigen(
    const ShellBlockData<ShellRange1>& iblk,
    const ShellData& jsh,
    const ShellBlockData<ShellRange2>& Xblk,
    Ref<TwoBodyThreeCenterInt>& ints,
    TwoBodyOper::type int_type
){
  auto rv = std::make_shared<ThreeCenterIntContainer>(
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

template <typename ShellRange>
CADFCLHF::ThreeCenterIntContainerPtr
CADFCLHF::ints_to_eigen(
    const ShellBlockData<ShellRange>& iblk,
    const ShellData& jsh,
    const ShellData& Xsh,
    Ref<TwoBodyThreeCenterInt>& ints,
    TwoBodyOper::type int_type
){
  auto rv = std::make_shared<ThreeCenterIntContainer>(
      iblk.nbf * jsh.nbf, Xsh.nbf
  );
  int block_offset = 0;
  for(auto ish: shell_range(iblk)) {
    const auto& ints_ptr = ints_to_eigen(ish, jsh, Xsh, ints, int_type);
    rv->block(
        block_offset * jsh.nbf, 0,
        ish.nbf*jsh.nbf, Xsh.nbf
    ) = *ints_ptr;
    block_offset += ish.nbf;
  }
  return rv;
}

template <typename ShellRange>
CADFCLHF::ThreeCenterIntContainerPtr
CADFCLHF::ints_to_eigen(
    const ShellData& ish,
    const ShellBlockData<ShellRange>& jblk,
    const ShellData& Xsh,
    Ref<TwoBodyThreeCenterInt>& ints,
    TwoBodyOper::type int_type
){
  auto rv = std::make_shared<ThreeCenterIntContainer>(
      ish.nbf * jblk.nbf, Xsh.nbf
  );
  int block_offset = 0;
  for(auto jsh : shell_range(jblk)) {
    const auto& ints_ptr = ints_to_eigen(ish, jsh, Xsh, ints, int_type);
    for(auto mu : function_range(ish)){
      rv->block(
          mu.bfoff_in_shell*jblk.nbf + block_offset, 0,
          jsh.nbf, Xsh.nbf
      ) = ints_ptr->block(mu.bfoff_in_shell*jsh.nbf, 0, jsh.nbf, Xsh.nbf);
    }
    block_offset += jsh.nbf;
  }
  return rv;
}

} // end namespace sc

#endif /* _chemistry_qc_scf_cadfclhf_h */
