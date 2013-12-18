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
#include <util/misc/property.h>

typedef typename boost::integer_range<int> int_range;


namespace sc {

// Forward declarations
class ShellData;
class BasisFunctionData;

class shell_iter_wrapper {
    Ref<GaussianBasisSet> basis_;
    Ref<GaussianBasisSet> dfbasis_;
    int first_shell_;
    int last_shell_;
  public:
    shell_iter_wrapper(
        const Ref<GaussianBasisSet>& basis,
        const Ref<GaussianBasisSet>& dfbasis = 0,
        int first_shell = 0,
        int last_shell = -1
    );
    ShellIter begin() const;
    ShellIter end() const;
};

class function_iter_wrapper {
    Ref<GaussianBasisSet> basis_;
    Ref<GaussianBasisSet> dfbasis_;
    int first_function_;
    int last_function_;
  public:
    function_iter_wrapper(
        const Ref<GaussianBasisSet>& basis,
        const Ref<GaussianBasisSet>& dfbasis = 0,
        int first_function = 0,
        int last_function = -1
    );
    BasisFunctionIter begin() const;
    BasisFunctionIter end() const;
};

class ShellIter {

  public:

    ShellIter() = delete;

    ShellIter(const Ref<GaussianBasisSet>& basis, int position);
    ShellIter(const Ref<GaussianBasisSet>& basis, const Ref<GaussianBasisSet>& dfbasis, int position);

    bool operator!=(const ShellIter& other) const;

    const ShellIter& operator++();

    const ShellData operator*() const;

  private:

    int pos_;
    Ref<GaussianBasisSet> basis_;
    Ref<GaussianBasisSet> dfbasis_;

};

class BasisFunctionIter {

  public:

    BasisFunctionIter() = delete;

    BasisFunctionIter(const Ref<GaussianBasisSet>& basis, int position);
    BasisFunctionIter(const Ref<GaussianBasisSet>& basis, const Ref<GaussianBasisSet>& dfbasis, int position);

    bool operator!=(const BasisFunctionIter& other) const;

    const BasisFunctionIter& operator++();

    const BasisFunctionData operator*() const;

  private:

    int pos_;
    Ref<GaussianBasisSet> basis_;
    Ref<GaussianBasisSet> dfbasis_;

};


#define ASSERT_SHELL_BOUNDS \
  assert(basis.nonnull() && "null basis"); \
  assert(index < basis->nshell() && "index larger than nshell()"); \
  assert(index >= 0 && "index less than 0")

struct ShellData {

    ShellData();

    ShellData(
        int idx,
        Ref<GaussianBasisSet> basis,
        Ref<GaussianBasisSet> dfbasis = 0
    );

    ShellData(const ShellData&);

    ShellData& operator=(const ShellData& other);

    int index;
    property<ShellData, int> bfoff;
    property<ShellData, int> nbf;
    property<ShellData, int> center;
    property<ShellData, int> atom_bfoff;
    property<ShellData, int> atom_shoff;
    property<ShellData, int> atom_nsh;
    property<ShellData, int> atom_nbf;
    property<ShellData, int> bfoff_in_atom;
    property<ShellData, int> shoff_in_atom;
    property<ShellData, int> atom_last_function;
    property<ShellData, int> atom_last_shell;
    property<ShellData, int> last_function;

    // Used when an auxiliary basis is set in the parent ShellIter.  Otherwise, set to -1
    property<ShellData, int> atom_dfshoff;
    property<ShellData, int> atom_dfbfoff;
    property<ShellData, int> atom_dfnbf;
    property<ShellData, int> atom_dfnsh;
    property<ShellData, int> atom_df_last_function;
    property<ShellData, int> atom_df_last_shell;

    operator int() { ASSERT_SHELL_BOUNDS; return index; }

    Ref<GaussianBasisSet> basis;
    Ref<GaussianBasisSet> dfbasis;

  private:
    inline int get_nbf() const { ASSERT_SHELL_BOUNDS; return basis->shell(index).nfunction(); }
    inline int get_bfoff() const { ASSERT_SHELL_BOUNDS; return basis->shell_to_function(index); }
    inline int get_center() const { ASSERT_SHELL_BOUNDS; return basis->shell_to_center(index); }
    inline int get_atom_bfoff() const { ASSERT_SHELL_BOUNDS; return basis->shell_to_function(atom_shoff); }
    inline int get_atom_shoff() const { ASSERT_SHELL_BOUNDS; return basis->shell_on_center(center, 0); }
    inline int get_atom_nsh() const { ASSERT_SHELL_BOUNDS; return basis->nshell_on_center(center); }
    inline int get_atom_nbf() const { ASSERT_SHELL_BOUNDS; return basis->nbasis_on_center(center); }
    inline int get_shoff_in_atom() const { ASSERT_SHELL_BOUNDS; return index - atom_shoff; }
    inline int get_bfoff_in_atom() const { ASSERT_SHELL_BOUNDS; return bfoff - atom_bfoff; }
    inline int get_atom_last_function() const { ASSERT_SHELL_BOUNDS; return atom_bfoff + atom_nbf - 1; }
    inline int get_last_function() const { ASSERT_SHELL_BOUNDS; return bfoff + nbf - 1; }
    inline int get_atom_last_shell() const { ASSERT_SHELL_BOUNDS; return atom_shoff + atom_nsh - 1; }
    inline int get_atom_dfshoff() const { assert(dfbasis.nonnull()); return dfbasis->shell_on_center(center, 0); }
    inline int get_atom_dfbfoff() const { assert(dfbasis.nonnull()); return dfbasis->shell_to_function(atom_dfshoff); }
    inline int get_atom_dfnsh() const { assert(dfbasis.nonnull()); return dfbasis->nshell_on_center(center); }
    inline int get_atom_dfnbf() const { assert(dfbasis.nonnull()); return dfbasis->nbasis_on_center(center); }
    inline int get_atom_df_last_function() const { assert(dfbasis.nonnull()); return atom_dfbfoff + atom_dfnbf - 1; }
    inline int get_atom_df_last_shell() const { assert(dfbasis.nonnull()); return atom_dfshoff + atom_dfnsh - 1; }
    inline int assert_not_initialized() const { assert(false && "ShellData object not initialized"); return -1; }

};


struct BasisFunctionData {

    BasisFunctionData(
        int idx,
        const Ref<GaussianBasisSet>& basis,
        const Ref<GaussianBasisSet>& dfbasis = 0
    );

    BasisFunctionData() { }

    int index;
    int shell_index;
    int shell_nbf;
    int shell_bfoff;
    int center;
    int atom_bfoff;
    int atom_shoff;
    int atom_nsh;
    int atom_nbf;
    int bfoff_in_atom;
    int bfoff_in_shell;
    int shoff_in_atom;
    int atom_last_function;
    int atom_last_shell;

    GaussianBasisSet::Shell* shell;

    // Used when an auxiliary basis is set in the parent BasisFunctionIter.  Otherwise, set to -1
    int atom_dfshoff = -1;
    int atom_dfbfoff = -1;
    int atom_dfnbf = -1;
    int atom_dfnsh = -1;
    int atom_df_last_function = -1;
    int atom_df_last_shell = -1;

    operator int() { return index; }

};

//============================================================================//
// shell_iterator

const shell_iter_wrapper
shell_iterator(const Ref<GaussianBasisSet>& basis, const Ref<GaussianBasisSet>& dfbasis = 0);

const shell_iter_wrapper
shell_iterator(const Ref<GaussianBasisSet>& basis, int last_shell);

const shell_iter_wrapper
shell_iterator(const Ref<GaussianBasisSet>& basis, int first_shell, int last_shell);

const shell_iter_wrapper
shell_iterator(const Ref<GaussianBasisSet>& basis, const Ref<GaussianBasisSet>& dfbasis, int last_shell);

const shell_iter_wrapper
shell_iterator(const Ref<GaussianBasisSet>& basis, const Ref<GaussianBasisSet>& dfbasis, int first_shell, int last_shell);

//============================================================================//
// function_iterator

const function_iter_wrapper
function_iterator(const Ref<GaussianBasisSet>& basis, const Ref<GaussianBasisSet>& dfbasis = 0);

const function_iter_wrapper
function_iterator(const Ref<GaussianBasisSet>& basis, int last_function, const Ref<GaussianBasisSet>& dfbasis = 0);

const function_iter_wrapper
function_iterator(const Ref<GaussianBasisSet>& basis, int first_function, int last_function, const Ref<GaussianBasisSet>& dfbasis = 0);

const function_iter_wrapper
function_iterator(const Ref<GaussianBasisSet>& basis, const Ref<GaussianBasisSet>& dfbasis, int first_function, int last_function);

const function_iter_wrapper
function_iterator(const Ref<GaussianBasisSet>& basis, const Ref<GaussianBasisSet>& dfbasis, int last_function);

const function_iter_wrapper
function_iterator(const ShellData& ish);

//============================================================================//

const shell_iter_wrapper
iter_shells_on_center(const Ref<GaussianBasisSet>& basis, int center, const Ref<GaussianBasisSet>& dfbasis = 0);

const function_iter_wrapper
iter_functions_on_center(const Ref<GaussianBasisSet>& basis, int center, const Ref<GaussianBasisSet>& dfbasis = 0);

const function_iter_wrapper
iter_functions_in_shell(const Ref<GaussianBasisSet>& basis, int shell, const Ref<GaussianBasisSet>& dfbasis = 0);

ShellData shell_data(Ref<GaussianBasisSet> basis, int ish, Ref<GaussianBasisSet> dfbasis = 0);

const BasisFunctionData function_data(const Ref<GaussianBasisSet>& basis, int ish, const Ref<GaussianBasisSet>& dfbasis = 0);

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
