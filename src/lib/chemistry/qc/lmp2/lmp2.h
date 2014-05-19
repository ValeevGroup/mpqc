
/*
 * Copyright 2009 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 *
 * This file is a part of the MPQC LMP2 library.
 *
 * The MPQC LMP2 library is free software: you can redistribute it
 * and/or modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _chemistry_qc_lmp2_lmp2_h
#define _chemistry_qc_lmp2_lmp2_h

#include <math.h>

#include <math/optimize/diis.h>
#include <math/scmat/repl.h>
#include <math/scmat/matrix.h>

#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/scf/scf.h>

#include <chemistry/qc/lmp2/sma.h>
#include <chemistry/qc/lmp2/domain.h>
#include <chemistry/qc/lmp2/util.h>
#include <chemistry/qc/lmp2/lcorr.h>
#include <chemistry/qc/lmp2/parallel.h>

namespace sc {

typedef std::vector<std::pair<int,int> > my_occ_pairs_t;

typedef std::set<std::pair<int,int> > k_2occ_local_pairs_t;

typedef std::set<std::pair<int,int> > k_3_4_occ_local_pairs_t;

typedef std::map<sma2::triplet<int,int,int>, domainmapvirbs_t > domainmap_triple;

/** \brief Computes the local second order perturbation theory energy.
 *
 * This class implements a massively parallel algorithm for the
 * quantum mechanical computation of energies of molecular systems
 * using local second-order M&oslash;ller-Plesset (LMP2) perturbation
 * theory. It is built upon the Massively Parallel Quantum Chemistry
 * program (MPQC). Both the storage requirement and the computational
 * time scale linearly with the molecular size.  High parallel
 * efficiency of the algorithm has been demonstrated for applications
 * employing up to 100 processors.
 *
 * The parallel algorithm is designed to be scalable, employing a
 * distributed data scheme for the two-electron integrals, avoiding
 * communication bottlenecks, and distributing tasks in all
 * computationally significant steps. A sparse data representation and
 * a set of generalized contraction routines have been developed to
 * allow efficient massively parallel implementation using distributed
 * sparse multidimensional arrays.
 *
 * The implementation of this software and its performance are
 * described in Nielsen, I. M. B.; Janssen, C. L.; <i>J. Chem. Theory
 * and Comput.</i>; 2007; pp. 71-79.
 *
*/
class LMP2: public LCorr {

    sc::Ref<sc::MessageGrp> msg_;

    sc::Ref<sc::SCF> ref_;
    
    int max_iter_;

    int nfzc_;
    
    double bound_;
    double S_threshold_;
    double integral_threshold_;
    double q1_threshold_;
    double q2_threshold_;
    double q3_threshold_;
    double q4_threshold_;
    double threshold_factor_;
    double distance_threshold_;
    double completeness_threshold_;

    // "canonical" or "pipek-mezey"
    std::string occ_orbitals_;
    // "canonical" or "projected_atomic"
    std::string vir_orbitals_;
    
    int parallel_transform123a_;
    int parallel_transform123b_;
    int parallel_transform4a_;
    int parallel_transform4b_;
    int parallel_iterations_;
    int completedomains_;
    bool singlevirb_;

    bool i_ge_j_;
    bool m_ge_n_;
    bool minimize_q2_;
    bool always_use_dist_t_;

    sc::Ref<sc::SelfConsistentExtrapolation> extrap_T_;

    /** The range for AO indices.  The block size is the number of basis
        functions in each shell. */
    sma2::Range ao_;
    /** The range for projected AO virtual indices.  The block size is
        the number of basis functions on each atom. */
    sma2::Range vir_;
    /** The range for active occupied indices.  The block size is one. */
    sma2::Range occ_act_;
    /** The projection matrix P(ao,vir) (AO -> projected AO).  This is the
        matrix P in [1].  P = I - D*S, where D is one half times the AO
        density matrix and S is the AO overlap matrix. */
    sma2::Array<2> P_;
    /** The AO to localized MO transformation matrix,  L(ao,occ_act). */
    sma2::Array<2> L_;
    /** The mapping from pairs of MO's to a set of virtual blocks
     * (virbs) for strong pairs. */
    domainmap_t domainmap_;
    /** The active occupied pairs stored on this node.  The
        orbitals numbers are offset by nfzc. */
    my_occ_pairs_t my_occ_pairs_;
    /** The vector containing for each active occupied mo, i, the set of
        mo's j for which the ij pair is included in T and K_2occ. */
    std::vector<std::set<int> > paired_occ_;

    // Maps occ pairs to node numbers.
    sc::Ref<sma2::PairMapping> pair_mapping_;

    /** K_2occ(i,j,r,s) = P(mu,r) * L(rho,i) * P(nu,s) * L(sig,j) *
       (mu,rho|nu,sig), where r and s are projected atomic orbitals and i
       and j are localized, occupied, active MO's.  Thus, the
       correspondence between indices and electrons 1 and 2 is (1,2,1,2).
       This arrangement is due to the contraction routines which require
       fixed indices to be first.  */
    sma2::Array<4> K_2occ_;

    /** The replicated doubles, T(occ_act,occ_act,vir,vir). This is only
        used if always_use_dist_t_ is false. */
    sma2::Array<4> T_;
    /** This stores another mapping of the replicated T.  It refers to the
        same data as T.  It allows the j index to be fixed. */
    sma2::Array<4> T_jirs_;
    /** These are the distributed T. */
    sma2::Array<4> T_local_;
    /// The number of elements in T, globally
    double T_n_element_;

    sma2::Array<2> F_occ_;
    sma2::Array<2> F_vir_;
    sma2::Array<2> S_;

    std::vector<double> F_diag_;

    /// Free storage for the above arrays.
    void clear();

    /// useful statistics
    bool analyze_occ_orbs_;
    std::vector<SCVector3> r_i_;      //!< "the-center-of-mass" of active occupied orbitals (replicated)
    std::vector<double> dist_ij_;     //!< spatial distance between domains of active orbitals i and j (replicated),
                                      //!< not the same as the distance between centers of mass of i and j
    std::vector<double> emp2_ij_;     //!< energy for (active) pair ij (replicated)
    void analyze_occ_orbs(const RefSCMatrix& scf_local);

    double compute_lmp2_energy();
    double compute_ecorr_lmp2();

    void rearrange_q2_all_ij(sma2::Array<4> &q2_K_oo);
    void rearrange_q2_i_ge_j(sma2::Array<4> &q2_K_oo);

    void compute_K_2occ(
        sc::RefSCMatrix &T_schwarz, sc::RefSCVector &T_schwarz_maxvec,
        double T_schwarz_maxval, sc::RefSCMatrix &Pmax,
        sc::RefSCVector &Pmaxvec,
        std::vector<std::vector<double> > &Lmax,
        sc::RefSCMatrix &Lshellmax,
        std::vector<std::vector<double> > &PPmax,
        const std::vector<std::vector<double> > &Dmax, double &Dmax_maxval, std::vector<double> &Lmaxvec,
        std::vector<std::multimap<double,int,std::greater<double> > > &L_map);

    void compute_Schwarz_screening_quantities(sc::Ref<sc::TwoBodyInt> &tbint,
                                              sc::RefSCMatrix &T_schwarz,
                                              sc::RefSCVector &T_schwarz_maxvec,
                                              double &T_schwarz_maxval,
                                              int nshell);
    
    void compute_P_screening_quantities(sc::RefSCMatrix &Pmax,
                                        sc::RefSCMatrix &Pmatrix,
                                        sc::RefSCVector &Pmaxvec,
                                        double &Pmax_maxval);

    void compute_L_and_D_screening_quantities(std::vector<std::vector<double> >  &Dmax,
                                              double &Dmax_maxval,
                                              std::vector<std::vector<double> > &Lmax,
                                              std::vector<double> &Lmaxvec,
                                              sc::RefSCMatrix &Lshellmax,
                                              std::vector<std::vector<int> > &L_blocks,
                                              std::vector<std::multimap<double,int,std::greater<double> > > &L_map,
                                              double &Lmax_maxval,
                                              double T_schwarz_maxval,
                                              double Pmax_maxval);

    void compute_PPmax_screening_matrix(std::vector<std::vector<double > > &PPmax,
                                        sc::RefSCMatrix &Pmax,
                                        const std::vector<std::set<int> >
                                           &virb_united_pair_domains,
                                        double T_schwarz_maxval,
                                        double Pmax_maxval,
                                        double Lmax_maxval);

    void compute_doubles_W();

    void compute_LMP2_residual(sma2::Array<4> &Res);
    void compute_LMP2_residual_SFT(sma2::Array<4> &R);
    void compute_LMP2_residual_SFTS(sma2::Array<4> &R);
    void compute_LMP2_residual_SFTS_i(sma2::Array<4> &R);
    void compute_LMP2_residual_SFTS_ij(sma2::Array<4> &R);

    double iterate_LMP2_equations(double energy_tolerance, double rms_tolerance);

    void compute_delta_T(sma2::Array<4> &Res,
                         sma2::Array<4> &delta_T);
    
  public:
    /** Construct an LMP2 object from KeyVal input.

        This reads the keywords in the table below. The <a
        href="./html/classLCorr.html">LCorr</a> input has additional
        keywords that are also read.

        <table border="1">

        <tr><td>%Keyword<td>Type<td>Default<td>Description

        <tr><td><tt>reference</tt><td>object<td>none<td> Gives an object that
        specializes OneBodyWavefunction which can be used as a reference LMP2
        wavefunction. This keyword is required and has no default.

        <tr><td><tt>max_iter</tt><td>int<td>100<td>The maximum number
        of iterations that will be used to converge the LMP2 amplitudes.

        <tr><td><tt>nfzc</tt><td>int or <tt>auto</tt><td>0<td>The number of
        frozen core orbitals. An integer can be specified or <tt>auto<tt> can
        be given if the number of frozen core orbitals is to be automatically
        determined.

        <tr><td><tt>distance_threshold</tt><td>double<td>15.0<td>
        The distance threshold used to create the domains, in bohr.

        <tr><td><tt>completeness_threshold</tt><td>double<td>0.02<td>
        The completeness threshold used to create the domains.

        </table>

        <p>In addition to the above options, are a variety of options that
        is useful to developers, listed in the following table.
        </p>

        <table border="1">

        <tr><td>Keyword<td>Type<td>Default<td>Description

        <tr><td><tt>S_threshold</tt><td>double<td>1.0e-6<td>
        The threshold used to determine which blocks of the overlap
        matrix are stored.

        <tr><td><tt>integral_threshold</tt><td>double<td>1.0e-8<td>
        A threshold used to prune the number of two electron integrals
        which are computed.

        <tr><td><tt>q1_threshold</tt><td>double<td>1.0e-8<td>
        A threshold used to prune the number of first quarter transformed
        integrals which are computed.

        <tr><td><tt>q2_threshold</tt><td>double<td>1.0e-8<td>
        A threshold used to prune the number of second quarter transformed
        integrals which are computed.

        <tr><td><tt>q3_threshold</tt><td>double<td>1.0e-8<td>
        A threshold used to prune the number of third quarter transformed
        integrals which are computed.

        <tr><td><tt>q4_threshold</tt><td>double<td>1.0e-8<td>
        A threshold used to prune the number of fourth quarter transformed
        integrals which are computed.

        <tr><td><tt>threshold_factor</tt><td>double<td>10.0<td>
        The desired energy tolerance is divided by this factor to obtain the
        threshold that determines which blocks of a variety of arrays are
        stored.

        <tr><td><tt>complete_domains</tt><td>bool<td>0<td>
        If true, all possible blocks are computed. The result in this
        case should be the same as the canonical MP2 energy.

        <tr><td><tt>single_virtual_block</tt><td>bool<td>0<td>
        All of the projected atomic virtual orbitals are grouped into a
        single, large block, rather than atom based blocks.

        <tr><td><tt>i_ge_j</tt><td>bool<td>1<td>
        If true, use i >= j in the formation of K_2occ.

        <tr><td><tt>m_ge_n</tt><td>bool<td>1<td>
        If true, use m >= n in the formation of K_2occ.

        <tr><td><tt>minimize_q2</tt><td>bool<td>1<td>
        If this and <tt>i_ge_j</tt> are true, use the i >= j symmetry when
        rearranging the second quarter tranformed integrals.

        <tr><td><tt>always_use_dist_t</tt><td>bool<td>0<td>
        If true, distribute the amplitudes, in which case an extra
        communication step is required each iteration.

        <tr><td><tt>extrap_T</tt><td>object<td>DIIS<td> This gives an object
        derived from SelfConsistentExtrapolation that will be used to
        extraplote the amplitudes.  By default DIIS(6,8,0.0,2,1) will be used.

        <tr><td><tt>occ_orbitals</tt><td>string<td><tt>pipek-mezey</tt><td>
        The method used to provide the occupied orbitals. This can
        be <tt>pipek-mezey</tt> or <tt>canonical</tt>.

        <tr><td><tt>vir_orbitals</tt><td>string<td><tt>projected_atomic</tt><td>
        The method used to provide the virtual orbitals. This can
        be <tt>projected_atomic</tt>, <tt>old_projected_atomic</tt>
        or <tt>canonical</tt>. If <tt>canonical</tt> then a single block
        is used for the virtual orbitals, rather than atom-based blocking.

        </table>

     */
    LMP2(const sc::Ref<sc::KeyVal> &);
    LMP2(sc::StateIn &);
    ~LMP2();
    void save_data_state(sc::StateOut &);
    void compute(void);
    int nelectron(void);
    sc::RefSymmSCMatrix density(void);
    double magnetic_moment() const;
    int value_implemented(void) const;

};

}

#endif
