
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

#ifndef _chemistry_qc_lmp2_lcorr_h
#define _chemistry_qc_lmp2_lcorr_h

#include <math.h>

#include <util/class/scexception.h>

#include <math/optimize/diis.h>
#include <math/scmat/repl.h>
#include <math/scmat/matrix.h>

#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/scf/scf.h>

#include <chemistry/qc/lmp2/sma.h>
#include <chemistry/qc/lmp2/domain.h>
#include <chemistry/qc/lmp2/util.h>

namespace sc {

/**  \brief A base class for local correlation methods.
*/
class LCorr: public sc::Wavefunction {
    /// The cutoff for S eigenvalues used to compute W.
    double W_eigval_threshold_;

    /** The set of basis functions that contribute to each virb.
     * The is needed for computing W. */
    std::vector<std::vector<int> > virb_to_bfns_;

    // compute_W fills these in:
    std::map<domainmapvirbs_t, sma2::Array<2> >      unique_W_;
    std::map<domainmapvirbs_t, std::vector<double> > unique_eigvals_;
    std::map<domainmapvirbs_t, sma2::Array<2> >      unique_F_tilde_;

    /// Places index_maps for converting W from SCMatrix form to Array<2> form
    /// in index_map1 and index_map2.  index_map2 is the map for the nonredundant
    /// index and index_map1 is for the redundant (PAO) index.
    void compute_W_index_maps(int blockdim, int blockdim_nonred,
                              std::set<int> &virbs,
                              std::vector<int> &index_map1,
                              std::vector<int> &index_map2,
                              const sma2::Range &vir);

  protected:

    int n_unique_W() const { return unique_W_.size(); }
    sma2::Array<2>            &unique_W(const domainmapvirbs_t &i);
    std::vector<double> &unique_eigvals(const domainmapvirbs_t &i);
    sma2::Array<2>      &unique_F_tilde(const domainmapvirbs_t &i);

    /// This must be called before compute_W is called.
    void init_virb_to_bfns(const sma2::Range &vir);

    void compute_W(domainmapvirbs_t &virset,
                   const sc::Ref<sc::GaussianBasisSet> &basis,
                   sc::RefSCMatrix &F_vir_mat, sc::RefSCMatrix &S_mat,
                   const sma2::Range &vir, int nocc_act,
                   double bound);

    /// Transform A to D using transformation matrices B,C: D = B^T*A*C
    void transform_array(sma2::Array<2> &A, sma2::Array<2> &B,
                         sma2::Array<2> &C, sma2::Array<2> &D,
                         const sc::Ref<sc::MessageGrp> &msg);

    /// Release stored data.
    void clear();

    /// Print input parameters.
    void print_parameters() const;
    
  public:
    /** Construct an LCorr object from KeyVal input.

        This reads the following input:
        <table border="1">

        <tr><td>Keyword<td>Type<td>Default<td>Description

        <tr><td><tt>W_eigval_threshold</tt><td>double<td>1.0e-6<td>
        The threshold for eigenvalues in the domain overlap
        matrices.

        </table>

     */
    LCorr(const sc::Ref<sc::KeyVal> &);
    LCorr(sc::StateIn &);
    ~LCorr();
    void save_data_state(sc::StateOut &);
};

}

#endif // _chemistry_qc_lmp2_lccsd_h
