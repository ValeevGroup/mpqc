//
// utils.h
//
// Copyright (C) 2005 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
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

#ifdef __GNUG__
#pragma interface
#endif

#include <vector>
#include <math/scmat/matrix.h>
#include <chemistry/qc/mbptr12/spin.h>
#include <chemistry/qc/mbptr12/pairiter.h>

#ifndef _chemistry_qc_mbptr12_utils_h
#define _chemistry_qc_mbptr12_utils_h

namespace sc {
  
  class MOIndexSpace;
  
  /** Takes the 4-index quantity <ij|A|kl> and returns,
      depending on the value of the PureSpinCase2 spin,
      a Singlet or Triplet spinadapted matrix.
    */
  template <PureSpinCase2 spin>
  RefSCMatrix spinadapt(const RefSCMatrix &A,
                          const Ref<MOIndexSpace> &bra,
                          const Ref<MOIndexSpace> &ket);
  /** Antisymmetrizes 4-index quantity <ij|A|kl> -> <ij|A|kl> - <ij|A|lk>
      and saves to Aanti. Row dimension of A has to be an integer multiple of
      bra->rank()*bra->rank(). Same for ket.
    */
  void antisymmetrize(RefSCMatrix& Aanti, const RefSCMatrix& A,
                      const Ref<MOIndexSpace>& bra,
                      const Ref<MOIndexSpace>& ket,
                      bool accumulate = false);
  /** Generalization of the above.
      Antisymmetrizes 4-index quantity <ij|A|kl>.
      antisymmetrize only makes sense
      if either bra1==bra2 or ket1==ket2.
      If bra1==bra2: <ij|A|kl> - <ji|A|kl>.
      If ket1==ket2: <ij|A|kl> - <ij|A|lk>.
      Row dimension of A has to be an integer multiple of
      bra1->rank()*bra2->rank(). Same for ket.
      The row dimension of Aanti: bra1==bra2 ? bra1->rank()*(bra1->rank()-1)/2.
      The col dimension of Aanti: ket1==ket2 ? ket1->rank()*(ket1->rank()-1)/2
    */
  template <bool accumulate>
    void antisymmetrize(RefSCMatrix& Aanti, const RefSCMatrix& A,
                        const Ref<MOIndexSpace>& bra1,
                        const Ref<MOIndexSpace>& bra2,
                        const Ref<MOIndexSpace>& ket1,
                        const Ref<MOIndexSpace>& ket2);
  /** Specialization of the above to symmetric matrices. bra2 must be equal to bra1, hence not needed.
      Antisymmetrizes 4-index quantity <ij|A|kl>.
      <ij|Aanti|kl> = <ij|A|kl> - <ij|A|lk>.
      Dimension of A has to be an integer multiple of
      bra1->rank()*bra1->rank().
      The dimension of Aanti: bra1->rank()*(bra1->rank()-1)/2.
    */
  template <bool accumulate>
    void antisymmetrize(RefSymmSCMatrix& Aanti, const RefSymmSCMatrix& A,
                        const Ref<MOIndexSpace>& bra1);
  /** Symmetrizes 4-index quantity <ij|A|kl> -> 1/2 * (<ij|A|kl> + <ji|A|lk>)
      and saves to Asymm. Row dimension has to be an integer multiple of
      bra->rank()*bra->rank(). Same for ket. Asymm and A can be the same matrix.
    */
  template <bool Accumulate>
    void symmetrize(RefSCMatrix& Asymm, const RefSCMatrix& A,
                    const Ref<MOIndexSpace>& bra,
                    const Ref<MOIndexSpace>& ket);

  /** Generalization of the above. Symmetrizes 4-index quantity <ij|A|kl> -> 1/2 * (<ij|A|kl> + <ji|A|lk>),
      where either bra or ket may have been antisymmetrized.
      If BraSymm==ASymm && KetSymm==ASymm: <ij|A|kl> -> 1/2 * (<ij|A|kl> + <ji|A|lk>)
      If BraSymm==ASymm && KetSymm==AntiSymm: <ij|A|kl> -> 1/2 * (<ij|A|kl> + <ji|A|lk>) = 1/2 * (<ij|A|kl> - <ji|A|kl>)
      If BraSymm==AntiSymm && KetSymm==ASymm: <ij|A|kl> -> 1/2 * (<ij|A|kl> + <ji|A|lk>) = 1/2 * (<ij|A|kl> - <ij|A|lk>)
      If BraSymm==ASymm && KetSymm==Symm: <ij|A|kl> -> 1/2 * (<ij|A|kl> + <ji|A|lk>) = 1/2 * (<ij|A|kl> + <ji|A|kl>)
      If BraSymm==Symm && KetSymm==ASymm: <ij|A|kl> -> 1/2 * (<ij|A|kl> + <ji|A|lk>) = 1/2 * (<ij|A|kl> + <ij|A|lk>)
      etc.
      and saves to Asymm. Row dimension of A has to be pairdim<BraSymm>(bra1->rank(),bra2->rank()).n().
      Same for ket. Asymm and A can be the same matrix.
    */
  template <bool Accumulate, sc::fastpairiter::PairSymm BraSymm, sc::fastpairiter::PairSymm KetSymm>
    void symmetrize12(RefSCMatrix& Asymm, const RefSCMatrix& A,
                      const Ref<MOIndexSpace>& bra1,
                      const Ref<MOIndexSpace>& bra2,
                      const Ref<MOIndexSpace>& ket1,
                      const Ref<MOIndexSpace>& ket2);

  /** Symmetrizes/antisymmetrizes bra and/or ket. Sanity is checked as fully as possible.
      If DstBraSymm!=ASymm, bra1 must equal bra2.
      If DstKetSymm!=ASymm, ket1 must equal ket2.
      Row dimension of A has to be an integer multiple of
      bra1->rank()*bra2->rank(). Same for ket.
      Row dimension of A has to be pairdim<SrcBraSymm>(bra1->rank(),bra2->rank()).n().
      Same for ket, etc.
      Asymm and A cannot be the same matrix.
    */
  template <bool Accumulate,
    sc::fastpairiter::PairSymm SrcBraSymm,
    sc::fastpairiter::PairSymm SrcKetSymm,
    sc::fastpairiter::PairSymm DstBraSymm,
    sc::fastpairiter::PairSymm DstKetSymm
    >
    void symmetrize(RefSCMatrix& Aanti, const RefSCMatrix& A,
                    const Ref<MOIndexSpace>& bra1,
                    const Ref<MOIndexSpace>& bra2,
                    const Ref<MOIndexSpace>& ket1,
                    const Ref<MOIndexSpace>& ket2);

  /** Converts RefDiagSCMatrix to std::vector<double>
  */
  std::vector<double> convert(const RefDiagSCMatrix& A);
  
  /// print out the Fortran-style matrix
  void print_f77_mat(const std::string& comment,
                     const double* A,
                     unsigned int nrow,
                     unsigned int ncol,
                     bool transpose = false);
  
  /// Returns the lower triangle of the matrix B (which should be symmetric)
  RefSymmSCMatrix to_lower_triangle(const RefSCMatrix& B);
  
}

#endif

