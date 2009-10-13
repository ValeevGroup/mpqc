//
// prediagon.cc: some functions related to efficient inversion of B-cX
//
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@qtp.ufl.edu>
// Maintainer: TS
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

#include <sstream>
#include <algorithm>
#include <cassert>
#include <chemistry/qc/mbptr12/blas.h>
#include <chemistry/qc/basis/orthog.h>
#include <chemistry/qc/ccr12/ccr12_info.h>
#include <chemistry/qc/ccr12/mtensor.h>

using namespace sc;
using namespace std;

void CCR12_Info::prediagon(RefDiagSCMatrix& eigvals, RefSCMatrix& eigvecs) {
	// CCR12_Info::B_ and CCR12_Info::X_ required.
	// B_ and X_ are RefSymmSCMatrix objects.
	assert(B_.nonnull() && X_.nonnull());

	OverlapOrthog xorthog(OverlapOrthog::Symmetric, X_, X_.kit(), OverlapOrthog::default_lindep_tol());
	RefSCMatrix mtilde = xorthog.basis_to_orthog_basis();

	// So far, I haven't thought about the reduced dimensional version...
	// Guessing occupied pairs can hardly be linearly dependent with each other.
	assert(mtilde.coldim() == X_.dim());

	RefSCMatrix btilde = mtilde.t() * B_ * mtilde;

	// TODO can I diagonalize RefSCMatrix using only one RefSCMatrix for eigen vectors?
	RefSCMatrix U(B_.dim(), B_.dim(), B_.kit());
	RefSCMatrix V(B_.dim(), B_.dim(), B_.kit());
	RefDiagSCMatrix beig(mtilde.coldim(), mtilde.kit());
	btilde.svd(U, beig, V);

	eigvals = beig;
	eigvecs = mtilde * U;

//#define LOCAL_DEBUG_PREDIAGON
#ifdef LOCAL_DEBUG_PREDIAGON
	RefDiagSCMatrix unit(mtilde.coldim(), mtilde.kit());
	unit.assign(1.0);

	(B_ - 3.0 * X_).gi().print();
    (eigvecs * (eigvals - 3.0 * unit).gi() * eigvecs.t()).print();
#endif

}


void CCR12_Info::denom_contraction(const Ref<Tensor>& in, Ref<Tensor>& out) {
  const size_t singles = maxtilesize() * maxtilesize();
  const size_t doubles = singles * singles;
  double* k_a0      = mem()->malloc_local_double(doubles);
  double* k_a1      = mem()->malloc_local_double(doubles);
  double* k_c       = mem()->malloc_local_double(doubles);
  double* k_c_sort  = mem()->malloc_local_double(singles);

  const int nocc_act = naoa();
  MTensor<4>::tile_ranges iiii(4, MTensor<4>::tile_range(0, noab()));
  MTensor<4>::element_ranges iiii_erange(4, MTensor<4>::element_range(0, nocc_act) );
  vector<long> amap;
  {
    vector<int> intmap = sc::map(*(r12evalinfo()->refinfo()->occ_act_sb(Alpha)), *corr_space(), false);
    amap.resize(intmap.size());
    std::copy(intmap.begin(), intmap.end(), amap.begin());
  }

  // this loop structure minimizes the O(o^8) operation...
  int count = 0;
  for (long h3b = 0L; h3b < noab(); ++h3b) {
    for (long h4b = h3b; h4b < noab(); ++h4b) {
      for (int h3 = 0; h3 < get_range(h3b); ++h3) {
        for (int h4 = 0; h4 < get_range(h4b); ++h4, ++count) {
          if (get_offset(h3b) + h3 >= get_offset(h4b) + h4) continue;
          if (count % mem()->n() != mem()->me() ) continue;

          // orbital energies
          const double eh3 = get_orb_energy(get_offset(h3b) + h3);
          const double eh4 = get_orb_energy(get_offset(h4b) + h4);

          // current denominator tensor
          stringstream ss;
          ss << "denom" << count;
          Ref<Tensor> denom = new Tensor(ss.str(), mem());
          offset_gt2(denom, false);
          MTensor<4> D(this, denom.pointer(), iiii);

          RefSymmSCMatrix refxminusb = X() * (eh3 + eh4) - B();
          RefSymmSCMatrix refinverse = refxminusb.gi();

          D.convert(refinverse, nocc_act, nocc_act, false, false,
                    amap, amap, amap, amap, &iiii_erange);

          for (long h1b = 0L; h1b < noab(); ++h1b) {
            for (long h2b = h1b; h2b < noab(); ++h2b) {

              if (!restricted() || get_spin(h3b)+get_spin(h4b)+get_spin(h1b)+get_spin(h2b) != 8L) {
                if (get_spin(h3b)+get_spin(h4b) == get_spin(h1b)+get_spin(h2b)) {
                  if ((get_sym(h3b)^(get_sym(h4b)^(get_sym(h1b)^get_sym(h2b)))) == irrep_t()) {

                    // target size
                    const long dimc_singles = get_range(h1b)*get_range(h2b);
                    std::fill(k_c_sort, k_c_sort+dimc_singles, 0.0);

                    for (long h5b = 0L; h5b < noab(); ++h5b) {
                      for (long h6b = h5b; h6b < noab(); ++h6b) {

                        if (get_spin(h3b)+get_spin(h4b) == get_spin(h5b)+get_spin(h6b)) {
                          if ((get_sym(h3b)^(get_sym(h4b)^(get_sym(h5b)^get_sym(h6b)))) == irrep_e()) {

                            long h3b_0, h4b_0, h5b_0, h6b_0;
                            restricted_4(h3b, h4b, h5b, h6b, h3b_0, h4b_0, h5b_0, h6b_0);
                            long h5b_1, h6b_1, h1b_1, h2b_1;
                            restricted_4(h5b, h6b, h1b, h2b, h5b_1, h6b_1, h1b_1, h2b_1);
                            const int dim_common = get_range(h5b) * get_range(h6b);
                            const int dima0_sort = get_range(h3b) * get_range(h4b);
                            const int dima1_sort = get_range(h1b) * get_range(h2b);

                            // read tilde V (redundant read is involved...)
                            in->get_block(h4b_0+noab()*(h3b_0+noab()*(h6b_0+noab()*(h5b_0))), k_a0);

                            // read denominator
                            denom->get_block(h6b_1+noab()*(h5b_1+noab()*(h2b_1+noab()*(h1b_1))), k_a1);
                            double factor = 1.0;
                            if (h5b == h6b) factor *= 0.5;
                            const int unit = 1;
                            const int h34 = h4 + get_range(h4b) * h3;
                            const int stride = get_range(h3b) * get_range(h4b);
                            const double one = 1.0;
                            F77_DGEMV("t", &dim_common, &dima1_sort,
                            		&factor, k_a1, &dim_common,
                            		k_a0 + h34, &stride,
                            		&one, k_c_sort, &unit);
                          }
                        }
                      }
                    }
                    const long dimc = get_range(h3b)*get_range(h4b)*get_range(h1b)*get_range(h2b);
                    std::fill(k_c, k_c+dimc, 0.0);
                    {
                      const size_t h34 = h4 + get_range(h4b) * h3;
                      const size_t stride = get_range(h3b) * get_range(h4b);
                      size_t iall = h34;
                      size_t i = 0;
                      for (int h1 = 0; h1 != get_range(h1b); ++h1) {
                        for (int h2 = 0; h2 != get_range(h2b); ++h2, iall += stride, ++i) {
                          k_c[iall] += k_c_sort[i];
                        }
                      }
                    }
                    // taking care of permutation symmetry.
                    if (h3b == h4b) {
                      const size_t h43 = h3 + get_range(h4b) * h4;
                      const size_t stride = get_range(h3b) * get_range(h4b);
                      size_t iall = h43;
                      size_t i = 0;
                      for (int h1 = 0; h1 != get_range(h1b); ++h1) {
                        for (int h2 = 0; h2 != get_range(h2b); ++h2, iall += stride, ++i) {
                          k_c[iall] -= k_c_sort[i];
                        }
                      }
                    }
                    out->add_block(h4b+noab()*(h3b+noab()*(h2b+noab()*(h1b))),k_c);
                  }
                }
              }
            }
          }
        }
      } // orbital loops
    }
  }
  mem()->free_local_double(k_a1);
  mem()->free_local_double(k_a0);
  mem()->free_local_double(k_c_sort);
  mem()->free_local_double(k_c);
  mem()->sync();
}

