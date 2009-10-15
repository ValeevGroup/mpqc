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

    Ref<SCMatrixKit> kit = SCMatrixKit::default_matrixkit();

	OverlapOrthog xorthog(OverlapOrthog::Symmetric, X_, kit, OverlapOrthog::default_lindep_tol());
	RefSCMatrix mtilde = xorthog.basis_to_orthog_basis();

	// So far, I haven't thought about the reduced dimensional version...
	// Guessing occupied pairs can hardly be linearly dependent with each other.
	if(mtilde.coldim() != X_.dim()) {
	  throw ProgrammingError("Currently occupied pairs are assumed to be linearly independent.", __FILE__, __LINE__);
	}

	RefSCMatrix btilde = mtilde.t() * B_ * mtilde;
	RefSymmSCMatrix btilde_symm(btilde.coldim(), kit);
	btilde_symm.assign_subblock(btilde, 0, btilde.nrow()-1, 0, btilde.ncol()-1);

	RefDiagSCMatrix beig(mtilde.coldim(), kit);
	RefSCMatrix U(B_.dim(), B_.dim(), kit);
	btilde_symm.diagonalize(beig, U);

	eigvals = beig;
	eigvecs = mtilde * U;

	// printing out some test cases if defined.
//#define LOCAL_DEBUG_PREDIAGON
#ifdef LOCAL_DEBUG_PREDIAGON
	RefDiagSCMatrix unit(mtilde.coldim(), mtilde.kit());
	unit.assign(1.0);

	const double coeff = 3.0;
	(B_ - coeff * X_).gi().print();
    (eigvecs * (eigvals - coeff * unit).gi() * eigvecs.t()).print();

    Ref<Tensor> in = new Tensor("testin", mem());
    offset_gt2(in, false);
    Ref<Tensor> out = new Tensor("testout", mem());
    offset_gt2(out, false);
	denom_contraction(in, out);
#endif

}

void CCR12_Info::denom_contraction(const Ref<Tensor>& in, Ref<Tensor>& out) {

#define USE_MATRIX_DECOMP
#ifdef USE_MATRIX_DECOMP

  // requires lmatrix_.
  assert(lmatrix_.nonnull());
  const int pair_size = lmatrix_.ncol();

  Ref<Tensor> ltensor1 = new Tensor("lternsor", mem());
  // first create a zero-cleared tensor with tile information.
  // pair basis is not treated tile-wise; therefore each block is the size of singles.
  {
    long size = 0L;
    for (long pair = 0; pair != pair_size; ++pair) { // not using blocks...
      for (long h1b = 0L; h1b < noab(); ++h1b) {
        for (long h2b = h1b; h2b < noab(); ++h2b) {
          ltensor1->input_offset(h2b + noab() * (h1b + noab() * pair), size);
          size += get_range(h1b) * get_range(h2b);
        }
      }
    }
    ltensor1->set_filesize(size);
    ltensor1->createfile();
  }

  Ref<Tensor> ltensor2 = ltensor1->clone();
  // then, fill in to this tensor (written by hand).
  // not sure about UHF reference
  {
    const size_t nocc_act = naoa();
    const size_t maxtile = maxtilesize() * maxtilesize();
    double* work = mem()->malloc_local_double(maxtile);
    double* work2 = mem()->malloc_local_double(maxtile);
    for (long pair = 0; pair != pair_size; ++pair) { // not using blocks...
      if (pair % mem()->n() != mem()->me()) continue;

      for (long h1b = 0L; h1b < noab(); ++h1b) {
        for (long h2b = h1b; h2b < noab(); ++h2b) {
          const size_t rh1b = get_range(h1b);
          const size_t rh2b = get_range(h2b);
          fill(work, work+rh1b*rh2b, 0.0);
          fill(work2, work2+rh1b*rh2b, 0.0);
          double* val = work;
          double* val2 = work2;
          for (int h1 = 0; h1 != rh1b; ++h1) {
            for (int h2 = 0; h2 != rh2b; ++h2, ++val, ++val2) {
              const size_t h1tot = h1 + get_offset(get_alpha(h1b));
              const size_t h2tot = h2 + get_offset(get_alpha(h2b));

              const size_t h21s = h2tot + nocc_act * h1tot;
              const size_t h21r = h1tot + nocc_act * h2tot;

              // Retrieving source data; not a efficient code since
              // the stride for read isn't small...
              const double src1 = lmatrix_(h21s, pair);
              const double src2 = lmatrix_(h21r, pair);

              // Anti-symmetrized with respect to h1 and h2.
              if (get_spin(h1b) == get_spin(h2b)) {
                *val = src1 - src2;
              } else {
                assert(get_spin(h1b) < get_spin(h2b));
                *val = src1;
              }
              *val2 = src1;
            }
          }
          ltensor1->put_block(h2b + noab() * (h1b + noab() * pair), work);
          ltensor2->put_block(h2b + noab() * (h1b + noab() * pair), work2);
        }
      }
    }
    mem()->free_local_double(work);
    mem()->free_local_double(work2);
  }

  Ref<Tensor> intermediate_tensor = ltensor1->clone();
  // Then evaluate tensor product L * in.
  // Denominator (e1+e2-eb) is applied also.
  {
    const size_t singles = maxtilesize() * maxtilesize();
    const size_t doubles = singles * singles;
    double* k_a0 = mem()->malloc_local_double(doubles);
    double* k_a1 = mem()->malloc_local_double(singles);
    double* k_c  = mem()->malloc_local_double(singles);

    for (int pair = 0; pair != pair_size; ++pair) {
      if (pair % mem()->n() != mem()->me()) continue;

      for (int h1b = 0; h1b < noab(); ++h1b) {
        for (int h2b = h1b; h2b < noab(); ++h2b) {

          // target size
          const size_t dimc_singles = get_range(h1b)*get_range(h2b);
          std::fill(k_c, k_c+dimc_singles, 0.0);

          for (int h5b = 0; h5b < noab(); ++h5b) {
            for (int h6b = h5b; h6b < noab(); ++h6b) {

              if (get_spin(h1b)+get_spin(h2b) == get_spin(h5b)+get_spin(h6b)) {
                if ((get_sym(h1b)^(get_sym(h2b)^(get_sym(h5b)^get_sym(h6b)))) == irrep_t()) {

                  long h5b_1, h6b_1, h1b_1, h2b_1;
                  restricted_4(h5b, h6b, h1b, h2b, h5b_1, h6b_1, h1b_1, h2b_1);
                  const int dim_common = get_range(h5b) * get_range(h6b);
                  const int dima1_sort = get_range(h1b) * get_range(h2b);

                  // read tilde from in
                  ltensor1->get_block(h6b+noab()*(h5b+noab()*pair), k_a1);

                  // read denominator
                  in->get_block(h2b_1+noab()*(h1b_1+noab()*(h6b_1+noab()*h5b_1)), k_a0);

                  const double factor = h5b == h6b ? 0.5 : 1.0;
                  const int unit = 1;
                  const double one = 1.0;

                  F77_DGEMV("n", &dima1_sort, &dim_common,
                          &factor, k_a0, &dima1_sort,
                          k_a1, &unit,
                          &one, k_c, &unit);
                }
              }
            }
          }
          // divide by the denominator
          const double diag_element = bdiag_->get_element(pair);
          int iall = 0;
          for (int h1 = 0; h1 != get_range(h1b); ++h1) {
            for (int h2 = 0; h2 != get_range(h2b); ++h2, ++iall) {
              const double eh1 = get_orb_energy(get_offset(h1b) + h1);
              const double eh2 = get_orb_energy(get_offset(h2b) + h2);
              k_c[iall] /= eh1 + eh2 - diag_element;
            }
          }
          intermediate_tensor->put_block(h2b+noab()*(h1b+noab()*pair), k_c);
        }
      }
    }
    mem()->free_local_double(k_a1);
    mem()->free_local_double(k_a0);
    mem()->free_local_double(k_c);
    mem()->sync();
  }

  // Finally creates the "out" tensor.
  {
    const size_t singles = maxtilesize() * maxtilesize();
    const size_t doubles = singles * singles;
    double* k_a0 = mem()->malloc_local_double(singles);
    double* k_a1 = mem()->malloc_local_double(singles);
    double* k_c  = mem()->malloc_local_double(doubles);
    for (int h1b = 0; h1b < noab(); ++h1b) {
      for (int h2b = h1b; h2b < noab(); ++h2b) {
        for (int h5b = 0; h5b < noab(); ++h5b) {
          for (int h6b = h5b; h6b < noab(); ++h6b) {
            const int rh1b = get_range(h1b);
            const int rh2b = get_range(h2b);
            const int rh5b = get_range(h5b);
            const int rh6b = get_range(h6b);

            const size_t tileoffset = h6b + noab() * (h5b + noab() * (h2b + noab() * h1b));
            if (!out->is_this_local(tileoffset)) continue;

            if (!restricted() || get_spin(h5b)+get_spin(h6b)+get_spin(h1b)+get_spin(h2b) != 8L) {
              if (get_spin(h5b)+get_spin(h6b) == get_spin(h1b)+get_spin(h2b)) {
                if ((get_sym(h5b)^(get_sym(h6b)^(get_sym(h1b)^get_sym(h2b)))) == irrep_t()) {
                  fill(k_c, k_c+rh1b*rh2b*rh5b*rh6b, 0.0);

                  // TODO we can perform this contraction with blocking of pairs, of course.
                  for (int pair = 0; pair != pair_size; ++pair) {
                    intermediate_tensor->get_block(h6b+noab()*(h5b+noab()*pair), k_a0);
                    ltensor2->get_block(h2b+noab()*(h1b+noab()*pair), k_a1);
                    const int r56 = rh5b * rh6b;
                    const int r12 = rh1b * rh2b;
                    const int unit = 1;
                    const double one = 1.0;

                    F77_DGEMM("n", "t", &r56, &r12, &unit, &one, k_a0, &r56, k_a1, &r12, &one, k_c, &r56);
                  }
                  out->add_block(tileoffset, k_c);
                }
              }
            }
          }
        }
      }
    }
    mem()->free_local_double(k_a1);
    mem()->free_local_double(k_a0);
    mem()->free_local_double(k_c);
    mem()->sync();
  }

#else

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

  RefDiagSCMatrix unit(lmatrix_.coldim(), lmatrix_.kit());
  unit.assign(1.0);

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

          // this does not help removing O(n^8) steps, as matrix-matrix multiplication
          // already involves that scaling.
#if 0
          RefSymmSCMatrix refxminusb = X() * (eh3 + eh4) - B();
          RefSymmSCMatrix refinverse = refxminusb.gi();
#else
          RefDiagSCMatrix kernel = (unit * (eh3 + eh4) - bdiag_).gi();
          RefSCMatrix geninverse = lmatrix_ * kernel * lmatrix_.t();
          RefSymmSCMatrix refinverse(X_.dim(), X_.kit());
          refinverse.assign_subblock(geninverse, 0, geninverse.nrow()-1, 0, geninverse.ncol()-1);
#endif
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
#endif

}

