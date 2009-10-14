//
// triples_utils.cc
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
#include <chemistry/qc/ccr12/ccr12_triples.h> 
#include <chemistry/qc/mbptr12/blas.h>
#include <chemistry/qc/ccr12/tensor.h>
#include <chemistry/qc/ccr12/mtensor.h>

using namespace sc;
using namespace std;


void CCR12_Triples::prediagon() {
    // CCR12_Triples::B_ and CCR12_Triples::X_ required.
    // B_ and X_ are RefSymmSCMatrix objects.
    assert(B_.nonnull() && X_.nonnull());

    ExEnv::out0() << setprecision(10);
    OverlapOrthog xorthog(OverlapOrthog::Symmetric, X_, X_.kit(), OverlapOrthog::default_lindep_tol());
    RefSCMatrix mtilde = xorthog.basis_to_orthog_basis();

    // So far, I haven't thought about the reduced dimensional version...
    // Guessing occupied pairs can hardly be linearly dependent with each other.
    if(mtilde.coldim() != X_.dim()) {
      throw ProgrammingError("Currently occupied pairs are assumed to be linearly independent.", __FILE__, __LINE__);
    }

    RefSCMatrix btilde = mtilde.t() * B_ * mtilde;
    RefSymmSCMatrix btilde_symm(btilde.coldim(), btilde.kit());
    btilde_symm.assign_subblock(btilde, 0, btilde.nrow()-1, 0, btilde.ncol()-1);

    RefDiagSCMatrix beig(mtilde.coldim(), mtilde.kit());
    RefSCMatrix U(B_.dim(), B_.dim(), B_.kit());
    btilde_symm.diagonalize(beig, U);

    bdiag_ = beig;
    lmatrix_ = mtilde * U;
}


void CCR12_Triples::fill_in_ltensors() {
  const long noab = z->noab();
  {
    ltensor1_ = new Tensor("ltensor1", z->mem());
    long size = 0L;

    for (long pair = 0; pair != pair_size_; ++pair) { // not using blocks...
      for (long h1b = 0L; h1b < noab; ++h1b) {
        for (long h2b = h1b; h2b < noab; ++h2b) {
          ltensor1_->input_offset(h2b + noab * (h1b + noab * pair), size);
          size += z->get_range(h1b) * z->get_range(h2b);
        }
      }
    }
    ltensor1_->set_filesize(size);
    ltensor1_->createfile();
    ltensor2_ = ltensor1_->clone();
  }

  {
    const size_t nocc_act = z->naoa();
    const size_t maxtile = z->maxtilesize() * z->maxtilesize();
    double* work = z->mem()->malloc_local_double(maxtile);
    double* work2 = z->mem()->malloc_local_double(maxtile);
    for (long pair = 0; pair != pair_size_; ++pair) { // not using blocks...
      if (pair % z->mem()->n() != z->mem()->me()) continue;

      for (long h1b = 0L; h1b < noab; ++h1b) {
        for (long h2b = h1b; h2b < noab; ++h2b) {
          const size_t rh1b = z->get_range(h1b);
          const size_t rh2b = z->get_range(h2b);
          fill(work, work+rh1b*rh2b, 0.0);
          fill(work2, work2+rh1b*rh2b, 0.0);
          double* val = work;
          double* val2 = work2;
          for (int h1 = 0; h1 != rh1b; ++h1) {
            for (int h2 = 0; h2 != rh2b; ++h2, ++val, ++val2) {
              const size_t h1tot = h1 + z->get_offset(z->get_alpha(h1b));
              const size_t h2tot = h2 + z->get_offset(z->get_alpha(h2b));

              const size_t h21s = h2tot + nocc_act * h1tot;
              const size_t h21r = h1tot + nocc_act * h2tot;

              // Retrieving source data; not a efficient code since
              // the stride for read isn't small...
              const double src1 = lmatrix_(h21s, pair);
              const double src2 = lmatrix_(h21r, pair);

              // Anti-symmetrized with respect to h1 and h2.
              if (z->get_spin(h1b) == z->get_spin(h2b)) {
                *val = src1 - src2;
              } else {
                assert(z->get_spin(h1b) < z->get_spin(h2b));
                *val = src1;
              }
              *val2 = src1;
            }
          }
          ltensor1_->put_block(h2b + noab * (h1b + noab * pair), work);
          ltensor2_->put_block(h2b + noab * (h1b + noab * pair), work2);
        }
      }
    }
    z->mem()->free_local_double(work);
    z->mem()->free_local_double(work2);
  }
}


void CCR12_Triples::offset_hhphhh(Ref<Tensor>& t) {
  long size=0L;
  for (long h4b=0L;h4b<z->noab();++h4b) { 
   for (long h5b=h4b;h5b<z->noab();++h5b) { 
    for (long p6b=z->noab();p6b<z->noab()+z->nvab();++p6b) { 
     for (long h1b=0L;h1b<z->noab();++h1b) { 
      for (long h2b=h1b;h2b<z->noab();++h2b) { 
       for (long h3b=h2b;h3b<z->noab();++h3b) { 
        if (!z->restricted() ||
        		z->get_spin(h4b)+z->get_spin(h5b)+z->get_spin(p6b)
        		+z->get_spin(h1b)+z->get_spin(h2b)+z->get_spin(h3b)!=12L) {
         if (z->get_spin(h4b)+z->get_spin(h5b)+z->get_spin(p6b)
        		 ==z->get_spin(h1b)+z->get_spin(h2b)+z->get_spin(h3b)) {
          if ((z->get_sym(h4b)^(z->get_sym(h5b)^(z->get_sym(p6b)^
        		  (z->get_sym(h1b)^(z->get_sym(h2b)^z->get_sym(h3b))))))==z->irrep_t()) {
           t->input_offset(h3b+z->noab()*(h2b+z->noab()*(h1b+z->noab()*(p6b-z->noab()+z->nvab()*(h5b+z->noab()*h4b)))),size);
           size+=z->get_range(h1b)*z->get_range(h2b)*z->get_range(h3b)*z->get_range(h4b)*z->get_range(h5b)*z->get_range(p6b);
          }
         }
        }
       }
      }
     }
    }
   }
  }
  t->set_filesize(size);
  t->createfile();
}


void CCR12_Triples::offset_bphhh(Ref<Tensor>& t) {
  long size=0L;
  const long noab = z->noab();
  const long nvab = z->nvab();
  for (long pair = 0L; pair != pair_size_; ++pair) {
    for (long p6b = noab; p6b < noab + nvab; ++p6b) {
      for (long h1b = 0L; h1b < noab; ++h1b) {
        for (long h2b = h1b; h2b < noab; ++h2b) {
          for (long h3b = h2b; h3b < noab; ++h3b) {
            t->input_offset(h3b+noab*(h2b+noab*(h1b+noab*(p6b-noab+nvab*pair))),size);
            size += z->get_range(h1b)*z->get_range(h2b)*z->get_range(h3b)*z->get_range(p6b);
          }
        }
      }
    }
  }
  t->set_filesize(size);
  t->createfile();
}


void CCR12_Triples::denom_contraction_new() {

  // ltensor1 and 2 must be filled.
  assert(ltensor1_.nonnull() && ltensor2_.nonnull());

  const long noab = z->noab();
  const long nvab = z->nvab();
  const size_t maxtile = z->maxtilesize();
  const size_t singles = maxtile * maxtile;
  const size_t doubles = singles * singles;
  const size_t triples = singles * doubles;
  Ref<Tensor> l_times_rhs = new Tensor("l_times_rhs", z->mem());

  {
    double* k_c  = z->mem()->malloc_local_double(doubles);
    double* k_a0 = z->mem()->malloc_local_double(triples);
    double* k_a1 = z->mem()->malloc_local_double(singles);

    offset_bphhh(l_times_rhs);

    for (long pair = 0; pair != pair_size_; ++pair) {
      if (pair % z->mem()->n() != z->mem()->me()) continue;

      for (long p6b = noab; p6b < noab+nvab; ++p6b) {
        for (long h1b = 0L; h1b < noab; ++h1b) {
          for (long h2b = h1b; h2b < noab; ++h2b) {
            for (long h3b = h2b; h3b < noab; ++h3b) {
              const size_t rh1b = z->get_range(h1b);
              const size_t rh2b = z->get_range(h2b);
              const size_t rh3b = z->get_range(h3b);
              const size_t rp6b = z->get_range(p6b);

              fill(k_c, k_c+rh1b*rh2b*rh3b*rp6b, 0.0);

              for (long h4b = 0L; h4b < noab; ++h4b) {
                for (long h5b = h4b; h5b < noab; ++h5b) {
                  if (z->get_spin(h1b)+z->get_spin(h2b)+z->get_spin(h3b)
                      == z->get_spin(h4b)+z->get_spin(h5b)+z->get_spin(p6b)) {
                    if ((z->get_sym(h1b)^(z->get_sym(h2b)^(z->get_sym(h3b)^
                        (z->get_sym(h4b)^(z->get_sym(h5b)^z->get_sym(p6b)))))) == z->irrep_t()) {
                      const size_t rh4b = z->get_range(h4b);
                      const size_t rh5b = z->get_range(h5b);

                      long h4b_1, h5b_1, p6b_1, h1b_1, h2b_1, h3b_1;
                      z->restricted_6(h4b, h5b, p6b, h1b, h2b, h3b,
                                      h4b_1, h5b_1, p6b_1, h1b_1, h2b_1, h3b_1);

                      doubles_intermediate_->get_block(h3b_1 + noab * (h2b_1 + noab *
                      (h1b_1 + noab * (p6b_1 - noab + nvab * (h5b_1 + noab * h4b_1)))), k_a0);

                      ltensor1_->get_block(h5b + noab * (h4b + noab * pair), k_a1);

                      const int unit = 1;
                      const double factor = h4b == h5b ? 0.5 : 1.0;
                      const double one = 1.0;
                      const int dim0 = rh1b * rh2b * rh3b * rp6b;
                      const int dim1 = rh4b * rh5b;
                      F77_DGEMV("n", &dim0, &dim1, &factor, k_a0, &dim0, k_a1, &unit, &one, k_c, &unit);
                    }
                  }
                }
              }
              const double diag_element = bdiag_->get_element(pair);
              int iall = 0;
              for (int p6 = 0; p6 != rp6b; ++p6) {
                const double ep6 = z->get_orb_energy(z->get_offset(p6b) + p6);
                for (int h1 = 0; h1 != rh1b; ++h1) {
                  const double eh1 = z->get_orb_energy(z->get_offset(h1b) + h1);
                  for (int h2 = 0; h2 != rh2b; ++h2) {
                    const double eh2 = z->get_orb_energy(z->get_offset(h2b) + h2);
                    for (int h3 = 0; h3 != rh3b; ++h3, ++iall) {
                      const double eh3 = z->get_orb_energy(z->get_offset(h3b) + h3);
                      k_c[iall] /= eh1 + eh2 + eh3 - ep6 - diag_element;
                    }
                  }
                }
              }
              l_times_rhs->put_block(h3b + noab * (h2b + noab * (h1b + noab * (p6b - noab + nvab * pair))), k_c);
            }
          }
        }
      }
    }

    z->mem()->free_local_double(k_c);
    z->mem()->free_local_double(k_a0);
    z->mem()->free_local_double(k_a1);
    z->mem()->sync();
  }

  // Now l_times_rhs is ready. Form rhs_intermediate!

  {
    double* k_c  = z->mem()->malloc_local_double(triples);
    double* k_a0 = z->mem()->malloc_local_double(doubles);
    double* k_a1 = z->mem()->malloc_local_double(singles);

    for (long h4b = 0L; h4b < noab; ++h4b) {
      for (long h5b = h4b; h5b < noab; ++h5b) {
        for (long p6b = noab; p6b < noab+nvab; ++p6b) {
          for (long h1b = 0L; h1b < noab; ++h1b) {
            for (long h2b = h1b; h2b < noab; ++h2b) {
              for (long h3b = h2b; h3b < noab; ++h3b) {
                const size_t tag = h3b + noab * (h2b + noab * (h1b + noab * (p6b - noab + nvab * (h5b + noab * h4b))));
                if (!rhs_intermediate_->is_this_local(tag)) continue;

                if (z->get_spin(h1b)+z->get_spin(h2b)+z->get_spin(h3b)
                    == z->get_spin(h4b)+z->get_spin(h5b)+z->get_spin(p6b)) {
                  if ((z->get_sym(h1b)^(z->get_sym(h2b)^(z->get_sym(h3b)^
                      (z->get_sym(h4b)^(z->get_sym(h5b)^z->get_sym(p6b)))))) == z->irrep_t()) {
                    if (z->get_spin(h1b)+z->get_spin(h2b)+z->get_spin(h3b)
                        +z->get_spin(h4b)+z->get_spin(h5b)+z->get_spin(p6b) != 12L) {

                      const int rh4b = z->get_range(h4b);
                      const int rh5b = z->get_range(h5b);
                      const int rp6b = z->get_range(p6b);
                      const int rh1b = z->get_range(h1b);
                      const int rh2b = z->get_range(h2b);
                      const int rh3b = z->get_range(h3b);

                      fill(k_c, k_c+rh4b*rh5b*rp6b*rh1b*rh2b*rh3b, 0.0);

                      const int unit = 1;
                      const int dim0 = rh1b * rh2b * rh3b * rp6b;
                      const int dim1 = rh4b * rh5b;
                      const double one = 1.0;

                      for (int pair = 0; pair != pair_size_; ++pair) {
                        l_times_rhs->get_block(h3b+noab*(h2b+noab*(h1b+noab*(p6b-noab+nvab*pair))), k_a0);
                        ltensor2_->get_block(h5b+noab*(h4b+noab*pair), k_a1);
                        F77_DGEMM("n", "t", &dim0, &dim1, &unit, &one, k_a0, &dim0, k_a1, &dim1, &one, k_c, &dim0);
                      }

                      rhs_intermediate_->put_block(tag, k_c);

                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    z->mem()->free_local_double(k_c);
    z->mem()->free_local_double(k_a0);
    z->mem()->free_local_double(k_a1);
    z->mem()->sync();
  }


}


void CCR12_Triples::denom_contraction(){ 

  const size_t maxtile = z->maxtilesize();
  const size_t singles = maxtile * maxtile;
  const size_t doubles = singles * singles;
  const size_t triples = singles * doubles;  
  double* k_c      = z->mem()->malloc_local_double(triples);
  double* k_c_sort = z->mem()->malloc_local_double(singles);
  double* k_a0     = z->mem()->malloc_local_double(doubles);
  double* k_a1     = z->mem()->malloc_local_double(triples);

  const long noab = z->noab();
  const long nvab = z->nvab();

  const int nocc_act = z->naoa();
  MTensor<4>::tile_ranges iiii(4, MTensor<4>::tile_range(0, noab));
  MTensor<4>::element_ranges iiii_erange(4, MTensor<4>::element_range(0, nocc_act) );
  vector<long> amap;
  {
    vector<int> intmap = sc::map(*(z->r12eval()->occ_act(Alpha)), *z->corr_space(), false);
    amap.resize(intmap.size());
    std::copy(intmap.begin(), intmap.end(), amap.begin());
  }

  // the loops are reordered to minimize the inversion...
  int count = 0;
  for (long p6b = noab; p6b < noab+nvab; ++p6b) {
   for (long h1b = 0L; h1b < noab; ++h1b) {
    for (long h2b = h1b; h2b < noab; ++h2b) { 
     for (long h3b = h2b; h3b < noab; ++h3b) {
      const int rh1b = z->get_range(h1b);
      const int rh2b = z->get_range(h2b);
      const int rh3b = z->get_range(h3b);
      const int rp6b = z->get_range(p6b);
      for (int p6 = 0; p6 < rp6b; ++p6) {
       for (int h1 = 0; h1 < rh1b; ++h1) {
        for (int h2 = 0; h2 < rh2b; ++h2) {
         for (int h3 = 0; h3 < rh3b; ++h3, ++count) {
          if (z->get_offset(h1b) + h1 == z->get_offset(h2b) + h2) continue;
          if (z->get_offset(h1b) + h1 == z->get_offset(h3b) + h3) continue;
          if (z->get_offset(h3b) + h3 == z->get_offset(h2b) + h2) continue;

          if (count % z->mem()->n() != z->mem()->me() ) continue; 

          const double ep6 = z->get_orb_energy(z->get_offset(p6b) + p6);
          const double eh1 = z->get_orb_energy(z->get_offset(h1b) + h1);
          const double eh2 = z->get_orb_energy(z->get_offset(h2b) + h2);
          const double eh3 = z->get_orb_energy(z->get_offset(h3b) + h3);

          // current denominator tensor
          stringstream ss;
          ss << "denom" << count;
          Ref<Tensor> denom = new Tensor(ss.str(), z->mem());
          z->offset_gt2(denom, false);
          MTensor<4> D(z, denom.pointer(), iiii);

          RefSymmSCMatrix refxminusb = X_ * (eh1 + eh2 + eh3 - ep6) - B_;
          RefSymmSCMatrix refinverse = refxminusb.gi();

          D.convert(refinverse, nocc_act, nocc_act, false, false,
                    amap, amap, amap, amap, &iiii_erange);

          const size_t h3216 = h3 + rh3b * (h2 + rh2b * (h1 + rh1b * p6));
          const int stride = rh3b * rh2b * rh1b * rp6b;

          for (long h4b = 0L; h4b < noab; ++h4b) { 
           for (long h5b = h4b; h5b < noab; ++h5b) {

        	const int rh4b = z->get_range(h4b);
        	const int rh5b = z->get_range(h5b);

            // spin & spatial symmetries... 
            if (!z->restricted() || 
                z->get_spin(h4b)+z->get_spin(h5b)+z->get_spin(p6b)+z->get_spin(h1b)+z->get_spin(h2b)+z->get_spin(h3b) != 12L) {
             if (z->get_spin(h4b)+z->get_spin(h5b)+z->get_spin(p6b) == z->get_spin(h1b)+z->get_spin(h2b)+z->get_spin(h3b)) {
              if ((z->get_sym(h4b)^(z->get_sym(h5b)^(z->get_sym(p6b)^(z->get_sym(h1b)^(z->get_sym(h2b)^z->get_sym(h3b))))))
                  == z->irrep_t()) {

               std::fill(k_c_sort, k_c_sort+rh4b*rh5b, 0.0);

               // summation loops
               for (long h7b = 0L; h7b < noab; ++h7b) { 
                for (long h8b = h7b; h8b < noab; ++h8b) {
                 const int rh7b = z->get_range(h7b);
                 const int rh8b = z->get_range(h8b);

                 if (z->get_spin(h4b)+z->get_spin(h5b) == z->get_spin(h7b)+z->get_spin(h8b)) { 
                  if ((z->get_sym(h4b)^(z->get_sym(h5b)^(z->get_sym(h7b)^z->get_sym(h8b)))) == z->irrep_t()) { 

                   long h4b_0, h5b_0, h7b_0, h8b_0; 
                   z->restricted_4(h4b, h5b, h7b, h8b,
                                   h4b_0, h5b_0, h7b_0, h8b_0); 
  
                   long h7b_1, h8b_1, p6b_1, h1b_1, h2b_1, h3b_1; 
                   z->restricted_6(h7b, h8b, p6b, h1b, h2b, h3b,
                                   h7b_1, h8b_1, p6b_1, h1b_1, h2b_1, h3b_1); 

                   const int dim_common = rh7b * rh8b;
                   const int dima0_sort = rh4b * rh5b;
  
                   // read the "denominator" tensor
                   denom->get_block(h8b_0+noab*(h7b_0+noab*(h5b_0+noab*h4b_0)), k_a0);

                   // read the "doubles" intermediate (right-hand-side numerator)
                   doubles_intermediate_->get_block(h3b_1+noab*(h2b_1+noab*(h1b_1+noab
                                                  *(p6b_1-noab+nvab*(h8b_1+noab*h7b_1)))), k_a1);

                   double factor = 1.0; 
                   if (h7b == h8b) factor *= 0.5; 
                   const int unit = 1;
                   const double one = 1.0;
                   F77_DGEMV("t", &dim_common, &dima0_sort,
                		     &factor, k_a0, &dim_common,
                             k_a1 + h3216, &stride,
                             &one, k_c_sort, &unit);
                  }
                 }
                }
               }
               const size_t dimc = rh4b*rh5b*rp6b*rh1b*rh2b*rh3b;
               std::fill(k_c, k_c+dimc, 0.0);
               {
                size_t iall = h3216;
                size_t i = 0;
                for (int h4 = 0; h4 != rh4b; ++h4) {
                 for (int h5 = 0; h5 != rh5b; ++h5, iall += stride, ++i) {
                   k_c[iall] = k_c_sort[i];
                 }
                }
               }
               rhs_intermediate_->add_block(h3b+noab*(h2b+noab*(h1b+noab
            		                  *(p6b-noab+nvab*(h5b+noab*h4b)))), k_c);
              }
             }
            }
           }
          }
         }
        }
       }
      }
     }
    }
   }
  }
  z->mem()->free_local_double(k_c); 
  z->mem()->free_local_double(k_c_sort); 
  z->mem()->free_local_double(k_a0); 
  z->mem()->free_local_double(k_a1); 
  z->mem()->sync(); 
} 

