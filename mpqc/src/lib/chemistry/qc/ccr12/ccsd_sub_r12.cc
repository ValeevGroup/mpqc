//
// ccsd_sub_r12.cc
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


#include <algorithm>
#include <chemistry/qc/ccr12/ccsd_sub_r12.h>
#include <chemistry/qc/mbptr12/blas.h>
#include <chemistry/qc/ccr12/tensor.h>
#include <chemistry/qc/ccr12/mtensor.h>

using namespace sc;
using namespace std;
  
  
void CCSD_Sub_R12::denom_contraction(){ 

  const size_t singles = z->maxtilesize() * z->maxtilesize();
  const size_t doubles = singles * singles;
  double* k_a0      = z->mem()->malloc_local_double(doubles); 
  double* k_a1      = z->mem()->malloc_local_double(doubles); 
  double* k_c       = z->mem()->malloc_local_double(doubles); 
  double* k_c_sort  = z->mem()->malloc_local_double(singles); 

  const int nocc_act = z->naoa();
  MTensor<4>::tile_ranges iiii(4, MTensor<4>::tile_range(0, z->noab()));
  MTensor<4>::element_ranges iiii_erange(4, MTensor<4>::element_range(0, nocc_act) );
  vector<long> amap;
  {
    vector<int> intmap = sc::map(*(z->r12evalinfo()->refinfo()->occ_act_sb(Alpha)), *z->corr_space(), false);
    amap.resize(intmap.size());
    std::copy(intmap.begin(), intmap.end(), amap.begin());
  }

  // this loop structure minimizes the O(o^8) operation...
  for (long h3b = 0L; h3b < z->noab(); ++h3b) { 
    for (long h4b = h3b; h4b < z->noab(); ++h4b) { 
      int count = 0;
      for (int h3 = 0; h3 < z->get_range(h3b); ++h3) { 
        for (int h4 = 0; h4 < z->get_range(h4b); ++h4, ++count) { 
          if (count % z->mem()->n() != z->mem()->me() ) continue; 

          // orbital energies
          const double eh3 = z->get_orb_energy(z->get_offset(h3b) + h3);
          const double eh4 = z->get_orb_energy(z->get_offset(h4b) + h4);

          // current denominator tensor
          Ref<Tensor> denom = new Tensor("denom", z->mem());
          z->offset_gt2(denom, false);
          MTensor<4> D(z, denom.pointer(), iiii);

          RefSymmSCMatrix refxminusb = z->X() * (eh3 + eh4) - z->B();
          RefSymmSCMatrix refinverse = refxminusb.gi();

          D.convert(refinverse, nocc_act, nocc_act, false, false,
                    amap, amap, amap, amap, &iiii_erange);

          for (long h1b = 0L; h1b < z->noab(); ++h1b) {
            for (long h2b = h1b; h2b < z->noab(); ++h2b) {

              if (!z->restricted() || z->get_spin(h3b)+z->get_spin(h4b)+z->get_spin(h1b)+z->get_spin(h2b) != 8L) { 
                if (z->get_spin(h3b)+z->get_spin(h4b) == z->get_spin(h1b)+z->get_spin(h2b)) { 
                  if ((z->get_sym(h3b)^(z->get_sym(h4b)^(z->get_sym(h1b)^z->get_sym(h2b)))) == z->irrep_t()) { 

                    // target size
                    const long dimc_singles = z->get_range(h1b)*z->get_range(h2b);
                    std::fill(k_c_sort, k_c_sort+dimc_singles, 0.0);

                    for (long h5b = 0L; h5b < z->noab(); ++h5b) {
                      for (long h6b = h5b; h6b < z->noab(); ++h6b) {

                        if (z->get_spin(h3b)+z->get_spin(h4b) == z->get_spin(h5b)+z->get_spin(h6b)) {
                          if ((z->get_sym(h3b)^(z->get_sym(h4b)^(z->get_sym(h5b)^z->get_sym(h6b)))) == z->irrep_e()) {

                            long h3b_0, h4b_0, h5b_0, h6b_0; 
                            z->restricted_4(h3b, h4b, h5b, h6b, h3b_0, h4b_0, h5b_0, h6b_0); 
                            long h5b_1, h6b_1, h1b_1, h2b_1; 
                            z->restricted_4(h5b, h6b, h1b, h2b, h5b_1, h6b_1, h1b_1, h2b_1); 
                            const int dim_common = z->get_range(h5b) * z->get_range(h6b); 
                            const int dima0_sort = z->get_range(h3b) * z->get_range(h4b); 
                            const int dima1_sort = z->get_range(h1b) * z->get_range(h2b); 

                            // read tilde V (redundant read is involved...)
                            tildeV_->get_block(h4b_0+z->noab()*(h3b_0+z->noab()*(h6b_0+z->noab()*(h5b_0))), k_a0);

                            // read denominator
                            denom->get_block(h6b_1+z->noab()*(h5b_1+z->noab()*(h2b_1+z->noab()*(h1b_1))), k_a1);
                            double factor = 1.0;
                            if (h5b == h6b) factor *= 0.5;
                            const int unit = 1;
                            const int h34 = h4 + z->get_range(h4b) * h3;
                            const int stride = z->get_range(h3b) * z->get_range(h4b);
                            const double one = 1.0;
                            F77_DGEMV("t", &dim_common, &dima1_sort, &factor, k_a1, &dim_common, k_a0 + h34, &stride, &one, k_c_sort, &unit);
                          }
                        } 
                      }
                    }
                    const long dimc = z->get_range(h3b)*z->get_range(h4b)*z->get_range(h1b)*z->get_range(h2b);
                    std::fill(k_c, k_c+dimc, 0.0);
                    {
                      const size_t h34 = h4 + z->get_range(h4b) * h3;
                      const size_t stride = z->get_range(h3b) * z->get_range(h4b);
                      size_t iall = 0L;
                      for (int h1 = 0; h1 != z->get_range(h1b); ++h1) {
                        for (int h2 = 0; h2 != z->get_range(h2b); ++h2, ++iall) {
                          k_c[iall * stride + h34] = k_c_sort[iall];
                        }
                      }
                    }
                    intermediate_->add_block(h4b+z->noab()*(h3b+z->noab()*(h2b+z->noab()*(h1b))),k_c); 
                  }
                }
              }
            } 
          } 
        }
      } // orbital loops
    } 
  } 
  z->mem()->free_local_double(k_a1); 
  z->mem()->free_local_double(k_a0); 
  z->mem()->free_local_double(k_c_sort); 
  z->mem()->free_local_double(k_c); 
  z->mem()->sync(); 
} 

