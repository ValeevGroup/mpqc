//
// triples_singles.cc
//
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@theochem.uni-stuttgart.de>
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
#include <chemistry/qc/ccr12/ccr12_triples.h> 
#include <chemistry/qc/ccr12/tensor.h>
using namespace sc;


void CCR12_Triples::singles_ig(Ref<Tensor>& out){

  const size_t maxtile = z->maxtilesize();
  const size_t singles = maxtile * maxtile;
  const size_t doubles = singles * singles;
  const size_t triples = singles * doubles;
  double* k_c      = z->mem()->malloc_local_double(triples);
  double* k_c_sort = z->mem()->malloc_local_double(triples);
  double* k_a0     = z->mem()->malloc_local_double(singles);
  double* k_a0_sort= z->mem()->malloc_local_double(singles);
  double* k_a1     = z->mem()->malloc_local_double(doubles);
  double* k_a1_sort= z->mem()->malloc_local_double(doubles);
  
  const long noab = z->noab();
  const long nvab = z->nvab();


  for (long h4b = 0L; h4b < noab; ++h4b) {
    for (long x5b = h4b; x5b < noab+nvab; ++x5b) {
      for (long p6b = noab; p6b < noab+nvab; ++p6b) {
        for (long h1b = 0L; h1b < noab; ++h1b) {
          for (long h2b = 0L; h2b < noab; ++h2b) {
            for (long h3b = h2b; h3b < noab; ++h3b) {
              long tileoffset;
              if (h1b < h2b) {
                tileoffset = h3b+noab*(h2b+noab*(h1b+noab*(p6b-noab+nvab*(x5b+(noab+nvab)*h4b))));
              } else if (h2b <= h1b && h1b < h3b) {
                tileoffset = h3b+noab*(h1b+noab*(h2b+noab*(p6b-noab+nvab*(x5b+(noab+nvab)*h4b))));
              } else if (h3b <= h1b) {
                tileoffset = h1b+noab*(h3b+noab*(h2b+noab*(p6b-noab+nvab*(x5b+(noab+nvab)*h4b))));
              }
              if (out->is_this_local(tileoffset)) {
                if (!z->restricted() || z->get_spin(h4b)+z->get_spin(x5b)+z->get_spin(p6b)
                                       +z->get_spin(h1b)+z->get_spin(h2b)+z->get_spin(h3b) != 12L) {
                  if (z->get_spin(h4b)+z->get_spin(x5b)+z->get_spin(p6b) ==
                      z->get_spin(h1b)+z->get_spin(h2b)+z->get_spin(h3b)) {
                    if ((z->get_sym(h4b)^(z->get_sym(x5b)^(z->get_sym(p6b)^
                        (z->get_sym(h1b)^(z->get_sym(h2b)^z->get_sym(h3b)))))) == (z->irrep_t()^z->irrep_e())) {

                      const size_t dimc = z->get_range(h4b)*z->get_range(x5b)*z->get_range(p6b)
                                         *z->get_range(h1b)*z->get_range(h2b)*z->get_range(h3b);
                      std::fill(k_c_sort, k_c_sort+dimc, 0.0);

                      if (z->get_spin(p6b) == z->get_spin(h1b)) {
                        if ((z->get_sym(p6b)^z->get_sym(h1b)) == z->irrep_t()) {

                          long p6b_0,h1b_0;
                          z->restricted_2(p6b,h1b,p6b_0,h1b_0);
                          long h4b_1,x5b_1,h2b_1,h3b_1;
                          z->restricted_4(h4b,x5b,h2b,h3b,h4b_1,x5b_1,h2b_1,h3b_1);

                          const long dim_common=1L;
                          const long dima0_sort=z->get_range(p6b)*z->get_range(h1b);
                          const long dima1_sort=z->get_range(h4b)*z->get_range(x5b)*z->get_range(h2b)*z->get_range(h3b);
                          z->t1()->get_block(h1b_0+noab*(p6b_0-noab), k_a0);
                          z->sort_indices2(k_a0, k_a0_sort,
                              z->get_range(p6b),z->get_range(h1b),
                              1,0, +1.0);
                          z->vd2_gen()->get_block(h3b_1+z->nab()*(h2b_1+z->nab()*(x5b_1+(noab+nvab)*h4b_1)),k_a1);
                          z->sort_indices4(k_a1,k_a1_sort,
                              z->get_range(h4b),z->get_range(x5b),
                              z->get_range(h2b),z->get_range(h3b),
                              3,2,1,0, +1.0);

                          double factor=1.0;
                          z->smith_dgemm(dima0_sort, dima1_sort, dim_common,
                              factor, k_a0_sort, dim_common,
                              k_a1_sort, dim_common,
                              1.0, k_c_sort, dima0_sort);
                        }
                      }
                      if (h2b >= h1b) {
                        z->sort_indices6(k_c_sort,k_c,
                            z->get_range(h3b),z->get_range(h2b),z->get_range(x5b),
                            z->get_range(h4b),z->get_range(h1b),z->get_range(p6b),
                            3,2,5,4,1,0, +1.0);
                        out->add_block(h3b+noab*(h2b+noab*(h1b+noab*(p6b-noab+nvab*(x5b+(noab+nvab)*h4b)))), k_c);
                      }
                      if (h3b >= h1b && h1b >= h2b) {
                        z->sort_indices6(k_c_sort,k_c,
                            z->get_range(h3b),z->get_range(h2b),z->get_range(x5b),
                            z->get_range(h4b),z->get_range(h1b),z->get_range(p6b),
                            3,2,5,1,4,0, -1.0);
                        out->add_block(h3b+noab*(h1b+noab*(h2b+noab*(p6b-noab+nvab*(x5b+(noab+nvab)*h4b)))), k_c);
                      }
                      if (h1b >= h3b) {
                        z->sort_indices6(k_c_sort,k_c,
                            z->get_range(h3b),z->get_range(h2b),z->get_range(x5b),
                            z->get_range(h4b),z->get_range(h1b),z->get_range(p6b),
                            3,2,5,1,0,4, +1.0);
                        out->add_block(h1b+noab*(h3b+noab*(h2b+noab*(p6b-noab+nvab*(x5b+(noab+nvab)*h4b)))), k_c);
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
  z->mem()->free_local_double(k_a0_sort);
  z->mem()->free_local_double(k_a1);
  z->mem()->free_local_double(k_a1_sort);
  z->mem()->sync();
}


void CCR12_Triples::singles(){ 

  Ref<Tensor> out = singles_intermediate_;

  const size_t maxtile = z->maxtilesize();
  const size_t singles = maxtile * maxtile;
  const size_t doubles = singles * singles;
  const size_t triples = singles * doubles;
  double* k_c      = z->mem()->malloc_local_double(triples);
  double* k_c_sort = z->mem()->malloc_local_double(triples);
  double* k_a0     = z->mem()->malloc_local_double(singles);
  double* k_a0_sort= z->mem()->malloc_local_double(singles);
  double* k_a1     = z->mem()->malloc_local_double(doubles);
  double* k_a1_sort= z->mem()->malloc_local_double(doubles);
  
  const long noab = z->noab();
  const long nvab = z->nvab();


  for (long h4b = 0L; h4b < noab; ++h4b) {
    for (long h5b = h4b; h5b < noab; ++h5b) {
      for (long p6b = noab; p6b < noab+nvab; ++p6b) {
        for (long h1b = 0L; h1b < noab; ++h1b) {
          for (long h2b = 0L; h2b < noab; ++h2b) {
            for (long h3b = h2b; h3b < noab; ++h3b) {
              long tileoffset;
              if (h1b < h2b) {
                tileoffset = h3b+noab*(h2b+noab*(h1b+noab*(p6b-noab+nvab*(h5b+noab*h4b))));
              } else if (h2b <= h1b && h1b < h3b) {
                tileoffset = h3b+noab*(h1b+noab*(h2b+noab*(p6b-noab+nvab*(h5b+noab*h4b))));
              } else if (h3b <= h1b) {
                tileoffset = h1b+noab*(h3b+noab*(h2b+noab*(p6b-noab+nvab*(h5b+noab*h4b))));
              }
              if (out->is_this_local(tileoffset)) {
                if (!z->restricted() || z->get_spin(h4b)+z->get_spin(h5b)+z->get_spin(p6b)
                                       +z->get_spin(h1b)+z->get_spin(h2b)+z->get_spin(h3b) != 12L) {
                  if (z->get_spin(h4b)+z->get_spin(h5b)+z->get_spin(p6b) ==
                      z->get_spin(h1b)+z->get_spin(h2b)+z->get_spin(h3b)) {
                    if ((z->get_sym(h4b)^(z->get_sym(h5b)^(z->get_sym(p6b)^
                        (z->get_sym(h1b)^(z->get_sym(h2b)^z->get_sym(h3b)))))) == (z->irrep_t()^z->irrep_e())) {

                      const size_t dimc = z->get_range(h4b)*z->get_range(h5b)*z->get_range(p6b)
                                         *z->get_range(h1b)*z->get_range(h2b)*z->get_range(h3b);
                      std::fill(k_c_sort, k_c_sort+dimc, 0.0);

                      if (z->get_spin(p6b) == z->get_spin(h1b)) {
                        if ((z->get_sym(p6b)^z->get_sym(h1b)) == z->irrep_t()) {

                          long p6b_0,h1b_0;
                          z->restricted_2(p6b,h1b,p6b_0,h1b_0);
                          long h4b_1,h5b_1,h2b_1,h3b_1;
                          z->restricted_4(h4b,h5b,h2b,h3b,h4b_1,h5b_1,h2b_1,h3b_1);

                          const long dim_common=1L;
                          const long dima0_sort=z->get_range(p6b)*z->get_range(h1b);
                          const long dima1_sort=z->get_range(h4b)*z->get_range(h5b)*z->get_range(h2b)*z->get_range(h3b);
                          z->t1()->get_block(h1b_0+noab*(p6b_0-noab), k_a0);
                          z->sort_indices2(k_a0, k_a0_sort,
                              z->get_range(p6b),z->get_range(h1b),
                              1,0, +1.0);
                          z->vd2()->get_block(h3b_1+z->nab()*(h2b_1+z->nab()*(h5b_1+noab*h4b_1)),k_a1);
                          z->sort_indices4(k_a1,k_a1_sort,
                              z->get_range(h4b),z->get_range(h5b),
                              z->get_range(h2b),z->get_range(h3b),
                              3,2,1,0, +1.0);

                          double factor=1.0;
                          z->smith_dgemm(dima0_sort, dima1_sort, dim_common,
                              factor, k_a0_sort, dim_common,
                              k_a1_sort, dim_common,
                              1.0, k_c_sort, dima0_sort);
                        }
                      }
                      if (h2b >= h1b) {
                        z->sort_indices6(k_c_sort,k_c,
                            z->get_range(h3b),z->get_range(h2b),z->get_range(h5b),
                            z->get_range(h4b),z->get_range(h1b),z->get_range(p6b),
                            3,2,5,4,1,0, +1.0);
                        out->add_block(h3b+noab*(h2b+noab*(h1b+noab*(p6b-noab+nvab*(h5b+noab*h4b)))), k_c);
                      }
                      if (h3b >= h1b && h1b >= h2b) {
                        z->sort_indices6(k_c_sort,k_c,
                            z->get_range(h3b),z->get_range(h2b),z->get_range(h5b),
                            z->get_range(h4b),z->get_range(h1b),z->get_range(p6b),
                            3,2,5,1,4,0, -1.0);
                        out->add_block(h3b+noab*(h1b+noab*(h2b+noab*(p6b-noab+nvab*(h5b+noab*h4b)))), k_c);
                      }
                      if (h1b >= h3b) {
                        z->sort_indices6(k_c_sort,k_c,
                            z->get_range(h3b),z->get_range(h2b),z->get_range(h5b),
                            z->get_range(h4b),z->get_range(h1b),z->get_range(p6b),
                            3,2,5,1,0,4, +1.0);
                        out->add_block(h1b+noab*(h3b+noab*(h2b+noab*(p6b-noab+nvab*(h5b+noab*h4b)))), k_c);
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
  z->mem()->free_local_double(k_a0_sort);
  z->mem()->free_local_double(k_a1);
  z->mem()->free_local_double(k_a1_sort);
  z->mem()->sync();
} 

