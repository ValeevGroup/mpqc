//
// tform.timpl.h
//
// Copyright (C) 2011 Edward Valeev
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

#ifndef _mpqc_src_lib_chemistry_qc_libint2_tformtimpl_h
#define _mpqc_src_lib_chemistry_qc_libint2_tformtimpl_h

#include <chemistry/qc/libint2/int1e.h>
#include <chemistry/qc/libint2/int2e.h>
#include <chemistry/qc/libint2/gto.h>

namespace sc {

  template <unsigned int ntypes> void
  Int1eLibint2::norm_contrcart_(double *data)
  {
    const GaussianShell& sh1 = *int_shell1_;
    const GaussianShell& sh2 = *int_shell2_;

    // if there are no cartesian functions, no need to re-normalize -- spherical harmonics are always normalized
    if (!sh1.has_cartesian() && !sh2.has_cartesian())
      return;

    const Ref<GTOInfo>& gto = GTOInfo::instance();

    for (int gc1=0; gc1<sh1.ncontraction(); gc1++) {
      const unsigned int am1 = sh1.am(gc1);
      const unsigned int n1 = sh1.nfunction(gc1);
      const int is_pure1 = sh1.is_pure(gc1);
      const double* ptr1 = is_pure1 ? gto->fp1() : gto->norm(am1);

      for (int gc2=0; gc2<sh2.ncontraction(); gc2++) {
        const unsigned int am2 = sh2.am(gc2);
        const unsigned int n2 = sh2.nfunction(gc2);
        const int is_pure2 = sh2.is_pure(gc2);
        if ( (am1 <= 1 && am2 <= 1) || (is_pure1 && is_pure2) ) {
          data += n1*n2*ntypes;
          continue;
        }
        const double* ptr2 = is_pure2 ? gto->fp1() : gto->norm(am2);

    /*------------------------------------------------------------------------------------
      Normalize contracted integrals - right here each cartesian component in the shell
      has the same normalization coefficient so that only components with radial parts of
      x^l, y^l, and z^l are normalized to unity. After this block of code all basis
      functions are normalized to unity. Needed this so that integrals in terms of
      puream i-functions were computed properly.
     ------------------------------------------------------------------------------------*/
        for (unsigned ii = 0; ii < n1; ii++) {
          const double norm_i = ptr1[ii];
          for (unsigned jj = 0; jj < n2; jj++) {
            const double norm_ij = norm_i * ptr2[jj];
            for(unsigned int t=0; t<ntypes; ++t)
              *(data++) *= norm_ij;
          }
        }

      }} // end of gen contr loops

  }

  template <unsigned int ntypes>
  void
  Int2eLibint2::norm_contrcart_(double *data)
  {
    const GaussianShell& sh1 = *int_shell1_;
    const GaussianShell& sh2 = *int_shell2_;
    const GaussianShell& sh3 = *int_shell3_;
    const GaussianShell& sh4 = *int_shell4_;

    const Ref<GTOInfo>& gto = GTOInfo::instance();

    for (int gc1=0; gc1<sh1.ncontraction(); gc1++) {
      const unsigned int am1 = sh1.am(gc1);
      const unsigned int n1 = sh1.nfunction(gc1);
      const int is_pure1 = sh1.is_pure(gc1);
      const double* ptr1 = is_pure1 ? gto->fp1() : gto->norm(am1);

      for (int gc2=0; gc2<sh2.ncontraction(); gc2++) {
        const unsigned int am2 = sh2.am(gc2);
        const unsigned int n2 = sh2.nfunction(gc2);
        const int is_pure2 = sh2.is_pure(gc2);
        const double* ptr2 = is_pure2 ? gto->fp1() : gto->norm(am2);

        for (int gc3=0; gc3<sh3.ncontraction(); gc3++) {
          const unsigned int am3 = sh3.am(gc3);
          const unsigned int n3 = sh3.nfunction(gc3);
          const int is_pure3 = sh3.is_pure(gc3);
          const double* ptr3 = is_pure3 ? gto->fp1() : gto->norm(am3);

          for (int gc4=0; gc4<sh4.ncontraction(); gc4++) {
            const unsigned int am4 = sh4.am(gc4);
            const unsigned int n4 = sh4.nfunction(gc4);
            const int is_pure4 = sh4.is_pure(gc4);
            const double* ptr4 = is_pure4 ? gto->fp1() : gto->norm(am4);

            if (is_pure1 && is_pure2 && is_pure3 && is_pure4) {
              data += n1*n2*n3*n4*ntypes;
              continue;
            }

    /*------------------------------------------------------------------------------------
      Normalize contracted integrals - right here each cartesian component in the shell
      has the same normalization coefficient so that only components with radial parts of
      x^l, y^l, and z^l are normalized to unity. After this block of code all basis
      functions are normalized to unity. Needed this so that integrals in terms of
      puream i-functions were computed properly.
     ------------------------------------------------------------------------------------*/
            for (unsigned ii = 0; ii < n1; ii++) {
              const double norm_i = ptr1[ii];
              for (unsigned jj = 0; jj < n2; jj++) {
                const double norm_ij = norm_i * ptr2[jj];
                for (unsigned kk = 0; kk < n3; kk++) {
                  const double norm_ijk = norm_ij * ptr3[kk];
                  for (unsigned ll = 0; ll < n4; ll++) {
                    const double norm_ijkl = norm_ijk * ptr4[ll];
                    for(unsigned int t=0; t<ntypes; ++t)
                      *(data++) *= norm_ijkl;
                  }
                }
              }
            }

          }}}} // end of gen contr loops

  }

} // end of namespace sc

#endif // end of header guard

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
