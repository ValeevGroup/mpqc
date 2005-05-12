//
// libint2_utils.h
//
// Copyright (C) 2001 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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

#ifndef _chemistry_qc_libint2_libint2utils_h
#define _chemistry_qc_libint2_libint2utils_h

namespace {
  typedef Libint_t prim_data;
  inline void assign_FjT(prim_data* Data, const int& jmax, const double* FjT, const double& scale) {
    switch (jmax) {
      case 12:
      Data->__ss_1_over_r_12_ss___up_12[0] = FjT[12] * scale;
      case 11:
      Data->__ss_1_over_r_12_ss___up_11[0] = FjT[11] * scale;
      case 10:
      Data->__ss_1_over_r_12_ss___up_10[0] = FjT[10] * scale;
      case 9:
      Data->__ss_1_over_r_12_ss___up_9[0] = FjT[9] * scale;
      case 8:
      Data->__ss_1_over_r_12_ss___up_8[0] = FjT[8] * scale;
      case 7:
      Data->__ss_1_over_r_12_ss___up_7[0] = FjT[7] * scale;
      case 6:
      Data->__ss_1_over_r_12_ss___up_6[0] = FjT[6] * scale;
      case 5:
      Data->__ss_1_over_r_12_ss___up_5[0] = FjT[5] * scale;
      case 4:
      Data->__ss_1_over_r_12_ss___up_4[0] = FjT[4] * scale;
      case 3:
      Data->__ss_1_over_r_12_ss___up_3[0] = FjT[3] * scale;
      case 2:
      Data->__ss_1_over_r_12_ss___up_2[0] = FjT[2] * scale;
      case 1:
      Data->__ss_1_over_r_12_ss___up_1[0] = FjT[1] * scale;
      case 0:
      Data->__ss_1_over_r_12_ss___up_0[0] = FjT[0] * scale;
      break;
      default:
      throw std::logic_error("assign_FjT() -- max_am exceeded");
    }
  }
};

#endif

