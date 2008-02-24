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
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(28))
      case 28:
      Data->LIBINT_T_SS_EREP_SS(28)[0] = FjT[28] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(27))
      case 27:
      Data->LIBINT_T_SS_EREP_SS(27)[0] = FjT[27] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(26))
      case 26:
      Data->LIBINT_T_SS_EREP_SS(26)[0] = FjT[26] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(25))
      case 25:
      Data->LIBINT_T_SS_EREP_SS(25)[0] = FjT[25] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(24))
      case 24:
      Data->LIBINT_T_SS_EREP_SS(24)[0] = FjT[24] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(23))
      case 23:
      Data->LIBINT_T_SS_EREP_SS(23)[0] = FjT[23] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(22))
      case 22:
      Data->LIBINT_T_SS_EREP_SS(22)[0] = FjT[22] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(21))
      case 21:
      Data->LIBINT_T_SS_EREP_SS(21)[0] = FjT[21] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(20))
      case 20:
      Data->LIBINT_T_SS_EREP_SS(20)[0] = FjT[20] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(19))
      case 19:
      Data->LIBINT_T_SS_EREP_SS(19)[0] = FjT[19] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(18))
      case 18:
      Data->LIBINT_T_SS_EREP_SS(18)[0] = FjT[18] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(17))
      case 17:
      Data->LIBINT_T_SS_EREP_SS(17)[0] = FjT[17] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(16))
      case 16:
      Data->LIBINT_T_SS_EREP_SS(16)[0] = FjT[16] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(15))
      case 15:
      Data->LIBINT_T_SS_EREP_SS(15)[0] = FjT[15] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(14))
      case 14:
      Data->LIBINT_T_SS_EREP_SS(14)[0] = FjT[14] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(13))
      case 13:
      Data->LIBINT_T_SS_EREP_SS(13)[0] = FjT[13] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(12))
      case 12:
      Data->LIBINT_T_SS_EREP_SS(12)[0] = FjT[12] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(11))
      case 11:
      Data->LIBINT_T_SS_EREP_SS(11)[0] = FjT[11] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(10))
      case 10:
      Data->LIBINT_T_SS_EREP_SS(10)[0] = FjT[10] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(9))
      case 9:
      Data->LIBINT_T_SS_EREP_SS(9)[0] = FjT[9] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(8))
      case 8:
      Data->LIBINT_T_SS_EREP_SS(8)[0] = FjT[8] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(7))
      case 7:
      Data->LIBINT_T_SS_EREP_SS(7)[0] = FjT[7] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(6))
      case 6:
      Data->LIBINT_T_SS_EREP_SS(6)[0] = FjT[6] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(5))
      case 5:
      Data->LIBINT_T_SS_EREP_SS(5)[0] = FjT[5] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(4))
      case 4:
      Data->LIBINT_T_SS_EREP_SS(4)[0] = FjT[4] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(3))
      case 3:
      Data->LIBINT_T_SS_EREP_SS(3)[0] = FjT[3] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(2))
      case 2:
      Data->LIBINT_T_SS_EREP_SS(2)[0] = FjT[2] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(1))
      case 1:
      Data->LIBINT_T_SS_EREP_SS(1)[0] = FjT[1] * scale;
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(0))
      case 0:
      Data->LIBINT_T_SS_EREP_SS(0)[0] = FjT[0] * scale;
#endif
      break;
      default:
      throw std::logic_error("assign_FjT() -- max_am exceeded");
    }
  }

  inline void assign_ss_r12m1g12_ss(prim_data* Data, const int& jmax, const double* ss_r12m1g12_ss, const double& scale) {
    switch (jmax) {
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(28))
      case 28:
      Data->LIBINT_T_SS_Km1G12_SS(28)[0] = ss_r12m1g12_ss[28] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(27))
      case 27:
      Data->LIBINT_T_SS_Km1G12_SS(27)[0] = ss_r12m1g12_ss[27] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(26))
      case 26:
      Data->LIBINT_T_SS_Km1G12_SS(26)[0] = ss_r12m1g12_ss[26] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(25))
      case 25:
      Data->LIBINT_T_SS_Km1G12_SS(25)[0] = ss_r12m1g12_ss[25] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(24))
      case 24:
      Data->LIBINT_T_SS_Km1G12_SS(24)[0] = ss_r12m1g12_ss[24] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(23))
      case 23:
      Data->LIBINT_T_SS_Km1G12_SS(23)[0] = ss_r12m1g12_ss[23] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(22))
      case 22:
      Data->LIBINT_T_SS_Km1G12_SS(22)[0] = ss_r12m1g12_ss[22] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(21))
      case 21:
      Data->LIBINT_T_SS_Km1G12_SS(21)[0] = ss_r12m1g12_ss[21] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(20))
      case 20:
      Data->LIBINT_T_SS_Km1G12_SS(20)[0] = ss_r12m1g12_ss[20] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(19))
      case 19:
      Data->LIBINT_T_SS_Km1G12_SS(19)[0] = ss_r12m1g12_ss[19] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(18))
      case 18:
      Data->LIBINT_T_SS_Km1G12_SS(18)[0] = ss_r12m1g12_ss[18] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(17))
      case 17:
      Data->LIBINT_T_SS_Km1G12_SS(17)[0] = ss_r12m1g12_ss[17] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(16))
      case 16:
      Data->LIBINT_T_SS_Km1G12_SS(16)[0] = ss_r12m1g12_ss[16] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(15))
      case 15:
      Data->LIBINT_T_SS_Km1G12_SS(15)[0] = ss_r12m1g12_ss[15] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(14))
      case 14:
      Data->LIBINT_T_SS_Km1G12_SS(14)[0] = ss_r12m1g12_ss[14] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(13))
      case 13:
      Data->LIBINT_T_SS_Km1G12_SS(13)[0] = ss_r12m1g12_ss[13] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(12))
      case 12:
      Data->LIBINT_T_SS_Km1G12_SS(12)[0] = ss_r12m1g12_ss[12] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(11))
      case 11:
      Data->LIBINT_T_SS_Km1G12_SS(11)[0] = ss_r12m1g12_ss[11] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(10))
      case 10:
      Data->LIBINT_T_SS_Km1G12_SS(10)[0] = ss_r12m1g12_ss[10] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(9))
      case 9:
      Data->LIBINT_T_SS_Km1G12_SS(9)[0] = ss_r12m1g12_ss[9] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(8))
      case 8:
      Data->LIBINT_T_SS_Km1G12_SS(8)[0] = ss_r12m1g12_ss[8] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(7))
      case 7:
      Data->LIBINT_T_SS_Km1G12_SS(7)[0] = ss_r12m1g12_ss[7] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(6))
      case 6:
      Data->LIBINT_T_SS_Km1G12_SS(6)[0] = ss_r12m1g12_ss[6] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(5))
      case 5:
      Data->LIBINT_T_SS_Km1G12_SS(5)[0] = ss_r12m1g12_ss[5] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(4))
      case 4:
      Data->LIBINT_T_SS_Km1G12_SS(4)[0] = ss_r12m1g12_ss[4] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(3))
      case 3:
      Data->LIBINT_T_SS_Km1G12_SS(3)[0] = ss_r12m1g12_ss[3] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(2))
      case 2:
      Data->LIBINT_T_SS_Km1G12_SS(2)[0] = ss_r12m1g12_ss[2] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(1))
      case 1:
      Data->LIBINT_T_SS_Km1G12_SS(1)[0] = ss_r12m1g12_ss[1] * scale;
#endif
#if LIBINT2_DEFINED(r12kg12,LIBINT_T_SS_Km1G12_SS(0))
      case 0:
      Data->LIBINT_T_SS_Km1G12_SS(0)[0] = ss_r12m1g12_ss[0] * scale;
#endif
      break;
      default:
      throw std::logic_error("assign_ss_r12m1g12_ss() -- max_am exceeded");
    }
  }
};

#endif

