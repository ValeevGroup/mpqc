//
// macros.h
//
// Copyright (C) 2014 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: Mar 5, 2014
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

#ifndef _chemistry_qc_scf_cadf_macros_h
#define _chemistry_qc_scf_cadf_macros_h

#define M_DUMP(M) std::cout << #M << " is " << M.rows() << " x " << M.cols() << std::endl;

#define M_ROW_ASSERT(M1, M2) \
  if(M1.rows() != M2.rows()) { \
    boost::lock_guard<boost::mutex> _tmp(debug_print_mutex); \
    cout << "assertion failed.  Rows not equal:  M1 => " << M1.rows() << " x " << M1.cols() << ", M2 => " << M2.rows() << " x " << M2.cols() << endl; \
    assert(false); \
  }
#define M_COL_ASSERT(M1, M2) \
  if(M1.cols() != M2.cols()) { \
    boost::lock_guard<boost::mutex> _tmp(debug_print_mutex); \
    cout << "assertion failed. Cols not equal:  M1 => " << M1.rows() << " x " << M1.cols() << ", M2 => " << M2.rows() << " x " << M2.cols() << endl; \
    assert(false); \
  }

#define M_PROD_CHECK(R, M1, M2) \
  if(R.rows() != M1.rows() || R.cols() != M2.cols() || M1.cols() != M2.rows()) { \
    boost::lock_guard<boost::mutex> _tmp(debug_print_mutex); \
    cout << "can't perform multiplication: (" << R.rows() << " x " << R.cols() << ") = (" << M1.rows() << " x " << M1.cols() << ") * (" << M2.rows() << " x " << M2.cols() << ")" << endl; \
    assert(false); \
  }

#define M_DOT_CHECK(M1, M2) \
  if(1 != M1.rows() || 1 != M2.cols() || M1.cols() != M2.rows()) { \
    boost::lock_guard<boost::mutex> _tmp(debug_print_mutex); \
    cout << "can't perform multiplication: (" << 1 << " x " << 1 << ") = (" << M1.rows() << " x " << M1.cols() << ") * (" << M2.rows() << " x " << M2.cols() << ")" << endl; \
    assert(false); \
  }

#define DECOMP_PRINT(D, M) \
  {\
    boost::lock_guard<boost::mutex> _tmp(debug_print_mutex); \
    cout << "Decomposition:  " << D.matrixQR().rows() << " x " << D.matrixQR().cols() << ", M => " << M.rows() << " x " << M.cols() << endl; \
  }

#define M_EQ_ASSERT(M1, M2) M_ROW_ASSERT(M1, M2); M_COL_ASSERT(M1, M2);

#define M_BLOCK_ASSERT(M, b1a, b1b, b2a, b2b) \
  if(b1a < 0) { \
    std::cout << "assertion 1 failed.  data: " << "(" << M.rows() << ", " << M.cols() << ", " << b1a << ", " << b1b << ", " << b2a << ", " << b2b << ")" << std::endl; \
  } \
  else if(b1b < 0) { \
    std::cout << "assertion 2 failed.  data: " << "(" << M.rows() << ", " << M.cols() << ", " << b1a << ", " << b1b << ", " << b2a << ", " << b2b << ")" << std::endl; \
  } \
  else if(b1a > M.rows() - b2a) { \
    std::cout << "assertion 3 failed.  data: " << "(" << M.rows() << ", " << M.cols() << ", " << b1a << ", " << b1b << ", " << b2a << ", " << b2b << ")" << std::endl; \
  } \
  else if(b1b > M.cols() - b2b) { \
    std::cout << "assertion 4 failed.  data: " << "(" << M.rows() << ", " << M.cols() << ", " << b1a << ", " << b1b << ", " << b2a << ", " << b2b << ")" << std::endl; \
  }

#define DUMP(expr) sc::ExEnv::out0() << #expr << " = " << (expr) << std::endl;
#define DUMP2(expr1, expr2) sc::ExEnv::out0() << #expr1 << " = " << (expr1) << ", " << #expr2 << " = " << (expr2) << std::endl;
#define DUMP3(expr1, expr2, expr3) sc::ExEnv::out0() << #expr1 << " = " << (expr1) << ", " << #expr2 << " = " << (expr2) << ", " << #expr3 << " = " << (expr3) << std::endl;
#define DUMP4(expr1, expr2, expr3, expr4) ExEnv::out0() << #expr1 << " = " << (expr1) << ", " << #expr2 << " = " << (expr2) << ", " << #expr3 << " = " << (expr3) << ", " << #expr4 << " = " << (expr4) << std::endl;
#define DUMP5(expr1, expr2, expr3, expr4, expr5) std::cout << #expr1 << " = " << (expr1) << ", " << #expr2 << " = " << (expr2) << ", " << #expr3 << " = " << (expr3) << ", " << #expr4 << " = " << (expr4) << ", " << #expr5 << " = " << (expr5) << std::endl;
#define out_assert(a, op, b) assert(a op b || ((ExEnv::out0() << "Failed assertion output: " << #a << " ( = " << a << ") " << #op << " " << #b <<  " ( = " << b << ")" << std::endl), false))


#endif /* _chemistry_qc_scf_cadf_macros_h */
