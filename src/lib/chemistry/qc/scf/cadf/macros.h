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

// TODO Move this to a more global, reasonable place
#define DECLARE_ATTRIBUTE_LESS_STRUCT_NAMED(struct_name, attr_name) \
  template<typename T> \
  struct struct_name \
  { \
    bool operator()(const T& a, const T& b) const { \
      return a.##attr_name < b.##attr_name; \
    } \
  }
#define DECLARE_ATTRIBUTE_LESS_STRUCT(attr_name) \
    DECLARE_ATTRIBUTE_LESS_STRUCT_NAMED(attr_name##_less, attr_name)


#define M_DUMP(M) std::cout << #M << " is " << M.rows() << " x " << M.cols() << std::endl;

#define M_ELEM_DUMP(M, r, c) \
    ExEnv::out0() << #M << " is " << M.rows() << " x " << M.cols() << ", trying to access element (" \
    << r << ", " << c << ")." << std::endl;

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
#define DUMPme(expr) sc::ExEnv::outn() << "  Node " << scf_grp_->me() << ": " << #expr << " = " << (expr) << std::endl;
#define DUMPn(expr) \
    sc::ExEnv::out0() << "Dumping expression " << #expr << " from all nodes:" << std::endl; \
    for(int __inode = 0; __inode < scf_grp_->n(); ++__inode) {\
      if(__inode == scf_grp_->me()) { \
        sc::ExEnv::outn() << "  Node " << scf_grp_->me() << ": " << #expr << " = " << expr << std::endl; \
      } \
      scf_grp_->sync(); \
    }
#define DUMP2(expr1, expr2) sc::ExEnv::out0() << #expr1 << " = " << (expr1) << ", " << #expr2 << " = " << (expr2) << std::endl;
#define DUMP2n(expr1, expr2) sc::ExEnv::outn() << "me = " << scf_grp_->me() << ", " << #expr1 << " = " << (expr1) << ", " << #expr2 << " = " << (expr2) << std::endl;
#define DUMP3(expr1, expr2, expr3) sc::ExEnv::out0() << #expr1 << " = " << (expr1) << ", " << #expr2 << " = " << (expr2) << ", " << #expr3 << " = " << (expr3) << std::endl;
#define DUMP4(expr1, expr2, expr3, expr4) ExEnv::out0() << #expr1 << " = " << (expr1) << ", " << #expr2 << " = " << (expr2) << ", " << #expr3 << " = " << (expr3) << ", " << #expr4 << " = " << (expr4) << std::endl;
#define DUMP5(expr1, expr2, expr3, expr4, expr5) std::cout << #expr1 << " = " << (expr1) << ", " << #expr2 << " = " << (expr2) << ", " << #expr3 << " = " << (expr3) << ", " << #expr4 << " = " << (expr4) << ", " << #expr5 << " = " << (expr5) << std::endl;
#define out_assert(a, op, b) assert(a op b || ((ExEnv::out0() << "Failed assertion output: " << #a << " ( = " << a << ") " << #op << " " << #b <<  " ( = " << b << ")" << std::endl), false))
#define DEBUG_DELETE_THIS

#define resize_and_zero_matrix(mat, nrows, ncols) mat.resize(nrows, ncols); mat = std::remove_reference<decltype(mat)>::type::Zero(nrows, ncols)
#define declare_and_zero_matrix(type, mat, nrows, ncols) type mat(nrows, ncols); mat = type::Zero(nrows, ncols)

#define MAKE_MATRIX(type, name, rows, cols) type name(rows, cols); auto __##name##_mem_holder = hold_memory(rows*cols*sizeof(Eigen::internal::traits<type>::Scalar) + sizeof(type));

//#define PRINT_AS_GRID(o, itervar, iterable, expr, fmt, nperline) \
//    { int __numdone = 0; \
//      for(auto&& itervar : iterable) { \
//        if(__numdone != 0 && __numdone % nperline == 0) { o << std::endl << sc::indent; } \
//        o << std::string(sc::scprintf(fmt, (expr))); \
//        ++__numdone; \
//      } \
//    }

#define PRINT_STRING_AS_GRID(o, itervar, iterable, str, nperline) \
    { int __numdone = 0; \
      for(auto&& itervar : iterable) { \
        if(__numdone != 0 && __numdone % nperline == 0) { o << std::endl << sc::indent; } \
        o << str << " "; \
        ++__numdone; \
      } \
    }

#define PRINT_AS_GRID(o, itervar, iterable, expr, fmt, nperline) PRINT_STRING_AS_GRID(o, itervar, iterable, (std::string(sc::scprintf(fmt, (expr)).str())), nperline)

#endif /* _chemistry_qc_scf_cadf_macros_h */
