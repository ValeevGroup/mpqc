//
// sharedptr.h
//
// Copyright (C) 2013 MPQC Developers
//
// Author: David Hollman <dhollman@vt.edu>
// Maintainer: DSH, EV
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

/*******************************************************************************
 * Temporary fix to Eclipse parsing issues with shared_ptr.  Also makes it
 * possible to use some other shared_ptr implementation as a std::shared_ptr
 * if C++11 is not available.
 *******************************************************************************/

#ifndef _lib_util_misc_sharedptr_h
#define _lib_util_misc_sharedptr_h
#include <memory>

#ifdef ECLIPSE_PARSER_ONLY

#include <cassert>
namespace std {
  template <typename T>
  class shared_ptr {
    public:
      // Act like a dumb pointer.  This should never actually get compiled.
      shared_ptr(T* ptr){ assert(false && "should never get here"); ptr_ = ptr; }
      T operator*(){ assert(false && "should never get here"); return *ptr_; }
      const T operator*() const { assert(false && "should never get here"); return *ptr_; }
      T* operator->(){ assert(false && "should never get here"); return ptr_; }
      const T operator->() const { assert(false && "should never get here"); return *ptr_; }
      T* get(){ assert(false && "should never get here"); return ptr_; }
      const T* get() const { assert(false && "should never get here"); return ptr_; }

    protected:
      T* ptr_;
  };
  template <typename T>
  class unique_ptr {
    public:
      // Act like a dumb pointer.  This should never actually get compiled.
      shared_ptr(T* ptr){ assert(false && "should never get here"); ptr_ = ptr; }
      T operator*(){ assert(false && "should never get here"); return *ptr_; }
      const T operator*() const { assert(false && "should never get here"); return *ptr_; }
      T* operator->(){ assert(false && "should never get here"); return ptr_; }
      const T operator->() const { assert(false && "should never get here"); return *ptr_; }
      T* get(){ assert(false && "should never get here"); return ptr_; }
      const T* get() const { assert(false && "should never get here"); return ptr_; }

    protected:
      T* ptr_;
  };

  template <typename T>
  struct make_shared {
      template<typename... Args>
      static shared_ptr<T>
      operator()(Args args...){
        return shared_ptr<T>(new T(args...));
      }
  };



  template <typename T> using weak_ptr = shared_ptr<T>;

}

#endif

#endif
