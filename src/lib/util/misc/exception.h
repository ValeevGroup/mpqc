//
// exception.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Joseph Kenny <jpkenny@sandia.gov>
// Maintainer: LPS
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

#ifndef _util_misc_exception_h
#define _util_misc_exception_h

#include <cstddef>
#include <exception>
#include <sstream>

#if __cplusplus > 199711L
// C++11
#define MPQC__NOEXCEPT noexcept
#else
#define MPQC__NOEXCEPT throw()
#endif

namespace sc {

  /** This is a std::exception specialization that records information
      about where an exception took place. It serves as the basis for all Exceptions thrown by MPQC.
   */
  class Exception: public std::exception {
      const char *description_;
      const char *file_;
      int line_;

    public:
      /** Create an Exception.

          @param description a description of the problem.
          @param file the file name where the problem occured.
          @param line the line number where the exception occured.

          It is suggested that the special macros __FILE__ and __LINE__ be
          given as the \p file and \p line arguments, respectively.
      */
      Exception(const char *description = 0,
                const char *file = 0,
                int line = 0) MPQC__NOEXCEPT;
      Exception(const Exception&) MPQC__NOEXCEPT;
      ~Exception() MPQC__NOEXCEPT;

      /** Reimplementation of std::exception::what().  The returned
          std::string is only valid for the lifetime of this object. */
      virtual const char* what() const MPQC__NOEXCEPT;

      /// Returns a description of what caused the exception.  May return
      /// null.
      const char *description() const MPQC__NOEXCEPT { return description_; }
      /// Returns the name of the file in which the exception was created.
      /// May return null.
      const char *file() const MPQC__NOEXCEPT { return file_; }
      /// Returns the line number where the exception was created.
      /// May return 0, if unknown.
      int line() const MPQC__NOEXCEPT { return line_; }

  };

}

#endif

