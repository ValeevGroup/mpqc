//
// scexception.h
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

#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_misc_scexception_h
#define _util_misc_scexception_h
#endif

#ifndef _util_class_class_h
#include <util/class/class.h>
#endif

#include <exception>
#include <sstream>
#include <vector>

namespace sc {

/** This is a std::exception specialization that records information
    about where an exception took place.
 */
class SCException: public std::exception {
    const char *description_;
    const char *file_;
    int line_;
    const ClassDesc* class_desc_;
    std::ostringstream *elaboration_;

  public:
    SCException(const char *description = 0,
                const char *file = 0,
                int line = 0,
                const ClassDesc *class_desc = 0) throw();
    SCException(const SCException&) throw();
    ~SCException() throw();

    /** Reimplementation of std::exception::what().  The returned
        std::string is only valid for the lifetime of this object. */
    const char* what() throw();

    const char *description() const throw() { return description_; }
    const char *file() const throw() { return file_; }
    int line() const throw() { return line_; }
    const ClassDesc *class_desc() const throw() { return class_desc_; }

    /** Returns a stream where addition information about the exception can
        be written.  This will throw if it is impossible to elaborate
        (possibly due to low memory), so it must be used in a try block. */
    std::ostream &elaborate();
};

/////////////////////////////////////////////////////////////////////////
// Algorithm Exceptions

/** This exception is thrown whenever a problem with an algorithm is
    encountered.  Usually, a class derived from this is thrown, such as
    MaxIterExceeded or ToleranceExeeded.
*/
class AlgorithmException: public SCException {

  public:
    AlgorithmException(const char *description = 0,
                       const char *file = 0,
                       int line = 0,
                       const ClassDesc *class_desc = 0) throw();
    AlgorithmException(const AlgorithmException&) throw();
    ~AlgorithmException() throw();
};

/** This is thrown when an iterative algorithm attempts to use more
    iterations than allowed.
 */
class MaxIterExceeded: public AlgorithmException {
    int max_iter_;

  public:
    MaxIterExceeded(const char *description = 0,
                    const char *file = 0,
                    int line = 0,
                    int maxiter = 0,
                    const ClassDesc *class_desc = 0) throw();
    MaxIterExceeded(const MaxIterExceeded&) throw();
    ~MaxIterExceeded() throw();

    int max_iter() const throw() { return max_iter_; }
};

/** This is thrown when when some tolerance is exceeded.
 */
class ToleranceExceeded: public AlgorithmException {
    const double tolerance_;
    const double value_;

public:
    ToleranceExceeded(const char *description = 0,
                      const char *file = 0,
                      int line = 0,
                      double tol=0,
                      double val=0,
                      const ClassDesc *class_desc = 0) throw();
    ToleranceExceeded(const ToleranceExceeded&) throw();
    ~ToleranceExceeded() throw();
    double tolerance() throw() { return tolerance_; }
    double value() throw() { return value_; }
};

}
