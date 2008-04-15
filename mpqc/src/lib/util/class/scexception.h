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

#ifndef _util_class_class_h
#include <util/class/class.h>
#endif

#include <stddef.h>
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
    const char *exception_type_;
    mutable char *elaboration_c_str_;
    std::ostringstream *elaboration_;

  public:
    /** Create an SCException.

        @param description a description of the problem.
        @param file the file name where the problem occured.
        @param line the line number where the exception occured.
        @param class_desc the ClassDesc for the object causing the
        exception.
        @param exception_type the classname of the SCException
        specialization. The default is "SCException".

        It is suggested that the special macros __FILE__ and __LINE__ be
        given as the \p file and \p line arguments, respectively.
    */
    SCException(const char *description = 0,
                const char *file = 0,
                int line = 0,
                const ClassDesc *class_desc = 0,
                const char *exception_type = "SCException") throw();
    SCException(const SCException&) throw();
    ~SCException() throw();

    /** Reimplementation of std::exception::what().  The returned
        std::string is only valid for the lifetime of this object. */
    const char* what() const throw();

    /// Returns a description of what caused the exception.  May return
    /// null.
    const char *description() const throw() { return description_; }
    /// Returns the name of the file in which the exception was created.
    /// May return null.
    const char *file() const throw() { return file_; }
    /// Returns the line number where the exception was created.
    /// May return 0, if unknown.
    int line() const throw() { return line_; }
    /// Returns the class descriptor of the object which generated the
    /// exception. May return null.
    const ClassDesc *class_desc() const throw() { return class_desc_; }
    /// Returns the classname of the exception.  May return null.
    const char *exception_type() const throw() { return exception_type_; }

    /** Returns a stream where addition information about the exception can
        be written.  This will throw if a valid stream cannot be returned
        (possibly due to low memory). */
    std::ostream &elaborate();
};

// ///////////////////////////////////////////////////////////////////////
// Programming Error Exceptions

/** This is thrown when a situations arises that should be impossible.
 */
class ProgrammingError: public SCException {

  public:
    /** Create a ProgrammingError exception.

        @param description a description of the problem.
        @param file the file name where the problem occured.
        @param line the line number where the exception occured.
        @param class_desc the ClassDesc for the object causing the
        exception.
        @param exception_type the classname of the SCException
        specialization. The default is "ProgrammingError".

        It is suggested that the special macros __FILE__ and __LINE__ be
        given as the \p file and \p line arguments, respectively.
    */
    ProgrammingError(const char *description = 0,
                     const char *file = 0,
                     int line = 0,
                     const ClassDesc *class_desc = 0,
                     const char *exception_type = "ProgrammingError") throw();
    ProgrammingError(const ProgrammingError&) throw();
    ~ProgrammingError() throw();
};

/** This is thrown when an attempt is made to use a feature that
    is not yet implemented.
 */
class FeatureNotImplemented: public ProgrammingError {

  public:
    /** Create a FeatureNotImplemented exception.

        @param description a description of the problem.
        @param file the file name where the problem occured.
        @param line the line number where the exception occured.
        @param class_desc the ClassDesc for the object causing the
        exception.
        @param exception_type the classname of the SCException
        specialization. The default is "FeatureNotImplemented".

        It is suggested that the special macros __FILE__ and __LINE__ be
        given as the \p file and \p line arguments, respectively.
    */
    FeatureNotImplemented(const char *description = 0,
                          const char *file = 0,
                          int line = 0,
                          const ClassDesc *class_desc = 0,
                          const char *exception_type = "FeatureNotImplemented")
        throw();
    FeatureNotImplemented(const FeatureNotImplemented&) throw();
    ~FeatureNotImplemented() throw();
};

// ///////////////////////////////////////////////////////////////////////
// Input Error Exceptions

/** This is thrown when invalid input is provided.  Note that sometimes
    input can be internally generated, so what logically would be a
    ProgrammingError could result in an InputError being thrown.
 */
class InputError: public SCException {
    const char *keyword_;
    char *value_;

  public:
    /** Create a InputError exception.

        @param description a description of the problem.
        @param file the file name where the problem occured.
        @param line the line number where the exception occured.
        @param keyword the keyword that was being read.
        @param value the value associated with the keyword
        @param class_desc the ClassDesc for the object causing the
        exception.
        @param exception_type the classname of the SCException
        specialization. The default is "InputError".

        It is suggested that the special macros __FILE__ and __LINE__ be
        given as the \p file and \p line arguments, respectively.
    */
    InputError(const char *description = 0,
               const char *file = 0,
               int line = 0,
               const char *keyword = 0,
               const char *value = 0,
               const ClassDesc *class_desc = 0,
               const char *exception_type = "InputError") throw();
    InputError(const InputError&) throw();
    ~InputError() throw();
    /// Return the keyword having an erroneous value.
    const char *keyword() const throw() { return keyword_; }
    /// Return the erroneous value which caused this exception to be
    /// thrown.
    const char *value() const throw() { return value_; }
};

// ///////////////////////////////////////////////////////////////////////
// System Exceptions

/** This is thrown when a system problem occurs.
 */
class SystemException: public SCException {

  public:
    /** Create a SystemException exception.

        @param description a description of the problem.
        @param file the file name where the problem occured.
        @param line the line number where the exception occured.
        @param class_desc the ClassDesc for the object causing the
        exception.
        @param exception_type the classname of the SCException
        specialization. The default is "SystemException".

        It is suggested that the special macros __FILE__ and __LINE__ be
        given as the \p file and \p line arguments, respectively.
    */
    SystemException(const char *description = 0,
                    const char *file = 0,
                    int line = 0,
                    const ClassDesc *class_desc = 0,
                    const char *exception_type = "SystemException") throw();
    SystemException(const SystemException&) throw();
    ~SystemException() throw();
};

/** This is thrown when a memory allocation fails.
 */
class MemAllocFailed: public SystemException {
    size_t nbyte_;

  public:
    /** Create a MemAllocFailed exception.

        @param description a description of the problem.
        @param file the file name where the problem occured.
        @param line the line number where the exception occured.
        @param nbyte the size of the attempted allocation.
        @param class_desc the ClassDesc for the object causing the
        exception.
        @param exception_type the classname of the SCException
        specialization. The default is "MemAllocFailed".

        It is suggested that the special macros __FILE__ and __LINE__ be
        given as the \p file and \p line arguments, respectively.
    */
    MemAllocFailed(const char *description = 0,
                   const char *file = 0,
                   int line = 0,
                   size_t nbyte = 0,
                   const ClassDesc *class_desc = 0,
                   const char *exception_type = "MemAllocFailed") throw();
    MemAllocFailed(const MemAllocFailed&) throw();
    ~MemAllocFailed() throw();

    /// Returns the number of bytes used in the failed allocation attempt.
    size_t nbyte() const throw() { return nbyte_; }
};

/** This is thrown when an operation on a file fails.
 */
class FileOperationFailed: public SystemException {
  public:
    enum FileOperation { Unknown, OpenR, OpenW, OpenRW,
                         Close, Read, Write, Corrupt, Other };

  private:
    const char *filename_;
    FileOperation operation_;

  public:
    /** Create a FileOperationFailure exception.

        @param description a description of the problem.
        @param file the file name where the problem occured.
        @param line the line number where the exception occured.
        @param filename the name of the file for which the operation failed.
        @param operation the type of fail operation that resulted in the
        failure.
        @param class_desc the ClassDesc for the object causing the
        exception.
        @param exception_type the classname of the SCException
        specialization. The default is "FileOperationFailure".

        It is suggested that the special macros __FILE__ and __LINE__ be
        given as the \p file and \p line arguments, respectively.
    */
    FileOperationFailed(const char *description = 0,
                   const char *source_file = 0,
                   int line = 0,
                   const char *filename = 0,
                   FileOperation operation = Unknown,
                   const ClassDesc *class_desc = 0,
                   const char *exception_type = "FileOperationFailed") throw();
    FileOperationFailed(const FileOperationFailed&) throw();
    ~FileOperationFailed() throw();

    /** Returns the file name of the file that caused the error, if known.
        Otherwise 0 is returned. */
    const char * filename() const throw() { return filename_; }
    /// Return the file operation that failed as a FileOperation enum.
    FileOperation operation() const throw() { return operation_; }
};

/** This is thrown when an system call fails with an errno.
 */
class SyscallFailed: public SystemException {
    const char *syscall_;
    int err_;

  public:
    /** Create a SyscallFailed exception.

        @param description a description of the problem.
        @param file the file name where the problem occured.
        @param line the line number where the exception occured.
        @param syscall the name of the syscall that failed.
        @param err the error returned by the failed syscall.
        @param class_desc the ClassDesc for the object causing the
        exception.
        @param exception_type the classname of the SCException
        specialization. The default is "SyscallFailed".

        It is suggested that the special macros __FILE__ and __LINE__ be
        given as the \p file and \p line arguments, respectively.
    */
    SyscallFailed(const char *description = 0,
                  const char *source_file = 0,
                  int line = 0,
                  const char *syscall = 0,
                  int err = 0, 
                  const ClassDesc *class_desc = 0,
                  const char *exception_type = "SyscallFailed") throw();
    SyscallFailed(const SyscallFailed&) throw();
    ~SyscallFailed() throw();

    /** Returns the file name of the file that caused the error, if known.
        Otherwise 0 is returned. */
    const char * syscall() const throw() { return syscall_; }
    /// Return the error code that the system call returned.
    int err() const throw() { return err_; }
};

// ///////////////////////////////////////////////////////////////////////
// Algorithm Exceptions

/** This exception is thrown whenever a problem with an algorithm is
    encountered.
*/
class AlgorithmException: public SCException {

  public:
    /** Create an AlgorithmException.

        @param description a description of the problem.
        @param file the file name where the problem occured.
        @param line the line number where the exception occured.
        @param class_desc the ClassDesc for the object causing the
        exception.
        @param exception_type the classname of the SCException
        specialization. The default is "AlgorithmException".

        It is suggested that the special macros __FILE__ and __LINE__ be
        given as the \p file and \p line arguments, respectively.
    */
    AlgorithmException(const char *description = 0,
                       const char *file = 0,
                       int line = 0,
                       const ClassDesc *class_desc = 0,
                       const char *exception_type = "AlgorithmException")
        throw();
    AlgorithmException(const AlgorithmException&) throw();
    ~AlgorithmException() throw();
};

/** This is thrown when an iterative algorithm attempts to use more
    iterations than allowed.
 */
class MaxIterExceeded: public AlgorithmException {
    int max_iter_;

  public:
    /** Create a MaxIterExceeded exception.

        @param description a description of the problem.
        @param file the file name where the problem occured.
        @param line the line number where the exception occured.
        @param maxiter the maximum number of iterations.
        @param class_desc the ClassDesc for the object causing the
        exception.
        @param exception_type the classname of the SCException
        specialization. The default is "MaxIterExceeded".

        It is suggested that the special macros __FILE__ and __LINE__ be
        given as the \p file and \p line arguments, respectively.
    */
    MaxIterExceeded(const char *description = 0,
                    const char *file = 0,
                    int line = 0,
                    int maxiter = 0,
                    const ClassDesc *class_desc = 0,
                    const char *exception_type = "MaxIterExceeded") throw();
    MaxIterExceeded(const MaxIterExceeded&) throw();
    ~MaxIterExceeded() throw();

    /// Return the maximum number of iterations.
    int max_iter() const throw() { return max_iter_; }
};

/** This is thrown when when some tolerance is exceeded.
 */
class ToleranceExceeded: public AlgorithmException {
    double tolerance_;
    double value_;

public:
    /** Create a ToleranceExceeded exception.

        @param description a description of the problem.
        @param file the file name where the problem occured.
        @param line the line number where the exception occured.
        @param tol the required tolerance.
        @param val the value which was obtained.
        @param class_desc the ClassDesc for the object causing the
        exception.
        @param exception_type the classname of the SCException
        specialization. The default is "ToleranceExceeded".

        It is suggested that the special macros __FILE__ and __LINE__ be
        given as the \p file and \p line arguments, respectively.
    */
    ToleranceExceeded(const char *description = 0,
                      const char *file = 0,
                      int line = 0,
                      double tol=0,
                      double val=0,
                      const ClassDesc *class_desc = 0,
                      const char *exception_type = "ToleranceExceeded") throw();
    ToleranceExceeded(const ToleranceExceeded&) throw();
    ~ToleranceExceeded() throw();
    /// Return the required tolerance.
    double tolerance() throw() { return tolerance_; }
    /// Return the value which was obtained.
    double value() throw() { return value_; }
};

// ///////////////////////////////////////////////////////////////////////
// Limit Exceeded Exceptions

/** This is thrown when a limit is exceeded.  It is more general than
    ToleranceExceeded.  For problems that are numerical in nature and use
    double types, then ToleranceExceeded should be used instead.
*/
template <class T>
class LimitExceeded: public SCException {
    T limit_;
    T value_;

public:
    /** Create a LimitExceeded exception.

        @param description a description of the problem.
        @param file the file name where the problem occured.
        @param line the line number where the exception occured.
        @param lim the limit.
        @param val the value which was obtained.
        @param class_desc the ClassDesc for the object causing the
        exception.
        @param exception_type the classname of the SCException
        specialization. The default is "LimitedExceeded".

        It is suggested that the special macros __FILE__ and __LINE__ be
        given as the \p file and \p line arguments, respectively.
    */
    LimitExceeded(const char *description,
                  const char *file,
                  int line,
                  T lim,
                  T val,
                  const ClassDesc *class_desc = 0,
                  const char *exception_type = "LimitExceeded") throw():
      SCException(description, file, line, class_desc, exception_type),
      limit_(lim), value_(val)
        {
          try {
              elaborate() << "value:       " << value_
                          << std::endl
                          << "limit:       " << limit_
                          << std::endl;
            }
          catch(...) {
            }
        }
    LimitExceeded(const LimitExceeded&ref) throw():
      SCException(ref),
      limit_(ref.limit_), value_(ref.value_)
        {
        }
    ~LimitExceeded() throw() {}
    /// The limit which was exceeded.
    T tolerance() throw() { return limit_; }
    /// The value which exceeded the limit.
    T value() throw() { return value_; }
};

}

#endif

