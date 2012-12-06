//
// scexception.cc
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

#include <errno.h>
#include <string.h>

#include <stdexcept>
#include <sstream>
#include <util/state/stateio.h>
#include <util/misc/scexception.h>

using namespace std;
using namespace sc;

////////////////////////////////////////////////////////////////////////
// SCException

SCException::SCException(const char *description,
                         const char *file,
                         int line,
                         const ClassDesc *class_desc,
                         const char *exception_type) throw():
  description_(description),
  file_(file),
  line_(line),
  class_desc_(class_desc),
  exception_type_(exception_type),
  elaboration_c_str_(0),
  backtrace_("=scbacktrace=: ")
{
  try {
      elaboration_ = new ostringstream;
      if (exception_type_) {
          elaborate() << "exception:   " << exception_type_
                      << std::endl;
        }
      if (description_) {
           elaborate() << "description: " << description
                      << std::endl;
        }
      if (file_) {
          elaborate() << "location:    " << file << ":" << line
                      << std::endl;
        }
      if (class_desc_) {
          elaborate() << "class:       " << class_desc_->name()
                      << std::endl;
        }
      const size_t nframes_to_skip = 1;
      elaborate() << "backtrace:" << std::endl << backtrace_.str(nframes_to_skip);
    }
  catch (...) {
      // info in the elaboration is incomplete, so delete it
      // and at least the basic description will be available
      delete elaboration_;
      elaboration_ = 0;
    }
}

SCException::SCException(const SCException& ref) throw(): 
  description_(ref.description_),
  file_(ref.file_),
  line_(ref.line_),
  class_desc_(ref.class_desc_),
  backtrace_(ref.backtrace_)
{
  elaboration_c_str_ = 0;
  elaboration_ = 0;
  if (ref.elaboration_) {
      try {
          elaboration_ = new ostringstream;
          elaborate() << ref.elaboration_->str();
        }
      catch (...) {
          delete elaboration_;
          elaboration_ = 0;
        }
    }
}

SCException::~SCException() throw()
{
  try{ ExEnv::out0().flush(); ExEnv::err0().flush(); }
  catch(...) {}
  delete elaboration_;
  delete elaboration_c_str_;
}

const char* 
SCException::what() const throw()
{
  try {
      if (elaboration_) {
          std::string elab(elaboration_->str());
          delete[] elaboration_c_str_;
          elaboration_c_str_ = 0;
          elaboration_c_str_ = new char[1+elab.size()];
          for (int i=0; i<elab.size(); i++) elaboration_c_str_[i] = elab[i];
          elaboration_c_str_[elab.size()] = '\0';
          return elaboration_c_str_;
        }
    }
  catch (...) {
      // Ignore the exception and return the next available string.
    }
  if (description_) {
      return description_;
    }
  return "No information available for SCException";
}

std::ostream &
SCException::elaborate()
{
  if (!elaboration_) {
      throw std::runtime_error("SCException::elaborate(): cannot elaborate");
    }
  return *elaboration_;
}

////////////////////////////////////////////////////////////////////////
// InputError

InputError::InputError(
    const char *description,
    const char *file,
    int line,
    const char *keyword,
    const char *value,
    const ClassDesc *class_desc,
    const char *exception_type) throw():
  SCException(description, file, line, class_desc, exception_type),
  keyword_(keyword)
{
  try {
      if (value) {
          value_ = new char[strlen(value)+1];
          if (value_) strcpy(value_, value);
        }
      else {
          value_ = 0;
        }
    }
  catch (...) {
    value_ = 0;
    }

  try {
      if (keyword_)
          elaborate() << "keyword:     " << keyword_ << std::endl;
      if (value_)
          elaborate() << "value:       " << value_ << std::endl;
    }
  catch (...) {
    }
}
  
InputError::InputError(const InputError& ref) throw():
  SCException(ref),
  keyword_(ref.keyword_)
{
}

InputError::~InputError() throw()
{
  delete[] value_;
}

////////////////////////////////////////////////////////////////////////
// ProgrammingError

ProgrammingError::ProgrammingError(
    const char *description,
    const char *file,
    int line,
    const ClassDesc *class_desc,
    const char *exception_type) throw():
  SCException(description, file, line, class_desc, exception_type)
{
}
  
ProgrammingError::ProgrammingError(const ProgrammingError& ref) throw():
  SCException(ref)
{
}

ProgrammingError::~ProgrammingError() throw()
{
}

////////////////////////////////////////////////////////////////////////
// FeatureNotImplemented

FeatureNotImplemented::FeatureNotImplemented(
    const char *description,
    const char *file,
    int line,
    const ClassDesc *class_desc,
    const char *exception_type) throw():
  ProgrammingError(description, file, line, class_desc, exception_type)
{
}
  
FeatureNotImplemented::FeatureNotImplemented(const FeatureNotImplemented& ref)
    throw():
  ProgrammingError(ref)
{
}

FeatureNotImplemented::~FeatureNotImplemented() throw()
{
}

////////////////////////////////////////////////////////////////////////
// SystemException

SystemException::SystemException(
    const char *description,
    const char *file,
    int line,
    const ClassDesc *class_desc,
    const char *exception_type) throw():
  SCException(description, file, line, class_desc, exception_type)
{
}
  
SystemException::SystemException(const SystemException& ref) throw():
  SCException(ref)
{
}

SystemException::~SystemException() throw()
{
}

////////////////////////////////////////////////////////////////////////
// Memory Allocation Failure

MemAllocFailed::MemAllocFailed(const char *description,
                               const char *file,
                               int line,
                               size_t nbyte,
                               const ClassDesc *class_desc,
                               const char *exception_type) throw():
  SystemException(description, file, line, class_desc, exception_type),
  nbyte_(nbyte)
{ 
  try {
      if (nbyte_) {
          elaborate() << "nbyte:       "
                      << nbyte
                      << std::endl;
        }
    }
  catch(...) {
    }
}

MemAllocFailed::MemAllocFailed(const MemAllocFailed& ref) throw():
  SystemException(ref), nbyte_(ref.nbyte_)
{ 
}

MemAllocFailed::~MemAllocFailed() throw()
{
}

////////////////////////////////////////////////////////////////////////
// File Operation Failure

FileOperationFailed::FileOperationFailed(const char *description,
                                         const char *file,
                                         int line,
                                         const char *filename,
                                         FileOperation op,
                                         const ClassDesc *class_desc,
                                         const char *exception_type) throw():
  SystemException(description, file, line, class_desc, exception_type),
  filename_(filename),
  operation_(op)
{ 
  try {
      if (filename_) {
          elaborate() << "file name:   "
                      << filename_
                      << std::endl;
        }
      elaborate() << "file op:     ";
      switch (operation_) {
        case Unknown:
            elaborate() << "Unknown";
            break;
        case OpenR:
            elaborate() << "OpenR";
            break;
        case OpenW:
            elaborate() << "OpenW";
            break;
        case OpenRW:
            elaborate() << "OpenRW";
            break;
        case Close:
            elaborate() << "Close";
            break;
        case Read:
            elaborate() << "Read";
            break;
        case Write:
            elaborate() << "Write";
            break;
        case Corrupt:
            elaborate() << "Corrupt";
            break;
        case Other:
            elaborate() << "Other";
            break;
        default:
            elaborate() << "Invalid";
        }
      elaborate() << std::endl;
    }
  catch(...) {
    }
}

FileOperationFailed::FileOperationFailed(const FileOperationFailed& ref) throw():
  SystemException(ref), filename_(ref.filename_), operation_(ref.operation_)
{ 
}

FileOperationFailed::~FileOperationFailed() throw()
{
}

////////////////////////////////////////////////////////////////////////
// Syscall Failure

SyscallFailed::SyscallFailed(const char *description,
                             const char *file,
                             int line,
                             const char *syscall,
                             int err,
                             const ClassDesc *class_desc,
                             const char *exception_type) throw():
  SystemException(description, file, line, class_desc, exception_type),
  syscall_(syscall),
  err_(err)
{ 
  try {
      if (err_ == 0) {
          err_ = errno;
        }
      if (syscall_) {
          elaborate() << "system call: "
                      << syscall_
                      << std::endl;
        }
      elaborate() << "error:       "
                  << strerror(err_)
                  << " (" << err_ << ")"
                  << std::endl;
    }
  catch(...) {
    }
}

SyscallFailed::SyscallFailed(const SyscallFailed& ref) throw():
  SystemException(ref), syscall_(ref.syscall_), err_(ref.err_)
{ 
}

SyscallFailed::~SyscallFailed() throw()
{
}

////////////////////////////////////////////////////////////////////////
// AlgorithmException

AlgorithmException::AlgorithmException(
    const char *description,
    const char *file,
    int line,
    const ClassDesc *class_desc,
    const char *exception_type) throw():
  SCException(description, file, line, class_desc, exception_type)
{
}
  
AlgorithmException::AlgorithmException(const AlgorithmException& ref) throw():
  SCException(ref)
{
}

AlgorithmException::~AlgorithmException() throw()
{
}

////////////////////////////////////////////////////////////////////////
// MaxIterExceeded

MaxIterExceeded::MaxIterExceeded(const char *description,
                                 const char *file,
                                 int line,
                                 int maxiter,
                                 const ClassDesc *class_desc,
                                 const char *exception_type) throw():
  AlgorithmException(description, file, line, class_desc, exception_type),
  max_iter_(maxiter)
{ 
  try {
      elaborate() << "max_iter:    "
                  << maxiter
                  << std::endl;
    }
  catch(...) {
    }
}

MaxIterExceeded::MaxIterExceeded(const MaxIterExceeded& ref) throw():
  AlgorithmException(ref), max_iter_(ref.max_iter_)
{ 
}

MaxIterExceeded::~MaxIterExceeded() throw()
{
}

//////////////////////////////////////////////////////////////////////
// ToleranceExceeded

ToleranceExceeded::ToleranceExceeded(const char *description,
                                     const char *file,
                                     int line,
                                     double tol,
                                     double val,
                                     const ClassDesc *class_desc,
                                     const char *exception_type) throw():
  AlgorithmException(description, file, line, class_desc, exception_type),
  tolerance_(tol), value_(val)
{
  try {
      elaborate() << "value:       " << value_
                  << std::endl
                  << "tolerance:   " << tolerance_
                  << std::endl;
    }
  catch(...) {
    }
}

ToleranceExceeded::ToleranceExceeded(const ToleranceExceeded& ref) throw():
  AlgorithmException(ref),
  tolerance_(ref.tolerance_), value_(ref.value_)
{
}

ToleranceExceeded::~ToleranceExceeded() throw()
{
}
