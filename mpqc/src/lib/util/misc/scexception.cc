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

#ifdef __GNUC__
#pragma implementation
#endif

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
                         const ClassDesc *class_desc) throw():
  description_(description),
  file_(file),
  line_(line),
  class_desc_(class_desc)
{
  try {
      elaboration_ = new ostringstream;
      elaborate() << "exception: " << description
                  << std::endl
                  << "location: " << file << ":" << line
                  << std::endl;
      if (class_desc_) {
          elaborate() << "throwing class: " << class_desc_->name()
                      << std::endl;
        }
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
  class_desc_(ref.class_desc_)
{
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
}

const char* 
SCException::what() throw()
{
  if (elaboration_) {
      return elaboration_->str().c_str();
    }
  return description_;
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
// AlgorithmException

AlgorithmException::AlgorithmException(
    const char *description,
    const char *file,
    int line,
    const ClassDesc *class_desc) throw():
  SCException(description, file, line, class_desc)
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
                                 const ClassDesc *class_desc) throw():
  AlgorithmException(description, file, line, class_desc),
  max_iter_(maxiter)
{ 
  try {
      elaborate() << "exceeded maximum number of iterations: max_iter = "
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
                                     const ClassDesc *class_desc) throw():
  AlgorithmException(description, file, line, class_desc),
  tolerance_(tol), value_(val)
{
  try {
      elaborate() << "exceeded tolerance: value = " << value_
                  << " tolerance = " << tolerance_
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
