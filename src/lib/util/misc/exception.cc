//
// exception.cc
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

#include <util/misc/exception.h>
#include <util/misc/exenv.h>

using namespace std;
using namespace sc;

////////////////////////////////////////////////////////////////////////
// Exception

Exception::Exception(const char *description,
                     const char *file,
                     int line) MPQC__NOEXCEPT:
  description_(description),
  file_(file),
  line_(line)
{
}

Exception::Exception(const Exception& ref) MPQC__NOEXCEPT:
    description_(ref.description_),
    file_(ref.file_),
    line_(ref.line_)
{
}

Exception::~Exception() MPQC__NOEXCEPT
{
  try{ ExEnv::out0().flush(); ExEnv::err0().flush(); }
  catch(...) {}
}

const char* 
Exception::what() const MPQC__NOEXCEPT
{
  try {
      std::ostringstream oss;
      if (description_) {
        oss << "Exception: " << description_ << std::endl;
      }
      if (file_) {
        oss   << "Exception: location = " << file_ << ":" << line_ << std::endl;
      }
      if (description_ || file_)
        return oss.str().c_str();
      else
        return "";
    }
  catch (...) {}

  return "No information available for Exception";
}
