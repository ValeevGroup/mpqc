//
// gaussbas_i.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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

#ifdef INLINE_FUNCTIONS
#define INLINE inline
#else
#define INLINE
#endif

INLINE const char*
GaussianBasisSet::name() const
{
  return name_;
}

INLINE int
GaussianBasisSet::nprimitive() const
{
  return nprim_;
}

INLINE int
GaussianBasisSet::nshell() const
{
  return nshell_;
}

INLINE int
GaussianBasisSet::nbasis() const
{
  return nbasis_;
}

INLINE int
GaussianBasisSet::shell_to_function(int i) const
{
  return shell_to_function_(i);
}

INLINE const GaussianShell&
GaussianBasisSet::operator()(int i) const
{
  return *shell_[i];
}

INLINE const GaussianShell&
GaussianBasisSet::operator[](int i) const
{
  return *shell_[i];
}

INLINE GaussianShell&
GaussianBasisSet::operator()(int i)
{
  return *shell_[i];
}

INLINE GaussianShell&
GaussianBasisSet::operator[](int i)
{
  return *shell_[i];
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
