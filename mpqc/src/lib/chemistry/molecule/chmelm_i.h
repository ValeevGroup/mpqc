//
// chmelm_i.h
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

QUERY_FUNCTION_IMPL(const char* ,name);
QUERY_FUNCTION_IMPL(const char* ,symbol);
QUERY_FUNCTION_IMPL(int   ,number);
QUERY_FUNCTION_IMPL(double,mass);
QUERY_FUNCTION_IMPL(int   ,family);
QUERY_FUNCTION_IMPL(int   ,row);
QUERY_FUNCTION_IMPL(int   ,valence);
QUERY_FUNCTION_IMPL(int   ,melting_pt);
QUERY_FUNCTION_IMPL(int   ,boiling_pt);
QUERY_FUNCTION_IMPL(int   ,first_ion_pt);
QUERY_FUNCTION_CONV_IMPL(double,bragg_radius,ANGSTROMS_TO_AU);
QUERY_FUNCTION_IMPL(double,electronegativity);
QUERY_FUNCTION_IMPL(double,specific_heat);
QUERY_FUNCTION_IMPL(double,density);
QUERY_FUNCTION_CONV_IMPL(double,atomic_radius,ANGSTROMS_TO_AU);
QUERY_FUNCTION_CONV_IMPL(double,vdw_radius,ANGSTROMS_TO_AU);

INLINE double
ChemicalElement::charge() const
{
  return (double) number();
}

#undef INLINE
