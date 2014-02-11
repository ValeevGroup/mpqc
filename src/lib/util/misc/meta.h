//
// meta.h
//
// Copyright (C) 2014 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: Feb 11, 2014
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

#ifndef _util_misc_meta_h
#define _util_misc_meta_h

// Boost includes
#include <boost/mpl/back.hpp>
#include <boost/mpl/pop_back.hpp>
#include <boost/mpl/empty.hpp>

namespace sc { namespace meta {

namespace mpl = boost::mpl;

////////////////////////////////////////////////////////////////////////////////
// splat:  a metafunction for "splatting" a boost::mpl::vector (or any other
//   boost::mpl Front Extensible Sequence) into variadic arguments of another
//   metafunction or type.


namespace { // anonymous namespace to hide splat_helper from outside world

  template<template<typename...> class T, bool empty, typename C, typename... Types>
  struct splat_helper {
      typedef typename mpl::pop_back<C>::type rest;
      typedef typename splat_helper<
          T, mpl::empty<rest>::value, rest, typename mpl::back<C>::type, Types...
      >::type type;
  };

  template<template<typename...> class T, typename C, typename... Types>
  struct splat_helper<T, true, C, Types...> {
      typedef T<Types...> type;
  };

  template<typename value_type, template<value_type...> class T, bool empty, typename C, value_type... values>
  struct splat_helper_values {
      typedef typename mpl::pop_back<C>::type rest;
      typedef typename splat_helper_values<
          value_type, T, mpl::empty<rest>::value, rest, mpl::back<C>::type::value, values...
      >::type type;
  };

  template<typename value_type, template<value_type...> class T, typename C, value_type... values>
  struct splat_helper_values<value_type, T, true, C, values...> {
      typedef T<values...> type;
  };

} // end anonymous namespace

template<template<typename...> class T, typename C>
struct splat_types
  : splat_helper<T, mpl::empty<C>::value, C>
{ };

// Container contains types with a value attribute that is an int
template<typename value_type, template<value_type...> class T, typename C>
struct splat_values
  : splat_helper_values<value_type, T, mpl::empty<C>::value, C>
{ };



////////////////////////////////////////////////////////////////////////////////

} } // end namespace sc::meta

#endif /* _util_misc_meta_h */
