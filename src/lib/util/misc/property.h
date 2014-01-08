//
// property.h
//
// Copyright (C) 2013 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: Dec 13, 2013
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

#ifndef _util_misc_property_h
#define _util_misc_property_h

#include <functional>
#include <utility>

namespace sc {

// Taken from:
//    http://www.daniweb.com/software-development/cpp/threads/443096/c-property-template

/**
 * Helper class to connect a 'property' in a
 * c++ class to getter/setter methods
 * \tparam ClassType   The type of the class in which the property is found.
 * \tparam MemberType  The type of the property value.
 */
template <typename ClassType, typename MemberType>
class property {
public:
  typedef property<ClassType, MemberType> self; // good to have a "self" type.
  typedef MemberType (ClassType::*getter_type)() const;
  typedef void (ClassType::*setter_type)(const MemberType&);

  /**
   * \brief ctor
   * NOTE: Avoid using leading underscores, technically, the C++ standard forbids them.
   */
  property(ClassType* aTarget,
           getter_type aGetter = NULL,
           setter_type aSetter = NULL) :
           target(aTarget),
           getter(aGetter),
           setter(aSetter) {
    if(!getter)
      throw std::logic_error("Cannot have a NULL getter function for a property!");
  };

  // Delete the default copy and move constructors, since they will make this refer to an old instance
  property(self&& other) : target(std::move(other.target)), getter(std::move(other.getter)), setter(std::move(other.setter)) { }
  self& operator=(self&) = delete;
  self& operator=(const self&) = delete;

  /**
   * assignment operator.  Provides setter
   *  parent.field = value
   * functionality
   */
  self& operator=(const MemberType& value) {
    if(setter)  // for non-nullity of setter pointer.
      (target->*setter)(value);
    return *this;
  };

  /**
   * conversion operator, provides 'getter'
   */
  operator MemberType() const {
    return (target->*getter)();
  };

  void set_target(ClassType* a) { target = a; }
  void set_getter(getter_type g){ getter = g; }
  void set_setter(setter_type s){ setter = s; }
  void set_target_getter_setter(ClassType* a, getter_type g, setter_type s = NULL) {
    target = a;
    getter = g;
    setter = s;
  }

private:
  ClassType* target;
  getter_type getter;
  setter_type setter;

};



} // end namespace sc

#endif /* end _util_misc_property_h */
