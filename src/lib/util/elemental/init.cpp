//
// init.cpp
//
// Copyright (C) 2013 Drew Lewis
//
// Authors: Drew Lewis
// Maintainer: Drew Lewis and Edward Valeev
//
// This file is part of the MPQC Toolkit.
//
// The MPQC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The MPQC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the MPQC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifdef __GNUG__
#pragma implementation
#endif

#include <mpqc_config.h>

#ifdef HAVE_ELEMENTAL
# ifndef GRID_INSTANTIATE_STATIC_TEMPLATES
# define GRID_INSTANTIATE_STATIC_TEMPLATES
# endif
#include <elemental.hpp>
#endif

#include <util/elemental/init.hpp>
#include <util/misc/exenv.h>
#include <util/misc/scexception.h>

bool mpqc::ELEMETNALRuntime::mpqc_initialized_elemental_ = false;

bool mpqc::ELEMETNALRuntime::initialized(){
  return elem::Initialized();
}

void mpqc::ELEMETNALRuntime::initialize() {
  if(not elem::Initialized()){
    elem::Initialize(sc::ExEnv::argc(), sc::ExEnv::argv());
    mpqc_initialized_elemental_ = true;
  } else {
    mpqc_initialized_elemental_ = false;
  }
}

void mpqc::ELEMETNALRuntime::finalize() {
  if(elem::Initialized()){
    if(mpqc_initialized_elemental_){
      elem::Finalize();
    }
    mpqc_initialized_elemental_ = false;
  }
}
