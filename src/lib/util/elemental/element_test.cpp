//
// element_test.cpp
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

#include <util/misc/regtime.h>
#include <iostream>
#include <util/elemental/grid.hpp>
#include <util/elemental/init.hpp>
#include <elemental.hpp>


using namespace sc;
using namespace mpqc;

int main(int argc, char** argv){

  ExEnv::init(argc, argv);
  ELEMETNALRuntime::initialize();
  mpqc::Grid* grid = new mpqc::Grid();
  elem::DistMatrix<double> mat(16,16, *grid->elemGrid());
  elem::Identity(mat, 16, 16);
  elem::Print(mat);

}

