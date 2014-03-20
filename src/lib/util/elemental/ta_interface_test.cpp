//
// ta_interface_test.cpp
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

#include <util/elemental/ta_interface.hpp>
#include <util/elemental/grid.hpp>
#include <util/elemental/init.hpp>
#include <util/madness/init.h>
#include <util/madness/world.h>
#include <tiled_array.h>

using namespace sc;
using namespace mpqc;

int main(int argc, char** argv){

  ExEnv::init(argc, argv);

  MADNESSRuntime::initialize();
  ELEMETNALRuntime::initialize();

  Grid* grid = new Grid();
  World* world_ = new World();

  std::vector<unsigned int> blocking;
  unsigned int block_size = 4;
  unsigned int matrix_size = 16;
  for(auto i = 0; i < matrix_size; i += block_size) {
    blocking.push_back(i);
  }
  std::vector<::TiledArray::TiledRange1> blocking2(
          2, ::TiledArray::TiledRange1(blocking.begin(), blocking.end()));
  ::TiledArray::TiledRange trange(blocking2.begin(), blocking2.end());

  ::TiledArray::Array<double, 2> array(*world_->madworld(), trange);
  array.set_all_local(2.0);
  if(world_->madworld()->rank()==0){
    std::cout << array << std::endl;
  }

  elem::DistMatrix<double> mat = array_to_distmat(array, *grid->elemGrid());
  elem::Print(mat,"TA Copy");
  elem::Identity(mat,matrix_size,matrix_size);
  Print(mat,"Idenitity");

  distmat_to_array(array, mat);

  world_->madworld()->gop.fence();

  if(!world_->madworld()->rank()){
    std::cout << "array = \n" << array << std::endl;
  }

  array.get_world().gop.fence();

  ELEMETNALRuntime::finalize();
  MADNESSRuntime::finalize();

  return 0;
}

