//
// eigensolver.cpp
//
// Copyright (C) 2014 Drew Lewis
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

#include <elemental.hpp>
#include <math/elemental/eigensolver.hpp>
#include <tiled_array.h>


namespace mpqc {
namespace TA {

  using TAMatrix = TiledArray::Array<double, 2, TiledArray::Tensor<double> >;
  using EMatrix = elem::DistMatrix<double>;

  void
  eigensolver_D(const TAMatrix &F,
                const TAMatrix &S,
                TAMatrix &D,
                int occ){

    // Get mats from TA
    EMatrix EF = TiledArray::array_to_elem(F, elem::DefaultGrid()); // Fock
    EMatrix ES = TiledArray::array_to_elem(S, elem::DefaultGrid()); // Overlap

    // Vec and value storage
    elem::DistMatrix<double , elem::VR, elem::STAR> vals(elem::DefaultGrid());
    EMatrix vecs(elem::DefaultGrid());

    // Compute vecs and values
    elem::mpi::Barrier(elem::mpi::COMM_WORLD);
    double t0 = madness::wall_time();
    // Compute only evals needed to form density.
    elem::HermitianGenDefiniteEig(elem::AXBX,elem::LOWER,EF,ES,vals,vecs,
                                  int(0), // Get evals between most negative
                                  occ-1, //  and the occupation.
                                  elem::ASCENDING);
    elem::mpi::Barrier(elem::mpi::COMM_WORLD);
    double t1 = madness::wall_time();

    // Timer
    if(F.get_world().rank()==0){
      std::cout << "\tEigen Solve took " << t1 - t0 << " s" << std::endl;
    }

    // Make Density Matrix
    EMatrix Density(elem::DefaultGrid());
    elem::Gemm(elem::NORMAL, elem::TRANSPOSE, 1.0, vecs, vecs, Density);

    // Copy back to TA
    TiledArray::elem_to_array(D,Density);
    F.get_world().gop.fence();
  }


  void
  eigensolver_Coeff(const TAMatrix &F,
                const TAMatrix &S,
                TAMatrix &C){

    // Get mats from TA
    EMatrix EF = TiledArray::array_to_elem(F, elem::DefaultGrid()); // Fock
    EMatrix ES = TiledArray::array_to_elem(S, elem::DefaultGrid()); // Overlap

    // Vec and value storage
    elem::DistMatrix<double , elem::VR, elem::STAR> vals(elem::DefaultGrid());
    EMatrix vecs(elem::DefaultGrid());

    // Compute vecs and values
    elem::mpi::Barrier(elem::mpi::COMM_WORLD);
    double t0 = madness::wall_time();
    // Compute only evals needed to form density.
    elem::HermitianGenDefiniteEig(elem::AXBX,elem::LOWER,EF,ES,vals,vecs,
                                  elem::ASCENDING);
    elem::mpi::Barrier(elem::mpi::COMM_WORLD);
    double t1 = madness::wall_time();

    // Timer
    if(F.get_world().rank()==0){
      std::cout << "\tEigen Solve took " << t1 - t0 << " s" << std::endl;
    }

    // Copy back to TA
    TiledArray::elem_to_array(C,vecs);
    F.get_world().gop.fence();
  }

} // namespace mpqc::TA
} // namespace mpqc


