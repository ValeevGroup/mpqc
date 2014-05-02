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

  // Wraps a pair where the eigenvalues are first and vectors are second.
  // Similar to Mathematica Notation
  using EigenSystem = std::pair<elem::DistMatrix<double, elem::VR, elem::STAR>,
                                EMatrix>;

  /*
   * HermitianGenEigensolver and get occupied vectors are two helper functions
   * to reduce code replication
   */
  // Wrapper for elem::HermitianGenDefiniteEig sorts vals ascending
  EigenSystem
  HermitianGenEigensolver(const TAMatrix &F, const TAMatrix &S){
    // Get mats from TA
    EMatrix EF = TiledArray::array_to_elem(F, elem::DefaultGrid()); // Fock
    EMatrix ES = TiledArray::array_to_elem(S, elem::DefaultGrid()); // Overlap

    // Vec and value storage
    elem::DistMatrix<double , elem::VR, elem::STAR> vals(elem::DefaultGrid());
    EMatrix vecs(elem::DefaultGrid());

    // Compute vecs and values
    elem::mpi::Barrier(elem::mpi::COMM_WORLD);
    double t0 = madness::wall_time();
    // Compute evals evals
    elem::HermitianGenDefiniteEig(elem::AXBX, elem::LOWER,EF,ES,vals,vecs,
                                  elem::ASCENDING);
    elem::mpi::Barrier(elem::mpi::COMM_WORLD);
    double t1 = madness::wall_time();

    // Timer
    if(F.get_world().rank()==0){
      std::cout << "\tCall to HermitianGenEigensolver took " << t1 - t0
                << " s" << std::endl;
    }

    return std::make_pair(vals, vecs);
  }

  // Function takes the full C matrix, a trange1 for the long dimensions, and
  // the number of occupied vectors.
  TAMatrix
  get_occupied_vectors(EMatrix &vecs, TiledArray::TiledRange1 trange1, int occ,
                       madness::World &world){

    // Grab the occupied vectors
    EMatrix C = elem::View(vecs, 0, 0,  vecs.Height(), occ);

    /*
     * Get number of tiles to pack low rank dimension into, rule of thumb will be
     * 1/5 the number of tiles of the full dimension
     * */
    auto ntile_guess = std::max(static_cast<int>(trange1.tiles().second/5), 1);
    auto size_guess = occ/ntile_guess; // How many elements per tile.

    // Get range for short side
    std::vector<std::size_t> short_side_blocksize;

    for(auto i = 0; i < occ; i += size_guess){
      short_side_blocksize.push_back(i);
    }
    short_side_blocksize.push_back(occ);

    // Get the TiledRange1 for the reduced rank dimension.
    TiledArray::TiledRange1 lowrank_t1(short_side_blocksize.begin(),
                                       short_side_blocksize.end());

    // Create array to initialized TiledRange for Coeff. matrix.
    // Will create N rows and occ cols
    std::array<TiledArray::TiledRange1,2> trange_array{{trange1, lowrank_t1}};

    // Make trange for C
    TiledArray::TiledRange trange_c(trange_array.begin(), trange_array.end());

    // Make C Array
    TAMatrix TA_C(world, trange_c);
    TA_C.set_all_local(0.0);
    world.gop.fence();

    // Copy back to TA
    TiledArray::elem_to_array(TA_C, C);
    world.gop.fence();

    return TA_C;
  }


  TAMatrix
  eigensolver_D(const TAMatrix &F, const TAMatrix &S, int occ){

    // Get eigensystem
    EigenSystem esys = HermitianGenEigensolver(F,S);

    // Grab occupied vectors
    TAMatrix C_occ = get_occupied_vectors(esys.second, F.trange().data()[0],
                                          occ, F.get_world());
    // Contract over occupied index i to form D_{AO}
    TAMatrix D = C_occ("mu,") * C_occ("nu,i");
    F.get_world().gop.fence(); // Fence to prevent data from going out of scope
    std::cout << "D = \n" << D << std::endl;

    return D;
  }

  TAMatrix
  eigensolver_occ_Coeff(const TAMatrix &F,
                const TAMatrix &S,
                int occ){
    // Get Eigensystem
    EigenSystem esys = HermitianGenEigensolver(F,S);

    return get_occupied_vectors(esys.second, F.trange().data()[0], occ,
                                F.get_world());
  }

  TAMatrix
  eigensolver_full_Coeff(const TAMatrix &F,
                const TAMatrix &S){

    // Create matrix to return
    TAMatrix C(F.get_world(), F.trange());

    // Call eigenvalue solver
    EigenSystem esys = HermitianGenEigensolver(F,S);

    // Copy back to TA
    TiledArray::elem_to_array(C,esys.second);
    F.get_world().gop.fence();
    return C;
  }



} // namespace mpqc::TA
} // namespace mpqc


