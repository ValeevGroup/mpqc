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
    elem::HermitianGenDefiniteEig(elem::AXBX,elem::LOWER,EF,ES,vals,vecs,
                                  elem::ASCENDING);
    elem::mpi::Barrier(elem::mpi::COMM_WORLD);
    double t1 = madness::wall_time();

    // Timer
    if(F.get_world().rank()==0){
      std::cout << "\tEigen Solve took " << t1 - t0 << " s" << std::endl;
    }

    // Make Density Matrix and get view of vecs corresponding to occupied
    EMatrix Density(elem::DefaultGrid());
    EMatrix C = elem::View(vecs, 0, 0, vecs.Height(), occ);
    elem::Gemm(elem::NORMAL, elem::TRANSPOSE, 1.0, C, C, Density);

    // Copy back to TA
    TiledArray::elem_to_array(D,Density);
    F.get_world().gop.fence();
  }

  TAMatrix
  eigensolver_occ_Coeff(const TAMatrix &F,
                const TAMatrix &S,
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
    elem::HermitianGenDefiniteEig(elem::AXBX,elem::LOWER,EF,ES,vals,vecs,
                                  elem::ASCENDING);
    elem::mpi::Barrier(elem::mpi::COMM_WORLD);
    double t1 = madness::wall_time();

    // Timer
    if(F.get_world().rank()==0){
      std::cout << "\tEigen Solve took " << t1 - t0 << " s" << std::endl;
    }

    EMatrix C = elem::View(vecs, 0, 0,  vecs.Height(), occ);

    EMatrix CT;
    elem::Transpose(C,CT);

    // Get range for tall side
    TiledArray::TiledRange1 height_t1 = F.trange().data()[0];

    // How many tiles over does the occupied vector sit.
    std::array<unsigned long, 2>  occ_pos{{0,static_cast<unsigned long>(occ-1)}};
    auto ntile_guess = F.trange().element_to_tile(occ_pos)[1] + 1;
    auto size_guess = occ/ntile_guess; // How big should tiles be

    // Get range for short side
    std::vector<std::size_t> short_side_blocksize;

    for(auto i = 0; i < occ; i += size_guess){
      short_side_blocksize.push_back(i);
    }
    short_side_blocksize.push_back(occ);

    TiledArray::TiledRange1 trange1_c(short_side_blocksize.begin(),
                                      short_side_blocksize.end());
    std::array<TiledArray::TiledRange1,2> trange_array{{trange1_c, height_t1}};

    // Make trange for C
    TiledArray::TiledRange trange_c(trange_array.begin(), trange_array.end());

    // Make C Array
    TAMatrix TA_C(F.get_world(), trange_c);
    TA_C.set_all_local(0.0);
    F.get_world().gop.fence();

    // Copy back to TA
    TiledArray::elem_to_array(TA_C, CT);
    F.get_world().gop.fence();
    return TA_C;
  }

  void
  eigensolver_full_Coeff(const TAMatrix &F,
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


