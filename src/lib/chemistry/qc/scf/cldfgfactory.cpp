//
// cldfgfactory.cpp
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

#include <chemistry/qc/scf/cldfgfactory.hpp>
#include <elemental.hpp>

using namespace mpqc;
using namespace mpqc::TA;
using namespace sc;

using TAMatrix = ClDfGFactory::TAMatrix;
using EMatrix = elem::DistMatrix<double>;

mpqc::TA::ClDfGFactory::ClDfGFactory(sc::Ref<sc::IntegralLibint2> integral,
                                     sc::Ref<TiledBasisSet> basis,
                                     sc::Ref<TiledBasisSet> dfbasis,
                                     const TAMatrix& density,
                                     sc::Ref<World> world) :
        integral_(integral), basis_(basis), dfbasis_(dfbasis),
        density_(density), world_(world)
{}

TiledArray::expressions::TensorExpression<TAMatrix::eval_type>
mpqc::TA::ClDfGFactory::operator ()( const std::string& v) {

  // Get the user input
  std::vector<std::string> input;
  std::istringstream is(v);
  std::string index;
  while(std::getline(is, index, ',')){ // must use ',' for char literal
    input.push_back(index);
  }

  if(input.size()!=2){
    std::cout << "input size 2 != " << input.size() << std::endl;
    throw std::invalid_argument("size of input to mpqc::TA::ClDfGFactory"
                                " was not equal to 2");
  }

  // Check to see if we computed ints yet.
  if (!df_ints_.is_initialized()) {
    compute_symetric_df_ints();
    world_->madworld()->gop.fence();
  }

 /*
  * Construct strings which will be used to generate the expressions comma for seperation
  * all sequence must start with either i or m hence they don't have commas
  */
  const std::string &i = input[0];
  const std::string &j = "," + input[1];
 /*
  * These strings use random sequence which are unlikely to have collisions.
  * This means that possibly if someone calls this function in an expression and
  * happens to use the same strings as one of the ones below something bad might
  * happen. This is unlikely due to the length and random nature of these strings.
  */
  const std::string m("qgmtkxexut");
  const std::string n(",ihgtpakqif");
  const std::string X(",ssaymsoids");

  auto expr = 2 * (df_ints_(i+j+X) * (density_(m+n) * df_ints_(m+n+X)))
                - (df_ints_(i+n+X) * (density_(m+n) * df_ints_(m+j+X)));

  return expr;
}

void mpqc::TA::ClDfGFactory::compute_symetric_df_ints() {

  // Get three center ints
  using eri3pool = IntegralEnginePool<sc::Ref<sc::TwoBodyThreeCenterInt> >;
  integral_->set_basis(basis_, basis_, dfbasis_);
  std::shared_ptr<eri3pool>
    eri3_ptr(new eri3pool(integral_->electron_repulsion3()->clone()));

  // Using the df_ints as temporary storage for the twobody three center ints
  df_ints_ =  Integrals(*world_->madworld(), eri3_ptr, basis_, dfbasis_);


  // Get two center ints
  using eri2pool = IntegralEnginePool<sc::Ref<sc::TwoBodyTwoCenterInt> >;
  integral_->set_basis(dfbasis_, dfbasis_);
  std::shared_ptr<eri2pool>
    eri2_ptr(new eri2pool(integral_->electron_repulsion2()->clone()));

  TAMatrix eri2_ints = Integrals(*world_->madworld(), eri2_ptr, dfbasis_);

  // Copy two body two center ints to elemental
  EMatrix eri2_elem =
      TiledArray::array_to_elem(eri2_ints, elem::DefaultGrid());
  world_->madworld()->gop.fence();

#if 0 // Inverse sqrt option disabled for now.
  // Perform the sqrt inverse of the twobody two center integrals
  // Eigen vectors and values
  EMatrix vectors(elem::DefaultGrid());
  elem::DistMatrix<double, elem::VR, elem::STAR> values(elem::DefaultGrid());

  // Eigensolver uses eri2_elem as storage so that matrix can't be trusted anymore
  elem::HermitianEig(elem::LOWER, eri2_elem, values, vectors);

  // Take the sqrt of the inverse of the eigenvalues
  auto local_size = values.LocalHeight() * values.LocalWidth();
  std::for_each(values.Buffer(), values.Buffer() + local_size,
                [](double& i){i = 1.0/(sqrt(i));});

  /*
   * Make a copy of the eigenvectors so that we can apply the matrix product
   * \mathbf{U}\mathbf{D} = \mathbf{E}. Then \mathbf{E}\mathbf{U}^T = \mathbf{M}^{-1/2}
   */
  auto E(vectors);
  elem::mpi::Barrier(elem::mpi::COMM_WORLD);
  elem::DiagonalScale(elem::RIGHT, elem::NORMAL, values, E);

  // Make sqrt inverse and put it in  eri2_elem
  elem::Gemm(elem::NORMAL, elem::TRANSPOSE, 1.0, E, vectors, eri2_elem);
#endif

  // Compute the cholesky inverse matrix.
  elem::Cholesky(elem::UPPER, eri2_elem);
  elem::MakeTriangular(elem::UPPER, eri2_elem);
  elem::TriangularInverse(elem::UPPER,elem::NON_UNIT,eri2_elem);

  // Copy back to TA
  ::TiledArray::elem_to_array(eri2_ints, eri2_elem);

  // Create df_ints_ tensor
  df_ints_ = df_ints_("i,j,P") * eri2_ints("P,X");

}

