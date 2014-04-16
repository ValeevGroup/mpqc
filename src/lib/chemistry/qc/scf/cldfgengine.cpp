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

#include <chemistry/qc/scf/cldfgengine.hpp>
#include <util/madness/world.h>
#include <chemistry/qc/basis/integralenginepool.hpp>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/libint2/libint2.h>
#include <chemistry/qc/basis/tiledbasisset.hpp>
#include <chemistry/qc/basis/taskintegrals.hpp>
#include <elemental.hpp>

using namespace mpqc;
using namespace mpqc::TA;
using namespace sc;

using TAMatrix = ClDFGEngine::TAMatrix;
using return_type = ClDFGEngine::return_type;
using EMatrix = elem::DistMatrix<double>;

sc::ClassDesc ClDFGEngine::class_desc_(
                typeid(mpqc::TA::ClDFGEngine),
                "TA.ClDFGEngine", 1, "public GEngineBase",
                0, 0, 0);

mpqc::TA::ClDFGEngine::ClDFGEngine(sc::Ref<sc::IntegralLibint2> integral,
                                     sc::Ref<TiledBasisSet> basis,
                                     sc::Ref<TiledBasisSet> dfbasis,
                                     TAMatrix *density,
                                     sc::Ref<World> world) :
        integral_(integral), basis_(basis), dfbasis_(dfbasis),
        density_(density), world_(world), density_set_(true)
{}

mpqc::TA::ClDFGEngine::ClDFGEngine(const sc::Ref<sc::KeyVal> &kv) : integral_(),
        basis_(), dfbasis_(), world_() {

  // Initialize everything
  world_ << kv->describedclassvalue("world");
  basis_ << kv->describedclassvalue("basis");
  dfbasis_ << kv->describedclassvalue("dfbasis");
  integral_ << kv->describedclassvalue("integrals");

  if(world_.null()){
    world_ = new World;
  }

  if(basis_.null()){ // Check that user didn't give guassianbasis
    sc::Ref<sc::GaussianBasisSet> bs;
    bs << kv->describedclassvalue("basis");
    if (bs.null())
      throw sc::InputError("missing \"basis\" keyword",
                           __FILE__, __LINE__, "basis", "", this->class_desc());
    else {
      basis_ = new TiledBasisSet(bs);
    }
  }

  if(dfbasis_.null()){ // First check that it wasn't a GBS
    sc::Ref<sc::GaussianBasisSet> bs;
    bs << kv->describedclassvalue("dfbasis");
    if(!bs.null()){ // If it was a GBS then convert it
      dfbasis_ = new TiledBasisSet(bs);
    }
    else { // If there was no dfbasis then make one.
      sc::Ref<sc::AssignedKeyVal> akv = new AssignedKeyVal;
      akv->assign("molecule", basis_->molecule().pointer());
      akv->assign("name", "cc-pVDZ/JKFIT");
      akv->assign("ntile", 2);
      dfbasis_ = new TiledBasisSet(static_cast<sc::Ref<sc::KeyVal> >(akv));
    }
  }

  if(integral_.null()){
    integral_ = new IntegralLibint2;
  }

}

return_type
mpqc::TA::ClDFGEngine::operator ()( const std::string& v) {

  // Get the user input
  std::vector<std::string> input; // Vector of input strings
  std::istringstream is(v);
  std::string index; // Temp string
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
    double int_start = madness::wall_time();
    compute_symetric_df_ints();
    world_->madworld()->gop.fence();
    double int_end = madness::wall_time();
    if(world_->madworld()->rank()==0){
      std::cout << "\tComputing the rank 3 symmetric tensor took " << int_end - int_start << " seconds" << std::endl;
    }
  }


  // Use coefficients if we have them
  if(coefficients_set()){
    return coefficient_contraction(input);
  }
  else if(densities_set()){
    return density_contraction(input);
  }
  else{
    std::cout << "Got to throw contraction" << std::endl;
    throw;
  }
}

// Do contraction with density
return_type
mpqc::TA::ClDFGEngine::density_contraction(const std::vector<std::string> &input){
   /*
  * Construct strings which will be used to generate the expressions comma for seperation
  * all sequence must start with either i or m hence they don't have commas
  */
  const std::string &i = input.at(0);
  const std::string &j = "," + input.at(1);
 /*
  * These strings are unlikely to have collisions.
  * This means that possibly if someone calls this function in an expression and
  * happens to use the same strings as one of the ones below something bad might
  * happen. This is unlikely due to the length and nature of these strings.
  */
  const std::string m("mpqc::TA::ClDfGFactory_m_density");
  const std::string n(",mpqc::TA::ClDfGFactory_n_density");
  const std::string X(",mpqc::TA::ClDfGFactory_X_density");

  // just for conveience
  const TAMatrix &dens = *density_;

  auto expr = 2 * (df_ints_(i+j+X) * (dens(m+n) * df_ints_(m+n+X)))
                - (df_ints_(i+n+X) * (dens(m+n) * df_ints_(m+j+X)));

  return expr;
}

// Do contraction with coefficients
return_type
mpqc::TA::ClDFGEngine::coefficient_contraction(
        const std::vector<std::string> &input){
  /*
  * Construct strings which will be used to generate the expressions comma for seperation
  * all sequence must start with either i or m hence they don't have commas
  */
  const std::string &i = input.at(0);
  const std::string &j = "," + input.at(1);
 /*
  * These strings are unlikely to have collisions.
  * This means that possibly if someone calls this function in an expression and
  * happens to use the same strings as one of the ones below something bad might
  * happen. This is unlikely due to the length and nature of these strings.
  */
  const std::string m("mpqc::TA::ClDfGFactory_m_coeff");
  const std::string n(",mpqc::TA::ClDfGFactory_n_coeff");
  const std::string X(",mpqc::TA::ClDfGFactory_X_coeff");
  const std::string Z(",mpqc::TA::ClDfGFactory_Z_coeff");

  // Terms where comma's need to be added or removed
  const std::string nE("mpqc::TA::ClDfGFactory_n_coeff");
  const std::string ZE("mpqc::TA::ClDfGFactory_Z_coeff");
  const std::string iE = "," + input.at(0);

  // just for conveience
  const TAMatrix &C = *coeff_;

  // Precompute Exch Term
  df_K_ = C(nE+Z) * df_ints_(nE+iE+X);

  if(densities_set()){
    // Update D
    TAMatrix &D = *density_;
    D("mu, nu") = C("mu,i") * C("nu,i");

    auto expr = 2 * (df_ints_(i+j+X) * (D(m+n) * df_ints_(m+n+X) ) )
                  - (df_K_(ZE+iE+X) * df_K_(ZE+j+X) );
    return expr;
  }
  else {
    auto expr = 2 * (df_ints_(i+j+X) * ( (C(n+Z)) * df_K_(ZE+n+X) ) )
                  - (df_K_(ZE+iE+X) * df_K_(ZE+j+X) );
    return expr;
  }
}

void
mpqc::TA::ClDFGEngine::set_densities(std::vector<TAMatrix *> densities) {
  MPQC_ASSERT(densities.size() == 1);
  density_ = densities.at(0);
  density_set_ = true;
}

bool
mpqc::TA::ClDFGEngine::densities_set() {
  return density_set_;
}

void mpqc::TA::ClDFGEngine::set_coefficients(std::vector<TAMatrix*> coeffs) {
  MPQC_ASSERT(coeffs.size() == 1);
  coeff_ = coeffs.at(0);
  coeff_set_ = true;
}

bool mpqc::TA::ClDFGEngine::coefficients_set() {
  return coeff_set_;
}

bool mpqc::TA::ClDFGEngine::using_coeff() {
  return coeff_set_;
}

// Only use density if coefficients have not been set.
bool mpqc::TA::ClDFGEngine::using_density() {
  return (!coeff_set_ && density_set_);
}

// Compute the integrals
void
mpqc::TA::ClDFGEngine::compute_symetric_df_ints() {

  // Set basis and grab a clone of the engine we need
  mutex::global::lock(); // <<< Begin Critical Section
    integral_->set_basis(basis_, basis_, dfbasis_);
    auto eri3_clone = integral_->electron_repulsion3()->clone();
  mutex::global::unlock(); // <<< End Critical Section

  // Make an integral engine pool out of our clone
  using eri3pool = IntegralEnginePool<sc::Ref<sc::TwoBodyThreeCenterInt> >;
  auto eri3_ptr = std::make_shared<eri3pool>(eri3_clone);

  // Using the df_ints as temporary storage for the twobody three center ints
  double computing_eri3_ints_time0 = madness::wall_time();
  df_ints_ =  Integrals(*world_->madworld(), eri3_ptr, basis_, dfbasis_);
  world_->madworld()->gop.fence();
  double computing_eri3_ints_time1 = madness::wall_time();

  // Get two center ints
  // Set basis and grab a clone of the engine we need
  mutex::global::lock(); // <<< Begin Critical Section
    integral_->set_basis(dfbasis_, dfbasis_);
    auto eri2_clone = integral_->electron_repulsion2()->clone();
  mutex::global::unlock(); // <<< End Critical Section

  using eri2pool = IntegralEnginePool<sc::Ref<sc::TwoBodyTwoCenterInt> >;
  auto eri2_ptr = std::make_shared<eri2pool>(eri2_clone);

  double computing_eri2_ints_time0 = madness::wall_time();
  TAMatrix eri2_ints = Integrals(*world_->madworld(), eri2_ptr, dfbasis_);
  world_->madworld()->gop.fence();
  double computing_eri2_ints_time1 = madness::wall_time();

  double int_time = (computing_eri2_ints_time1 - computing_eri2_ints_time0) +
                    (computing_eri3_ints_time1 - computing_eri3_ints_time0);

  if(world_->madworld()->rank()==0){
    std::cout << "\tTook " << int_time << " s" <<
            " to compute eri3 and eri2 intiegrals." << std::endl;
  }

  // Copy two body two center ints to elemental
  double inverse_time0 = madness::wall_time();
  EMatrix eri2_elem =
      TiledArray::array_to_elem(eri2_ints, elem::DefaultGrid());
  world_->madworld()->gop.fence(); // makesure we finish copy


  // Compute the cholesky inverse matrix.
  elem::Cholesky(elem::LOWER, eri2_elem);
  elem::TriangularInverse(elem::LOWER, elem::NON_UNIT, eri2_elem);
  elem::MakeTriangular(elem::LOWER, eri2_elem);
  elem::mpi::Barrier(elem::mpi::COMM_WORLD);
  double inverse_time1 = madness::wall_time();

  // Copy back to TA
  ::TiledArray::elem_to_array(eri2_ints, eri2_elem);
  world_->madworld()->gop.fence(); // makesure we finish copy
  if(world_->madworld()->rank()==0){
    std::cout << "\tTook " << inverse_time1 - inverse_time0 << " s" <<
            " to compute eri2 inverse." << std::endl;
  }

  double Final_contraction0 = madness::wall_time();
  // Create df_ints_ tensor from eri3(i,j,P) * U_{eri2}^{-1}(P,X)
  df_ints_ = df_ints_("i,j,P") * eri2_ints("X,P");
  world_->madworld()->gop.fence(); // so eri2_ints doesn't go out of scope.
  double Final_contraction1 = madness::wall_time();
  if(world_->madworld()->rank()==0){
    std::cout << "\tTook " << Final_contraction1 - Final_contraction0 << " s" <<
            " to do big contraction." << std::endl;
  }

}

