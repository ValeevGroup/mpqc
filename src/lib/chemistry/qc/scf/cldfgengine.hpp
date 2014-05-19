//
// cldfgfactory.hpp
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

#ifndef MPQC_CHEMISTRY_QC_SCF_CLDFGFACTORY_HPP
#define MPQC_CHEMISTRY_QC_SCF_CLDFGFACTORY_HPP

#include <chemistry/qc/scf/gengine_base.hpp>
#include <vector>
#include <tiledarray.h>
#include <util/class/class.h>

// Foward Declarations
namespace sc {
  class IntegralLibint2;
} // namespace sc

namespace mpqc{
  class World;

  namespace TA{
    class TiledBasisSet;
  } // namespace TA
} // namespace mpqc



namespace mpqc {
  namespace TA {

    class ClDFGEngine: public GEngineBase {
    public:
      typedef GEngineBase::TAMatrix TAMatrix;

      ClDFGEngine(sc::Ref<sc::IntegralLibint2> integral,
                   sc::Ref<TiledBasisSet> basis,
                   sc::Ref<TiledBasisSet> dfbasis,
                   TAMatrix *density,
                   sc::Ref<World> world);

      ClDFGEngine(const sc::Ref<sc::KeyVal> &kv);

      // No defaults copies or assignments.
      ClDFGEngine() = delete;
      ClDFGEngine operator=(const ClDFGEngine&) = delete;
      ClDFGEngine(const ClDFGEngine&) = delete;

      virtual ~ClDFGEngine() override = default;

      virtual TAMatrix
      operator()(const std::string v) override final;


      virtual void
      set_densities(std::vector<TAMatrix *>) override final;

      virtual bool
      densities_set() override final;

      virtual void
      set_coefficients(std::vector<TAMatrix*>) override final;

      virtual bool
      coefficients_set() override final;

      virtual bool
      using_coeff() override final;

      virtual bool
      using_density() override final;

    private:
      // Form G using coefficents
      TAMatrix
      coefficient_contraction(const std::vector<std::string> v);

      // Form G using density
      TAMatrix
      density_contraction(const std::vector<std::string> v);

      // Compute integrals needed for contraction
      void compute_symetric_df_ints();
      // integral object
      sc::Ref<sc::IntegralLibint2> integral_;

      // Basis sets
      sc::Ref<TiledBasisSet> basis_;
      sc::Ref<TiledBasisSet> dfbasis_;

      // Ptrs to density matrix or coeff matrix and our world object
      TAMatrix *density_ = nullptr;
      TAMatrix *coeff_ = nullptr;
      sc::Ref<World> world_;

      // Bools for which method we are using
      bool density_set_ = false;
      bool coeff_set_ = false;

      // Tensor that holds the integrals which have been combined with the
      // sqrt inverse of the two body two center integrals.
      TiledArray::Array<double, 3> df_ints_;
      TiledArray::Array<double, 3> df_K_; // Holds exchange intermediate

      static sc::ClassDesc class_desc_;
    };

  } // namespace TA
} // namespace mpqc

#endif /* MPQC_CHEMISTRY_QC_SCF_CLDFGFACTORY_HPP */
