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

#include <chemistry/qc/scf/gfactory.hpp>
#include <util/madness/world.h>
#include <chemistry/qc/libint2/libint2.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/taskintegrals.hpp>
#include <chemistry/qc/basis/integralenginepool.hpp>

namespace mpqc {
  namespace TA {

    class ClDfGFactory : public GFactory {
    public:
      typedef GFactory::TAMatrix TAMatrix;
      ClDfGFactory(sc::Ref<sc::IntegralLibint2> integral,
                   sc::Ref<TiledBasisSet> basis,
                   sc::Ref<TiledBasisSet> dfbasis,
                   const TAMatrix &density,
                   sc::Ref<World> world);
      virtual ~ClDfGFactory() override = default;

      virtual TiledArray::expressions::TensorExpression<TAMatrix::eval_type>
      operator()(const std::string &v) override;

    private:
      void compute_symetric_df_ints();
      // integral object
      sc::Ref<sc::IntegralLibint2> integral_;

      // Basis sets
      sc::Ref<TiledBasisSet> basis_;
      sc::Ref<TiledBasisSet> dfbasis_;

      // Reference to density matrix and our world object
      const TAMatrix &density_;
      sc::Ref<World> world_;

      // Tensor that holds the integrals which have been combined with the
      // sqrt inverse of the two body two center integrals.
      TiledArray::Array<double, 3> df_ints_;

      static sc::ClassDesc class_desc_;
    };

  } // namespace TA
} // namespace mpqc

#endif /* MPQC_CHEMISTRY_QC_SCF_CLDFGFACTORY_HPP */
