//
// tiledbasisset.hpp
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

#ifndef CHEMISTRY_QC_BASIS_TILEDBASISSET_HPP
#define CHEMISTRY_QC_BASIS_TILEDBASISSET_HPP

#include <tiled_array.h>

#include "basis.h"
#include "gaussbas.h"
#include <chemistry/molecule/atom.h>
#include <chemistry/molecule/molecule.h>

#include <mpqc/utility/foreach.hpp>
#include "shellorder.hpp"
#include "kcluster.hpp"

namespace mpqc {
namespace TA{

  // TODO add reference for k clustering
  /**
   * TiledBasisSet is a GaussianBasisSet in which basis functions/shells are grouped into tiles ("blocks").
   * This implementation uses the k-means clustering algorithm.
   */
  class TiledBasisSet: public sc::GaussianBasisSet {

    public:
      typedef ShellOrder::Shell Shell;
      typedef sc::Ref<sc::GaussianBasisSet> Basis;
      typedef ShellOrder::ShellRange ShellRange;

      /**
       * Constructs a TiledBasisSet from a sc::Keyval object.
       * The full list of keywords
       * that are accepted is below.
       *  <table border="1">
       *
       *  <tr><td><b>%Keyword</b><td><b>Type</b><td><b>Default</b><td><b>Description</b>
       *
       *  <tr><td><tt>basis</tt><td>GaussianBasisSet<td>none<td>mother basis set object. If not given, will try to instantiate
       *                        this object using keywords from the current scope.
       *
       *  <tr><td><tt>ntiles</tt><td>integer<td>1<td>number of tiles to produce. Need to define some heuristics.
       *
       *  </table>
       *
       * @param[in] keyval is a sc::Ref<sc::Keyval> object which contains
       *     information about the molecule and basis set desired.
       */
      TiledBasisSet(const sc::Ref<sc::KeyVal> &keyval);

      /**
       * Constructs a TiledBasisSet from a sc::GaussianBasisSet object.
       */
      TiledBasisSet(const sc::Ref<sc::GaussianBasisSet>& bs,
                    size_t ntiles = 1);

      /**
       * Returns the tiled range object describing the tiling of basis functions in this basis.
       * @return TiledArray::TiledRange1
       */
      TiledArray::TiledRange1 trange1() const;

      TiledBasisSet(sc::StateIn& s);

      void save_data_state(sc::StateOut& s);

      /// Print a detailed description of the basis set.
      virtual void print(std::ostream& = sc::ExEnv::out0()) const override;

    private:
      std::size_t ntiles_;
      ShellRange SRange_;

      static std::string converted_name(const std::string& name);

  }; // class TiledBasisSet

} // namespace TA
} // namespace mpqc

#endif /* CHEMISTRY_QC_BASIS_TILEDBASISSET_HPP */
