//
// linkage.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifndef _chemistry_molecule_linkage_h
#define _chemistry_molecule_linkage_h

#include <chemistry/molecule/coor.h>
#include <chemistry/molecule/taylor.h>
#include <chemistry/molecule/molfreq.h>
#include <chemistry/molecule/molrender.h>
#include <chemistry/molecule/molshape.h>
#include <chemistry/molecule/findisp.h>
#include <chemistry/molecule/frag.h>
#include <chemistry/molecule/molden.h>

#include <util/render/linkage.h>
#include <math/scmat/linkage.h>
#include <math/optimize/linkage.h>

namespace sc {

ForceLink<RedundMolecularCoor> molecule_force_link_a_;
ForceLink<CartMolecularCoor> molecule_force_link_b_;
ForceLink<SymmMolecularCoor> molecule_force_link_c_;
ForceLink<TaylorMolecularEnergy> molecule_force_link_d_;
ForceLink<MolecularFrequencies> molecule_force_link_e_;
ForceLink<RenderedStickMolecule> molecule_force_link_f_;
ForceLink<RenderedBallMolecule> molecule_force_link_g_;
ForceLink<RenderedMolecularSurface> molecule_force_link_h_;
ForceLink<VDWShape> molecule_force_link_i_;
ForceLink<DiscreteConnollyShape> molecule_force_link_j_;
ForceLink<ConnollyShape> molecule_force_link_k_;
ForceLink<FinDispMolecularHessian> molecule_force_link_l_;
ForceLink<FinDispMolecularGradient> molecule_force_link_m_;
ForceLink<MolecularFragment> molecule_force_link_n_;
ForceLink<WriteMolden> molecule_force_link_o_;

}

#endif
