//
// taylor.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
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

#ifndef _chemistry_molecule_taylor_h
#define _chemistry_molecule_taylor_h

#include <chemistry/molecule/energy.h>
#include <chemistry/molecule/coor.h>

namespace sc {

  /// @addtogroup ChemistryMolecule
  /// @{

// the molecular energy as a taylor expansion
class TaylorMolecularEnergy: public MolecularEnergy {
  private:
    // the coordinates
    Ref<SetIntCoor> coordinates_;

    // The force constants (only the unique ones are given) to arbitrary
    // order.  If nonunique force constants are put here, then the answer
    // will be wrong
    std::vector<std::vector<int> > force_constant_index_;
    std::vector<double> force_constant_value_;
    
    // the dimension of coordinates_;
    RefSCDimension dim_;

    // the expansion point
    RefSCVector expansion_point_;

    // the energy at the expansion point
    double e0_;

    // the maximum order derivative that can be computed
    int maxorder_;

  protected:
    bool analytic_gradient_implemented() const;
    bool analytic_hessian_implemented() const;

  public:
    TaylorMolecularEnergy(const Ref<KeyVal>&);
    TaylorMolecularEnergy(StateIn&);
    ~TaylorMolecularEnergy();
    void save_data_state(StateOut&);
    void print(std::ostream& = ExEnv::out0()) const;
    void compute();
    int value_implemented() const;
};

/// @}
// end of addtogroup ChemistryMolecule

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
