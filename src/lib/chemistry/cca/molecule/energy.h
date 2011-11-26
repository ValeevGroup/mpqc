//
// energy.h
//
// Copyright (C) 2007 Sandia National Laboratories
//
// Author: Joseph Kenny <jpkenny@sandia.gov>
// Maintainer: Joseph Kenny
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

#ifndef _chemistry_cca_molecule_energy_h
#define _chemistry_cca_molecule_energy_h

#include <iostream>

#include <util/misc/ccaenv.h>
#include <math/optimize/function.h>
#include <chemistry/molecule/energy.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/deriv.h>

#include <Chemistry_QC_ModelFactoryInterface.hxx>
#include <ChemistryOpt_CoordinateModelInterface.hxx>

namespace sc {

/** The MolecularEnergyCCA concrete class inherits from the MolecularEnergy class.
It computes the energy of the molecule as a function of the geometry using CCA components.
The coordinate system used (internal or cartesian) is determined by the underlying coordinate
model component.  */
class MolecularEnergyCCA: public MolecularEnergy {

  private:
    std::string factory_name_;
    std::string coor_model_name_;
    Chemistry::MoleculeInterface cca_molecule_;
    ChemistryOpt::CoordinateModelInterface coor_model_;
    Chemistry::QC::ModelFactoryInterface factory_;
    gov::cca::ComponentID cm_id_;
    gov::cca::ComponentID fac_id_;
    gov::cca::ConnectionID con1_;
    gov::cca::ConnectionID con2_;
    sidl::array<double> sidlx_;
    sidl::array<double> sidlg_;
    RefSCDimension scdim_;
    
  protected:

  public:
    MolecularEnergyCCA(const MolecularEnergyCCA&);
    MolecularEnergyCCA(const Ref<KeyVal>&);
    MolecularEnergyCCA(StateIn&);
    ~MolecularEnergyCCA();

    void save_data_state(StateOut&);
    int value_implemented() const { return 1; }
    void init_model();
    void compute();
    void set_x(const RefSCVector& v);
};

}

#endif
