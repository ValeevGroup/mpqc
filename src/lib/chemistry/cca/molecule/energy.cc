//
// energy.cc
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

#include <stdlib.h>
#include <math.h>
#include <stdexcept>
#include <util/misc/string.h>

#include <util/misc/scexception.h>
#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <math/scmat/local.h>
#include <util/keyval/keyval.h>
#include <chemistry/cca/molecule/energy.h>
#include <gov_cca.hxx>
#include <ChemistryCXX_Molecule.hxx>

using namespace std;
using namespace sc;

static ClassDesc MolecularEnergyCCA_cd(
  typeid(MolecularEnergyCCA),"MolecularEnergyCCA",1,"public MolecularEnergy",
  0, create<MolecularEnergyCCA>, create<MolecularEnergyCCA>);

MolecularEnergyCCA::MolecularEnergyCCA(const MolecularEnergyCCA& mole):
  MolecularEnergy(mole)
{

}

MolecularEnergyCCA::MolecularEnergyCCA(const Ref<KeyVal>& keyval):
  MolecularEnergy(keyval)
{
  
  coor_model_name_ = keyval->stringvalue("coordinate_model");
  factory_name_ = keyval->stringvalue("factory");
  init_model(); 

}

MolecularEnergyCCA::~MolecularEnergyCCA()
{

}

MolecularEnergyCCA::MolecularEnergyCCA(StateIn&s):
  MolecularEnergy(s)
{

}

void
MolecularEnergyCCA::save_data_state(StateOut&s)
{

}

void
MolecularEnergyCCA::init_model()
{
  // grab cca environment
  gov::cca::Services &services = *CCAEnv::get_services();
  gov::cca::ports::BuilderService &bs = *CCAEnv::get_builder_service();
  gov::cca::TypeMap &type_map = *CCAEnv::get_type_map();
  gov::cca::ComponentID &my_id = *CCAEnv::get_component_id();

  try{
    // coor_model already exists and registered in external framework case
    coor_model_ = sidl::babel_cast<ChemistryOpt::CoordinateModelInterface>(
      services.getPort("CoordinateModelInterface") );
  }
  catch(...) {
    // no coor model yet, must be embedded framework
  }

  if( !coor_model_ ) {

    // this happens in the embedded framework case
    services.registerUsesPort("CoorModelPort", "ChemsitryOpt.CoordinateModelInterface", type_map );
   
    cca_molecule_ = ChemistryCXX::Molecule::_create();
    cca_molecule_.initialize(molecule()->natom(),0,"bohr");
    for( int i=0; i<molecule()->natom(); ++i) {
      cca_molecule_.set_atomic_number(i,molecule()->Z(i));
      for( int j=0; j<3; ++j)
        cca_molecule_.set_cart_coor( i, j, molecule()->r(i,j) );
    }

    Chemistry::QC::ModelFactoryInterface fac;
    gov::cca::ConnectionID tcon;
    services.registerUsesPort( "ModelFactoryPort", "Chemistry.QC.ModelFactoryInterface", type_map );
    fac_id_ = bs.createInstance( "Factory", factory_name_, type_map ); 
    tcon = bs.connect( my_id, "ModelFactoryPort", fac_id_, "ModelFactoryInterface" );
    fac = sidl::babel_cast<Chemistry::QC::ModelFactoryInterface> (
      services.getPort("ModelFactoryPort") );
    fac.set_molecule( cca_molecule_ );
    services.releasePort( "ModelFactoryPort" );
    bs.disconnect( tcon, 0 );
    services.unregisterUsesPort( "ModelFactoryPort" );

    cm_id_ = bs.createInstance( "CoorModel", coor_model_name_, type_map );
    con1_ = bs.connect( my_id, "CoorModelPort", cm_id_,"CoordinateModelInterface" );
    coor_model_ = sidl::babel_cast<ChemistryOpt::CoordinateModelInterface> (
      services.getPort("CoorModelPort") );
    con2_ = bs.connect( cm_id_, "ModelFactoryInterface", fac_id_, "ModelFactoryInterface" );
    coor_model_.initialize();
  }

  scdim_ = new SCDimension( coor_model_.get_n_coor() );
  sidlx_ = coor_model_.get_coor();
  sidlg_ = sidl::array<double>::create1d(scdim_.n());
  RefSCVector v( scdim_, matrixkit_ );
  for( int i=0; i<scdim_.n(); ++i )
    v[i] = sidlx_.get(i);
  MolecularEnergy::set_x( v );

}

void
MolecularEnergyCCA::set_x(const RefSCVector& v)
{
  MolecularEnergy::set_x( v );

  for( int i=0; i<coor_model_.get_n_coor(); ++i )
    sidlx_.set(i,v[i]); 

}

void
MolecularEnergyCCA::compute()
{
  double energy;
  RefSCVector g( scdim_, matrixkit_ );
  coor_model_.get_energy_and_gradient(sidlx_,energy,sidlg_);
  set_energy(energy);
  for( int i=0; i<scdim_.n(); ++i )
    g[i] = sidlg_.get(i);
  set_gradient( g );
}


