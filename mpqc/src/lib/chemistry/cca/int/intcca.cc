//
// intcca.cc
//
// Copyright (C) 2004 Sandia National Laboratories
//
// Author: Joe Kenny <jpkenny@sandia.gov>
// Maintainer: Joe Kenny
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

#include <Chemistry_QC_GaussianBasis_DerivCentersInterface.hxx>

#include <sstream>
#include <iostream>
#include <vector>
#include <set>
#include <util/state/stateio.h>
#include <util/misc/ccaenv.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/cca/int/intcca.h>
#include <chemistry/cca/int/obintcca.h>
#include <chemistry/cca/int/tbintcca.h>
#include <util/class/scexception.h>
#include <Chemistry_QC_GaussianBasis_IntegralDescrInterface.hxx>
#include <ChemistryIntegralDescrCXX_CompositeIntegralDescr.hxx>
#include <ChemistryIntegralDescrCXX_OverlapIntegralDescr.hxx>
#include <ChemistryIntegralDescrCXX_PdotNuclearPIntegralDescr.hxx>
#include <ChemistryIntegralDescrCXX_PcrossNuclearPIntegralDescr.hxx>
#include <ChemistryIntegralDescrCXX_KineticIntegralDescr.hxx>
#include <ChemistryIntegralDescrCXX_NuclearIntegralDescr.hxx>
#include <ChemistryIntegralDescrCXX_HCoreIntegralDescr.hxx>
#include <ChemistryIntegralDescrCXX_DipoleIntegralDescr.hxx>
#include <ChemistryIntegralDescrCXX_QuadrupoleIntegralDescr.hxx>
#include <ChemistryIntegralDescrCXX_Eri4IntegralDescr.hxx>
#include <ChemistryIntegralDescrCXX_R12IntegralDescr.hxx>
#include <ChemistryIntegralDescrCXX_R12T1IntegralDescr.hxx>
#include <ChemistryIntegralDescrCXX_R12T2IntegralDescr.hxx>
#include <ChemistryIntegralDescrCXX_DerivCenters.hxx>
#include <ChemistryIntegralDescrCXX_DipoleData.hxx>
#include <chemistry/cca/int/cartit.h>
#include <chemistry/cca/int/tform.h>
#include <MPQC_IntV3EvaluatorFactory.hxx>

using namespace std;
using namespace sc;
using namespace Chemistry::QC::GaussianBasis;
using namespace ChemistryIntegralDescrCXX;

namespace auxintv3 {
  CartesianIter* new_cartesian_iter(int l);
  RedundantCartesianIter* new_redundant_cartesian_iter(int l);
  RedundantCartesianSubIter* new_redundant_cartesian_sub_iter(int l);
}

static int factory_instance_number=0;
static int sfactory_instance_number=0;
static int intv3port_instance_number=0;

static ClassDesc IntegralCCA_cd(
  typeid(IntegralCCA),"IntegralCCA",1,"public Integral",
  0, create<IntegralCCA>, create<IntegralCCA>);

extern Ref<Integral> default_integral;

/* may need to add optional "eval_factory" argument to this method in
integral class to get this capability
Integral*
Integral::get_default_integral()
{
  if (default_integral.null())
    default_integral = new IntegralCCA();

  return default_integral;
}
*/

IntegralCCA::IntegralCCA( const Ref<GaussianBasisSet> &b1,
                          const Ref<GaussianBasisSet> &b2,
                          const Ref<GaussianBasisSet> &b3,
                          const Ref<GaussianBasisSet> &b4,
                          std::string default_sf,
                          std::string factory_type,
                          sidl::array<std::string> types,
                          sidl::array<std::string> derivs,
                          sidl::array<std::string> sfacs,
                          bool intv3_order
                         ):
  Integral(b1,b2,b3,b4),
  default_subfactory_(default_sf), factory_type_(factory_type),
  types_(types), derivs_(derivs), sfacs_(sfacs), intv3_order_(intv3_order)
{
  superfactory_type_ = "Chemistry.IntegralSuperFactory";

  eval_req_ = CompositeIntegralDescr::_create();

  initialize_transforms();

  init_factory();

  init_generators();
}

IntegralCCA::IntegralCCA(const Ref<KeyVal> &keyval):
  Integral(keyval)
{
  superfactory_type_ = "Chemistry.IntegralSuperFactory";

  eval_req_ = CompositeIntegralDescr::_create();

  //------------
  // parse input
  //------------

  // for debugging/benchmarking
  intv3_order_ = false;
  bool tempbool = keyval->booleanvalue("intv3_order");
  if ( keyval->error() == KeyVal::OK ) {
    intv3_order_ = tempbool;
    if( intv3_order_ )
      ExEnv::out0() << indent
                    << "Using intv3 integral ordering by user request\n";
  }
  
  // get evaluator factory type (default to SuperFactory)
  factory_type_ = keyval->stringvalue("evaluator_factory");
  if ( keyval->error() != KeyVal::OK ) {
    factory_type_ = superfactory_type_;
  }

  // get subfactory configuration
  if ( factory_type_ == superfactory_type_ ) {

    // get a default
    default_subfactory_ = keyval->stringvalue("default_subfactory");
    if ( keyval->error() == KeyVal::OK ) 
      ExEnv::out0() << indent << "Default subfactory: " 
                    << default_subfactory_ << std::endl;
    else {
      ExEnv::out0() << indent 
                    << "Using MPQC.IntV3EvaluatorFactory for default subfactory"
	            << std::endl;
      default_subfactory_ = "MPQC.IntV3EvaluatorFactory";
    }

    // grab explicit factory associations
    int ntype = keyval->count("type");
    int nderiv = keyval->count("deriv");
    int nsfac = keyval->count("subfactory");
    if( ntype != nsfac ) 
      throw InputError("ntype != nsfac",__FILE__,__LINE__);
    if( nderiv > 0 && nderiv != ntype )
      throw InputError("nderiv != ntype",__FILE__,__LINE__);
    types_ = sidl::array<string>::create1d(ntype);
    derivs_ = sidl::array<string>::create1d(ntype);
    sfacs_ = sidl::array<string>::create1d(ntype);
    for( int i=0; i<ntype; ++i) {
      std::string tp( keyval->stringvalue("type",i) );
      types_.set( i, tp );
      std_types_.push_back( tp );
      std::string sf( keyval->stringvalue("subfactory",i) );
      sfacs_.set( i, sf );
      std_sfacs_.push_back( sf );
      std::string dv;
      if( nderiv )
        dv = keyval->stringvalue("deriv",i);
      else
        dv = "n";
      derivs_.set( i, dv );
      std_derivs_.push_back( dv );
    }
   
  }

  initialize_transforms();

  //-----------------
  // get eval factory
  //-----------------

  init_factory();

  //-------------------------------------------------
  // set up function objects for evaluator generation
  //-------------------------------------------------

  init_generators();

}

void
IntegralCCA::save_data_state(StateOut& s)
{
  Integral::save_data_state(s);
  s.put(superfactory_type_);
  s.put(buffer_);
  s.put(factory_type_);
  s.put(default_subfactory_);
  s.put(std_types_);
  s.put(std_sfacs_);
  s.put(std_derivs_);
}

IntegralCCA::IntegralCCA(StateIn& s): 
  Integral(s)
{

  initialize_transforms();

  eval_req_ = CompositeIntegralDescr::_create();

  s.get(superfactory_type_);
  s.get(buffer_);
  s.get(factory_type_);
  s.get(default_subfactory_);
  s.get(std_types_);
  s.get(std_sfacs_);
  s.get(std_derivs_);
  int ntype = std_types_.size();
  types_ = sidl::array<string>::create1d(ntype);
  derivs_ = sidl::array<string>::create1d(ntype);
  sfacs_ = sidl::array<string>::create1d(ntype);
  for( int i=0; i<ntype; ++i) {
    types_.set( i, std_types_[i] );
    sfacs_.set( i, std_sfacs_[i] );
    derivs_.set( i, std_derivs_[i] );
  }

  init_factory();
  init_generators();
}

void
IntegralCCA::init_factory()
{
  // grab cca environment
  gov::cca::Services &services = *CCAEnv::get_services();
  gov::cca::ports::BuilderService &bs = *CCAEnv::get_builder_service();
  gov::cca::TypeMap &type_map = *CCAEnv::get_type_map();
  gov::cca::ComponentID &my_id = *CCAEnv::get_component_id();

  // get a (super) evaluator factory
  int nsubfac;
  IntegralSuperFactoryInterface superfac;
  if( factory_type_ == superfactory_type_ ) {

    ostringstream fname;
      fname << "evaluator_factory" << factory_instance_number;
    ostringstream fname_port;
      fname_port << "IntegralSuperFactory" << factory_instance_number;
    ++factory_instance_number;

    // get the super factory
    //fac_id_ = bs.createInstance(fname.str(),factory_type_,type_map);
    fac_id_ = bs.createInstance(fname.str(),"Chemistry.IntegralSuperFactory",type_map);
    services.registerUsesPort(
      fname_port.str(),
      "Chemistry.QC.GaussianBasis.IntegralSuperFactoryInterface",
                              type_map);
    fac_con_ = bs.connect(my_id,fname_port.str(),
                          fac_id_,"IntegralSuperFactoryInterface");
    eval_factory_ = babel_cast<
      Chemistry::QC::GaussianBasis::IntegralSuperFactoryInterface> (
        services.getPort(fname_port.str()) );
    superfac = babel_cast<
      Chemistry::QC::GaussianBasis::IntegralSuperFactoryInterface> (
        services.getPort(fname_port.str()) );

    eval_factory_.set_default_subfactory( default_subfactory_ );
    eval_factory_.set_subfactory_config( types_, derivs_, sfacs_ );

    // get sub factories
    set<string> subfac_set;
    map<string,gov::cca::ComponentID> subfac_name_to_id;
    nsubfac = sfacs_.length();
    for( int i=0; i<nsubfac; ++i)
      subfac_set.insert( sfacs_.get(i) );
    if( !subfac_set.count(default_subfactory_) )
      subfac_set.insert( default_subfactory_ );
    nsubfac = subfac_set.size();
    set<string>::iterator iter;
    for (iter = subfac_set.begin(); iter != subfac_set.end(); iter++) {
      if( (*iter != "MPQC.IntV3EvaluatorFactory") && intv3_order_ )
        throw InputError( 
          "intv3_order can only be used with MPQC.IntV3EvaluatorFactory",
          __FILE__, __LINE__ );
      ExEnv::out0() << indent << "Instantiating: " << *iter << std::endl;
      ostringstream sfname;
      sfname << "subfactory" << sfactory_instance_number;
      ++sfactory_instance_number;
      subfac_name_to_id[sfname.str()] =
        bs.createInstance(sfname.str(),
                          *iter,
                          type_map);
      if( (*iter == "MPQC.IntV3EvaluatorFactory") && intv3_order_ ) {
        ostringstream intv3_port;
        intv3_port << "IntV3Port" << intv3port_instance_number;
        ++intv3port_instance_number;
        services.registerUsesPort(intv3_port.str(),
                                  "MPQC.IntV3EvaluatorFactory",
                                  type_map);
        gov::cca::ConnectionID conid = 
          bs.connect( my_id,intv3_port.str(),
                      subfac_name_to_id[sfname.str()],
                      "IntV3EvaluatorFactory");
        MPQC::IntV3EvaluatorFactory fac = 
          babel_cast<MPQC::IntV3EvaluatorFactory>(
            services.getPort(intv3_port.str()) );
        fac.set_reorder(false);
      }    
    }
    // connect factories with super factory
    sidl::array<string> sfac_portnames = superfac.add_uses_ports(nsubfac);
    map<string,gov::cca::ComponentID>::iterator miter;
    vector<gov::cca::ConnectionID> subfac_conids;
    int portname_iter=-1;
    for( miter = subfac_name_to_id.begin();
         miter != subfac_name_to_id.end(); miter++) {
      subfac_conids.push_back(
          bs.connect( fac_id_, sfac_portnames.get(++portname_iter),
                      (*miter).second, "IntegralEvaluatorFactoryInterface") );
    }
  }
  else { 
    //do a straightforward hook up
  }

  eval_factory_.set_storage(0);
}

void
IntegralCCA::init_generators()
{
  obgen_ = onebody_generator( this, eval_factory_, !intv3_order_ );
  obgen_.set_basis( bs1_, bs2_ );
  sc_eval_factory< OneBodyInt, onebody_generator>
    ob( obgen_ );
  get_onebody = ob;

  obdgen_ = onebody_deriv_generator( this, eval_factory_, !intv3_order_ );
  obdgen_.set_basis( bs1_, bs2_ );
  sc_eval_factory< OneBodyDerivInt, onebody_deriv_generator >
    obd( obdgen_ );
  get_onebody_deriv = obd;

  tbgen_ = twobody_generator( this, eval_factory_ );
  tbgen_.set_basis( bs1_, bs2_, bs3_, bs4_ );
  sc_eval_factory< TwoBodyInt, twobody_generator >
    tb( tbgen_ );
  get_twobody = tb;

  tbdgen_ = twobody_deriv_generator( this, eval_factory_ );
  tbdgen_.set_basis( bs1_, bs2_, bs3_, bs4_ );
  sc_eval_factory< TwoBodyDerivInt, twobody_deriv_generator >
    tbd( tbdgen_ );
  get_twobody_deriv = tbd;
}

IntegralCCA::~IntegralCCA()
{
  free_transforms();
}

Integral*
IntegralCCA::clone()
{
  return new IntegralCCA( bs1_, bs2_, bs3_, bs4_,
                          default_subfactory_, factory_type_,
                          types_, derivs_, sfacs_, intv3_order_ );
}

Ref<OneBodyInt>
IntegralCCA::overlap()
{
  eval_req_.clear();
  OverlapIntegralDescr desc =
    OverlapIntegralDescr::_create();
  cca_dcs_.push_back( ChemistryIntegralDescrCXX::DerivCenters::_create() );
  desc.set_deriv_centers( cca_dcs_.back() );
  descs_.push_back( desc );
  eval_req_.add_descr( desc );

  return get_onebody( eval_req_ );
}

Ref<OneBodyInt>
IntegralCCA::p_dot_nuclear_p()
{
  eval_req_.clear();
  PdotNuclearPIntegralDescr desc =
    PdotNuclearPIntegralDescr::_create();
  cca_dcs_.push_back( ChemistryIntegralDescrCXX::DerivCenters::_create() );
  desc.set_deriv_centers( cca_dcs_.back() );
  descs_.push_back( desc );
  eval_req_.add_descr( desc );

  return get_onebody( eval_req_ );
}

Ref<OneBodyInt>
IntegralCCA::p_cross_nuclear_p()
{
  eval_req_.clear();
  PcrossNuclearPIntegralDescr desc =
    PcrossNuclearPIntegralDescr::_create();
  cca_dcs_.push_back( ChemistryIntegralDescrCXX::DerivCenters::_create() );
  desc.set_deriv_centers( cca_dcs_.back() );
  descs_.push_back( desc );
  eval_req_.add_descr( desc );

  // need to make the the (x,y,z) components are correctly ordered
  throw FeatureNotImplemented("p_cross_nuclear_p needs more work",
                              __FILE__,__LINE__);

  return get_onebody( eval_req_ );
}

Ref<OneBodyInt>
IntegralCCA::kinetic()
{
  eval_req_.clear();
  KineticIntegralDescr desc =
    KineticIntegralDescr::_create();
  cca_dcs_.push_back( ChemistryIntegralDescrCXX::DerivCenters::_create() );
  desc.set_deriv_centers( cca_dcs_.back() );
  descs_.push_back( desc );
  eval_req_.add_descr( desc );
  
  return get_onebody( eval_req_ );
}

Ref<OneBodyInt>
IntegralCCA::nuclear()
{
  eval_req_.clear();
  NuclearIntegralDescr desc =
    NuclearIntegralDescr::_create();
  cca_dcs_.push_back( ChemistryIntegralDescrCXX::DerivCenters::_create() );
  desc.set_deriv_centers( cca_dcs_.back() );
  descs_.push_back( desc );
  eval_req_.add_descr( desc );
  
  return get_onebody( eval_req_ );
}

Ref<OneBodyInt>
IntegralCCA::hcore()
{
  eval_req_.clear();
  HCoreIntegralDescr desc =
    HCoreIntegralDescr::_create();
  cca_dcs_.push_back( ChemistryIntegralDescrCXX::DerivCenters::_create() );
  desc.set_deriv_centers( cca_dcs_.back() );
  descs_.push_back( desc );
  eval_req_.add_descr( desc );
  
  return get_onebody( eval_req_ );
}

Ref<OneBodyInt>
IntegralCCA::point_charge(const Ref<PointChargeData>& dat)
{
  throw FeatureNotImplemented("point_charge not implemented",
                              __FILE__,__LINE__);
}

Ref<OneBodyInt>
IntegralCCA::efield_dot_vector(const Ref<EfieldDotVectorData>&dat)
{
  throw FeatureNotImplemented("efield_dot_vector not iplemented",
                              __FILE__,__LINE__);
}

Ref<OneBodyInt>
IntegralCCA::dipole(const Ref<DipoleData>& dat)
{
  sidl::array<double> sidl_origin =  sidl::array<double>::create1d(3);

  if( dat.pointer() )
    for(int i=0; i<3; ++i)
      sidl_origin.set(i,dat->origin[i]);
  else
    for(int i=0; i<3; ++i)
      sidl_origin.set(i,0.0);

  ChemistryIntegralDescrCXX::DipoleData cca_dat
    = ChemistryIntegralDescrCXX::DipoleData::_create();
  cca_dat.set_origin( sidl_origin );

  eval_req_.clear();
  DipoleIntegralDescr desc 
    = DipoleIntegralDescr::_create();
  desc.set_dipole_data( cca_dat );
  desc.set_deriv_lvl(0);
  cca_dcs_.push_back( ChemistryIntegralDescrCXX::DerivCenters::_create() );
  desc.set_deriv_centers( cca_dcs_.back() );
  descs_.push_back( desc );
  eval_req_.add_descr( desc );

  return get_onebody( eval_req_ ); 
}

Ref<OneBodyInt>
IntegralCCA::quadrupole(const Ref<DipoleData>& dat)
{
  sidl::array<double> sidl_origin =  sidl::array<double>::create1d(3);

  if( dat.pointer() )
    for(int i=0; i<3; ++i)
      sidl_origin.set(i,dat->origin[i]);
  else
    for(int i=0; i<3; ++i)
      sidl_origin.set(i,0.0);

  ChemistryIntegralDescrCXX::DipoleData cca_dat
    = ChemistryIntegralDescrCXX::DipoleData::_create();
  cca_dat.set_origin( sidl_origin );

  eval_req_.clear();
  QuadrupoleIntegralDescr desc
    = QuadrupoleIntegralDescr::_create();
  desc.set_dipole_data( cca_dat );
  desc.set_deriv_lvl(0);
  cca_dcs_.push_back( ChemistryIntegralDescrCXX::DerivCenters::_create() );
  desc.set_deriv_centers( cca_dcs_.back() );
  descs_.push_back( desc );
  eval_req_.add_descr( desc );

  return get_onebody( eval_req_ );
}

Ref<OneBodyDerivInt>
IntegralCCA::overlap_deriv()
{
  eval_req_.clear();
  OverlapIntegralDescr desc =
    OverlapIntegralDescr::_create();
  desc.set_deriv_lvl(1);
  cca_dcs_.push_back( ChemistryIntegralDescrCXX::DerivCenters::_create() );
  desc.set_deriv_centers( cca_dcs_.back() );
  descs_.push_back( desc );
  eval_req_.add_descr( desc );
  
  return get_onebody_deriv( eval_req_ );
}

Ref<OneBodyDerivInt>
IntegralCCA::kinetic_deriv()
{
  eval_req_.clear();
  KineticIntegralDescr desc =
    KineticIntegralDescr::_create();
  desc.set_deriv_lvl(1);
  cca_dcs_.push_back( ChemistryIntegralDescrCXX::DerivCenters::_create() );
  desc.set_deriv_centers( cca_dcs_.back() );
  descs_.push_back( desc );
  eval_req_.add_descr( desc );
  
  return get_onebody_deriv( eval_req_ );
}

Ref<OneBodyDerivInt>
IntegralCCA::nuclear_deriv()
{
  eval_req_.clear();
  NuclearIntegralDescr desc =
    NuclearIntegralDescr::_create();
  desc.set_deriv_lvl(1);
  cca_dcs_.push_back( ChemistryIntegralDescrCXX::DerivCenters::_create() );
  desc.set_deriv_centers( cca_dcs_.back() );
  descs_.push_back( desc );
  eval_req_.add_descr( desc );
  
  return get_onebody_deriv( eval_req_ );
}

Ref<OneBodyDerivInt>
IntegralCCA::hcore_deriv()
{
  eval_req_.clear();
  HCoreIntegralDescr desc =
    HCoreIntegralDescr::_create();
  desc.set_deriv_lvl(1);
  cca_dcs_.push_back( ChemistryIntegralDescrCXX::DerivCenters::_create() );
  desc.set_deriv_centers( cca_dcs_.back() );
  descs_.push_back( desc );
  eval_req_.add_descr( desc );
  
  return get_onebody_deriv( eval_req_ );
}

Ref<TwoBodyInt>
IntegralCCA::electron_repulsion()
{
  eval_req_.clear();
  Eri4IntegralDescr desc =
    Eri4IntegralDescr::_create();
  cca_dcs_.push_back( ChemistryIntegralDescrCXX::DerivCenters::_create() );
  desc.set_deriv_centers( cca_dcs_.back() );
  descs_.push_back( desc );
  eval_req_.add_descr( desc );
  
  return get_twobody( eval_req_ );
}

Ref<TwoBodyDerivInt>
IntegralCCA::electron_repulsion_deriv()
{
  eval_req_.clear();
  Eri4IntegralDescr desc =
    Eri4IntegralDescr::_create();
  desc.set_deriv_lvl(1);
  cca_dcs_.push_back( ChemistryIntegralDescrCXX::DerivCenters::_create() );
  desc.set_deriv_centers( cca_dcs_.back() );
  descs_.push_back( desc );
  eval_req_.add_descr( desc );
  
  return get_twobody_deriv( eval_req_ );
}

Ref<TwoBodyInt>
IntegralCCA::grt()
{
  eval_req_.clear();
  Eri4IntegralDescr desc =
    Eri4IntegralDescr::_create();
  desc.set_deriv_lvl(0);
  cca_dcs_.push_back( ChemistryIntegralDescrCXX::DerivCenters::_create() );
  desc.set_deriv_centers( cca_dcs_.back() );
  descs_.push_back( desc );
  eval_req_.add_descr( desc );

  R12IntegralDescr desc2 = 
    R12IntegralDescr::_create();
  desc2.set_deriv_lvl(0);
  cca_dcs_.push_back( ChemistryIntegralDescrCXX::DerivCenters::_create() );
  desc2.set_deriv_centers( cca_dcs_.back() );
  descs_.push_back( desc2 );
  eval_req_.add_descr( desc2 );

  R12T1IntegralDescr desc3 = 
    R12T1IntegralDescr::_create();
  desc3.set_deriv_lvl(0);
  cca_dcs_.push_back( ChemistryIntegralDescrCXX::DerivCenters::_create() );
  desc3.set_deriv_centers( cca_dcs_.back() );
  descs_.push_back( desc3 );
  eval_req_.add_descr( desc3 );

  R12T2IntegralDescr desc4 = 
    R12T2IntegralDescr::_create();
  desc4.set_deriv_lvl(0);
  cca_dcs_.push_back( ChemistryIntegralDescrCXX::DerivCenters::_create() );
  desc4.set_deriv_centers( cca_dcs_.back() );
  descs_.push_back( desc4 );
  eval_req_.add_descr( desc4 );

  return get_twobody( eval_req_ );
}

void
IntegralCCA::set_basis(const Ref<GaussianBasisSet> &b1,
                       const Ref<GaussianBasisSet> &b2,
                       const Ref<GaussianBasisSet> &b3,
                       const Ref<GaussianBasisSet> &b4)
{
  free_transforms();
  Integral::set_basis(b1,b2,b3,b4);
  initialize_transforms();

  // delete &get_onebody; invalid???
  obgen_.set_basis( bs1_, bs2_ ); 
  sc_eval_factory< OneBodyInt, onebody_generator >
    ob( obgen_ );
  get_onebody = ob;

  // delete &get_onebody_deriv; invalid???
  obdgen_.set_basis( bs1_, bs2_ );
  sc_eval_factory< OneBodyDerivInt, onebody_deriv_generator >
    obd( obdgen_ );
  get_onebody_deriv = obd;

  // delete &get_twobody; invalid???
  tbgen_.set_basis( bs1_, bs2_, bs3_, bs4_ );
  sc_eval_factory< TwoBodyInt, twobody_generator >
    tb( tbgen_ ); 
  get_twobody = tb;

  // delete &get_twobody_deriv; invalid???
  tbdgen_.set_basis( bs1_, bs2_, bs3_, bs4_ );
  sc_eval_factory< TwoBodyDerivInt, twobody_deriv_generator >
    tbd( tbdgen_ );
  get_twobody_deriv = tbd;
}

CartesianIter *
IntegralCCA::new_cartesian_iter(int l)
{
  if( intv3_order_ )
    return new_cartesian_iterV3(l);

  return new CartesianIterCCA(l);
}

RedundantCartesianIter *
IntegralCCA::new_redundant_cartesian_iter(int l)
{
  if( intv3_order_ )
    return new_redundant_cartesian_iterV3(l);

  return new RedundantCartesianIterCCA(l);
}

RedundantCartesianSubIter *
IntegralCCA::new_redundant_cartesian_sub_iter(int l)
{
  if( intv3_order_ )
    return new_redundant_cartesian_sub_iterV3(l);

  return new RedundantCartesianSubIterCCA(l);
}

SphericalTransformIter *
IntegralCCA::new_spherical_transform_iter(int l, int inv, int subl)
{
  if( intv3_order_ )
    return new_spherical_transform_iterV3(l,inv,subl);
 
  if (l>maxl_ || l<0)
      throw ProgrammingError("new_spherical_transform_iter: bad l",
                             __FILE__,__LINE__);
  if (subl == -1) subl = l;
  if (subl < 0 || subl > l || (l-subl)%2 != 0)
      throw ProgrammingError("new_spherical_transform_iter: bad subl",
                             __FILE__,__LINE__);
  if (inv)
      return new SphericalTransformIter(ist_[l][(l-subl)/2]);
  return new SphericalTransformIter(st_[l][(l-subl)/2]);
}

const SphericalTransform *
IntegralCCA::spherical_transform(int l, int inv, int subl)
{
  if( intv3_order_ )
    return spherical_transformV3(l,inv,subl);   

  if (l>maxl_ || l<0)
      throw ProgrammingError("spherical_transform_iter: bad l",
                             __FILE__,__LINE__);
  if (subl == -1) subl = l;
  if (subl < 0 || subl > l || (l-subl)%2 != 0)
      throw ProgrammingError("spherical_transform_iter: bad subl",
                             __FILE__,__LINE__);
  if (inv)
      return ist_[l][(l-subl)/2];
  return st_[l][(l-subl)/2];
}

void
IntegralCCA::free_transforms()
{
  if( intv3_order_ )
    free_transformsV3();

  else {
    int i,j;
    for (i=0; i<=maxl_; i++) {
      for (j=0; j<=i/2; j++) {
        delete st_[i][j];
        delete ist_[i][j];
      }
      delete[] st_[i];
      delete[] ist_[i];
    }
    if (maxl_ >= 0) {
      delete[] st_;
      delete[] ist_;
    }
    st_ = NULL;
    ist_ = NULL;
  }

}

void
IntegralCCA::initialize_transforms()
{
  if( intv3_order_ )
    initialize_transformsV3();  

  else {
    maxl_ = -1;
    int maxam;
    maxam = bs1_.nonnull()?bs1_->max_angular_momentum():-1;
    if (maxl_ < maxam) maxl_ = maxam;
    maxam = bs2_.nonnull()?bs2_->max_angular_momentum():-1;
    if (maxl_ < maxam) maxl_ = maxam;
    maxam = bs3_.nonnull()?bs3_->max_angular_momentum():-1;
    if (maxl_ < maxam) maxl_ = maxam;
    maxam = bs4_.nonnull()?bs4_->max_angular_momentum():-1;
    if (maxl_ < maxam) maxl_ = maxam;
 
    if (maxl_ >= 0) {
      st_ = new SphericalTransform**[maxl_+1];
      ist_ = new ISphericalTransform**[maxl_+1];
      int i,j;
      for (i=0; i<=maxl_; i++) {
        st_[i] = new SphericalTransform*[i/2+1];
        ist_[i] = new ISphericalTransform*[i/2+1];
        for (j=0; j<=i/2; j++) {
          st_[i][j] = new SphericalTransformCCA(i,i-2*j);
          ist_[i][j] = new ISphericalTransformCCA(i,i-2*j);
        }
      }
    }
    else {
      st_ = NULL;
      ist_ = NULL;
    }
  }

}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
