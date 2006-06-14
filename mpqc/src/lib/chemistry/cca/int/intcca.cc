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

#include <Chemistry_QC_GaussianBasis_DerivCenters.hh>

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
#include <Chemistry_QC_GaussianBasis_IntegralDescr.hh>
#include <Chemistry_CompositeIntegralDescr.hh>
#include <Chemistry_BaseIntegralDescr.hh>
#include <Chemistry_OverlapIntegralDescr.hh>
#include <Chemistry_PdotNuclearPIntegralDescr.hh>
#include <Chemistry_PcrossNuclearPIntegralDescr.hh>
#include <Chemistry_KineticIntegralDescr.hh>
#include <Chemistry_NuclearIntegralDescr.hh>
#include <Chemistry_HCoreIntegralDescr.hh>
#include <Chemistry_DipoleIntegralDescr.hh>
#include <Chemistry_QuadrupoleIntegralDescr.hh>
#include <Chemistry_Eri4IntegralDescr.hh>
#include <Chemistry_R12IntegralDescr.hh>
#include <Chemistry_R12T1IntegralDescr.hh>
#include <Chemistry_R12T2IntegralDescr.hh>
#include <Chemistry_DerivCenters.hh>
#include <Chemistry_DipoleData.hh>
#ifdef INTV3_ORDER
  #include <chemistry/qc/intv3/cartitv3.h>
  #include <chemistry/qc/intv3/tformv3.h>
#else
  #include <chemistry/cca/int/cartit.h>
  #include <chemistry/cca/int/tform.h>
#endif

using namespace std;
using namespace sc;
using namespace Chemistry::QC::GaussianBasis;

static int factory_instance_number=0;
static int sfactory_instance_number=0;

static ClassDesc IntegralCCA_cd(
  typeid(IntegralCCA),"IntegralCCA",1,"public Integral",
  0, create<IntegralCCA>, 0);

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

/* It'll take a fair bit of work to make superfactory work in ccafe
IntegralCCA::IntegralCCA( IntegralEvaluatorFactory eval_factory, 
			  bool use_opaque,
			  const Ref<GaussianBasisSet> &b1,
			  const Ref<GaussianBasisSet> &b2,
			  const Ref<GaussianBasisSet> &b3,
			  const Ref<GaussianBasisSet> &b4 ):
  Integral(b1,b2,b3,b4), eval_factory_(eval_factory), use_opaque_(use_opaque)
{
  initialize_transforms();
}
*/

IntegralCCA::IntegralCCA( bool use_opaque,
                          const Ref<GaussianBasisSet> &b1,
                          const Ref<GaussianBasisSet> &b2,
                          const Ref<GaussianBasisSet> &b3,
                          const Ref<GaussianBasisSet> &b4,
                          std::string default_sf,
                          std::string factory_type,
                          sidl::array<std::string> types,
                          sidl::array<std::string> derivs,
                          sidl::array<std::string> sfacs
                         ):
  Integral(b1,b2,b3,b4), use_opaque_(use_opaque),
  default_subfactory_(default_sf), factory_type_(factory_type),
  types_(types), derivs_(derivs), sfacs_(sfacs)
{
  superfactory_type_ = "Chemistry.IntegralSuperFactory";

  cca_dc_ = Chemistry::DerivCenters::_create();
  eval_req_ = Chemistry::CompositeIntegralDescr::_create();

  initialize_transforms();

  init_factory();

  init_generators();
}


IntegralCCA::IntegralCCA(const Ref<KeyVal> &keyval):
  Integral(keyval)
{

  superfactory_type_ = "Chemistry.IntegralSuperFactory";

  cca_dc_ = Chemistry::DerivCenters::_create();
  eval_req_ = Chemistry::CompositeIntegralDescr::_create();

  initialize_transforms();

  //------------
  // parse input
  //------------

  // get integral buffer type
  string buffer = keyval->stringvalue("integral_buffer");
  if ( keyval->error() != KeyVal::OK ) buffer = "opaque";
  if ( buffer == "opaque" ) use_opaque_ = true;
  else if ( buffer == "array" ) use_opaque_ = false;
  else 
    throw InputError( "integral_buffer must be either opaque or array",
		      __FILE__, __LINE__,"integral_buffer",
		      buffer.c_str(),class_desc() );
  
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
      types_.set( i, keyval->stringvalue("type",i) );
      sfacs_.set( i, keyval->stringvalue("subfactory",i) );
      if( nderiv )
        derivs_.set( i, keyval->stringvalue("deriv",i) );
      else
        derivs_.set( i, "n" );
    }
  }

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
IntegralCCA::init_factory()
{
  // grab cca environment
  gov::cca::Services &services = *CCAEnv::get_services();
  gov::cca::ports::BuilderService &bs = *CCAEnv::get_builder_service();
  gov::cca::TypeMap &type_map = *CCAEnv::get_type_map();
  gov::cca::ComponentID &my_id = *CCAEnv::get_component_id();

  // get a (super) evaluator factory
  if( factory_type_ == superfactory_type_ ) {

    ostringstream fname;
      fname << "evaluator_factory" << factory_instance_number;
    ostringstream fname_port;
      fname_port << "IntegralSuperFactory" << factory_instance_number;
    ++factory_instance_number;

    // get the super factory
    fac_id_ = bs.createInstance(fname.str(),factory_type_,type_map);
    services.registerUsesPort(fname_port.str(),
                              "Chemistry.QC.GaussianBasis.IntegralSuperFactory",
                              type_map);
    fac_con_ = bs.connect(my_id,fname_port.str(),
                          fac_id_,"IntegralSuperFactory");
    eval_factory_ = services.getPort(fname_port.str());
    IntegralSuperFactory superfac = services.getPort(fname_port.str());

    eval_factory_.set_default_subfactory( default_subfactory_ );
    eval_factory_.set_subfactory_config( types_, derivs_, sfacs_ );

    // get sub factories
    set<string> subfac_set;
    map<string,gov::cca::ComponentID> subfac_name_to_id;
    int nsubfac = sfacs_.length();
    for( int i=0; i<nsubfac; ++i)
      subfac_set.insert( sfacs_.get(i) );
    if( !subfac_set.count(default_subfactory_) )
      subfac_set.insert( default_subfactory_ );
    nsubfac = subfac_set.size();
    set<string>::iterator iter;
    for (iter = subfac_set.begin(); iter != subfac_set.end(); iter++) {
      std::cerr << "Instantiating: " << *iter << std::endl;
      ostringstream sfname;
      sfname << "subfactory" << sfactory_instance_number;
      ++sfactory_instance_number;
      subfac_name_to_id[sfname.str()] =
        bs.createInstance(sfname.str(),
                          *iter,
                          type_map);
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
                      (*miter).second, "IntegralEvaluatorFactory") );
    }
  }
  else { 
    //do a straightforward hook up
  }
}

void
IntegralCCA::init_generators()
{
  obgen_ = onebody_generator( this, eval_factory_, use_opaque_ );
  obgen_.set_basis( bs1_, bs2_ );
  sc_eval_factory< OneBodyInt, onebody_generator>
    ob( obgen_ );
  get_onebody = ob;

  obdgen_ = onebody_deriv_generator( this, eval_factory_, use_opaque_ );
  obdgen_.set_basis( bs1_, bs2_ );
  sc_eval_factory< OneBodyDerivInt, onebody_deriv_generator >
    obd( obdgen_ );
  get_onebody_deriv = obd;

  tbgen_ = twobody_generator( this, eval_factory_, use_opaque_ );
  tbgen_.set_basis( bs1_, bs2_, bs3_, bs4_ );
  sc_eval_factory< TwoBodyInt, twobody_generator >
    tb( tbgen_ );
  get_twobody = tb;

  tbdgen_ = twobody_deriv_generator( this, eval_factory_, use_opaque_ );
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
  return new IntegralCCA( use_opaque_, bs1_, bs2_, bs3_, bs4_,
                          default_subfactory_, factory_type_,
                          types_, derivs_, sfacs_ );
}

CartesianIter *
IntegralCCA::new_cartesian_iter(int l)
{
#ifdef INTV3_ORDER
  return new CartesianIterV3(l);
#else
  return new CartesianIterCCA(l);
#endif
}

RedundantCartesianIter *
IntegralCCA::new_redundant_cartesian_iter(int l)
{
#ifdef INTV3_ORDER
  return new RedundantCartesianIterV3(l);
#else
  return new RedundantCartesianIterCCA(l);
#endif
}

RedundantCartesianSubIter *
IntegralCCA::new_redundant_cartesian_sub_iter(int l)
{
#ifdef INTV3_ORDER
  return new RedundantCartesianSubIterV3(l);
#else
  return new RedundantCartesianSubIterCCA(l);
#endif
}

SphericalTransformIter *
IntegralCCA::new_spherical_transform_iter(int l, int inv, int subl)
{
#ifdef INTV3_ORDER

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

#else

  // CINTS version
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

#endif

}

const SphericalTransform *
IntegralCCA::spherical_transform(int l, int inv, int subl)
{
#ifdef INTV3_ORDER

  // INTV3 version
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

#else

  // CINTS version
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

#endif

}

Ref<OneBodyInt>
IntegralCCA::overlap()
{
  eval_req_.clear();
  Chemistry::OverlapIntegralDescr desc =
    Chemistry::OverlapIntegralDescr::_create();
  desc.set_deriv_centers( cca_dc_ );
  eval_req_.add_descr( desc );

  return get_onebody( eval_req_ );
}

Ref<OneBodyInt>
IntegralCCA::p_dot_nuclear_p()
{
  eval_req_.clear();
  Chemistry::PdotNuclearPIntegralDescr desc =
    Chemistry::PdotNuclearPIntegralDescr::_create();
  desc.set_deriv_centers( cca_dc_ );
  eval_req_.add_descr( desc );

  return get_onebody( eval_req_ );
}

Ref<OneBodyInt>
IntegralCCA::p_cross_nuclear_p()
{
  eval_req_.clear();
  Chemistry::PcrossNuclearPIntegralDescr desc =
    Chemistry::PcrossNuclearPIntegralDescr::_create();
  desc.set_deriv_centers( cca_dc_ );
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
  Chemistry::KineticIntegralDescr desc =
    Chemistry::KineticIntegralDescr::_create();
  desc.set_deriv_centers( cca_dc_ );
  eval_req_.add_descr( desc );
  
  return get_onebody( eval_req_ );
}

Ref<OneBodyInt>
IntegralCCA::nuclear()
{
  eval_req_.clear();
  Chemistry::NuclearIntegralDescr desc =
    Chemistry::NuclearIntegralDescr::_create();
  desc.set_deriv_centers( cca_dc_ );
  eval_req_.add_descr( desc );
  
  return get_onebody( eval_req_ );
}

Ref<OneBodyInt>
IntegralCCA::hcore()
{
  eval_req_.clear();
  Chemistry::HCoreIntegralDescr desc =
    Chemistry::HCoreIntegralDescr::_create();
  desc.set_deriv_centers( cca_dc_ );
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

  Chemistry::DipoleData cca_dat
    = Chemistry::DipoleData::_create();
  cca_dat.set_origin( sidl_origin );

  eval_req_.clear();
  Chemistry::DipoleIntegralDescr desc 
    = Chemistry::DipoleIntegralDescr::_create();
  desc.set_dipole_data( cca_dat );
  desc.set_deriv_lvl(0);
  desc.set_deriv_centers( cca_dc_ );
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

  Chemistry::DipoleData cca_dat
    = Chemistry::DipoleData::_create();
  cca_dat.set_origin( sidl_origin );

  eval_req_.clear();
  Chemistry::QuadrupoleIntegralDescr desc
    = Chemistry::QuadrupoleIntegralDescr::_create();
  desc.set_dipole_data( cca_dat );
  desc.set_deriv_lvl(0);
  desc.set_deriv_centers( cca_dc_ );
  eval_req_.add_descr( desc );

  return get_onebody( eval_req_ );
}

Ref<OneBodyDerivInt>
IntegralCCA::overlap_deriv()
{
  eval_req_.clear();
  Chemistry::OverlapIntegralDescr desc =
    Chemistry::OverlapIntegralDescr::_create();
  desc.set_deriv_lvl(1);
  desc.set_deriv_centers( cca_dc_ );
  eval_req_.add_descr( desc );
  
  return get_onebody_deriv( eval_req_ );
}

Ref<OneBodyDerivInt>
IntegralCCA::kinetic_deriv()
{
  eval_req_.clear();
  Chemistry::KineticIntegralDescr desc =
    Chemistry::KineticIntegralDescr::_create();
  desc.set_deriv_lvl(1);
  desc.set_deriv_centers( cca_dc_ );
  eval_req_.add_descr( desc );
  
  return get_onebody_deriv( eval_req_ );
}

Ref<OneBodyDerivInt>
IntegralCCA::nuclear_deriv()
{
  eval_req_.clear();
  Chemistry::NuclearIntegralDescr desc =
    Chemistry::NuclearIntegralDescr::_create();
  desc.set_deriv_lvl(1);
  desc.set_deriv_centers( cca_dc_ );
  eval_req_.add_descr( desc );
  
  return get_onebody_deriv( eval_req_ );
}

Ref<OneBodyDerivInt>
IntegralCCA::hcore_deriv()
{
  eval_req_.clear();
  Chemistry::HCoreIntegralDescr desc =
    Chemistry::HCoreIntegralDescr::_create();
  desc.set_deriv_lvl(1);
  desc.set_deriv_centers( cca_dc_ );
  eval_req_.add_descr( desc );
  
  return get_onebody_deriv( eval_req_ );
}

Ref<TwoBodyInt>
IntegralCCA::electron_repulsion()
{
  eval_req_.clear();
  Chemistry::Eri4IntegralDescr desc =
    Chemistry::Eri4IntegralDescr::_create();
  desc.set_deriv_centers( cca_dc_ );
  eval_req_.add_descr( desc );
  
  return get_twobody( eval_req_ );
}

Ref<TwoBodyDerivInt>
IntegralCCA::electron_repulsion_deriv()
{
  eval_req_.clear();
  Chemistry::Eri4IntegralDescr desc =
    Chemistry::Eri4IntegralDescr::_create();
  desc.set_deriv_lvl(1);
  desc.set_deriv_centers( cca_dc_ );
  eval_req_.add_descr( desc );
  
  return get_twobody_deriv( eval_req_ );
}

Ref<TwoBodyInt>
IntegralCCA::grt()
{
  eval_req_.clear();
  Chemistry::QC::GaussianBasis::IntegralDescr desc =
    Chemistry::Eri4IntegralDescr::_create();
  desc.set_deriv_lvl(0);
  desc.set_deriv_centers( cca_dc_ );
  eval_req_.add_descr( desc );

  desc = Chemistry::R12IntegralDescr::_create();
  desc.set_deriv_lvl(0);
  desc.set_deriv_centers( cca_dc_ );
  eval_req_.add_descr( desc );

  desc = Chemistry::R12T1IntegralDescr::_create();
  desc.set_deriv_lvl(0);
  desc.set_deriv_centers( cca_dc_ );
  eval_req_.add_descr( desc );

  desc = Chemistry::R12T2IntegralDescr::_create();
  desc.set_deriv_lvl(0);
  desc.set_deriv_centers( cca_dc_ );
  eval_req_.add_descr( desc );

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

void
IntegralCCA::free_transforms()
{
#ifdef INTV3_ORDER

  // INTV3 version
  int i,j;
  for (i=0; i<=maxl_; i++) {
      for (j=0; j<=i/2; j++) {
          delete st_[i][j];
          delete ist_[i][j];
        }
      delete[] st_[i];
      delete[] ist_[i];
    }
  delete[] st_;
  delete[] ist_;
  st_ = 0;
  ist_ = 0;

#else

  // CINTS version
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

#endif

}

 void
 IntegralCCA::initialize_transforms()
 {
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
 
#ifdef INTV3_ORDER

   // INTV3 version
   st_ = new SphericalTransform**[maxl_+1];
   ist_ = new ISphericalTransform**[maxl_+1];
   int i,j;
   for (i=0; i<=maxl_; i++) {
       st_[i] = new SphericalTransform*[i/2+1];
       ist_[i] = new ISphericalTransform*[i/2+1];
       for (j=0; j<=i/2; j++) {
           st_[i][j] = new SphericalTransformV3(i,i-2*j);
           ist_[i][j] = new ISphericalTransformV3(i,i-2*j);
         }
     }
#else

  // CINTS version
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

#endif

}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
