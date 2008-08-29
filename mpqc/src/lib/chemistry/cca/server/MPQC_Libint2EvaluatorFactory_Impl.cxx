// 
// File:          MPQC_Libint2EvaluatorFactory_Impl.cxx
// Symbol:        MPQC.Libint2EvaluatorFactory-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.Libint2EvaluatorFactory
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 
#include "MPQC_Libint2EvaluatorFactory_Impl.hxx"

// 
// Includes for all method dependencies.
// 
#ifndef included_Chemistry_QC_GaussianBasis_CompositeDescrInterface_hxx
#include "Chemistry_QC_GaussianBasis_CompositeDescrInterface.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_DescrInterface_hxx
#include "Chemistry_QC_GaussianBasis_DescrInterface.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_IntegralEvaluator1Interface_hxx
#include "Chemistry_QC_GaussianBasis_IntegralEvaluator1Interface.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_IntegralEvaluator2Interface_hxx
#include "Chemistry_QC_GaussianBasis_IntegralEvaluator2Interface.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_IntegralEvaluator3Interface_hxx
#include "Chemistry_QC_GaussianBasis_IntegralEvaluator3Interface.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_IntegralEvaluator4Interface_hxx
#include "Chemistry_QC_GaussianBasis_IntegralEvaluator4Interface.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_MolecularInterface_hxx
#include "Chemistry_QC_GaussianBasis_MolecularInterface.hxx"
#endif
#ifndef included_gov_cca_CCAException_hxx
#include "gov_cca_CCAException.hxx"
#endif
#ifndef included_gov_cca_Services_hxx
#include "gov_cca_Services.hxx"
#endif
#ifndef included_sidl_BaseInterface_hxx
#include "sidl_BaseInterface.hxx"
#endif
#ifndef included_sidl_ClassInfo_hxx
#include "sidl_ClassInfo.hxx"
#endif
#ifndef included_sidl_RuntimeException_hxx
#include "sidl_RuntimeException.hxx"
#endif
#ifndef included_sidl_NotImplementedException_hxx
#include "sidl_NotImplementedException.hxx"
#endif
// DO-NOT-DELETE splicer.begin(MPQC.Libint2EvaluatorFactory._includes)

#include <string>
#include <stdexcept>
#include <vector>
#include "basis_cca_to_sc.h"

#include "MPQC_IntegralEvaluator1.hxx"
#include "MPQC_IntegralEvaluator2.hxx"
#include "MPQC_IntegralEvaluator3.hxx"
#include "MPQC_IntegralEvaluator4.hxx"

#include "Chemistry_QC_GaussianBasis_DescrInterface.hxx"
#include "Chemistry_QC_GaussianBasis_MultipoleDescrInterface.hxx"

#include "ChemistryDescrCXX_CompositeDescr.hxx"
#include "ChemistryDescrCXX_OverlapDescr.hxx"
#include "ChemistryDescrCXX_KineticDescr.hxx"
#include "ChemistryDescrCXX_NuclearDescr.hxx"
#include "ChemistryDescrCXX_HCoreDescr.hxx"
#include "ChemistryDescrCXX_PointChargeDescr.hxx"
#include "ChemistryDescrCXX_EfieldDotVectorDescr.hxx"
#include "ChemistryDescrCXX_DipoleDescr.hxx"
#include "ChemistryDescrCXX_QuadrupoleDescr.hxx"
#include "ChemistryDescrCXX_Eri2Descr.hxx"
#include "ChemistryDescrCXX_Eri3Descr.hxx"
#include "ChemistryDescrCXX_Eri4Descr.hxx"

using namespace sc;
using namespace Chemistry;
using namespace Chemistry::QC::GaussianBasis;
using namespace ChemistryDescrCXX;

// DO-NOT-DELETE splicer.end(MPQC.Libint2EvaluatorFactory._includes)

// speical constructor, used for data wrapping(required).  Do not put code here unless you really know what you're doing!
MPQC::Libint2EvaluatorFactory_impl::Libint2EvaluatorFactory_impl() : StubBase(
  reinterpret_cast< void*>(::MPQC::Libint2EvaluatorFactory::_wrapObj(
  reinterpret_cast< void*>(this))),false) , _wrapped(true){ 
  // DO-NOT-DELETE splicer.begin(MPQC.Libint2EvaluatorFactory._ctor2)
  // Insert-Code-Here {MPQC.Libint2EvaluatorFactory._ctor2} (ctor2)
  // DO-NOT-DELETE splicer.end(MPQC.Libint2EvaluatorFactory._ctor2)
}

// user defined constructor
void MPQC::Libint2EvaluatorFactory_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.Libint2EvaluatorFactory._ctor)

  integral_ = new sc::IntegralLibint2();
  DescrInterface desc;

  //
  // Create registry of integral types which can be computed 
  //
  cdesc_ = ChemistryDescrCXX::CompositeDescr::_create();

  desc = OverlapDescr::_create();
  desc.set_deriv_lvl(0);
  cdesc_.add_descr( desc );

  desc = KineticDescr::_create();
  desc.set_deriv_lvl(0);
  cdesc_.add_descr( desc );

  desc = NuclearDescr::_create();
  desc.set_deriv_lvl(0);
  cdesc_.add_descr( desc );

  desc = HCoreDescr::_create();
  desc.set_deriv_lvl(0);
  cdesc_.add_descr( desc );

  desc = DipoleDescr::_create();
  desc.set_deriv_lvl(0);
  cdesc_.add_descr( desc );

  desc = QuadrupoleDescr::_create();
  desc.set_deriv_lvl(0);
  cdesc_.add_descr( desc );

  desc = Eri4Descr::_create();
  desc.set_deriv_lvl(0);
  cdesc_.add_descr( desc );

  //
  // Create registry of integral types for which derivatives cannot be computed 
  //
  cdesc_no_deriv_ = CompositeDescr::_create();

  desc = OverlapDescr::_create();
  cdesc_no_deriv_.add_descr( desc );

  desc = KineticDescr::_create();
  cdesc_no_deriv_.add_descr( desc );

  desc = NuclearDescr::_create();
  cdesc_no_deriv_.add_descr( desc );

  desc = HCoreDescr::_create();
  cdesc_no_deriv_.add_descr( desc );

  desc = DipoleDescr::_create();
  cdesc_no_deriv_.add_descr( desc );

  desc = QuadrupoleDescr::_create();
  cdesc_no_deriv_.add_descr( desc );

  desc = Eri4Descr::_create();
  cdesc_no_deriv_.add_descr( desc );

  // DO-NOT-DELETE splicer.end(MPQC.Libint2EvaluatorFactory._ctor)
}

// user defined destructor
void MPQC::Libint2EvaluatorFactory_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.Libint2EvaluatorFactory._dtor)
  // Insert-Code-Here {MPQC.Libint2EvaluatorFactory._dtor} (destructor)
  // DO-NOT-DELETE splicer.end(MPQC.Libint2EvaluatorFactory._dtor)
}

// static class initializer
void MPQC::Libint2EvaluatorFactory_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.Libint2EvaluatorFactory._load)
  // Insert-Code-Here {MPQC.Libint2EvaluatorFactory._load} (class initialization)
  // DO-NOT-DELETE splicer.end(MPQC.Libint2EvaluatorFactory._load)
}

// user defined static methods: (none)

// user defined non-static methods:
/**
 * Method:  get_name[]
 */
::std::string
MPQC::Libint2EvaluatorFactory_impl::get_name_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Libint2EvaluatorFactory.get_name)

  return std::string("MPQC.Libint2EvaluatorFactory");

  // DO-NOT-DELETE splicer.end(MPQC.Libint2EvaluatorFactory.get_name)
}

/**
 * Method:  get_descriptor[]
 */
::Chemistry::QC::GaussianBasis::CompositeDescrInterface
MPQC::Libint2EvaluatorFactory_impl::get_descriptor_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Libint2EvaluatorFactory.get_descriptor)
  // Insert-Code-Here {MPQC.Libint2EvaluatorFactory.get_descriptor} (get_descriptor method)

  return cdesc_no_deriv_;

  // DO-NOT-DELETE splicer.end(MPQC.Libint2EvaluatorFactory.get_descriptor)
}

/**
 * Method:  is_supported[]
 */
bool
MPQC::Libint2EvaluatorFactory_impl::is_supported_impl (
  /* in */::Chemistry::QC::GaussianBasis::DescrInterface desc ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Libint2EvaluatorFactory.is_supported)

  const std::string type = desc.get_type();
  int dlevel = desc.get_deriv_lvl();
  for( int i=0; i<cdesc_.get_n_descr(); ++i ) {
    if( type == cdesc_.get_descr(i).get_type() &&
        dlevel == cdesc_.get_descr(i).get_deriv_lvl() )
      return true;
  }

  return false;

  // DO-NOT-DELETE splicer.end(MPQC.Libint2EvaluatorFactory.is_supported)
}

/**
 *  Set available storage
 * @param storage Available storage in bytes 
 */
void
MPQC::Libint2EvaluatorFactory_impl::set_storage_impl (
  /* in */int64_t storage ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Libint2EvaluatorFactory.set_storage)

  storage_ = storage;
  integral_->set_storage( storage_ );

  // DO-NOT-DELETE splicer.end(MPQC.Libint2EvaluatorFactory.set_storage)
}

/**
 *  Get a 1-center integral evaluator
 * @param desc Integral set descriptor
 * @return 1-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator1Interface
MPQC::Libint2EvaluatorFactory_impl::get_evaluator1_impl (
  /* in */::Chemistry::QC::GaussianBasis::CompositeDescrInterface desc,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs1 ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Libint2EvaluatorFactory.get_evaluator1)

  throw std::runtime_error("Libint2EvaluatorFactory: unsupported evaluator1 type");
 
  // DO-NOT-DELETE splicer.end(MPQC.Libint2EvaluatorFactory.get_evaluator1)
}

/**
 *  Get a 2-center integral evaluator
 * @param desc Integral set descriptor
 * @return 2-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator2Interface
MPQC::Libint2EvaluatorFactory_impl::get_evaluator2_impl (
  /* in */::Chemistry::QC::GaussianBasis::CompositeDescrInterface desc,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs1,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs2 ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Libint2EvaluatorFactory.get_evaluator2)

  //
  // Copied over from CintsEvaluatorFactory
  //

  using std::vector;
  vector< Ref<OneBodyInt> > obint_vec;
  vector< Ref<OneBodyDerivInt> > obderivint_vec;
  vector< Ref<TwoBodyTwoCenterInt> > tbtcint_vec;

  bool is_ob, is_obderiv, is_tbtc;

  MPQC::IntegralEvaluator2 eval = MPQC::IntegralEvaluator2::_create();
  eval.set_basis( bs1, bs2 );

  sc::Ref<sc::GaussianBasisSet> sc_bs1 = basis_cca_to_sc( bs1 );
  sc::Ref<sc::GaussianBasisSet> sc_bs2;
  if( bs1.isSame(bs2) ) sc_bs2.assign_pointer(sc_bs1.pointer());
  else sc_bs2 = basis_cca_to_sc( bs2 );
  integral_->set_basis(sc_bs1,sc_bs2);

  for( int i=0; i<desc.get_n_descr(); ++i ) {

    is_ob = is_obderiv = is_tbtc = false;

    DescrInterface idesc = desc.get_descr(i);
    const int ideriv = idesc.get_deriv_lvl();
    std::string itype = idesc.get_type();

    if( itype == "overlap"  && ideriv == 0 ) {
      obint_vec.push_back( integral_->overlap() );
      is_ob = true;
    }
    else if( itype == "kinetic" && ideriv == 0 ) {
      obint_vec.push_back( integral_->kinetic() );
      is_ob = true;
    }
    else if( itype == "nuclear" && ideriv == 0 ) {
      obint_vec.push_back( integral_->nuclear() );
      is_ob = true;
    }
    else if( itype == "hcore" && ideriv == 0 ) {
      obint_vec.push_back( integral_->hcore() );
      is_ob = true;
    }
    else if( itype == "dipole" && ideriv == 0 ) {
      Chemistry::QC::GaussianBasis::MultipoleDescrInterface ddesc;
      ddesc = 
        sidl::babel_cast<Chemistry::QC::GaussianBasis::MultipoleDescrInterface>(
          idesc );
      sidl::array<double> origin = ddesc.get_origin();
      double sc_origin[3];
      sc_origin[0] = origin.get(0);
      sc_origin[1] = origin.get(1);
      sc_origin[2] = origin.get(2);
      dipole_data_ = new sc::DipoleData(sc_origin);
      obint_vec.push_back( integral_->dipole(dipole_data_) );
      is_ob = true;
    }
    else if( itype == "quadrupole" && ideriv == 0 ) {
      Chemistry::QC::GaussianBasis::MultipoleDescrInterface ddesc;
      ddesc = 
        sidl::babel_cast<
          Chemistry::QC::GaussianBasis::MultipoleDescrInterface>( 
            idesc);
      sidl::array<double> origin = ddesc.get_origin();
      double sc_origin[3];
      sc_origin[0] = origin.get(0);
      sc_origin[1] = origin.get(1);
      sc_origin[2] = origin.get(2);
      quad_data_ = new sc::DipoleData(sc_origin);
      obint_vec.push_back( integral_->quadrupole(quad_data_) );
      is_ob = true;
    }
    else
      throw std::runtime_error("Libint2EvaluatorFactory: unsupported evaluator2 type");

    // multiple types could be a problem here
    if( obint_vec.size() && is_ob  ) {
      eval.add_evaluator( (void*) obint_vec.back().pointer(),
                          idesc );
      obint_vec_.push_back( obint_vec );
    }
    else
      throw std::runtime_error("Libint2EvaluatorFactory: unsupported evaluator2 type");
      
  }

  return eval;

  // DO-NOT-DELETE splicer.end(MPQC.Libint2EvaluatorFactory.get_evaluator2)
}

/**
 *  Get a 3-center integral evaluator
 * @param desc Integral set descriptor
 * @return 3-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator3Interface
MPQC::Libint2EvaluatorFactory_impl::get_evaluator3_impl (
  /* in */::Chemistry::QC::GaussianBasis::CompositeDescrInterface desc,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs1,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs2,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs3 ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Libint2EvaluatorFactory.get_evaluator3)

  throw std::runtime_error("Libint2EvaluatorFactory: unsupported evaluator3 type");
 
  // DO-NOT-DELETE splicer.end(MPQC.Libint2EvaluatorFactory.get_evaluator3)
}

/**
 *  Get a 4-center integral evaluator
 * @param desc Integral set descriptor
 * @return 4-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator4Interface
MPQC::Libint2EvaluatorFactory_impl::get_evaluator4_impl (
  /* in */::Chemistry::QC::GaussianBasis::CompositeDescrInterface desc,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs1,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs2,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs3,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs4 ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Libint2EvaluatorFactory.get_evaluator4)

  //
  // Copied over from CintsEvaluatorFactory and modified
  //
  
  bool is_tb, is_tbderiv;
  using std::vector;
  vector< Ref<TwoBodyInt> > tbint_vec;
  vector< Ref<TwoBodyDerivInt> > tbderivint_vec;

  //
  // Decide if the requested integrals can computed usign a composite evaluator
  //
  CompositeDescrInterface comp
    = CompositeDescr::_create();
  bool need_comp = false;
  for( int i=0; i<desc.get_n_descr(); ++i) {
    std::string tp = desc.get_descr(i).get_type();
    // insert the logic here, perhaps it needs to be encapsulated in a separate Libint2-specific class
    //if( tp != "eri")
    //  need_comp = true;
  }

  MPQC::IntegralEvaluator4 eval = MPQC::IntegralEvaluator4::_create();
  eval.set_basis( bs1, bs2, bs3, bs4 );

  sc::Ref<sc::GaussianBasisSet> sc_bs1 = basis_cca_to_sc( bs1 );
  sc::Ref<sc::GaussianBasisSet> sc_bs2, sc_bs3, sc_bs4;
  if( bs1.isSame(bs2) ) sc_bs2.assign_pointer(sc_bs1.pointer());
  else sc_bs2 = basis_cca_to_sc( bs2 );
  if( bs2.isSame(bs3) ) sc_bs3.assign_pointer(sc_bs2.pointer());
  else sc_bs3 = basis_cca_to_sc( bs3 );
  if( bs3.isSame(bs4) ) sc_bs4.assign_pointer(sc_bs3.pointer());
  else sc_bs4 = basis_cca_to_sc( bs4);
  integral_->set_basis(sc_bs1,sc_bs2,sc_bs3,sc_bs4);

  for( int i=0; i<desc.get_n_descr(); ++i ) {

    is_tb = is_tbderiv = false;

    DescrInterface idesc = desc.get_descr(i);
    const int ideriv = idesc.get_deriv_lvl();
    std::string itype = idesc.get_type();

    if( itype == "eri4"  && ideriv == 0 ) {
      tbint_vec.push_back( integral_->electron_repulsion() );
      is_tb = true;
    }
    else
      throw std::runtime_error("Libint2EvaluatorFactory: unsupported evaluator4 type");

    if( tbint_vec.size() && is_tb ) {
      // turning off integral caching for now -- need update to interface?
      //tbint_vec.back().pointer()->set_integral_storage(storage_);
      tbint_vec.back().pointer()->set_integral_storage(0);
      eval.add_evaluator( (void*) tbint_vec.back().pointer(), idesc );
      tbint_vec_.push_back( tbint_vec );
    }
    else
      throw std::runtime_error("Libint2EvaluatorFactory: unsupported evaluator4 type");
  }

  return eval;

  // DO-NOT-DELETE splicer.end(MPQC.Libint2EvaluatorFactory.get_evaluator4)
}

/**
 *  This should be called when the object is no longer needed.
 * No other members may be called after finalize. 
 */
int32_t
MPQC::Libint2EvaluatorFactory_impl::finalize_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Libint2EvaluatorFactory.finalize)

  // don't think need to do anything here
    
  // DO-NOT-DELETE splicer.end(MPQC.Libint2EvaluatorFactory.finalize)
}

/**
 *  Starts up a component presence in the calling framework.
 * @param services the component instance's handle on the framework world.
 * Contracts concerning services and setServices:
 * 
 * The component interaction with the CCA framework
 * and Ports begins on the call to setServices by the framework.
 * 
 * This function is called exactly once for each instance created
 * by the framework.
 * 
 * The argument services will never be nil/null.
 * 
 * Those uses ports which are automatically connected by the framework
 * (so-called service-ports) may be obtained via getPort during
 * setServices.
 */
void
MPQC::Libint2EvaluatorFactory_impl::setServices_impl (
  /* in */::gov::cca::Services services ) 
// throws:
//     ::gov::cca::CCAException
//     ::sidl::RuntimeException
{
  // DO-NOT-DELETE splicer.begin(MPQC.Libint2EvaluatorFactory.setServices)

  //
  // Copied over from CintsEvaluatorFactory
  //
  
  services_ = services;
  if (services_._is_nil()) return;

  try {
      services_.addProvidesPort( 
        *this, "IntegralEvaluatorFactoryInterface",
        "Chemistry.QC.GaussianBasis.IntegralEvaluatorFactoryInterface", 0);
  }
  catch (gov::cca::CCAException e) {
    std::cout << "Error using services: "
              << e.getNote() << std::endl;
  }

  // DO-NOT-DELETE splicer.end(MPQC.Libint2EvaluatorFactory.setServices)
}


// DO-NOT-DELETE splicer.begin(MPQC.Libint2EvaluatorFactory._misc)
// Insert-Code-Here {MPQC.Libint2EvaluatorFactory._misc} (miscellaneous code)
// DO-NOT-DELETE splicer.end(MPQC.Libint2EvaluatorFactory._misc)

