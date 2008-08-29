// 
// File:          MPQC_IntV3EvaluatorFactory_Impl.cxx
// Symbol:        MPQC.IntV3EvaluatorFactory-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.IntV3EvaluatorFactory
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 
#include "MPQC_IntV3EvaluatorFactory_Impl.hxx"

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
// DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory._includes)
#include "Chemistry_QC_GaussianBasis_DescrInterface.hxx"

#include "MPQC_IntegralEvaluator1.hxx"
#include "MPQC_IntegralEvaluator2.hxx"
#include "MPQC_IntegralEvaluator3.hxx"
#include "MPQC_IntegralEvaluator4.hxx"
#include <sstream>
#include <stdexcept>
#include "basis_cca_to_sc.h"

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
#include "ChemistryDescrCXX_BaseDescr.hxx"

#include "Chemistry_QC_GaussianBasis_MultipoleDescrInterface.hxx"

using namespace std;
using namespace sc;
using namespace Chemistry;
using namespace Chemistry::QC::GaussianBasis;
using namespace ChemistryDescrCXX;

// DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory._includes)

// speical constructor, used for data wrapping(required).  Do not put code here unless you really know what you're doing!
MPQC::IntV3EvaluatorFactory_impl::IntV3EvaluatorFactory_impl() : StubBase(
  reinterpret_cast< void*>(::MPQC::IntV3EvaluatorFactory::_wrapObj(
  reinterpret_cast< void*>(this))),false) , _wrapped(true){ 
  // DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory._ctor2)
  // Insert-Code-Here {MPQC.IntV3EvaluatorFactory._ctor2} (ctor2)
  // DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory._ctor2)
}

// user defined constructor
void MPQC::IntV3EvaluatorFactory_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory._ctor)

  reorder_ = true;

  integral_ = new sc::IntegralV3();
  DescrInterface desc;

  cdesc_ = CompositeDescr::_create();

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

  desc = Eri4Descr::_create();
  desc.set_deriv_lvl(0);
  cdesc_.add_descr( desc );

  desc = OverlapDescr::_create();
  desc.set_deriv_lvl(1);
  cdesc_.add_descr(desc);

  desc = KineticDescr::_create();
  desc.set_deriv_lvl(1);
  cdesc_.add_descr(desc);

  desc = NuclearDescr::_create();
  desc.set_deriv_lvl(1);
  cdesc_.add_descr(desc);

  desc = HCoreDescr::_create();
  desc.set_deriv_lvl(1);
  cdesc_.add_descr(desc);

  desc = Eri4Descr::_create();
  desc.set_deriv_lvl(1);
  cdesc_.add_descr(desc);


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

  desc = Eri4Descr::_create();
  cdesc_no_deriv_.add_descr( desc );

  desc = OverlapDescr::_create();
  cdesc_no_deriv_.add_descr(desc);

  desc = KineticDescr::_create();
  cdesc_no_deriv_.add_descr(desc);

  desc = NuclearDescr::_create();
  cdesc_no_deriv_.add_descr(desc);

  desc = HCoreDescr::_create();
  cdesc_no_deriv_.add_descr(desc);

  desc = Eri4Descr::_create();
  cdesc_no_deriv_.add_descr(desc);

  // DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory._ctor)
}

// user defined destructor
void MPQC::IntV3EvaluatorFactory_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory._dtor)
  // DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory._dtor)
}

// static class initializer
void MPQC::IntV3EvaluatorFactory_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory._load)
  // Insert-Code-Here {MPQC.IntV3EvaluatorFactory._load} (class initialization)
  // DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory._load)
}

// user defined static methods: (none)

// user defined non-static methods:
/**
 * Method:  set_reorder[]
 */
void
MPQC::IntV3EvaluatorFactory_impl::set_reorder_impl (
  /* in */bool reorder ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory.set_reorder)

  reorder_ = reorder;

  // DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory.set_reorder)
}

/**
 * Method:  get_name[]
 */
::std::string
MPQC::IntV3EvaluatorFactory_impl::get_name_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory.get_name)

  return string("MPQC.IntV3EvaluatorFactory");

  // DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory.get_name)
}

/**
 * Method:  get_descriptor[]
 */
::Chemistry::QC::GaussianBasis::CompositeDescrInterface
MPQC::IntV3EvaluatorFactory_impl::get_descriptor_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory.get_descriptor)

  return cdesc_no_deriv_;

  // DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory.get_descriptor)
}

/**
 * Method:  is_supported[]
 */
bool
MPQC::IntV3EvaluatorFactory_impl::is_supported_impl (
  /* in */::Chemistry::QC::GaussianBasis::DescrInterface desc ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory.is_supported)

  string type = desc.get_type();
  int dlevel = desc.get_deriv_lvl();
  for( int i=0; i<cdesc_.get_n_descr(); ++i ) {
    if( type == cdesc_.get_descr(i).get_type() &&
        dlevel == cdesc_.get_descr(i).get_deriv_lvl() )
      return true;
  }

  return false;

  // DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory.is_supported)
}

/**
 *  Set available storage
 * @param storage Available storage in bytes 
 */
void
MPQC::IntV3EvaluatorFactory_impl::set_storage_impl (
  /* in */int64_t storage ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory.set_storage)

  storage_ = storage;
  integral_->set_storage( storage_ );

  // DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory.set_storage)
}

/**
 *  Get a 1-center integral evaluator
 * @param desc Integral set descriptor
 * @return 1-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator1Interface
MPQC::IntV3EvaluatorFactory_impl::get_evaluator1_impl (
  /* in */::Chemistry::QC::GaussianBasis::CompositeDescrInterface desc,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs1 ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory.get_evaluator1)

  //obocint_vec_.clear();
  bool is_oboc;

  vector< Ref<OneBodyOneCenterInt> > obocint_vec;

  MPQC::IntegralEvaluator1 eval = MPQC::IntegralEvaluator1::_create();
  eval.set_basis( bs1 );

  for( int i=0; i<desc.get_n_descr(); ++i ) {
  
    is_oboc = false;
    
    DescrInterface idesc = desc.get_descr(i);
    int ideriv = idesc.get_deriv_lvl();
    string itype = idesc.get_type();

    // this needs additional data
    //if( itype == "pointcharge"  && ideriv == 0 )
    //  obint = integral_->point_charge1( ??? );
    /*else*/
    throw runtime_error("IntV3EvaluatorFactory: unsupported integral type");

    if( obocint_vec.size() && is_oboc ) {
      eval.add_evaluator( (void*) obocint_vec.back().pointer(), idesc );
      obocint_vec_.push_back( obocint_vec );
    }
  }

  if( reorder_ )
    eval.init_reorder();
  
  return eval;
  
  // DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory.get_evaluator1)
}

/**
 *  Get a 2-center integral evaluator
 * @param desc Integral set descriptor
 * @return 2-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator2Interface
MPQC::IntV3EvaluatorFactory_impl::get_evaluator2_impl (
  /* in */::Chemistry::QC::GaussianBasis::CompositeDescrInterface desc,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs1,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs2 ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory.get_evaluator2)

  //obint_vec_.clear();
  //obderivint_vec_.clear();
  //tbtcint_vec_.clear();
  bool is_ob, is_obderiv, is_tbtc;

  vector< Ref<OneBodyInt> > obint_vec;
  vector< Ref<OneBodyDerivInt> > obderivint_vec;
  vector< Ref<TwoBodyTwoCenterInt> > tbtcint_vec;

  MPQC::IntegralEvaluator2 eval = MPQC::IntegralEvaluator2::_create();
  eval.set_basis( bs1, bs2 );

  sc::Ref<sc::GaussianBasisSet> sc_bs1 = basis_cca_to_sc( bs1 );
  sc::Ref<sc::GaussianBasisSet> sc_bs2;
  if( bs1.isSame(bs2) ) 
   sc_bs2.assign_pointer(sc_bs1.pointer());
  else sc_bs2 = basis_cca_to_sc( bs2 );
  integral_->set_basis(sc_bs1,sc_bs2);
  
  for( int i=0; i<desc.get_n_descr(); ++i ) {
   
    is_ob = is_obderiv = is_tbtc = false;
    
    DescrInterface idesc = desc.get_descr(i);
    int ideriv = idesc.get_deriv_lvl();
    string itype = idesc.get_type();

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
          idesc);
      sidl::array<double> origin = ddesc.get_origin();
      double sc_origin[3];
      sc_origin[0] = origin.get(0);
      sc_origin[1] = origin.get(1);
      sc_origin[2] = origin.get(2);
      dipole_data_ = new sc::DipoleData(sc_origin);
      obint_vec.push_back( integral_->dipole(dipole_data_) );
      is_ob = true;
    }
    else if( itype == "pcrossnuclearp" && ideriv == 0 ) {
      obint_vec.push_back( integral_->p_cross_nuclear_p() );
      is_ob = true;
    }
    else if( itype == "pdotnuclearp" && ideriv == 0 ) {
      obint_vec.push_back( integral_->p_dot_nuclear_p() );
      is_ob = true;
    }
    // these need additional data
    //else if( itype == "pointcharge2" && ideriv == 0 )
    //  obint_ = integral_->point_charge( ??? );
    //else if( itype == "efield_dot_vector" && ideriv == 0 )
    //  obint_ = integral_->efield_dot_vector( ??? );
    //else if( itype == "dipole" && ideriv == 0 )
    //  obint_ = integral_->dipole( ??? );
    //else if( itype == "quadrupole" && ideriv == 0 )
    //  obint_ = integal_->quadrupole( ??? );
    else if( itype == "overlap" && ideriv == 1 ) {
      obderivint_vec.push_back( integral_->overlap_deriv() );
      is_obderiv = true;
    }
    else if( itype == "kinetic" && ideriv == 1 ) {
      obderivint_vec.push_back( integral_->kinetic_deriv() );
      is_obderiv = true;
    }
    else if( itype == "nuclear" && ideriv == 1 ) {
      obderivint_vec.push_back( integral_->nuclear_deriv() );
      is_obderiv = true;
    }
    else if( itype == "hcore" && ideriv == 1 ) {
      obderivint_vec.push_back( integral_->hcore_deriv() );
      is_obderiv = true;
    }
    else if( itype == "eri2" && ideriv == 0 ) {
      tbtcint_vec.push_back( integral_->electron_repulsion2() );
      is_obderiv = true;
    }
    else 
      throw runtime_error("IntV3EvaluatorFactory: unsupported integral set");

    if( obint_vec.size() && is_ob ) {
      eval.add_evaluator( (void*) obint_vec.back().pointer(), 
                          idesc );
      obint_vec_.push_back( obint_vec );
    }
    else if( obderivint_vec.size() && is_obderiv ) {
      eval.add_evaluator( (void*) obderivint_vec.back().pointer(), 
                          idesc );
      obderivint_vec_.push_back( obderivint_vec );
    }
    else if( tbtcint_vec.size() && is_tbtc ) {
      // turning off integral caching for now -- need update to interface?
      //tbtcint_vec.back().pointer()->set_integral_storage(storage_);
      tbtcint_vec.back().pointer()->set_integral_storage(0);
      eval.add_evaluator( (void*) tbtcint_vec.back().pointer(), 
                          idesc );
      tbtcint_vec_.push_back( tbtcint_vec );
    }
  }

  if( reorder_ )
    eval.init_reorder();
  
  return eval;

  // DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory.get_evaluator2)
}

/**
 *  Get a 3-center integral evaluator
 * @param desc Integral set descriptor
 * @return 3-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator3Interface
MPQC::IntV3EvaluatorFactory_impl::get_evaluator3_impl (
  /* in */::Chemistry::QC::GaussianBasis::CompositeDescrInterface desc,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs1,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs2,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs3 ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory.get_evaluator3)

  //tb3cint_vec_.clear();
  bool is_tb3c;

  vector< Ref<TwoBodyThreeCenterInt> > tb3cint_vec;

  MPQC::IntegralEvaluator3 eval = MPQC::IntegralEvaluator3::_create();
  eval.set_basis( bs1, bs2, bs3 );

  sc::Ref<sc::GaussianBasisSet> sc_bs1 = basis_cca_to_sc( bs1 );
  sc::Ref<sc::GaussianBasisSet> sc_bs2, sc_bs3;
  if( bs1.isSame(bs2) ) sc_bs2.assign_pointer(sc_bs1.pointer());
  else sc_bs2 = basis_cca_to_sc( bs2 );
  if( bs2.isSame(bs3) ) sc_bs3.assign_pointer(sc_bs2.pointer());
  else sc_bs3 = basis_cca_to_sc( bs3 );
  integral_->set_basis(sc_bs1,sc_bs2,sc_bs3);


  for( int i=0; i<desc.get_n_descr(); ++i ) {
 
    is_tb3c = false;
    
    DescrInterface idesc = desc.get_descr(i);
    int ideriv = idesc.get_deriv_lvl();
    string itype = idesc.get_type();

    if( itype == "eri3"  && ideriv == 0 ) {
      tb3cint_vec.push_back( integral_->electron_repulsion3() );
      is_tb3c = true;
    }
    else 
      throw runtime_error("IntV3EvaluatorFactory: unsupported integral set");

    if( tb3cint_vec.size() && is_tb3c ) {
      // turning off integral caching for now -- need update to interface?
      // tb3cint_vec.back().pointer()->set_integral_storage(storage_);
      tb3cint_vec.back().pointer()->set_integral_storage(0);
      eval.add_evaluator( (void*) tb3cint_vec.back().pointer(), idesc );
      tb3cint_vec_.push_back( tb3cint_vec );
    }
  }

  if( reorder_ )
    eval.init_reorder();

  return eval;

  // DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory.get_evaluator3)
}

/**
 *  Get a 4-center integral evaluator
 * @param desc Integral set descriptor
 * @return 4-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator4Interface
MPQC::IntV3EvaluatorFactory_impl::get_evaluator4_impl (
  /* in */::Chemistry::QC::GaussianBasis::CompositeDescrInterface desc,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs1,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs2,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs3,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs4 ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory.get_evaluator4)

  //tbint_vec_.clear();
  //tbderivint_vec_.clear();
  bool is_tb, is_tbderiv;

  vector< Ref<TwoBodyInt> > tbint_vec;
  vector< Ref<TwoBodyDerivInt> > tbderivint_vec;

  MPQC::IntegralEvaluator4 eval = MPQC::IntegralEvaluator4::_create();
  eval.set_basis( bs1, bs2, bs3, bs4 );

  sc::Ref<sc::GaussianBasisSet> sc_bs1 = basis_cca_to_sc( bs1 );
  sc::Ref<sc::GaussianBasisSet> sc_bs2, sc_bs3, sc_bs4;
  if( bs1.isSame(bs2) ) sc_bs2.assign_pointer(sc_bs1.pointer());
  else sc_bs2 = basis_cca_to_sc( bs2 );
  if( bs2.isSame(bs3) ) sc_bs3.assign_pointer(sc_bs2.pointer());
  else sc_bs3 = basis_cca_to_sc( bs3 );
  if( bs3.isSame(bs3) ) sc_bs4.assign_pointer(sc_bs3.pointer());
  else sc_bs4 = basis_cca_to_sc( bs4);
  integral_->set_basis(sc_bs1,sc_bs2,sc_bs3,sc_bs4);

  for( int i=0; i<desc.get_n_descr(); ++i ) {
 
    is_tb = is_tbderiv = false;
    
    DescrInterface idesc = desc.get_descr(i);
    int ideriv = idesc.get_deriv_lvl();
    string itype = idesc.get_type();

    if( itype == "eri4"  && ideriv == 0 ) {
      tbint_vec.push_back( integral_->electron_repulsion() );
      is_tb = true;
    }
    else if( itype == "eri4" && ideriv == 1 ) {
      tbderivint_vec.push_back( integral_->electron_repulsion_deriv() );
      is_tbderiv = true;
    }
    else 
      throw runtime_error("IntV3EvaluatorFactory: unsupported integral set");

    if( tbint_vec.size() && is_tb ) {
      // turning off integral caching for now -- need update to interface?
      //tbint_vec.back().pointer()->set_integral_storage(storage_);
      tbint_vec.back().pointer()->set_integral_storage(0);
      eval.add_evaluator( (void*) tbint_vec.back().pointer(), idesc );
      tbint_vec_.push_back( tbint_vec );
    }
    else if( tbderivint_vec.size() && is_tbderiv ) {
      eval.add_evaluator( (void*) tbderivint_vec.back().pointer(), idesc );
      tbderivint_vec_.push_back( tbderivint_vec );
    }
  }

  if( reorder_ )
    eval.init_reorder();

  return eval;

  // DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory.get_evaluator4)
}

/**
 *  This should be called when the object is no longer needed.
 * No other members may be called after finalize. 
 */
int32_t
MPQC::IntV3EvaluatorFactory_impl::finalize_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory.finalize)
  // Insert-Code-Here {MPQC.IntV3EvaluatorFactory.finalize} (finalize method)

  // do something ? 

  // DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory.finalize)
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
MPQC::IntV3EvaluatorFactory_impl::setServices_impl (
  /* in */::gov::cca::Services services ) 
// throws:
//     ::gov::cca::CCAException
//     ::sidl::RuntimeException
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory.setServices)

  services_ = services;
  if (services_._is_nil()) return;

  try {
      services_.addProvidesPort( 
        *this, "IntegralEvaluatorFactoryInterface",
        "Chemistry.QC.GaussianBasis.IntegralEvaluatorFactoryInterface", 0);
      services_.addProvidesPort( *this, "IntV3EvaluatorFactory",
                    "MPQC.IntV3EvaluatorFactory",0);
  }
  catch (gov::cca::CCAException e) {
    std::cout << "Error using services: "
              << e.getNote() << std::endl;
  }

  // DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory.setServices)
}


// DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory._misc)
// DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory._misc)

