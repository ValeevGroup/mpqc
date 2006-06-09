// 
// File:          MPQC_IntV3EvaluatorFactory_Impl.cc
// Symbol:        MPQC.IntV3EvaluatorFactory-v0.2
// Symbol Type:   class
// Babel Version: 0.10.12
// Description:   Server-side implementation for MPQC.IntV3EvaluatorFactory
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// 
#include "MPQC_IntV3EvaluatorFactory_Impl.hh"

// DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory._includes)
#include "Chemistry_QC_GaussianBasis_IntegralDescr.hh"

#include "MPQC_IntegralEvaluator1.hh"
#include "MPQC_IntegralEvaluator2.hh"
#include "MPQC_IntegralEvaluator3.hh"
#include "MPQC_IntegralEvaluator4.hh"
#include <sstream>
#include <stdexcept>

#include "Chemistry_CompositeIntegralDescr.hh"
#include "Chemistry_OverlapIntegralDescr.hh"
#include "Chemistry_KineticIntegralDescr.hh"
#include "Chemistry_NuclearIntegralDescr.hh"
#include "Chemistry_HCoreIntegralDescr.hh"
#include "Chemistry_PointChargeIntegralDescr.hh"
#include "Chemistry_EfieldDotVectorIntegralDescr.hh"
#include "Chemistry_DipoleIntegralDescr.hh"
#include "Chemistry_QuadrupoleIntegralDescr.hh"
#include "Chemistry_Eri2IntegralDescr.hh"
#include "Chemistry_Eri3IntegralDescr.hh"
#include "Chemistry_Eri4IntegralDescr.hh"
#include "Chemistry_BaseIntegralDescr.hh"

#include "Chemistry_QC_GaussianBasis_DipoleIntegralDescr.hh"

sc::Ref<sc::GaussianBasisSet> 
basis_cca_to_sc( Chemistry::QC::GaussianBasis::Molecular& );

// DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory._includes)

// user-defined constructor.
void MPQC::IntV3EvaluatorFactory_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory._ctor)

  integral_ = new sc::IntegralV3();
  IntegralDescr desc;

  cdesc_ = Chemistry::CompositeIntegralDescr::_create();

  desc = OverlapIntegralDescr::_create();
  desc.set_deriv_lvl(0);
  cdesc_.add_descr( desc );

  desc = KineticIntegralDescr::_create();
  desc.set_deriv_lvl(0);
  cdesc_.add_descr( desc );

  desc = NuclearIntegralDescr::_create();
  desc.set_deriv_lvl(0);
  cdesc_.add_descr( desc );

  desc = HCoreIntegralDescr::_create(); 
  desc.set_deriv_lvl(0);
  cdesc_.add_descr( desc );

  desc = Chemistry::DipoleIntegralDescr::_create();
  desc.set_deriv_lvl(0);
  cdesc_.add_descr( desc );

  desc = Eri4IntegralDescr::_create();
  desc.set_deriv_lvl(0);
  cdesc_.add_descr( desc );

  desc = OverlapIntegralDescr::_create();
  desc.set_deriv_lvl(1);
  cdesc_.add_descr(desc);

  desc = KineticIntegralDescr::_create();
  desc.set_deriv_lvl(1);
  cdesc_.add_descr(desc);

  desc = NuclearIntegralDescr::_create();
  desc.set_deriv_lvl(1);
  cdesc_.add_descr(desc);

  desc = HCoreIntegralDescr::_create();
  desc.set_deriv_lvl(1);
  cdesc_.add_descr(desc);

  desc = Eri4IntegralDescr::_create();
  desc.set_deriv_lvl(1);
  cdesc_.add_descr(desc);


  cdesc_no_deriv_ = Chemistry::CompositeIntegralDescr::_create();

  desc = OverlapIntegralDescr::_create();
  cdesc_no_deriv_.add_descr( desc );

  desc = KineticIntegralDescr::_create();
  cdesc_no_deriv_.add_descr( desc );

  desc = NuclearIntegralDescr::_create();
  cdesc_no_deriv_.add_descr( desc );

  desc = HCoreIntegralDescr::_create();
  cdesc_no_deriv_.add_descr( desc );

  desc = Chemistry::DipoleIntegralDescr::_create();
  cdesc_no_deriv_.add_descr( desc );

  desc = Eri4IntegralDescr::_create();
  cdesc_no_deriv_.add_descr( desc );

  desc = OverlapIntegralDescr::_create();
  cdesc_no_deriv_.add_descr(desc);

  desc = KineticIntegralDescr::_create();
  cdesc_no_deriv_.add_descr(desc);

  desc = NuclearIntegralDescr::_create();
  cdesc_no_deriv_.add_descr(desc);

  desc = HCoreIntegralDescr::_create();
  cdesc_no_deriv_.add_descr(desc);

  desc = Eri4IntegralDescr::_create();
  cdesc_no_deriv_.add_descr(desc);

  // DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory._ctor)
}

// user-defined destructor.
void MPQC::IntV3EvaluatorFactory_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory._dtor)
  // DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory._dtor)
}

// static class initializer.
void MPQC::IntV3EvaluatorFactory_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory._load)
  // Insert-Code-Here {MPQC.IntV3EvaluatorFactory._load} (class initialization)
  // DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
/**
 * Method:  get_name[]
 */
::std::string
MPQC::IntV3EvaluatorFactory_impl::get_name ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory.get_name)

  return string("MPQC.IntV3EvaluatorFactory");

  // DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory.get_name)
}

/**
 * Method:  get_descriptor[]
 */
::Chemistry::QC::GaussianBasis::CompositeIntegralDescr
MPQC::IntV3EvaluatorFactory_impl::get_descriptor ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory.get_descriptor)

  return cdesc_no_deriv_;

  // DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory.get_descriptor)
}

/**
 * Method:  is_supported[]
 */
bool
MPQC::IntV3EvaluatorFactory_impl::is_supported (
  /* in */ ::Chemistry::QC::GaussianBasis::IntegralDescr desc ) 
throw ( 
  ::sidl::BaseException
){
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
 * Set available storage
 * @param storage Available storage in bytes 
 */
void
MPQC::IntV3EvaluatorFactory_impl::set_storage (
  /* in */ int64_t storage ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory.set_storage)

  storage_ = storage;

  // DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory.set_storage)
}

/**
 * Get a 1-center integral evaluator
 * @param desc Integral set descriptor
 * @return 1-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator1
MPQC::IntV3EvaluatorFactory_impl::get_evaluator1 (
  /* in */ ::Chemistry::QC::GaussianBasis::CompositeIntegralDescr desc,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1 ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory.get_evaluator1)

  //obocint_vec_.clear();
  bool is_oboc;

  vector< Ref<OneBodyOneCenterInt> > obocint_vec;

  integral_->set_storage(storage_);
  MPQC::IntegralEvaluator1 eval = MPQC::IntegralEvaluator1::_create();

  for( int i=0; i<desc.get_n_descr(); ++i ) {
  
    is_oboc = false;
    
    IntegralDescr idesc = desc.get_descr(i);
    int ideriv = idesc.get_deriv_lvl();
    string itype = idesc.get_type();

    std::cerr << "MPQC.IntV3EvaluatorFactory: creating " << itype
              << " evaluator\n";
 
    // this needs additional data
    //if( itype == "pointcharge"  && ideriv == 0 )
    //  obint = integral_->point_charge1( ??? );
    /*else*/
    throw runtime_error("IntV3EvaluatorFactory: unsupported integral type");

    if( obocint_vec.size() && is_oboc ) {
      eval.add_evaluator( (void*) obocint_vec.back().pointer(), desc );
      obocint_vec_.push_back( obocint_vec );
    }
  }

  eval.set_basis( bs1 );
  eval.init_reorder();
  
  return eval;
  
  // DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory.get_evaluator1)
}

/**
 * Get a 2-center integral evaluator
 * @param desc Integral set descriptor
 * @return 2-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator2
MPQC::IntV3EvaluatorFactory_impl::get_evaluator2 (
  /* in */ ::Chemistry::QC::GaussianBasis::CompositeIntegralDescr desc,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs2 ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory.get_evaluator2)

  //obint_vec_.clear();
  //obderivint_vec_.clear();
  //tbtcint_vec_.clear();
  bool is_ob, is_obderiv, is_tbtc;

  vector< Ref<OneBodyInt> > obint_vec;
  vector< Ref<OneBodyDerivInt> > obderivint_vec;
  vector< Ref<TwoBodyTwoCenterInt> > tbtcint_vec;

  integral_->set_storage(storage_);
  MPQC::IntegralEvaluator2 eval = MPQC::IntegralEvaluator2::_create();

  sc::Ref<sc::GaussianBasisSet> sc_bs1 = basis_cca_to_sc( bs1 );
  sc::Ref<sc::GaussianBasisSet> sc_bs2;
  if( bs1.isSame(bs2) ) 
   sc_bs2.assign_pointer(sc_bs1.pointer());
  else sc_bs2 = basis_cca_to_sc( bs2 );
  integral_->set_basis(sc_bs1,sc_bs2);
  
  for( int i=0; i<desc.get_n_descr(); ++i ) {
   
    is_ob = is_obderiv = is_tbtc = false;
    
    IntegralDescr idesc = desc.get_descr(i);
    int ideriv = idesc.get_deriv_lvl();
    string itype = idesc.get_type();

    std::cerr << "MPQC.IntV3EvaluatorFactory: creating " << itype 
              << " evaluator\n";

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

      sidl::SIDLException ex = sidl::SIDLException::_create();
      try {
        ex.setNote("IntV3 dipole doesn't work via CCA yet (reordering)");
        ex.add(__FILE__, __LINE__,"");
      }
      catch(...) { }
      throw ex;

      Chemistry::QC::GaussianBasis::DipoleIntegralDescr ddesc;
      ddesc = idesc;
      Chemistry::QC::GaussianBasis::DipoleData cca_data;
      cca_data = ddesc.get_dipole_data();
      sidl::array<double> origin = cca_data.get_origin();
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
      eval.add_evaluator( (void*) tbtcint_vec.back().pointer(), 
                          idesc );
      tbtcint_vec_.push_back( tbtcint_vec );
    }
  }

  eval.set_basis( bs1, bs2 );
  eval.init_reorder();
  
  return eval;

  // DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory.get_evaluator2)
}

/**
 * Get a 3-center integral evaluator
 * @param desc Integral set descriptor
 * @return 3-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator3
MPQC::IntV3EvaluatorFactory_impl::get_evaluator3 (
  /* in */ ::Chemistry::QC::GaussianBasis::CompositeIntegralDescr desc,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs2,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs3 ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory.get_evaluator3)

  //tb3cint_vec_.clear();
  bool is_tb3c;

  vector< Ref<TwoBodyThreeCenterInt> > tb3cint_vec;

  integral_->set_storage(storage_);
  MPQC::IntegralEvaluator3 eval = MPQC::IntegralEvaluator3::_create();

  sc::Ref<sc::GaussianBasisSet> sc_bs1 = basis_cca_to_sc( bs1 );
  sc::Ref<sc::GaussianBasisSet> sc_bs2, sc_bs3;
  if( bs1.isSame(bs2) ) sc_bs2.assign_pointer(sc_bs1.pointer());
  else sc_bs2 = basis_cca_to_sc( bs2 );
  if( bs2.isSame(bs3) ) sc_bs3.assign_pointer(sc_bs2.pointer());
  else sc_bs3 = basis_cca_to_sc( bs3 );
  integral_->set_basis(sc_bs1,sc_bs2,sc_bs3);


  for( int i=0; i<desc.get_n_descr(); ++i ) {
 
    is_tb3c = false;
    
    IntegralDescr idesc = desc.get_descr(i);
    int ideriv = idesc.get_deriv_lvl();
    string itype = idesc.get_type();

    std::cerr << "MPQC.IntV3EvaluatorFactory: creating " << itype
              << " evaluator\n";

    if( itype == "eri3"  && ideriv == 0 ) {
      tb3cint_vec.push_back( integral_->electron_repulsion3() );
      is_tb3c = true;
    }
    else 
      throw runtime_error("IntV3EvaluatorFactory: unsupported integral set");

    if( tb3cint_vec.size() && is_tb3c ) {
      eval.add_evaluator( (void*) tb3cint_vec.back().pointer(), desc );
      tb3cint_vec_.push_back( tb3cint_vec );
    }
  }

  eval.set_basis( bs1, bs2, bs3 );
  eval.init_reorder();

  return eval;

  // DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory.get_evaluator3)
}

/**
 * Get a 4-center integral evaluator
 * @param desc Integral set descriptor
 * @return 4-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator4
MPQC::IntV3EvaluatorFactory_impl::get_evaluator4 (
  /* in */ ::Chemistry::QC::GaussianBasis::CompositeIntegralDescr desc,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs2,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs3,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs4 ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory.get_evaluator4)

  //tbint_vec_.clear();
  //tbderivint_vec_.clear();
  bool is_tb, is_tbderiv;

  vector< Ref<TwoBodyInt> > tbint_vec;
  vector< Ref<TwoBodyDerivInt> > tbderivint_vec;

  integral_->set_storage(storage_);
  MPQC::IntegralEvaluator4 eval = MPQC::IntegralEvaluator4::_create();

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
    
    IntegralDescr idesc = desc.get_descr(i);
    int ideriv = idesc.get_deriv_lvl();
    string itype = idesc.get_type();

    std::cerr << "MPQC.IntV3EvaluatorFactory: creating " << itype
              << " evaluator\n";

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
      eval.add_evaluator( (void*) tbint_vec.back().pointer(), idesc );
      tbint_vec_.push_back( tbint_vec );
    }
    else if( tbderivint_vec.size() && is_tbderiv ) {
      eval.add_evaluator( (void*) tbderivint_vec.back().pointer(), idesc );
      tbderivint_vec_.push_back( tbderivint_vec );
    }
  }

  eval.set_basis( bs1, bs2, bs3, bs4 );
  eval.init_reorder();

  return eval;

  // DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory.get_evaluator4)
}

/**
 * Starts up a component presence in the calling framework.
 * @param services the component instance's handle on the framework world.
 * Contracts concerning Svc and setServices:
 * 
 * The component interaction with the CCA framework
 * and Ports begins on the call to setServices by the framework.
 * 
 * This function is called exactly once for each instance created
 * by the framework.
 * 
 * The argument Svc will never be nil/null.
 * 
 * Those uses ports which are automatically connected by the framework
 * (so-called service-ports) may be obtained via getPort during
 * setServices.
 */
void
MPQC::IntV3EvaluatorFactory_impl::setServices (
  /* in */ ::gov::cca::Services services ) 
throw ( 
  ::gov::cca::CCAException
){
  // DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory.setServices)

  services_ = services;
  if (services_._is_nil()) return;

  try {
      services_.addProvidesPort(self, "IntegralEvaluatorFactory",
                    "Chemistry.QC.GaussianBasis.IntegralEvaluatorFactory", 0);
  }
  catch (gov::cca::CCAException e) {
    std::cout << "Error using services: "
              << e.getNote() << std::endl;
  }

  // DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory.setServices)
}


// DO-NOT-DELETE splicer.begin(MPQC.IntV3EvaluatorFactory._misc)
// DO-NOT-DELETE splicer.end(MPQC.IntV3EvaluatorFactory._misc)

