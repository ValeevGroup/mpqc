// 
// File:          MPQC_CintsEvaluatorFactory_Impl.cc
// Symbol:        MPQC.CintsEvaluatorFactory-v0.2
// Symbol Type:   class
// Babel Version: 0.10.12
// Description:   Server-side implementation for MPQC.CintsEvaluatorFactory
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// 
#include "MPQC_CintsEvaluatorFactory_Impl.hh"

// DO-NOT-DELETE splicer.begin(MPQC.CintsEvaluatorFactory._includes)
// Insert-Code-Here {MPQC.CintsEvaluatorFactory._includes} (additional includes or code)
// DO-NOT-DELETE splicer.end(MPQC.CintsEvaluatorFactory._includes)

// user-defined constructor.
void MPQC::CintsEvaluatorFactory_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.CintsEvaluatorFactory._ctor)
  // Insert-Code-Here {MPQC.CintsEvaluatorFactory._ctor} (constructor)
  // DO-NOT-DELETE splicer.end(MPQC.CintsEvaluatorFactory._ctor)
}

// user-defined destructor.
void MPQC::CintsEvaluatorFactory_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.CintsEvaluatorFactory._dtor)
  // Insert-Code-Here {MPQC.CintsEvaluatorFactory._dtor} (destructor)
  // DO-NOT-DELETE splicer.end(MPQC.CintsEvaluatorFactory._dtor)
}

// static class initializer.
void MPQC::CintsEvaluatorFactory_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.CintsEvaluatorFactory._load)
  // Insert-Code-Here {MPQC.CintsEvaluatorFactory._load} (class initialization)
  // DO-NOT-DELETE splicer.end(MPQC.CintsEvaluatorFactory._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
/**
 * Method:  get_name[]
 */
::std::string
MPQC::CintsEvaluatorFactory_impl::get_name ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(MPQC.CintsEvaluatorFactory.get_name)
  // Insert-Code-Here {MPQC.CintsEvaluatorFactory.get_name} (get_name method)
  // DO-NOT-DELETE splicer.end(MPQC.CintsEvaluatorFactory.get_name)
}

/**
 * Method:  get_descriptor[]
 */
::Chemistry::QC::GaussianBasis::CompositeIntegralDescr
MPQC::CintsEvaluatorFactory_impl::get_descriptor ()
throw ( 
  ::sidl::BaseException
)
{
  // DO-NOT-DELETE splicer.begin(MPQC.CintsEvaluatorFactory.get_descriptor)
  // Insert-Code-Here {MPQC.CintsEvaluatorFactory.get_descriptor} (get_descriptor method)
  // DO-NOT-DELETE splicer.end(MPQC.CintsEvaluatorFactory.get_descriptor)
}

/**
 * Method:  get_max_deriv_lvls[]
 */
::sidl::array<int32_t>
MPQC::CintsEvaluatorFactory_impl::get_max_deriv_lvls (
  /* in */ ::Chemistry::QC::GaussianBasis::CompositeIntegralDescr desc ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(MPQC.CintsEvaluatorFactory.get_max_deriv_lvls)
  // Insert-Code-Here {MPQC.CintsEvaluatorFactory.get_max_deriv_lvls} (get_max_deriv_lvls method)
  // DO-NOT-DELETE splicer.end(MPQC.CintsEvaluatorFactory.get_max_deriv_lvls)
}

/**
 * Set available storage
 * @param storage Available storage in bytes 
 */
void
MPQC::CintsEvaluatorFactory_impl::set_storage (
  /* in */ int64_t storage ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(MPQC.CintsEvaluatorFactory.set_storage)
  // Insert-Code-Here {MPQC.CintsEvaluatorFactory.set_storage} (set_storage method)
  // DO-NOT-DELETE splicer.end(MPQC.CintsEvaluatorFactory.set_storage)
}

/**
 * Get a 1-center integral evaluator
 * @param desc Integral set descriptor
 * @return 1-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator1
MPQC::CintsEvaluatorFactory_impl::get_evaluator1 (
  /* in */ ::Chemistry::QC::GaussianBasis::CompositeIntegralDescr desc,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1 ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(MPQC.CintsEvaluatorFactory.get_evaluator1)
  // Insert-Code-Here {MPQC.CintsEvaluatorFactory.get_evaluator1} (get_evaluator1 method)
  // DO-NOT-DELETE splicer.end(MPQC.CintsEvaluatorFactory.get_evaluator1)
}

/**
 * Get a 2-center integral evaluator
 * @param desc Integral set descriptor
 * @return 2-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator2
MPQC::CintsEvaluatorFactory_impl::get_evaluator2 (
  /* in */ ::Chemistry::QC::GaussianBasis::CompositeIntegralDescr desc,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs2 ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(MPQC.CintsEvaluatorFactory.get_evaluator2)
  // Insert-Code-Here {MPQC.CintsEvaluatorFactory.get_evaluator2} (get_evaluator2 method)
  // DO-NOT-DELETE splicer.end(MPQC.CintsEvaluatorFactory.get_evaluator2)
}

/**
 * Get a 3-center integral evaluator
 * @param desc Integral set descriptor
 * @return 3-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator3
MPQC::CintsEvaluatorFactory_impl::get_evaluator3 (
  /* in */ ::Chemistry::QC::GaussianBasis::CompositeIntegralDescr desc,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs2,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs3 ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(MPQC.CintsEvaluatorFactory.get_evaluator3)
  // Insert-Code-Here {MPQC.CintsEvaluatorFactory.get_evaluator3} (get_evaluator3 method)
  // DO-NOT-DELETE splicer.end(MPQC.CintsEvaluatorFactory.get_evaluator3)
}

/**
 * Get a 4-center integral evaluator
 * @param desc Integral set descriptor
 * @return 4-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator4
MPQC::CintsEvaluatorFactory_impl::get_evaluator4 (
  /* in */ ::Chemistry::QC::GaussianBasis::CompositeIntegralDescr desc,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs1,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs2,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs3,
  /* in */ ::Chemistry::QC::GaussianBasis::Molecular bs4 ) 
throw ( 
  ::sidl::BaseException
){
  // DO-NOT-DELETE splicer.begin(MPQC.CintsEvaluatorFactory.get_evaluator4)
  // Insert-Code-Here {MPQC.CintsEvaluatorFactory.get_evaluator4} (get_evaluator4 method)
  // DO-NOT-DELETE splicer.end(MPQC.CintsEvaluatorFactory.get_evaluator4)
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
MPQC::CintsEvaluatorFactory_impl::setServices (
  /* in */ ::gov::cca::Services services ) 
throw ( 
  ::gov::cca::CCAException
){
  // DO-NOT-DELETE splicer.begin(MPQC.CintsEvaluatorFactory.setServices)
  // Insert-Code-Here {MPQC.CintsEvaluatorFactory.setServices} (setServices method)
  // DO-NOT-DELETE splicer.end(MPQC.CintsEvaluatorFactory.setServices)
}


// DO-NOT-DELETE splicer.begin(MPQC.CintsEvaluatorFactory._misc)
// Insert-Code-Here {MPQC.CintsEvaluatorFactory._misc} (miscellaneous code)
// DO-NOT-DELETE splicer.end(MPQC.CintsEvaluatorFactory._misc)

