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
#ifndef included_Chemistry_QC_GaussianBasis_CompositeIntegralDescrInterface_hxx
#include "Chemistry_QC_GaussianBasis_CompositeIntegralDescrInterface.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_IntegralDescrInterface_hxx
#include "Chemistry_QC_GaussianBasis_IntegralDescrInterface.hxx"
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
  // Insert-Code-Here {MPQC.Libint2EvaluatorFactory._ctor} (constructor)
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
::Chemistry::QC::GaussianBasis::CompositeIntegralDescrInterface
MPQC::Libint2EvaluatorFactory_impl::get_descriptor_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Libint2EvaluatorFactory.get_descriptor)
  // Insert-Code-Here {MPQC.Libint2EvaluatorFactory.get_descriptor} (get_descriptor method)
    
    // DO-DELETE-WHEN-IMPLEMENTING exception.begin()
    /*
     * This method has not been implemented
     */
    ::sidl::NotImplementedException ex = ::sidl::NotImplementedException::_create();
    ex.setNote("This method has not been implemented");
    ex.add(__FILE__, __LINE__, "get_descriptor");
    throw ex;
    // DO-DELETE-WHEN-IMPLEMENTING exception.end()
    
  // DO-NOT-DELETE splicer.end(MPQC.Libint2EvaluatorFactory.get_descriptor)
}

/**
 * Method:  is_supported[]
 */
bool
MPQC::Libint2EvaluatorFactory_impl::is_supported_impl (
  /* in */::Chemistry::QC::GaussianBasis::IntegralDescrInterface desc ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Libint2EvaluatorFactory.is_supported)
  // Insert-Code-Here {MPQC.Libint2EvaluatorFactory.is_supported} (is_supported method)
    
    // DO-DELETE-WHEN-IMPLEMENTING exception.begin()
    /*
     * This method has not been implemented
     */
    ::sidl::NotImplementedException ex = ::sidl::NotImplementedException::_create();
    ex.setNote("This method has not been implemented");
    ex.add(__FILE__, __LINE__, "is_supported");
    throw ex;
    // DO-DELETE-WHEN-IMPLEMENTING exception.end()
    
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
  // Insert-Code-Here {MPQC.Libint2EvaluatorFactory.set_storage} (set_storage method)
    
    // DO-DELETE-WHEN-IMPLEMENTING exception.begin()
    /*
     * This method has not been implemented
     */
    ::sidl::NotImplementedException ex = ::sidl::NotImplementedException::_create();
    ex.setNote("This method has not been implemented");
    ex.add(__FILE__, __LINE__, "set_storage");
    throw ex;
    // DO-DELETE-WHEN-IMPLEMENTING exception.end()
    
  // DO-NOT-DELETE splicer.end(MPQC.Libint2EvaluatorFactory.set_storage)
}

/**
 *  Get a 1-center integral evaluator
 * @param desc Integral set descriptor
 * @return 1-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator1Interface
MPQC::Libint2EvaluatorFactory_impl::get_evaluator1_impl (
  /* in */::Chemistry::QC::GaussianBasis::CompositeIntegralDescrInterface desc,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs1 ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Libint2EvaluatorFactory.get_evaluator1)
  // Insert-Code-Here {MPQC.Libint2EvaluatorFactory.get_evaluator1} (get_evaluator1 method)
    
    // DO-DELETE-WHEN-IMPLEMENTING exception.begin()
    /*
     * This method has not been implemented
     */
    ::sidl::NotImplementedException ex = ::sidl::NotImplementedException::_create();
    ex.setNote("This method has not been implemented");
    ex.add(__FILE__, __LINE__, "get_evaluator1");
    throw ex;
    // DO-DELETE-WHEN-IMPLEMENTING exception.end()
    
  // DO-NOT-DELETE splicer.end(MPQC.Libint2EvaluatorFactory.get_evaluator1)
}

/**
 *  Get a 2-center integral evaluator
 * @param desc Integral set descriptor
 * @return 2-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator2Interface
MPQC::Libint2EvaluatorFactory_impl::get_evaluator2_impl (
  /* in */::Chemistry::QC::GaussianBasis::CompositeIntegralDescrInterface desc,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs1,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs2 ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Libint2EvaluatorFactory.get_evaluator2)
  // Insert-Code-Here {MPQC.Libint2EvaluatorFactory.get_evaluator2} (get_evaluator2 method)
    
    // DO-DELETE-WHEN-IMPLEMENTING exception.begin()
    /*
     * This method has not been implemented
     */
    ::sidl::NotImplementedException ex = ::sidl::NotImplementedException::_create();
    ex.setNote("This method has not been implemented");
    ex.add(__FILE__, __LINE__, "get_evaluator2");
    throw ex;
    // DO-DELETE-WHEN-IMPLEMENTING exception.end()
    
  // DO-NOT-DELETE splicer.end(MPQC.Libint2EvaluatorFactory.get_evaluator2)
}

/**
 *  Get a 3-center integral evaluator
 * @param desc Integral set descriptor
 * @return 3-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator3Interface
MPQC::Libint2EvaluatorFactory_impl::get_evaluator3_impl (
  /* in */::Chemistry::QC::GaussianBasis::CompositeIntegralDescrInterface desc,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs1,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs2,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs3 ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Libint2EvaluatorFactory.get_evaluator3)
  // Insert-Code-Here {MPQC.Libint2EvaluatorFactory.get_evaluator3} (get_evaluator3 method)
    
    // DO-DELETE-WHEN-IMPLEMENTING exception.begin()
    /*
     * This method has not been implemented
     */
    ::sidl::NotImplementedException ex = ::sidl::NotImplementedException::_create();
    ex.setNote("This method has not been implemented");
    ex.add(__FILE__, __LINE__, "get_evaluator3");
    throw ex;
    // DO-DELETE-WHEN-IMPLEMENTING exception.end()
    
  // DO-NOT-DELETE splicer.end(MPQC.Libint2EvaluatorFactory.get_evaluator3)
}

/**
 *  Get a 4-center integral evaluator
 * @param desc Integral set descriptor
 * @return 4-center integral evaluator 
 */
::Chemistry::QC::GaussianBasis::IntegralEvaluator4Interface
MPQC::Libint2EvaluatorFactory_impl::get_evaluator4_impl (
  /* in */::Chemistry::QC::GaussianBasis::CompositeIntegralDescrInterface desc,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs1,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs2,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs3,
  /* in */::Chemistry::QC::GaussianBasis::MolecularInterface bs4 ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Libint2EvaluatorFactory.get_evaluator4)
  // Insert-Code-Here {MPQC.Libint2EvaluatorFactory.get_evaluator4} (get_evaluator4 method)
    
    // DO-DELETE-WHEN-IMPLEMENTING exception.begin()
    /*
     * This method has not been implemented
     */
    ::sidl::NotImplementedException ex = ::sidl::NotImplementedException::_create();
    ex.setNote("This method has not been implemented");
    ex.add(__FILE__, __LINE__, "get_evaluator4");
    throw ex;
    // DO-DELETE-WHEN-IMPLEMENTING exception.end()
    
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
  // Insert-Code-Here {MPQC.Libint2EvaluatorFactory.finalize} (finalize method)
    
    // DO-DELETE-WHEN-IMPLEMENTING exception.begin()
    /*
     * This method has not been implemented
     */
    ::sidl::NotImplementedException ex = ::sidl::NotImplementedException::_create();
    ex.setNote("This method has not been implemented");
    ex.add(__FILE__, __LINE__, "finalize");
    throw ex;
    // DO-DELETE-WHEN-IMPLEMENTING exception.end()
    
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
  // Insert-Code-Here {MPQC.Libint2EvaluatorFactory.setServices} (setServices method)
    
    // DO-DELETE-WHEN-IMPLEMENTING exception.begin()
    /*
     * This method has not been implemented
     */
    ::sidl::NotImplementedException ex = ::sidl::NotImplementedException::_create();
    ex.setNote("This method has not been implemented");
    ex.add(__FILE__, __LINE__, "setServices");
    throw ex;
    // DO-DELETE-WHEN-IMPLEMENTING exception.end()
    
  // DO-NOT-DELETE splicer.end(MPQC.Libint2EvaluatorFactory.setServices)
}


// DO-NOT-DELETE splicer.begin(MPQC.Libint2EvaluatorFactory._misc)
// Insert-Code-Here {MPQC.Libint2EvaluatorFactory._misc} (miscellaneous code)
// DO-NOT-DELETE splicer.end(MPQC.Libint2EvaluatorFactory._misc)

