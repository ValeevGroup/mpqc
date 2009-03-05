// 
// File:          MPQC_ComponentFactory_Impl.cxx
// Symbol:        MPQC.ComponentFactory-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.ComponentFactory
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 
#include "MPQC_ComponentFactory_Impl.hxx"

// 
// Includes for all method dependencies.
// 
#ifndef included_gov_cca_CCAException_hxx
#include "gov_cca_CCAException.hxx"
#endif
#ifndef included_gov_cca_Component_hxx
#include "gov_cca_Component.hxx"
#endif
#ifndef included_gov_cca_ComponentClassDescription_hxx
#include "gov_cca_ComponentClassDescription.hxx"
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
// DO-NOT-DELETE splicer.begin(MPQC.ComponentFactory._includes)
#include "AllBabelCCA.hxx"
#include "MPQC_ComponentClassDescription.hxx"
#include "ChemistryCXX_IntegralSuperFactory.hxx"
#include "MPQC_IntV3EvaluatorFactory.hxx"
#ifdef HAVE_CINTS
  #include "MPQC_CintsEvaluatorFactory.hxx"
#endif
#ifdef HAVE_LIBINT2
  #include "MPQC_Libint2EvaluatorFactory.hxx"
#endif
// DO-NOT-DELETE splicer.end(MPQC.ComponentFactory._includes)

// speical constructor, used for data wrapping(required).  Do not put code here unless you really know what you're doing!
MPQC::ComponentFactory_impl::ComponentFactory_impl() : StubBase(
  reinterpret_cast< void*>(::MPQC::ComponentFactory::_wrapObj(reinterpret_cast< 
  void*>(this))),false) , _wrapped(true){ 
  // DO-NOT-DELETE splicer.begin(MPQC.ComponentFactory._ctor2)
  // Insert-Code-Here {MPQC.ComponentFactory._ctor2} (ctor2)
  // DO-NOT-DELETE splicer.end(MPQC.ComponentFactory._ctor2)
}

// user defined constructor
void MPQC::ComponentFactory_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.ComponentFactory._ctor)
  //addDescription( "MPQC.IntegralEvaluatorFactory", "IntegralEvaluatorFactory");
  // DO-NOT-DELETE splicer.end(MPQC.ComponentFactory._ctor)
}

// user defined destructor
void MPQC::ComponentFactory_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.ComponentFactory._dtor)
  // add destruction details here
  // DO-NOT-DELETE splicer.end(MPQC.ComponentFactory._dtor)
}

// static class initializer
void MPQC::ComponentFactory_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.ComponentFactory._load)
  // guaranteed to be called at most once before any other method in this class
  // DO-NOT-DELETE splicer.end(MPQC.ComponentFactory._load)
}

// user defined static methods: (none)

// user defined non-static methods:
/**
 * Method:  addDescription[]
 */
void
MPQC::ComponentFactory_impl::addDescription_impl (
  /* in */const ::std::string& className,
  /* in */const ::std::string& classAlias ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.ComponentFactory.addDescription)

  MPQC::ComponentClassDescription cccd;
  cccd = MPQC::ComponentClassDescription::_create();
  cccd.initialize(className, classAlias);
  gov::cca::ComponentClassDescription gcccd = cccd;
  descriptions.push_back(gcccd);

  // DO-NOT-DELETE splicer.end(MPQC.ComponentFactory.addDescription)
}

/**
 *  
 * Collect the currently obtainable class name strings from
 * factories known to the builder and the from the
 * already instantiated components.
 * @return The list of class description, which may be empty, that are
 * known a priori to contain valid values for the className
 * argument of createInstance. 
 * @throws CCAException in the event of error.
 */
::sidl::array< ::gov::cca::ComponentClassDescription>
MPQC::ComponentFactory_impl::getAvailableComponentClasses_impl () 
// throws:
//     ::gov::cca::CCAException
//     ::sidl::RuntimeException

{
  // DO-NOT-DELETE splicer.begin(MPQC.ComponentFactory.getAvailableComponentClasses)

  size_t nd = descriptions.size();
  ::sidl::array< ::gov::cca::ComponentClassDescription> descArray  =
  ::sidl::array< ::gov::cca::ComponentClassDescription>::create1d(nd);
  for (size_t i = 0; i < nd; i++) {
    descArray.set(i, descriptions[i]);
  }
  return descArray;

  // DO-NOT-DELETE splicer.end(MPQC.ComponentFactory.getAvailableComponentClasses)
}

/**
 *  the component instance returned is nil if the name is unknown
 * to the factory. The component is raw: it has been constructed
 * but not initialized via setServices.
 */
::gov::cca::Component
MPQC::ComponentFactory_impl::createComponentInstance_impl (
  /* in */const ::std::string& className ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.ComponentFactory.createComponentInstance)

  if (className == "Chemistry.IntegralSuperFactory") { 
    ChemistryCXX::IntegralSuperFactory x =
      ChemistryCXX::IntegralSuperFactory::_create();
    gov::cca::Component c = x;
    return c;
  }

  else if (className == "MPQC.IntV3EvaluatorFactory") {
    MPQC::IntV3EvaluatorFactory x =
      MPQC::IntV3EvaluatorFactory::_create();
    gov::cca::Component c = x;
    return c;
  }

#ifdef HAVE_CINTS
  else if (className == "MPQC.CintsEvaluatorFactory") {
    MPQC::CintsEvaluatorFactory x =
      MPQC::CintsEvaluatorFactory::_create();
    gov::cca::Component c = x;
    return c;
  }
#endif

#ifdef HAVE_LIBINT2
  else if (className == "MPQC.Libint2EvaluatorFactory") {
    MPQC::Libint2EvaluatorFactory x =
      MPQC::Libint2EvaluatorFactory::_create();
    gov::cca::Component c = x;
    return c;
  }
#endif

  gov::cca::Component dummy;
  return dummy;

  // DO-NOT-DELETE splicer.end(MPQC.ComponentFactory.createComponentInstance)
}

/**
 *  reclaim any resources the factory may have associated with
 * the port it is using. This will occur after the
 * normal component shutdown  (ala componentrelease) is finished. 
 */
void
MPQC::ComponentFactory_impl::destroyComponentInstance_impl (
  /* in */const ::std::string& className,
  /* in */::gov::cca::Component c ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.ComponentFactory.destroyComponentInstance)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.ComponentFactory.destroyComponentInstance)
}


// DO-NOT-DELETE splicer.begin(MPQC.ComponentFactory._misc)
// Put miscellaneous code here
// DO-NOT-DELETE splicer.end(MPQC.ComponentFactory._misc)

