// 
// File:          MPQC_ComponentFactory_Impl.cc
// Symbol:        MPQC.ComponentFactory-v0.2
// Symbol Type:   class
// Babel Version: 0.10.12
// Description:   Server-side implementation for MPQC.ComponentFactory
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// xml-url       = /home/jpkenny/src/mpqc-libint2.build-shared/src/lib/chemistry/cca/../../../../lib/cca/repo/MPQC.ComponentFactory-v0.2.xml
// 
#include "MPQC_ComponentFactory_Impl.hh"

// DO-NOT-DELETE splicer.begin(MPQC.ComponentFactory._includes)
#include "dc/babel/babel-cca/AllBabelCCA.hh"
#include "MPQC_ComponentClassDescription.hh"
#include "MPQC_IntV3EvaluatorFactory.hh"
// DO-NOT-DELETE splicer.end(MPQC.ComponentFactory._includes)

// user-defined constructor.
void MPQC::ComponentFactory_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.ComponentFactory._ctor)
  //addDescription( "MPQC.IntegralEvaluatorFactory", "IntegralEvaluatorFactory");
  // DO-NOT-DELETE splicer.end(MPQC.ComponentFactory._ctor)
}

// user-defined destructor.
void MPQC::ComponentFactory_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.ComponentFactory._dtor)
  // add destruction details here
  // DO-NOT-DELETE splicer.end(MPQC.ComponentFactory._dtor)
}

// static class initializer.
void MPQC::ComponentFactory_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.ComponentFactory._load)
  // guaranteed to be called at most once before any other method in this class
  // DO-NOT-DELETE splicer.end(MPQC.ComponentFactory._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
/**
 * Method:  addDescription[]
 */
void
MPQC::ComponentFactory_impl::addDescription (
  /* in */ const ::std::string& className,
  /* in */ const ::std::string& classAlias ) 
throw () 
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
 *  Collect the currently obtainable class name strings from
 *  factories known to the builder and the from the
 *  already instantiated components.
 *  @return The list of class description, which may be empty, that are
 *   known a priori to contain valid values for the className
 *  argument of createInstance. 
 *  @throws CCAException in the event of error.
 */
::sidl::array< ::gov::cca::ComponentClassDescription>
MPQC::ComponentFactory_impl::getAvailableComponentClasses ()
throw ( 
  ::gov::cca::CCAException
)
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
 * the component instance returned is nil if the name is unknown
 * to the factory. The component is raw: it has been constructed
 * but not initialized via setServices.
 */
::gov::cca::Component
MPQC::ComponentFactory_impl::createComponentInstance (
  /* in */ const ::std::string& className ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.ComponentFactory.createComponentInstance)

  if (className == "MPQC.IntV3EvaluatorFactory") {
    MPQC::IntV3EvaluatorFactory x =
      MPQC::IntV3EvaluatorFactory::_create();
    gov::cca::Component c = x;
    return c;
  }

  gov::cca::Component dummy;
  return dummy;

  // DO-NOT-DELETE splicer.end(MPQC.ComponentFactory.createComponentInstance)
}

/**
 * reclaim any resources the factory may have associated with
 * the port it is using. This will occur after the
 * normal component shutdown  (ala componentrelease) is finished. 
 */
void
MPQC::ComponentFactory_impl::destroyComponentInstance (
  /* in */ const ::std::string& className,
  /* in */ ::gov::cca::Component c ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.ComponentFactory.destroyComponentInstance)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.ComponentFactory.destroyComponentInstance)
}


// DO-NOT-DELETE splicer.begin(MPQC.ComponentFactory._misc)
// Put miscellaneous code here
// DO-NOT-DELETE splicer.end(MPQC.ComponentFactory._misc)

