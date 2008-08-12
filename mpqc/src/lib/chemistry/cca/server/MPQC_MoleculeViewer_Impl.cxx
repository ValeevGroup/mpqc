// 
// File:          MPQC_MoleculeViewer_Impl.cxx
// Symbol:        MPQC.MoleculeViewer-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.MoleculeViewer
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 
#include "MPQC_MoleculeViewer_Impl.hxx"

// 
// Includes for all method dependencies.
// 
#ifndef included_Chemistry_MoleculeInterface_hxx
#include "Chemistry_MoleculeInterface.hxx"
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
// DO-NOT-DELETE splicer.begin(MPQC.MoleculeViewer._includes)

#include <sstream>
#include <iostream>

// DO-NOT-DELETE splicer.end(MPQC.MoleculeViewer._includes)

// speical constructor, used for data wrapping(required).  Do not put code here unless you really know what you're doing!
MPQC::MoleculeViewer_impl::MoleculeViewer_impl() : StubBase(reinterpret_cast< 
  void*>(::MPQC::MoleculeViewer::_wrapObj(reinterpret_cast< void*>(this))),
  false) , _wrapped(true){ 
  // DO-NOT-DELETE splicer.begin(MPQC.MoleculeViewer._ctor2)
  // Insert-Code-Here {MPQC.MoleculeViewer._ctor2} (ctor2)
  // DO-NOT-DELETE splicer.end(MPQC.MoleculeViewer._ctor2)
}

// user defined constructor
void MPQC::MoleculeViewer_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.MoleculeViewer._ctor)

  is_updated=0;

  // DO-NOT-DELETE splicer.end(MPQC.MoleculeViewer._ctor)
}

// user defined destructor
void MPQC::MoleculeViewer_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.MoleculeViewer._dtor)

#if USE_SOCKET
  try {
      socket_.close();
    }
  catch (std::exception &e) {
      std::cout << "NOTE: could not close viewer connection: "
                << e.what()
                << std::endl;
    }
#endif // USE_SOCKET

  // DO-NOT-DELETE splicer.end(MPQC.MoleculeViewer._dtor)
}

// static class initializer
void MPQC::MoleculeViewer_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.MoleculeViewer._load)
  // Insert-Code-Here {MPQC.MoleculeViewer._load} (class initialization)
  // DO-NOT-DELETE splicer.end(MPQC.MoleculeViewer._load)
}

// user defined static methods: (none)

// user defined non-static methods:
/**
 * Method:  set_molecule[]
 */
void
MPQC::MoleculeViewer_impl::set_molecule_impl (
  /* in */::Chemistry::MoleculeInterface molecule ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.MoleculeViewer.set_molecule)

  molecule_ = molecule;
  is_updated=1;

  // DO-NOT-DELETE splicer.end(MPQC.MoleculeViewer.set_molecule)
}

/**
 * Method:  set_coor[]
 */
void
MPQC::MoleculeViewer_impl::set_coor_impl (
  /* in */const ::std::string& coords ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.MoleculeViewer.set_coor)
  // Insert-Code-Here {MPQC.MoleculeViewer.set_coor} (set_coor method)
  // 
  // This method has not been implemented
  // 
    ::sidl::NotImplementedException ex = ::sidl::NotImplementedException::_create();
    ex.setNote("This method has not been implemented");
    ex.add(__FILE__, __LINE__, "set_coor");
    throw ex;
  // DO-NOT-DELETE splicer.end(MPQC.MoleculeViewer.set_coor)
}

/**
 * Method:  run_gui[]
 */
void
MPQC::MoleculeViewer_impl::run_gui_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.MoleculeViewer.run_gui)
  // Insert-Code-Here {MPQC.MoleculeViewer.run_gui} (run_gui method)
  // 
  // This method has not been implemented
  // 
    ::sidl::NotImplementedException ex = ::sidl::NotImplementedException::_create();
    ex.setNote("This method has not been implemented");
    ex.add(__FILE__, __LINE__, "run_gui");
    throw ex;
  // DO-NOT-DELETE splicer.end(MPQC.MoleculeViewer.run_gui)
}

/**
 * Method:  draw[]
 */
void
MPQC::MoleculeViewer_impl::draw_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.MoleculeViewer.draw)

  if (molecule_._is_nil()) {
      return;
    }

  std::cout << "drawing" << std::endl;

  for (int i=0; i<molecule_.get_n_atom(); i++) {
      std::cout << "====="
                << " " << molecule_.get_cart_coor(i,0)
                << " " << molecule_.get_cart_coor(i,1)
                << " " << molecule_.get_cart_coor(i,2)
                << std::endl;
    }

#if USE_SOCKET
  try {

      if (!socket_.initialized()) socket_.create();
      if (!socket_.bound()) socket_.bind(10002,12000);
      if (!socket_.connected()) socket_.connect("localhost", 10001);

      std::ostringstream str;
      const double f = 1.0/(1.0e-10/5.29177249e-11); // m
      for (int i=0; i<molecule_.get_n_atom(); i++) {
          str << " " << f*molecule_.get_cart_coor(i,0)
              << " " << f*molecule_.get_cart_coor(i,1)
              << " " << f*molecule_.get_cart_coor(i,2);
        }
      str << " ";
      socket_.write(str.str().c_str(), str.str().size());
      std::cout << "sent to socket: " << str.str() << std::endl;
    }
  catch (std::exception &e) {
      std::cout << "NOTE: Could not draw to viewer: "
                << e.what()
                << std::endl;
    }
#endif // USE_SOCKET

  // DO-NOT-DELETE splicer.end(MPQC.MoleculeViewer.draw)
}

/**
 *  This should be called when the object is no longer needed.
 * No other members may be called after finalize. 
 */
int32_t
MPQC::MoleculeViewer_impl::finalize_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.MoleculeViewer.finalize)
  // Insert-Code-Here {MPQC.MoleculeViewer.finalize} (finalize method)
  // 
  // This method has not been implemented
  // 
    ::sidl::NotImplementedException ex = ::sidl::NotImplementedException::_create();
    ex.setNote("This method has not been implemented");
    ex.add(__FILE__, __LINE__, "finalize");
    throw ex;
  // DO-NOT-DELETE splicer.end(MPQC.MoleculeViewer.finalize)
}

/**
 *  Starts up a component presence in the calling framework.
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
MPQC::MoleculeViewer_impl::setServices_impl (
  /* in */::gov::cca::Services services ) 
// throws:
//     ::gov::cca::CCAException
//     ::sidl::RuntimeException
{
  // DO-NOT-DELETE splicer.begin(MPQC.MoleculeViewer.setServices)

  services.addProvidesPort(*this,"MoleculeViewer","Chemistry.MoleculeViewer",
                           services.createTypeMap());

  // DO-NOT-DELETE splicer.end(MPQC.MoleculeViewer.setServices)
}


// DO-NOT-DELETE splicer.begin(MPQC.MoleculeViewer._misc)
// Insert-Code-Here {MPQC.MoleculeViewer._misc} (miscellaneous code)
// DO-NOT-DELETE splicer.end(MPQC.MoleculeViewer._misc)

