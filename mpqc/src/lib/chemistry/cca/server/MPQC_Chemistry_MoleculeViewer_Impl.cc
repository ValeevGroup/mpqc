// 
// File:          MPQC_Chemistry_MoleculeViewer_Impl.cc
// Symbol:        MPQC.Chemistry_MoleculeViewer-v0.2
// Symbol Type:   class
// Babel Version: 0.10.12
// Description:   Server-side implementation for MPQC.Chemistry_MoleculeViewer
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// 
#include "MPQC_Chemistry_MoleculeViewer_Impl.hh"

// DO-NOT-DELETE splicer.begin(MPQC.Chemistry_MoleculeViewer._includes)
#include <sstream>
#include <iostream>
// DO-NOT-DELETE splicer.end(MPQC.Chemistry_MoleculeViewer._includes)

// user-defined constructor.
void MPQC::Chemistry_MoleculeViewer_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_MoleculeViewer._ctor)
  is_updated=0;
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_MoleculeViewer._ctor)
}

// user-defined destructor.
void MPQC::Chemistry_MoleculeViewer_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_MoleculeViewer._dtor)

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

  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_MoleculeViewer._dtor)
}

// static class initializer.
void MPQC::Chemistry_MoleculeViewer_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_MoleculeViewer._load)
  // guaranteed to be called at most once before any other method in this class
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_MoleculeViewer._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
/**
 * Method:  set_molecule[]
 */
void
MPQC::Chemistry_MoleculeViewer_impl::set_molecule (
  /* in */ ::Chemistry::Molecule molecule ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_MoleculeViewer.set_molecule)
  molecule_ = molecule;
  is_updated=1;
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_MoleculeViewer.set_molecule)
}

/**
 * Method:  set_coor[]
 */
void
MPQC::Chemistry_MoleculeViewer_impl::set_coor (
  /* in */ const ::std::string& coords ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_MoleculeViewer.set_coor)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_MoleculeViewer.set_coor)
}

/**
 * Method:  run_gui[]
 */
void
MPQC::Chemistry_MoleculeViewer_impl::run_gui ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_MoleculeViewer.run_gui)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_MoleculeViewer.run_gui)
}

/**
 * Method:  draw[]
 */
void
MPQC::Chemistry_MoleculeViewer_impl::draw ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_MoleculeViewer.draw)
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
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_MoleculeViewer.draw)
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
MPQC::Chemistry_MoleculeViewer_impl::setServices (
  /* in */ ::gov::cca::Services services ) 
throw ( 
  ::gov::cca::CCAException
){
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_MoleculeViewer.setServices)
  services.addProvidesPort(self,"MoleculeViewer","Chemistry.MoleculeViewer",
                           services.createTypeMap());
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_MoleculeViewer.setServices)
}


// DO-NOT-DELETE splicer.begin(MPQC.Chemistry_MoleculeViewer._misc)
// Put miscellaneous code here
// DO-NOT-DELETE splicer.end(MPQC.Chemistry_MoleculeViewer._misc)

