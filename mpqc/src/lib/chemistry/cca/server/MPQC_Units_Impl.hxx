// 
// File:          MPQC_Units_Impl.hxx
// Symbol:        MPQC.Units-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.Units
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 

#ifndef included_MPQC_Units_Impl_hxx
#define included_MPQC_Units_Impl_hxx

#ifndef included_sidl_cxx_hxx
#include "sidl_cxx.hxx"
#endif
#ifndef included_MPQC_Units_IOR_h
#include "MPQC_Units_IOR.h"
#endif
#ifndef included_MPQC_Units_hxx
#include "MPQC_Units.hxx"
#endif
#ifndef included_Physics_UnitsInterface_hxx
#include "Physics_UnitsInterface.hxx"
#endif
#ifndef included_sidl_BaseClass_hxx
#include "sidl_BaseClass.hxx"
#endif
#ifndef included_sidl_BaseInterface_hxx
#include "sidl_BaseInterface.hxx"
#endif
#ifndef included_sidl_ClassInfo_hxx
#include "sidl_ClassInfo.hxx"
#endif


// DO-NOT-DELETE splicer.begin(MPQC.Units._hincludes)
#include <util/misc/units.h>
// DO-NOT-DELETE splicer.end(MPQC.Units._hincludes)

namespace MPQC { 

  /**
   * Symbol "MPQC.Units" (version 0.2)
   */
  class Units_impl : public virtual ::MPQC::Units 
  // DO-NOT-DELETE splicer.begin(MPQC.Units._inherits)
  // Insert-Code-Here {MPQC.Units._inherits} (optional inheritance here)
  // DO-NOT-DELETE splicer.end(MPQC.Units._inherits)

  {

  // All data marked protected will be accessable by 
  // descendant Impl classes
  protected:

    bool _wrapped;

    // DO-NOT-DELETE splicer.begin(MPQC.Units._implementation)
      sc::Ref<sc::Units> units;
    public:
      void set_units(const sc::Ref<sc::Units> &);
    // DO-NOT-DELETE splicer.end(MPQC.Units._implementation)

  public:
    // default constructor, used for data wrapping(required)
    Units_impl();
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
      Units_impl( struct MPQC_Units__object * ior ) : StubBase(ior,true), 
    ::Physics::UnitsInterface((ior==NULL) ? NULL : &((
      *ior).d_physics_unitsinterface)) , _wrapped(false) {
      ior->d_data = this;
      _ctor();
    }


    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~Units_impl() { _dtor(); }

    // user defined destruction
    void _dtor();

    // true if this object was created by a user newing the impl
    inline bool _isWrapped() {return _wrapped;}

    // static class initializer
    static void _load();

  public:


    /**
     *  Initializes the units as a human readable string
     * options are "angstroms" or "bohr" 
     */
    void
    initialize_impl (
      /* in */const ::std::string& unitname
    )
    ;


    /**
     *  Returns the units as a human readable string. 
     */
    ::std::string
    get_unit_name_impl() ;

    /**
     *  Converts from self's units to the given unit name. 
     */
    double
    convert_to_impl (
      /* in */const ::std::string& unitname
    )
    ;


    /**
     *  Converts to self's units from the given unit name. 
     */
    double
    convert_from_impl (
      /* in */const ::std::string& unitname
    )
    ;

  };  // end class Units_impl

} // end namespace MPQC

// DO-NOT-DELETE splicer.begin(MPQC.Units._hmisc)
// Insert-Code-Here {MPQC.Units._misc} (miscellaneous things)
// DO-NOT-DELETE splicer.end(MPQC.Units._hmisc)

#endif
