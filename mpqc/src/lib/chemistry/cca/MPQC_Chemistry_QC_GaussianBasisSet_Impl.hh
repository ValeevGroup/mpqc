// 
// File:          MPQC_Chemistry_QC_GaussianBasisSet_Impl.hh
// Symbol:        MPQC.Chemistry_QC_GaussianBasisSet-v0.2
// Symbol Type:   class
// Babel Version: 0.9.8
// Description:   Server-side implementation for MPQC.Chemistry_QC_GaussianBasisSet
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.9.8
// 

#ifndef included_MPQC_Chemistry_QC_GaussianBasisSet_Impl_hh
#define included_MPQC_Chemistry_QC_GaussianBasisSet_Impl_hh

#ifndef included_sidl_cxx_hh
#include "sidl_cxx.hh"
#endif
#ifndef included_MPQC_Chemistry_QC_GaussianBasisSet_IOR_h
#include "MPQC_Chemistry_QC_GaussianBasisSet_IOR.h"
#endif
// 
// Includes for all method dependencies.
// 
#ifndef included_Chemistry_Molecule_hh
#include "Chemistry_Molecule.hh"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_AngularType_hh
#include "Chemistry_QC_GaussianBasis_AngularType.hh"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_Atomic_hh
#include "Chemistry_QC_GaussianBasis_Atomic.hh"
#endif
#ifndef included_MPQC_Chemistry_QC_GaussianBasisSet_hh
#include "MPQC_Chemistry_QC_GaussianBasisSet.hh"
#endif
#ifndef included_sidl_BaseInterface_hh
#include "sidl_BaseInterface.hh"
#endif
#ifndef included_sidl_ClassInfo_hh
#include "sidl_ClassInfo.hh"
#endif


// DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_GaussianBasisSet._includes)
// Put additional includes or other arbitrary code here...
// DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_GaussianBasisSet._includes)

namespace MPQC { 

  /**
   * Symbol "MPQC.Chemistry_QC_GaussianBasisSet" (version 0.2)
   */
  class Chemistry_QC_GaussianBasisSet_impl
  // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_GaussianBasisSet._inherits)
  // Put additional inheritance here...
  // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_GaussianBasisSet._inherits)
  {

  private:
    // Pointer back to IOR.
    // Use this to dispatch back through IOR vtable.
    Chemistry_QC_GaussianBasisSet self;

    // DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_GaussianBasisSet._implementation)
    // Put additional implementation details here...
    // DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_GaussianBasisSet._implementation)

  private:
    // private default constructor (required)
    Chemistry_QC_GaussianBasisSet_impl() {} 

  public:
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    Chemistry_QC_GaussianBasisSet_impl( struct 
      MPQC_Chemistry_QC_GaussianBasisSet__object * s ) : self(s,
      true) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~Chemistry_QC_GaussianBasisSet_impl() { _dtor(); }

    // user defined destruction
    void _dtor();

  public:

    /**
     * user defined non-static method.
     */
    ::std::string
    get_label() throw () 
    ;
    /**
     * user defined non-static method.
     */
    int64_t
    get_n_basis() throw () 
    ;
    /**
     * user defined non-static method.
     */
    int64_t
    get_n_shell() throw () 
    ;
    /**
     * user defined non-static method.
     */
    int64_t
    get_max_angular_momentum() throw () 
    ;
    /**
     * user defined non-static method.
     */
    ::Chemistry::QC::GaussianBasis::AngularType
    get_angular_type() throw () 
    ;
    /**
     * user defined non-static method.
     */
    ::Chemistry::QC::GaussianBasis::Atomic
    get_atomic (
      /*in*/ int64_t atomnum
    )
    throw () 
    ;

    /**
     * user defined non-static method.
     */
    ::Chemistry::Molecule
    get_molecule() throw () 
    ;
  };  // end class Chemistry_QC_GaussianBasisSet_impl

} // end namespace MPQC

// DO-NOT-DELETE splicer.begin(MPQC.Chemistry_QC_GaussianBasisSet._misc)
// Put miscellaneous things here...
// DO-NOT-DELETE splicer.end(MPQC.Chemistry_QC_GaussianBasisSet._misc)

#endif
