// 
// File:          MPQC_GaussianBasisMolecular_Impl.hxx
// Symbol:        MPQC.GaussianBasisMolecular-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.GaussianBasisMolecular
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 

#ifndef included_MPQC_GaussianBasisMolecular_Impl_hxx
#define included_MPQC_GaussianBasisMolecular_Impl_hxx

#ifndef included_sidl_cxx_hxx
#include "sidl_cxx.hxx"
#endif
#ifndef included_MPQC_GaussianBasisMolecular_IOR_h
#include "MPQC_GaussianBasisMolecular_IOR.h"
#endif
#ifndef included_Chemistry_MoleculeInterface_hxx
#include "Chemistry_MoleculeInterface.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_AngularType_hxx
#include "Chemistry_QC_GaussianBasis_AngularType.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_AtomicInterface_hxx
#include "Chemistry_QC_GaussianBasis_AtomicInterface.hxx"
#endif
#ifndef included_Chemistry_QC_GaussianBasis_MolecularInterface_hxx
#include "Chemistry_QC_GaussianBasis_MolecularInterface.hxx"
#endif
#ifndef included_MPQC_GaussianBasisMolecular_hxx
#include "MPQC_GaussianBasisMolecular.hxx"
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


// DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisMolecular._includes)

#include <chemistry/qc/basis/basis.h>
#include <ChemistryCXX_Molecule.hxx>
#include <MPQC_GaussianBasisAtomic_Impl.hxx>

// DO-NOT-DELETE splicer.end(MPQC.GaussianBasisMolecular._includes)

namespace MPQC { 

  /**
   * Symbol "MPQC.GaussianBasisMolecular" (version 0.2)
   */
  class GaussianBasisMolecular_impl : public virtual 
    ::MPQC::GaussianBasisMolecular 
  // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisMolecular._inherits)
  // Insert-Code-Here {MPQC.GaussianBasisMolecular._inherits} (optional inheritance here)
  // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisMolecular._inherits)
  {

  // All data marked protected will be accessable by 
  // descendant Impl classes
  protected:

    bool _wrapped;

    // DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisMolecular._implementation)

    sc::GaussianBasisSet *gbs_ptr_;
    sc::Ref<sc::GaussianBasisSet> sc_gbs_;
    std::string label_;
    Chemistry::QC::GaussianBasis::AngularType angular_type_;
    ChemistryCXX::Molecule molecule_;
    MPQC::GaussianBasisAtomic *atomic_array_;
    int natom_;

    // DO-NOT-DELETE splicer.end(MPQC.GaussianBasisMolecular._implementation)

  public:
    // default constructor, used for data wrapping(required)
    GaussianBasisMolecular_impl();
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    GaussianBasisMolecular_impl( struct MPQC_GaussianBasisMolecular__object * s 
      ) : StubBase(s,true), _wrapped(false) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~GaussianBasisMolecular_impl() { _dtor(); }

    // user defined destruction
    void _dtor();

    // true if this object was created by a user newing the impl
    inline bool _isWrapped() {return _wrapped;}

    // static class initializer
    static void _load();

  public:

    /**
     * user defined non-static method.
     */
    void
    initialize_impl (
      /* in */void* scbasis,
      /* in */const ::std::string& label
    )
    ;

    /**
     * user defined non-static method.
     */
    void*
    sc_gbs_pointer_impl() ;

    /**
     *  Get the user specified name.
     * @return name 
     */
    ::std::string
    get_label_impl() ;

    /**
     *  Get the number of basis functions.
     * @return number of functions 
     */
    int64_t
    get_n_basis_impl() ;

    /**
     *  Get the number of shells.
     * @return number of shells 
     */
    int64_t
    get_n_shell_impl() ;

    /**
     *  Get the max angular momentum for any contraction in the 
     * basis set.
     * @return max angular momentum 
     */
    int32_t
    get_max_angular_momentum_impl() ;

    /**
     *  Get the angular type.
     * @return enum AngularType 
     */
    ::Chemistry::QC::GaussianBasis::AngularType
    get_angular_type_impl() ;

    /**
     *  Get an atomic basis set.
     * @param atomnum atom number 
     * @return Atomic 
     */
    ::Chemistry::QC::GaussianBasis::AtomicInterface
    get_atomic_impl (
      /* in */int64_t atomnum
    )
    ;


    /**
     *  Get the molecule.
     * @return Molecule 
     */
    ::Chemistry::MoleculeInterface
    get_molecule_impl() ;

    /**
     *  Print the molecular basis data. 
     */
    void
    print_molecular_impl() ;
  };  // end class GaussianBasisMolecular_impl

} // end namespace MPQC

// DO-NOT-DELETE splicer.begin(MPQC.GaussianBasisMolecular._misc)
// Insert-Code-Here {MPQC.GaussianBasisMolecular._misc} (miscellaneous things)
// DO-NOT-DELETE splicer.end(MPQC.GaussianBasisMolecular._misc)

#endif
