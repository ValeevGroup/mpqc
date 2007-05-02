// 
// File:          MPQC_CoordinateModel_Impl.hxx
// Symbol:        MPQC.CoordinateModel-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.CoordinateModel
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 

#ifndef included_MPQC_CoordinateModel_Impl_hxx
#define included_MPQC_CoordinateModel_Impl_hxx

#ifndef included_sidl_cxx_hxx
#include "sidl_cxx.hxx"
#endif
#ifndef included_MPQC_CoordinateModel_IOR_h
#include "MPQC_CoordinateModel_IOR.h"
#endif
#ifndef included_Chemistry_QC_ModelInterface_hxx
#include "Chemistry_QC_ModelInterface.hxx"
#endif
#ifndef included_ChemistryOpt_CoordinateModelInterface_hxx
#include "ChemistryOpt_CoordinateModelInterface.hxx"
#endif
#ifndef included_MPQC_CoordinateModel_hxx
#include "MPQC_CoordinateModel.hxx"
#endif
#ifndef included_gov_cca_CCAException_hxx
#include "gov_cca_CCAException.hxx"
#endif
#ifndef included_gov_cca_Component_hxx
#include "gov_cca_Component.hxx"
#endif
#ifndef included_gov_cca_Port_hxx
#include "gov_cca_Port.hxx"
#endif
#ifndef included_gov_cca_Services_hxx
#include "gov_cca_Services.hxx"
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
#ifndef included_sidl_RuntimeException_hxx
#include "sidl_RuntimeException.hxx"
#endif


// DO-NOT-DELETE splicer.begin(MPQC.CoordinateModel._includes)

#include <Chemistry_MoleculeViewerInterface.hxx>
#include <Chemistry_QC_ModelFactoryInterface.hxx>
#include <chemistry/molecule/coor.h>
#include "CoordinateModel.hxx"
#include "Chemistry_MoleculeInterface.hxx"
#include <gov_cca_ports_ParameterPortFactory.hxx>
#include <gov_cca_ports_ParameterPort.hxx>

// DO-NOT-DELETE splicer.end(MPQC.CoordinateModel._includes)

namespace MPQC { 

  /**
   * Symbol "MPQC.CoordinateModel" (version 0.2)
   */
  class CoordinateModel_impl : public virtual ::MPQC::CoordinateModel 
  // DO-NOT-DELETE splicer.begin(MPQC.CoordinateModel._inherits)
  ,public CcaChemGeneric::CoordinateModel
  // DO-NOT-DELETE splicer.end(MPQC.CoordinateModel._inherits)
  {

  // All data marked protected will be accessable by 
  // descendant Impl classes
  protected:

    bool _wrapped;

    // DO-NOT-DELETE splicer.begin(MPQC.CoordinateModel._implementation)

    gov::cca::Services services_;
    gov::cca::TypeMap tm_;
    gov::cca::ports::ParameterPortFactory ppf_;
    gov::cca::ports::ParameterPort pp_;

    CcaChemGeneric::CoordinateModel genericModel_;
    Chemistry::QC::ModelInterface model_;
    Chemistry::MoleculeInterface molecule_;
    sc::Ref<sc::MolecularCoor> scCoor_;
    sc::Ref<sc::Molecule> scMol_;
    sc::Ref<sc::SCMatrixKit> kit_;
    sc::Ref<sc::SCMatrixKit> rkit_;
    sc::RefSymmSCMatrix ihess_;
    double grad_rms_, grad_max_, disp_rms_, disp_max_, diagonal_scale_;
    bool multiple_guess_h_, use_current_geom_, use_diagonal_h_;
    std::string coordinates_, extra_bonds_;
    double convFrom_;
    bool have_guess_h_;
    enum {cart,symm,redund};
    int coorType_;
    int numCoor_;
    int natom3_;

    void draw();

    // DO-NOT-DELETE splicer.end(MPQC.CoordinateModel._implementation)

  public:
    // default constructor, used for data wrapping(required)
    CoordinateModel_impl();
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    CoordinateModel_impl( struct MPQC_CoordinateModel__object * s ) : StubBase(
      s,true), _wrapped(false) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~CoordinateModel_impl() { _dtor(); }

    // user defined destruction
    void _dtor();

    // true if this object was created by a user newing the impl
    inline bool _isWrapped() {return _wrapped;}

    // static class initializer
    static void _load();

  public:


    /**
     *  Registers and gets ports, and requests Model object(s) from the 
     * ModelFactory component(s). This must be the first method called 
     * following instantiation.
     */
    int32_t
    initialize_impl() ;

    /**
     *  Releases and unregisters ports.  This should be called when the
     * CoordinateModel object is no longer needed.
     */
    int32_t
    finalize_impl() ;

    /**
     *  Sets the contained chemistry Model object (currently unused as the
     * chemistry Model object is normally obtained from a ModelFactory 
     * during initialization).
     * @param model The chemistry model object.
     */
    void
    set_model_impl (
      /* in */::Chemistry::QC::ModelInterface model
    )
    ;


    /**
     *  Returns the contained chemistry Model object.
     * @return The chemistry Model object.
     */
    ::Chemistry::QC::ModelInterface
    get_model_impl() ;

    /**
     *  Returns the number of coordinates.
     * @return The number of coordinates. 
     */
    int32_t
    get_n_coor_impl() ;

    /**
     *  Returns the array of (cartesian or internal) coordinates which are 
     * being  optimized.
     * @return The array of coordinates which are being optimized.
     */
    ::sidl::array<double>
    get_coor_impl() ;

    /**
     *  Returns the energy of the currently contained model with the values
     * of the optimization coordinates given in x.  This requires
     * that the CoordinateModel updates the cartesian coordinates of a 
     * contained Molecule object (possibly requiring transformation) and set 
     * this Molecule object on a contained Model object, prior to calling
     * get_energy() on the Model object.
     * @param x The optimization coordinate values.
     * @return The energy of the chemistry model at x.
     */
    double
    get_energy_impl (
      /* in array<double> */::sidl::array<double> x
    )
    ;


    /**
     *  Returns the energy gradient of the currently contained model with 
     * the values of the optimization coordinates given in x.  This requires
     * that the CoordinateModel updates the cartesian coordinates of a
     * contained Molecule object (possibly requiring transformation) and set
     * this Molecule object on a contained Model object, prior to calling
     * get_gradient() on the Model object.  If the optimization coordinate
     * system is not cartesian, the gradient is transformed.
     * @param x The optimization coordinate values.
     * @return The energy gradient of the chemistry model at x.
     */
    ::sidl::array<double>
    get_gradient_impl (
      /* in array<double> */::sidl::array<double> x
    )
    ;


    /**
     *  Returns the energy Hessian of the currently contained model with
     * the values of the optimization coordinates given in x.  This requires
     * that the CoordinateModel updates the cartesian coordinates of a
     * contained Molecule object (possibly requiring transformation) and set
     * this Molecule object on a contained Model object, prior to calling
     * get_hessian() on the Model object.  If the optimization coordinate
     * system is not cartesian, the Hessian is transformed.
     * @param x The optimization coordinate values.
     * @return The energy Hessian of the chemistry model at x.
     */
    ::sidl::array<double>
    get_hessian_impl (
      /* in array<double> */::sidl::array<double> x
    )
    ;


    /**
     *  Sets f and g to the energy and energy gradient, respectively,
     * of the chemistry model at x.  This is similar to calling
     * get_energy() and get_gradient() separately, but set_molecule()
     * must be called on the Model object only once.  This is necessary
     * for some model implementations, as a second molecule update
     * would invalidate results from an energy computation.  An alternative
     * would be to always return the energy as well when get_gradient() is 
     * called.
     * @param x The optimization coordinate values.
     * @param f Variable that energy will be assigned to.
     * @param g Array that the gradient will be assigned to.
     */
    void
    get_energy_and_gradient_impl (
      /* in array<double> */::sidl::array<double> x,
      /* out */double& f,
      /* in array<double> */::sidl::array<double> g
    )
    ;


    /**
     *  Returns the product of the guess hessian inverse and an effective
     * gradient.  Probably unique to TAO's limited memory variable metric
     * algorithm, which uses this method to accomodate dense guess hessians.
     * "first_geom_ptr" provides the Cartesian coordinates for which the
     * guess Hessian should be computed (first_geom_ptr=0 for current
     * geometry).
     * @param effective_grad An effective gradient.
     * @param effective_step Array that effective step is assigned to.
     * @param first_geom     Pointer to array of Cartesians 
     */
    void
    guess_hessian_solve_impl (
      /* in array<double> */::sidl::array<double> effective_grad,
      /* in array<double> */::sidl::array<double> effective_step,
      /* in */void* first_geom
    )
    ;


    /**
     *  Determines if the optimization has converged, flag is set to 1
     * if convergence has been achieved and 0 otherwise.
     * @param flag Variable that convergence value is assigned to.
     */
    void
    checkConvergence_impl (
      /* inout */int32_t& flag
    )
    ;


    /**
     *  For visualization, possibly unused (?).  CoordinateModel objects
     * may callback to viewers that implement the Chemistry.MoleculeViewer 
     * interface, such as the cca-chem python GUI, making this method 
     * unnecessary.
     */
    void
    monitor_impl() ;

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
    setServices_impl (
      /* in */::gov::cca::Services services
    )
    // throws:
    //     ::gov::cca::CCAException
    //     ::sidl::RuntimeException
    ;

  };  // end class CoordinateModel_impl

} // end namespace MPQC

// DO-NOT-DELETE splicer.begin(MPQC.CoordinateModel._misc)
// Insert-Code-Here {MPQC.CoordinateModel._misc} (miscellaneous things)
// DO-NOT-DELETE splicer.end(MPQC.CoordinateModel._misc)

#endif
