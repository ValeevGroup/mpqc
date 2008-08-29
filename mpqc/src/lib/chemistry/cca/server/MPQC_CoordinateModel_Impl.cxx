// 
// File:          MPQC_CoordinateModel_Impl.cxx
// Symbol:        MPQC.CoordinateModel-v0.2
// Symbol Type:   class
// Description:   Server-side implementation for MPQC.CoordinateModel
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// 
#include "MPQC_CoordinateModel_Impl.hxx"

// 
// Includes for all method dependencies.
// 
#ifndef included_Chemistry_QC_ModelInterface_hxx
#include "Chemistry_QC_ModelInterface.hxx"
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
// DO-NOT-DELETE splicer.begin(MPQC.CoordinateModel._includes)

#include <iostream>
#include <sstream>
#include "util/keyval/keyvalval.h"
#include "math/scmat/matrix.h"
#include "math/scmat/local.h"
#include "math/scmat/repl.h"

#define DEFAULT_COORTYPE symm

sidl::array<double>
vector_to_array(const sc::RefSCVector &v);

sc::RefSCVector
array_to_vector(sidl::array<double>, const sc::RefSCVector &v);

sidl::array<double>
matrix_to_array(const sc::RefSymmSCMatrix &v);

sc::RefSymmSCMatrix
array_to_matrix(sidl::array<double>, const sc::RefSymmSCMatrix &v);

// DO-NOT-DELETE splicer.end(MPQC.CoordinateModel._includes)

// speical constructor, used for data wrapping(required).  Do not put code here unless you really know what you're doing!
MPQC::CoordinateModel_impl::CoordinateModel_impl() : StubBase(reinterpret_cast< 
  void*>(::MPQC::CoordinateModel::_wrapObj(reinterpret_cast< void*>(this))),
  false) , _wrapped(true){ 
  // DO-NOT-DELETE splicer.begin(MPQC.CoordinateModel._ctor2)
  // Insert-Code-Here {MPQC.CoordinateModel._ctor2} (ctor2)
  // DO-NOT-DELETE splicer.end(MPQC.CoordinateModel._ctor2)
}

// user defined constructor
void MPQC::CoordinateModel_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.CoordinateModel._ctor)
  // DO-NOT-DELETE splicer.end(MPQC.CoordinateModel._ctor)
}

// user defined destructor
void MPQC::CoordinateModel_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.CoordinateModel._dtor)
  // Insert-Code-Here {MPQC.CoordinateModel._dtor} (destructor)
  // DO-NOT-DELETE splicer.end(MPQC.CoordinateModel._dtor)
}

// static class initializer
void MPQC::CoordinateModel_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.CoordinateModel._load)
  // Insert-Code-Here {MPQC.CoordinateModel._load} (class initialization)
  // DO-NOT-DELETE splicer.end(MPQC.CoordinateModel._load)
}

// user defined static methods: (none)

// user defined non-static methods:
/**
 *  Gets ports, and requests Model object(s) from the 
 * ModelFactory component(s). This must be the first method called 
 * following instantiation.
 */
int32_t
MPQC::CoordinateModel_impl::initialize_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.CoordinateModel.initialize)
 
  int i;
  have_guess_h_ = false;

  std::cout << "\nInitializing MPQC::ChemistryOpt_CoordinateModel\n";

  gov::cca::TypeMap tm = pp_.readConfigurationMap();

  grad_rms_ = double( tm.getDouble("grad_rms",0.0030) );
  grad_max_ = double( tm.getDouble("grad_max",0.0045) );
  disp_rms_ = double( tm.getDouble("disp_rms",0.0012) );
  disp_max_ = double( tm.getDouble("disp_max",0.0018) );
  diagonal_scale_ = double( tm.getDouble("diagonal_scale",1.0) );

  multiple_guess_h_ = bool( tm.getBool("multiple_guess_h",1) );
  use_current_geom_ = bool( tm.getBool("use_current_geom",0) );
  use_diagonal_h_ = bool( tm.getBool("use_diagonal_h",0) );

  coordinates_ = std::string( tm.getString("coordinate_type",
                                           "failed coordinate_type fetch") );
  extra_bonds_ = std::string( tm.getString("extra_bonds",
                                           "failed extra_bonds fetch") );

  CcaChemGeneric::CoordinateModel::set_tolerances(grad_rms_, grad_max_,
						  disp_rms_, disp_max_);
  CcaChemGeneric::CoordinateModel::initialize(services_);

  //get matrix kits
  kit_ = new sc::LocalSCMatrixKit; 
  rkit_ = new sc::ReplSCMatrixKit;

  //get coordinate type
  std::cout << "  Using coordinate type: " << coordinates_ << std::endl;
  if(coordinates_ == "cartesian") coorType_ = cart;
  else if(coordinates_ == "symmetrized") coorType_ = symm;
  else if(coordinates_ == "redundant") coorType_ = redund;
  else {
    std::cout << "  Unrecognized coordinate type, using default\n";
    coorType_ = DEFAULT_COORTYPE;
  }  

  //get model and molecule
  model_ = CcaChemGeneric::CoordinateModel::get_model();
  molecule_ = model_.get_molecule();
  double conv = molecule_.get_units().convert_to("bohr");
  convFrom_ = molecule_.get_units().convert_from("bohr");
  int natom = molecule_.get_n_atom();

  std::cout << "\n  CoordinateModel: setting up coordinates\n";

  //create input strings for MPQC classes
  std::ostringstream input;
  input
    << " molecule<Molecule>:(\n"
    << "   symmetry = auto\n"
    << "   unit = bohr\n"
    << "   {  n atoms geometry }={\n";
  for(i=0;i<natom;++i) {
    input
    << "\t" << i << "\t" << molecule_.get_atomic_number(i)
    << "\t[  " << molecule_.get_cart_coor(i,0)*conv 
    << "  " << molecule_.get_cart_coor(i,1)*conv
    << "  " << molecule_.get_cart_coor(i,2)*conv << "  ]\n";
  }
  input
    << "   }\n )\n";

  std::cout << input.str();

  switch(coorType_) {
  case cart:
    input << " coor<CartMolecularCoor>: (\n"
          << "   molecule = $..:molecule\n"
          << "   )\n";
    break;
  case symm:
    input << " coor<SymmMolecularCoor>: (\n"  
          << "   molecule = $..:molecule\n"     
          << "   update_bmat = 1\n"
	  << "   cartesian_tolerance = 1e-9\n"
          << "   generator<IntCoorGen>: (\n"     
          << "     molecule = $..:..:molecule\n";
    if( extra_bonds_.size() != 0 )
      input << "     extra_bonds = [" << extra_bonds_ << "]\n";
    input << "   )\n )\n";
    break;
  case redund:
    input << " coor<RedundMolecularCoor>: (\n" 
          << "   molecule = $..:molecule\n"
          << "   update_bmat = 1\n"
	  << "   cartesian_tolerance = 1e-9\n"
          << "   generator<IntCoorGen>: (\n"
          << "     molecule = $..:..:molecule\n";
    if( extra_bonds_.size() != 0 )
      input << "     extra_bonds = [" << extra_bonds_ << "]\n";
    input << "   )\n )\n";  
    break;
  }

  std::cout << input.str();

  //create the MPQC classes
  sc::Ref<sc::ParsedKeyVal> kv = new sc::ParsedKeyVal();
  kv->parse_string(input.str().c_str());
  sc::Ref<sc::DescribedClass> dccoor = kv->describedclassvalue("coor");
  sc::Ref<sc::DescribedClass> dcmol = kv->describedclassvalue("molecule");
  scMol_ = dynamic_cast<sc::Molecule*>(dcmol.pointer());
  switch(coorType_) {
  case cart:
    scCoor_ = dynamic_cast<sc::CartMolecularCoor*>(dccoor.pointer());
    break;
  case symm:
    scCoor_ = dynamic_cast<sc::SymmMolecularCoor*>(dccoor.pointer());
    break;
  case redund:
    scCoor_ = dynamic_cast<sc::RedundMolecularCoor*>(dccoor.pointer());
    break;
  }
 
  scCoor_->print();
  numCoor_ = scCoor_->dim().n();
  natom3_ = scCoor_->dim_natom3().n();
  std::cout << "\n";

  // convergence checking needs this method invoked
  CcaChemGeneric::CoordinateModel::get_n_coor();

  return 0;

  // DO-NOT-DELETE splicer.end(MPQC.CoordinateModel.initialize)
}

/**
 *  Releases and unregisters ports.  This should be called when the
 * CoordinateModel object is no longer needed.
 */
int32_t
MPQC::CoordinateModel_impl::finalize_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.CoordinateModel.finalize)

  return CcaChemGeneric::CoordinateModel::finalize();

  // DO-NOT-DELETE splicer.end(MPQC.CoordinateModel.finalize)
}

/**
 *  Sets the contained chemistry Model object (currently unused as the
 * chemistry Model object is normally obtained from a ModelFactory 
 * during initialization).
 * @param model The chemistry model object.
 */
void
MPQC::CoordinateModel_impl::set_model_impl (
  /* in */::Chemistry::QC::ModelInterface model ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.CoordinateModel.set_model)

  CcaChemGeneric::CoordinateModel::set_model( model );

  // DO-NOT-DELETE splicer.end(MPQC.CoordinateModel.set_model)
}

/**
 *  Returns the contained chemistry Model object.
 * @return The chemistry Model object.
 */
::Chemistry::QC::ModelInterface
MPQC::CoordinateModel_impl::get_model_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.CoordinateModel.get_model)

  return CcaChemGeneric::CoordinateModel::get_model();

  // DO-NOT-DELETE splicer.end(MPQC.CoordinateModel.get_model)
}

/**
 *  Returns the number of coordinates.
 * @return The number of coordinates. 
 */
int32_t
MPQC::CoordinateModel_impl::get_n_coor_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.CoordinateModel.get_n_coor)

  return numCoor_;

  // DO-NOT-DELETE splicer.end(MPQC.CoordinateModel.get_n_coor)
}

/**
 *  Returns the array of (cartesian or internal) coordinates which are 
 * being  optimized.
 * @return The array of coordinates which are being optimized.
 */
::sidl::array<double>
MPQC::CoordinateModel_impl::get_coor_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.CoordinateModel.get_coor)

  ::sidl::array<double> sidlCoor;

  switch(coorType_) {
  case cart:
    sidlCoor = molecule_.get_coor();
    break;
  default:
    sc::Ref<sc::SCMatrixKit> kit = new sc::LocalSCMatrixKit;
    sc::RefSCVector scCoor = kit->vector(scCoor_->dim());
    scCoor_->to_internal(scCoor);
    sidlCoor = vector_to_array(scCoor);
    break;
  } 

  return sidlCoor;

  // DO-NOT-DELETE splicer.end(MPQC.CoordinateModel.get_coor)
}

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
MPQC::CoordinateModel_impl::get_energy_impl (
  /* in array<double> */::sidl::array<double> x ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.CoordinateModel.get_energy)

  std::cout << "***** MPQC ChemistryOpt Calculate Energy *****\n";
  double f;

  //get energy, transform coordinates if neeeded
  switch(coorType_) {
  case cart:
    f = CcaChemGeneric::CoordinateModel::get_energy(x);
    break;
  default:
    sc::RefSCVector scCoor = kit_->vector(scCoor_->dim());
    sidl::array<double> cartx = 
      sidl::array<double>::create1d(natom3_);
    array_to_vector(x,scCoor);
    scCoor_->to_cartesian(scMol_,scCoor);
    for(int i=0; i<(natom3_/3); ++i) 
      for( int j=0; j<3; ++j) 
	cartx.set(i*3+j, scMol_->r(i,j)*convFrom_);
    f = CcaChemGeneric::CoordinateModel::get_energy(cartx);
    scCoor_->print();
    break;
  }

  return f; 

  // DO-NOT-DELETE splicer.end(MPQC.CoordinateModel.get_energy)
}

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
MPQC::CoordinateModel_impl::get_gradient_impl (
  /* in array<double> */::sidl::array<double> x ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.CoordinateModel.get_gradient)

  std::cout << "***** MPQC ChemistryOpt Calculate Gradient *****\n";

  ::sidl::array<double> cartg = ::sidl::array<double>::create1d(natom3_);

  //get gradient, transform coordinates if needed
  switch(coorType_) {
  case cart:
    cartg.copy(CcaChemGeneric::CoordinateModel::get_gradient(x));
    break;	
  default:
    sc::RefSCVector scCoor = kit_->vector(scCoor_->dim());
    array_to_vector(x,scCoor);
    scCoor_->to_cartesian(scMol_,scCoor);
    sidl::array<double> cartx = 
      sidl::array<double>::create1d(natom3_);
    for(int i=0; i<(natom3_/3); ++i) 
      for( int j=0; j<3; ++j) 
	cartx.set(i*3+j, scMol_->r(i,j)*convFrom_);
    cartg.copy(CcaChemGeneric::CoordinateModel::get_gradient(cartx));
    scCoor_->print();
    break;
  }

  //transform gradient if using internals
  ::sidl::array<double> g;
  switch(coorType_) {
  case cart:
    g = cartg;
    break;
  default:
    sc::RefSCVector scCartGrad = rkit_->vector(scCoor_->dim_natom3());
    sc::RefSCVector scGrad = rkit_->vector(scCoor_->dim());
    array_to_vector(cartg,scCartGrad);
    scCoor_->to_internal(scGrad,scCartGrad);
    g.copy(vector_to_array(scGrad));
    break;
  }

  return g;

  // DO-NOT-DELETE splicer.end(MPQC.CoordinateModel.get_gradient)
}

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
MPQC::CoordinateModel_impl::get_hessian_impl (
  /* in array<double> */::sidl::array<double> x ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.CoordinateModel.get_hessian)

  std::cout << "***** MPQC ChemistryOpt Calculate Hessian *****\n";
  std::cout << "  WARNING: this method has not been tested yet\n";

  ::sidl::array<double> cartH =
      ::sidl::array<double>::create2dRow(natom3_,natom3_);

  //get hessian, transform coordinates if needed
  switch(coorType_) {
  case cart:
    cartH.copy(CcaChemGeneric::CoordinateModel::get_hessian(x));
    break;
  default:
    sc::RefSCVector scCoor = kit_->vector(scCoor_->dim());
    array_to_vector(x,scCoor);
    scCoor_->to_cartesian(scMol_,scCoor);
    sidl::array<double> cartx = 
      sidl::array<double>::create1d(natom3_);
    for(int i=0; i<(natom3_/3); ++i) 
      for( int j=0; j<3; ++j) 
	cartx.set(i*3+j, scMol_->r(i,j)*convFrom_);
    cartH.copy(CcaChemGeneric::CoordinateModel::get_hessian(cartx));
    scCoor_->print();
    break;
  }

  //transform Hessian if using internals
  ::sidl::array<double> H;
  switch(coorType_) {
  case cart:
    H = cartH;
    break;
  default:
    sc::RefSymmSCMatrix scCartH = 
      rkit_->symmmatrix(scCoor_->dim_natom3());
    sc::RefSymmSCMatrix scH = 
      rkit_->symmmatrix(scCoor_->dim());
    array_to_matrix(cartH,scCartH);
    scCoor_->to_internal(scH,scCartH);
    H = matrix_to_array(scH);
    break;
  }

  return H;

  // DO-NOT-DELETE splicer.end(MPQC.CoordinateModel.get_hessian)
}

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
MPQC::CoordinateModel_impl::get_energy_and_gradient_impl (
  /* in array<double> */::sidl::array<double> x,
  /* out */double& f,
  /* in array<double> */::sidl::array<double> g ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.CoordinateModel.get_energy_and_gradient)

  std::cout << "***** MPQC ChemistryOpt Calculate Energy and Gradient *****\n";

  ::sidl::array<double> cartg = ::sidl::array<double>::create1d(natom3_);

  //get gradient, transform coordinates if needed
  switch(coorType_) {
  case cart:
    CcaChemGeneric::CoordinateModel::get_energy_and_gradient(x,&f,g);
    break;
  default:
    sc::RefSCVector scCoor = kit_->vector(scCoor_->dim());
    array_to_vector(x,scCoor);
    scCoor_->to_cartesian(scMol_,scCoor);
    sidl::array<double> cartx =
      sidl::array<double>::create1d(natom3_);
    for(int i=0; i<(natom3_/3); ++i)
      for( int j=0; j<3; ++j)
	cartx.set(i*3+j, scMol_->r(i,j)*convFrom_);
    CcaChemGeneric::CoordinateModel::get_energy_and_gradient(cartx,&f,cartg);

    scCoor_->print();
    
    break;
  }

  //transform gradient if using internals
  switch(coorType_) {
  case cart:
    break;
  default:
    sc::RefSCVector scCartGrad = rkit_->vector(scCoor_->dim_natom3());
    sc::RefSCVector scGrad = rkit_->vector(scCoor_->dim());
    array_to_vector(cartg,scCartGrad);
    scCoor_->to_internal(scGrad,scCartGrad);
    g.copy(vector_to_array(scGrad));
    break;
  }

  // DO-NOT-DELETE splicer.end(MPQC.CoordinateModel.get_energy_and_gradient)
}

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
MPQC::CoordinateModel_impl::guess_hessian_solve_impl (
  /* in array<double> */::sidl::array<double> effective_grad,
  /* in array<double> */::sidl::array<double> effective_step,
  /* in */void* first_geom ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.CoordinateModel.guess_hessian_solve)

  sidl::array<double> *sidl_geom_ptr = 
    static_cast< sidl::array<double>* >( first_geom );
  
  if( !use_diagonal_h_ ) {
    if(multiple_guess_h_ && !use_current_geom_ && sidl_geom_ptr){
      std::cout << "Using geometry for first correction pair\n";
      sc::RefSCVector scCoor = kit_->vector(scCoor_->dim());
      array_to_vector( *sidl_geom_ptr, scCoor );
      scCoor_->to_cartesian(scMol_,scCoor);    
    }      

    if(  multiple_guess_h_ || !have_guess_h_  ) {
      std::cout << "Determining approximate Hessian\n";
      sc::RefSymmSCMatrix hess = rkit_->symmmatrix(scCoor_->dim());
      scCoor_->guess_hessian(hess);
      ihess_ = scCoor_->inverse_hessian(hess);
      have_guess_h_ = 1;
    }
  }
  else {
    std::cout << "Using diagonal Hessian scaled by " << diagonal_scale_ 
              << std::endl;
    ihess_ = rkit_->symmmatrix(scCoor_->dim());
    ihess_->unit();
    ihess_->scale(1.0/diagonal_scale_);
  }

  std::cout << "Solving approximate Hessian system\n";
  sc::RefSCVector scV = rkit_->vector(scCoor_->dim());
  array_to_vector(effective_grad, scV);
  sc::RefSCVector result = ihess_*scV;
  effective_step.copy(vector_to_array( result ));

  // DO-NOT-DELETE splicer.end(MPQC.CoordinateModel.guess_hessian_solve)
}

/**
 *  Determines if the optimization has converged, flag is set to 1
 * if convergence has been achieved and 0 otherwise.
 * @param flag Variable that convergence value is assigned to.
 */
void
MPQC::CoordinateModel_impl::checkConvergence_impl (
  /* inout */int32_t& flag ) 
{
  // DO-NOT-DELETE splicer.begin(MPQC.CoordinateModel.checkConvergence)

  CcaChemGeneric::CoordinateModel::checkConvergence( flag );

  // DO-NOT-DELETE splicer.end(MPQC.CoordinateModel.checkConvergence)
}

/**
 *  For visualization, possibly unused (?).  CoordinateModel objects
 * may callback to viewers that implement the Chemistry.MoleculeViewer 
 * interface, such as the cca-chem python GUI, making this method 
 * unnecessary.
 */
void
MPQC::CoordinateModel_impl::monitor_impl () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.CoordinateModel.monitor)
  // Insert-Code-Here {MPQC.CoordinateModel.monitor} (monitor method)
  // 
  // This method has not been implemented
  // 
  // DO-NOT-DELETE splicer.end(MPQC.CoordinateModel.monitor)
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
MPQC::CoordinateModel_impl::setServices_impl (
  /* in */::gov::cca::Services services ) 
// throws:
//     ::gov::cca::CCAException
//     ::sidl::RuntimeException
{
  // DO-NOT-DELETE splicer.begin(MPQC.CoordinateModel.setServices)

  services_ = services;
  if (services_._is_nil()) return;

  try {
      services_.addProvidesPort( *this,
                                "CoordinateModelInterface",
                                "ChemistryOpt.CoordinateModelInterface",
                                0);
      services_.registerUsesPort("ModelFactoryInterface",
                                 "Chemistry.QC.ModelFactoryInterface",
                                 0);
      services_.registerUsesPort("BackupModelFactoryInterface",
                                 "Chemistry.QC.ModelFactoryInterface",
                                 0);
      services_.registerUsesPort("MoleculeViewerInterface",
                                 "Chemistry.MoleculeViewerInterface",
                                 0);
      services_.registerUsesPort("ppf",
                                 "gov.cca.ports.ParameterPortFactory", 0);
  }
      catch (gov::cca::CCAException e) {
      std::cout << "Error using services: "
                << e.getNote() << std::endl;
  }

  try {

    tm_ = services_.createTypeMap();
    if(tm_._is_nil()) {
      std::cerr << "TypeMap is nill\n";
      abort();
    }
    ppf_ = sidl::babel_cast<gov::cca::ports::ParameterPortFactory>(
             services_.getPort("ppf") );
    ppf_.initParameterData(tm_, "CONFIG");
    ppf_.setBatchTitle(tm_,"MPQC CoordinateModel");
    ppf_.setGroupName(tm_,"options");
    ppf_.addRequestString(tm_, 
           "coordinate_type", 
           "coordinate type: cartesian, symmetrized, or redundant",
           "coordinate_type", "symmetrized");
    ppf_.addRequestDouble(tm_, 
           "grad_rms",
           "RMS gradient convergence tolerance",
           "grad_rms", 0.00030,0,1000000);
    ppf_.addRequestDouble(tm_,
           "grad_max",
           "Max gradient convergence tolerance",
           "grad_max", 0.00045,0,1000000);
    ppf_.addRequestDouble(tm_,
           "disp_rms",
           "RMS displacement convergence tolerance",
           "disp_rms", 0.00120,0,1000000);
    ppf_.addRequestDouble(tm_,
           "disp_max",
           "Max displacement convergence tolerance",
           "disp_max", 0.00180,0,1000000);
    ppf_.addRequestDouble(tm_,
           "diagonal_scale",
           "Diagonal guess Hessian scale factor",
           "diagonal_scale",1.0,0,1000000);
    ppf_.addRequestBoolean(tm_, 
           "multiple_guess_h",
           "Guess H at every guess_hessian_solve",
           "multiple_guess_h",1);
    ppf_.addRequestBoolean(tm_,
           "use_current_geom",
           "Guess Hessian at current geometry",
           "use_current_geom",0);
    ppf_.addRequestBoolean(tm_,
           "use_diagonal_h",
           "Use diagonal guess Hessian",
           "use_diagonal_h",0);
    ppf_.addRequestString(tm_,
           "coordinate_type",
           "coordinate type: cartesian, symmetrized, or redundant",
           "coordinate_type", "symmetrized");
    ppf_.addRequestString(tm_,
           "extra_bonds", 
           "extra_bonds vector",
           "extra_bonds", "");
    ppf_.addParameterPort(tm_, services_);
    services_.releasePort("ppf");

    pp_ = sidl::babel_cast<gov::cca::ports::ParameterPort>(
            services_.getPort("CONFIG") );
    if (pp_._is_nil()) {
      std::cerr << "getport failed\n";
      abort();
    }

  }
  catch(std::exception& e) {
    std::cerr << "Error in parameter port setup: " << e.what() << std::endl;
  }

  // DO-NOT-DELETE splicer.end(MPQC.CoordinateModel.setServices)
}


// DO-NOT-DELETE splicer.begin(MPQC.CoordinateModel._misc)
// Insert-Code-Here {MPQC.CoordinateModel._misc} (miscellaneous code)
// DO-NOT-DELETE splicer.end(MPQC.CoordinateModel._misc)

