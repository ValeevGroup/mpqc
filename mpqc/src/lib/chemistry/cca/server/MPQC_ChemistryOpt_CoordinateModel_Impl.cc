// 
// File:          MPQC_ChemistryOpt_CoordinateModel_Impl.cc
// Symbol:        MPQC.ChemistryOpt_CoordinateModel-v0.2
// Symbol Type:   class
// Babel Version: 0.10.12
// Description:   Server-side implementation for MPQC.ChemistryOpt_CoordinateModel
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// xml-url       = /home/jpkenny/src/mpqc-libint2.build-shared/src/lib/chemistry/cca/server/../../../../../lib/cca/repo/MPQC.ChemistryOpt_CoordinateModel-v0.2.xml
// 
#include "MPQC_ChemistryOpt_CoordinateModel_Impl.hh"

// DO-NOT-DELETE splicer.begin(MPQC.ChemistryOpt_CoordinateModel._includes)
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

// DO-NOT-DELETE splicer.end(MPQC.ChemistryOpt_CoordinateModel._includes)

// user-defined constructor.
void MPQC::ChemistryOpt_CoordinateModel_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(MPQC.ChemistryOpt_CoordinateModel._ctor)
  have_guess_h_ = 0;
  // DO-NOT-DELETE splicer.end(MPQC.ChemistryOpt_CoordinateModel._ctor)
}

// user-defined destructor.
void MPQC::ChemistryOpt_CoordinateModel_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(MPQC.ChemistryOpt_CoordinateModel._dtor)
  // add destruction details here
  // DO-NOT-DELETE splicer.end(MPQC.ChemistryOpt_CoordinateModel._dtor)
}

// static class initializer.
void MPQC::ChemistryOpt_CoordinateModel_impl::_load() {
  // DO-NOT-DELETE splicer.begin(MPQC.ChemistryOpt_CoordinateModel._load)
  // guaranteed to be called at most once before any other method in this class
  // DO-NOT-DELETE splicer.end(MPQC.ChemistryOpt_CoordinateModel._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
/**
 * Registers and gets ports, and requests Model object(s) from the 
 * ModelFactory component(s). This must be the first method called 
 * following instantiation.
 */
int32_t
MPQC::ChemistryOpt_CoordinateModel_impl::initialize ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.ChemistryOpt_CoordinateModel.initialize)
  
  int i;

  std::cout << "\nInitializing MPQC::ChemistryOpt_CoordinateModel\n";

  CcaChemGeneric::CoordinateModel::set_tolerances(grad_rms_->value,
						  grad_max_->value,
						  disp_rms_->value,
						  disp_max_->value);
  CcaChemGeneric::CoordinateModel::initialize(services_);

  //get matrix kits
  kit_ = new sc::LocalSCMatrixKit; 
  rkit_ = new sc::ReplSCMatrixKit;

  //get coordinate type
  std::string coorString;
  coorString = std::string( coordinates_->getValueString() );
  std::cout << "  Using coordinate type: " << coorString << std::endl;
  if(coorString == "cartesian") coorType_ = cart;
  else if(coorString == "symmetrized") coorType_ = symm;
  else if(coorString == "redundant") coorType_ = redund;
  else {
    std::cout << "  Unrecognized coordinate type, using default\n";
    coorType_ = DEFAULT_COORTYPE;
  }  
  services_.releasePort("CoordinateType");

  //get extra_bonds
  std::string bondsString( extra_bonds_->getValueString() );

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
    if( bondsString.size() != 0 )
      input << "     extra_bonds = [" << bondsString << "]\n";
    input << "   )\n )\n";
    break;
  case redund:
    input << " coor<RedundMolecularCoor>: (\n" 
          << "   molecule = $..:molecule\n"
          << "   update_bmat = 1\n"
	  << "   cartesian_tolerance = 1e-9\n"
          << "   generator<IntCoorGen>: (\n"
          << "     molecule = $..:..:molecule\n";
    if( bondsString.size() != 0 )
      input << "     extra_bonds = [" << bondsString << "]\n";
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
  // DO-NOT-DELETE splicer.end(MPQC.ChemistryOpt_CoordinateModel.initialize)
}

/**
 * Releases and unregisters ports.  This should be called when the
 * CoordinateModel object is no longer needed.
 */
int32_t
MPQC::ChemistryOpt_CoordinateModel_impl::finalize ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.ChemistryOpt_CoordinateModel.finalize)
  return CcaChemGeneric::CoordinateModel::finalize();
  // DO-NOT-DELETE splicer.end(MPQC.ChemistryOpt_CoordinateModel.finalize)
}

/**
 * Sets the contained chemistry Model object (currently unused as the
 * chemistry Model object is normally obtained from a ModelFactory 
 * during initialization).
 * @param model The chemistry model object.
 */
void
MPQC::ChemistryOpt_CoordinateModel_impl::set_model (
  /* in */ ::Chemistry::QC::Model model ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.ChemistryOpt_CoordinateModel.set_model)
  CcaChemGeneric::CoordinateModel::set_model( model );
  // DO-NOT-DELETE splicer.end(MPQC.ChemistryOpt_CoordinateModel.set_model)
}

/**
 * Returns the contained chemistry Model object.
 * @return The chemistry Model object.
 */
::Chemistry::QC::Model
MPQC::ChemistryOpt_CoordinateModel_impl::get_model ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.ChemistryOpt_CoordinateModel.get_model)
  return CcaChemGeneric::CoordinateModel::get_model();
  // DO-NOT-DELETE splicer.end(MPQC.ChemistryOpt_CoordinateModel.get_model)
}

/**
 * Returns the number of coordinates.
 * @return The number of coordinates. 
 */
int32_t
MPQC::ChemistryOpt_CoordinateModel_impl::get_n_coor ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.ChemistryOpt_CoordinateModel.get_n_coor)
  return numCoor_;
  // DO-NOT-DELETE splicer.end(MPQC.ChemistryOpt_CoordinateModel.get_n_coor)
}

/**
 * Returns the array of (cartesian or internal) coordinates which are 
 * being  optimized.
 * @return The array of coordinates which are being optimized.
 */
::sidl::array<double>
MPQC::ChemistryOpt_CoordinateModel_impl::get_coor ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.ChemistryOpt_CoordinateModel.get_coor)

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
  // DO-NOT-DELETE splicer.end(MPQC.ChemistryOpt_CoordinateModel.get_coor)
}

/**
 * Returns the energy of the currently contained model with the values
 * of the optimization coordinates given in x.  This requires
 * that the CoordinateModel updates the cartesian coordinates of a 
 * contained Molecule object (possibly requiring transformation) and set 
 * this Molecule object on a contained Model object, prior to calling
 * get_energy() on the Model object.
 * @param x The optimization coordinate values.
 * @return The energy of the chemistry model at x.
 */
double
MPQC::ChemistryOpt_CoordinateModel_impl::get_energy (
  /* in */ ::sidl::array<double> x ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.ChemistryOpt_CoordinateModel.get_energy)

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

  // DO-NOT-DELETE splicer.end(MPQC.ChemistryOpt_CoordinateModel.get_energy)
}

/**
 * Returns the energy gradient of the currently contained model with 
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
MPQC::ChemistryOpt_CoordinateModel_impl::get_gradient (
  /* in */ ::sidl::array<double> x ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.ChemistryOpt_CoordinateModel.get_gradient)

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

  // DO-NOT-DELETE splicer.end(MPQC.ChemistryOpt_CoordinateModel.get_gradient)
}

/**
 * Returns the energy Hessian of the currently contained model with
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
MPQC::ChemistryOpt_CoordinateModel_impl::get_hessian (
  /* in */ ::sidl::array<double> x ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.ChemistryOpt_CoordinateModel.get_hessian)

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

  // DO-NOT-DELETE splicer.end(MPQC.ChemistryOpt_CoordinateModel.get_hessian)
}

/**
 * Sets f and g to the energy and energy gradient, respectively,
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
MPQC::ChemistryOpt_CoordinateModel_impl::get_energy_and_gradient (
  /* in */ ::sidl::array<double> x,
  /* out */ double& f,
  /* in */ ::sidl::array<double> g ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.ChemistryOpt_CoordinateModel.get_energy_and_gradient)
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

  // DO-NOT-DELETE splicer.end(MPQC.ChemistryOpt_CoordinateModel.get_energy_and_gradient)
}

/**
 * Returns the product of the guess hessian inverse and an effective
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
MPQC::ChemistryOpt_CoordinateModel_impl::guess_hessian_solve (
  /* in */ ::sidl::array<double> effective_grad,
  /* in */ ::sidl::array<double> effective_step,
  /* in */ void* first_geom ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.ChemistryOpt_CoordinateModel.guess_hessian_solve)
  
  sidl::array<double> *sidl_geom_ptr = 
    static_cast< sidl::array<double>* >( first_geom );
  
  if(multiple_guess_h_->value && !use_current_geom_->value && sidl_geom_ptr){
    std::cout << "Using geometry for first correction pair\n";
    sc::RefSCVector scCoor = kit_->vector(scCoor_->dim());
    array_to_vector( *sidl_geom_ptr, scCoor );
    scCoor_->to_cartesian(scMol_,scCoor);    
  }      
  
  if(  multiple_guess_h_->value || !have_guess_h_  ) {
    std::cout << "Determining approximate Hessian\n";
    sc::RefSymmSCMatrix hess = rkit_->symmmatrix(scCoor_->dim());
    scCoor_->guess_hessian(hess);
    ihess_ = scCoor_->inverse_hessian(hess);
    have_guess_h_ = 1;
  }
  
  std::cout << "Solving approximate Hessian system\n";
  sc::RefSCVector scV = rkit_->vector(scCoor_->dim());
  array_to_vector(effective_grad, scV);
  sc::RefSCVector result = ihess_*scV;
  effective_step.copy(vector_to_array( result ));
  
  // DO-NOT-DELETE splicer.end(MPQC.ChemistryOpt_CoordinateModel.guess_hessian_solve)
}

/**
 * Determines if the optimization has converged, flag is set to 1
 * if convergence has been achieved and 0 otherwise.
 * @param flag Variable that convergence value is assigned to.
 */
void
MPQC::ChemistryOpt_CoordinateModel_impl::checkConvergence (
  /* inout */ int32_t& flag ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(MPQC.ChemistryOpt_CoordinateModel.checkConvergence)
  CcaChemGeneric::CoordinateModel::checkConvergence( flag );
  // DO-NOT-DELETE splicer.end(MPQC.ChemistryOpt_CoordinateModel.checkConvergence)
}

/**
 * For visualization, possibly unused (?).  CoordinateModel objects
 * may callback to viewers that implement the Chemistry.MoleculeViewer 
 * interface, such as the cca-chem python GUI, making this method 
 * unnecessary.
 */
void
MPQC::ChemistryOpt_CoordinateModel_impl::monitor ()
throw () 

{
  // DO-NOT-DELETE splicer.begin(MPQC.ChemistryOpt_CoordinateModel.monitor)
  // insert implementation here
  // DO-NOT-DELETE splicer.end(MPQC.ChemistryOpt_CoordinateModel.monitor)
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
MPQC::ChemistryOpt_CoordinateModel_impl::setServices (
  /* in */ ::gov::cca::Services services ) 
throw ( 
  ::gov::cca::CCAException
){
  // DO-NOT-DELETE splicer.begin(MPQC.ChemistryOpt_CoordinateModel.setServices)

  services_ = services;
  if (services_._is_nil()) return;

  try {
      services_.addProvidesPort(self,
                                "CoordinateModel",
                                "ChemistryOpt.CoordinateModel",
                                0);
      services_.registerUsesPort("ModelFactory",
                                 "Chemistry.QC.ModelFactory",
                                 0);
      services_.registerUsesPort("BackupModelFactory",
                                 "Chemistry.QC.ModelFactory",
                                 0);
      services_.registerUsesPort("MoleculeViewer",
                                 "Chemistry.MoleculeViewer",
                                 0);
      services_.registerUsesPort("CoordinateType",
				 "Util.StringProvider",
				 0);
  }
      catch (gov::cca::CCAException e) {
      std::cout << "Error using services: "
                << e.getNote() << std::endl;
  }

    // setup parameters
  try {
    
    if (services_._not_nil()) {
      gov::cca::TypeMap tm = services_.createTypeMap();
      ::gov::cca::Port self_port = self;
      services_.addProvidesPort(self_port,
				"string",
				"Util.StringProvider",tm);
      
      services_.registerUsesPort("classicParam",
				 "gov.cca.ParameterPortFactoryService",tm);
      gov::cca::Port p = services_.getPort("classicParam");
      ccaffeine::ports::PortTranslator portX = p;
      if(portX._not_nil()) {
	classic::gov::cca::Port *cp
	  =static_cast<classic::gov::cca::Port*>(portX.getClassicPort());
	if(!cp) {
	  std::cout << "Couldn't get classic port" << std::endl;
	  return;
	}
	ConfigurableParameterFactory *cpf
	  = dynamic_cast<ConfigurableParameterFactory *>(cp);
	ConfigurableParameterPort *pp = setup_parameters(cpf);
	classic::gov::cca::Port *clscp
	  = dynamic_cast<classic::gov::cca::Port*>(pp);
	if (!clscp) {
	  std::cout << "Couldn't cast to classic::gov::cca::Port"
		    << std::endl;
	}
	void *vp = static_cast<void*>(clscp);
	ccaffeine::ports::PortTranslator provideX
	  = ccaffeine::ports::PortTranslator::createFromClassic(vp);
	
	services_.addProvidesPort(provideX,
				  "configure", "ParameterPort", tm);
	
	services_.releasePort("classicParam");
	services_.unregisterUsesPort("classicParam");
      }
    }
    
  }
  catch(std::exception& e) {
    std::cout << "Exception caught: " << e.what() << std::endl;
  }

  // DO-NOT-DELETE splicer.end(MPQC.ChemistryOpt_CoordinateModel.setServices)
}


// DO-NOT-DELETE splicer.begin(MPQC.ChemistryOpt_CoordinateModel._misc)

ConfigurableParameterPort*
MPQC::ChemistryOpt_CoordinateModel_impl::setup_parameters(ConfigurableParameterFactory *cpf)
{
  ConfigurableParameterPort * pp = cpf->createConfigurableParameterPort();

  pp->setBatchTitle("PortTranslatorStarter Configuration");
  pp->setGroupName("Coordinate Model Input");
  grad_rms_ = new DoubleParameter("grad_rms",
				  "RMS gradient convergence tolerance",
				  "grad_rms", 0.00030,0,1000000);
  grad_max_ = new DoubleParameter("grad_max",
				  "Max gradient convergence tolerance",
				  "grad_max", 0.00045,0,1000000);
  disp_rms_ = new DoubleParameter("disp_rms",
				  "RMS displacement convergence tolerance",
				  "disp_rms", 0.00120,0,1000000);
  disp_max_ = new DoubleParameter("disp_max",
				  "Max displacement convergence tolerance",
				  "disp_max", 0.00180,0,1000000);
  multiple_guess_h_ = new BoolParameter("multiple_guess_h",
                                        "Guess H at every guess_hessian_solve",
					"multiple_guess_h",1);
  use_current_geom_ = new BoolParameter("use_current_geom",
					"Guess Hessian at current geometry",
					"use_current_geom",0);
  coordinates_ = new StringParameter("coordinate_type",
                       "Coordinate type: cartesian, symmetrized, or redundant",
                                     "coordinate_type", "symmetrized");
  extra_bonds_ = new StringParameter("extra_bonds", "extra_bonds vector",
				     "extra_bonds", "");
  pp->addRequest(grad_rms_);
  pp->addRequest(grad_max_);
  pp->addRequest(disp_rms_);  
  pp->addRequest(disp_max_);
  pp->addRequest(multiple_guess_h_);
  pp->addRequest(use_current_geom_);
  pp->addRequest(coordinates_);
  pp->addRequest(extra_bonds_);

  return pp;
}

// DO-NOT-DELETE splicer.end(MPQC.ChemistryOpt_CoordinateModel._misc)

