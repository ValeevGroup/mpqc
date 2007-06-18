//
// mbptr12.cc
//
// Copyright (C) 2001 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
// Maintainer: EV
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifdef __GNUC__
#pragma implementation
#endif

#include <stdexcept>
#include <sstream>
#include <util/misc/string.h>
#include <util/class/scexception.h>
#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <util/state/stateio.h>
#include <math/scmat/blocked.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/scf/clhf.h>
#include <chemistry/qc/scf/uhf.h>
#if HAVE_INTEGRALCINTS
#  include <chemistry/qc/cints/cints.h>
#endif
#if HAVE_INTEGRALLIBINT2
#  include <chemistry/qc/libint2/libint2.h>
#endif
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/transform_factory.h>
#include <chemistry/qc/mbptr12/gaussianfit.h>
#include <chemistry/qc/mbptr12/gaussianfit.timpl.h>
#include <chemistry/qc/mbptr12/linearr12.timpl.h>
#include <chemistry/qc/mbptr12/print.h>

// set to 1 to include the cusp region in the fit, other use Tew&Klopper's strategy
#define INCLUDE_CUSP_IN_GTG_FIT 0

using namespace std;
using namespace sc;
using namespace sc::LinearR12;

/*--------------------------------
  MBPT2_R12
 --------------------------------*/

static ClassDesc MBPT2_R12_cd(
  typeid(MBPT2_R12),"MBPT2_R12",8,"public MBPT2",
  0, create<MBPT2_R12>, create<MBPT2_R12>);

MBPT2_R12::MBPT2_R12(StateIn& s):
  MBPT2(s)
{
  r12eval_ << SavableState::restore_state(s);
  r12a_energy_ << SavableState::restore_state(s);
  r12ap_energy_ << SavableState::restore_state(s);
  r12app_energy_ << SavableState::restore_state(s);
  r12b_energy_ << SavableState::restore_state(s);
  r12c_energy_ << SavableState::restore_state(s);
  aux_basis_ << SavableState::restore_state(s);
  vir_basis_ << SavableState::restore_state(s);

  gbc_ = true;  ebc_ = true;
  if (s.version(::class_desc<MBPT2_R12>()) >= 3) {
    int gbc; s.get(gbc); gbc_ = (bool)gbc;
    int ebc; s.get(ebc); ebc_ = (bool)ebc;
  }
  if (s.version(::class_desc<MBPT2_R12>()) >= 6) {
    int ks_ebcfree; s.get(ks_ebcfree); ks_ebcfree_ = (bool)ks_ebcfree;
    int omit_P; s.get(omit_P); omit_P_ = (bool)omit_P;
  }
  
  abs_method_ = LinearR12::ABS_ABS;
  if (s.version(::class_desc<MBPT2_R12>()) >= 2) {
    int absmethod; s.get(absmethod); abs_method_ = (LinearR12::ABSMethod)absmethod;
  }
  int stdapprox; s.get(stdapprox); stdapprox_ = (LinearR12::StandardApproximation) stdapprox;
  if (s.version(::class_desc<MBPT2_R12>()) >= 7) {
    ansatz_ << SavableState::restore_state(s);
  }

  maxnabs_ = 2;
  if (s.version(::class_desc<MBPT2_R12>()) >= 5) {
    s.get(maxnabs_);
  }

  int spinadapted; s.get(spinadapted); spinadapted_ = (bool)spinadapted;

  include_mp1_ = false;
  if (s.version(::class_desc<MBPT2_R12>()) >= 4) {
    int include_mp1; s.get(include_mp1); include_mp1_ = static_cast<bool>(include_mp1);
  }

  int r12ints_method; s.get(r12ints_method);
    r12ints_method_ = static_cast<R12IntEvalInfo::StoreMethod::type>(r12ints_method);
  s.get(r12ints_file_);

  // didn't do the checks in the old code
  safety_check_ = false;
  if (s.version(::class_desc<MBPT2_R12>()) >= 8) {
      int safety_check; s.get(safety_check);
      safety_check_ = static_cast<bool>(safety_check);
  }
  // didn't ensure positive definite B in the old code
  posdef_B_ = LinearR12::PositiveDefiniteB_no;
  if (s.version(::class_desc<MBPT2_R12>()) >= 8) {
      int posdef_B; s.get(posdef_B);
      posdef_B_ = static_cast<LinearR12::PositiveDefiniteB>(posdef_B);
  }

  s.get(mp2_corr_energy_);
  s.get(r12_corr_energy_);

  twopdm_grid_ << SavableState::restore_state(s);
  s.get(plot_pair_function_[0]);
  s.get(plot_pair_function_[1]);
}

MBPT2_R12::MBPT2_R12(const Ref<KeyVal>& keyval):
  MBPT2(keyval)
{
  // Make sure can use the integral factory for linear R12
  check_integral_factory_();
  
  // Verify that this is a closed-shell or high-spin open-shell system
  CLHF* clhfref = dynamic_cast<CLHF*>(reference_.pointer());
  UHF* uhfref = dynamic_cast<UHF*>(reference_.pointer());
  if (clhfref == 0 && uhfref == 0) {
    throw FeatureNotImplemented("MBPT2_R12::MBPT2_R12: reference wavefunction is neither CLHF nor UHF",
    __FILE__,__LINE__);
  }
  const bool closedshell = (clhfref != 0);

  aux_basis_ = require_dynamic_cast<GaussianBasisSet*>(
    keyval->describedclassvalue("aux_basis").pointer(),
    "MBPT2_R12::MBPT2_R12\n"
    );
  if (aux_basis_.pointer() == NULL)
    aux_basis_ = basis();

  const bool abs_eq_obs = aux_basis_->equiv(basis());

  vir_basis_ = require_dynamic_cast<GaussianBasisSet*>(
    keyval->describedclassvalue("vir_basis").pointer(),
    "MBPT2_R12::MBPT2_R12\n"
    );
  if (vir_basis_.pointer() == NULL)
    vir_basis_ = basis();

  const bool vbs_eq_obs = vir_basis_->equiv(basis());

  // Default is to use R12 factor
  std::string corrfactor = keyval->stringvalue("corr_factor", KeyValValuestring("r12"));

  // Default method is MBPT2-R12/A'
  char *sa_string = keyval->pcharvalue("stdapprox",KeyValValuepchar("A'"));
  if ( !strcmp(sa_string,"A") ||
       !strcmp(sa_string,"a") ) {
    stdapprox_ = LinearR12::StdApprox_A;
  }
  else if ( !strcmp(sa_string,"Ap") ||
	    !strcmp(sa_string,"ap") ||
	    !strcmp(sa_string,"A'") ||
	    !strcmp(sa_string,"a'") ) {
    stdapprox_ = LinearR12::StdApprox_Ap;
  }
  else if ( !strcmp(sa_string,"App") ||
	    !strcmp(sa_string,"app") ||
	    !strcmp(sa_string,"A''") ||
	    !strcmp(sa_string,"a''") ) {
    stdapprox_ = LinearR12::StdApprox_App;
  }
  else if ( !strcmp(sa_string,"B") ||
	    !strcmp(sa_string,"b") ) {
    stdapprox_ = LinearR12::StdApprox_B;
  }
  else if ( !strcmp(sa_string,"C") ||
	    !strcmp(sa_string,"c") ) {
    stdapprox_ = LinearR12::StdApprox_C;
  }
  else {
    delete[] sa_string;
    throw std::runtime_error("MBPT2_R12::MBPT2_R12() -- unrecognized value for stdapprox");
  }

  // if no explicit correlation then set stdapprox to A
  if (sa_string == "none" || sa_string == "NONE") {
      stdapprox_ = LinearR12::StdApprox_A;
  }
  
  //
  // r12 correlation factor?
  //
  if (corrfactor == "r12" ||
      corrfactor == "R12") {
    corrfactor_ = new LinearR12::R12CorrelationFactor();
  }
  //
  // g12 correlation factor?
  //
  else if (corrfactor == "g12" || corrfactor == "G12") {
    if (keyval->exists("corr_param")) {
      typedef LinearR12::G12CorrelationFactor::CorrelationParameters CorrParams;
      CorrParams params;
      const int num_f12 = keyval->count("corr_param");
      if (num_f12 != 0) {
        // Do I have contracted functions?
        bool contracted = (keyval->count("corr_param",0) != 0);
        if (!contracted) {
          // Primitive functions only
          for(int f=0; f<num_f12; f++) {
            double exponent = keyval->doublevalue("corr_param", f);
            LinearR12::G12CorrelationFactor::ContractedGeminal vtmp;
            vtmp.push_back(std::make_pair(exponent,1.0));
            params.push_back(vtmp);
          }
        }
        else {
          // Contracted functions
          for(int f=0; f<num_f12; f++) {
            const int nprims = keyval->count("corr_param", f);
            if (nprims == 0)
              throw InputError("Contracted and primitive geminals cannot be mixed in the input", __FILE__, __LINE__);
            LinearR12::G12CorrelationFactor::ContractedGeminal vtmp;
            for(int p=0; p<nprims; p++) {
              if (keyval->count("corr_param", f, p) != 2)
                throw InputError("Invalid contracted geminal specification",__FILE__,__LINE__);
              double exponent = keyval->Va_doublevalue("corr_param", 3, f, p, 0);
              double coef = keyval->Va_doublevalue("corr_param", 3, f, p, 1);
              vtmp.push_back(std::make_pair(exponent,coef));
            }
            params.push_back(vtmp);
          }
        }
      }
      else {
        double exponent = keyval->doublevalue("corr_param");
        std::vector< std::pair<double,double> > vtmp;  vtmp.push_back(std::make_pair(exponent,1.0));
        params.push_back(vtmp);
      }
      // If stdapprox_ == C, no need for commutators
      if (stdapprox_ == LinearR12::StdApprox_C)
	  corrfactor_ = new LinearR12::G12NCCorrelationFactor(params);
      else
	  corrfactor_ = new LinearR12::G12CorrelationFactor(params);
    }
    else
      throw ProgrammingError("MBPT2_R12::MBPT2_R12() -- corr_param keyword must be given when corr_factor=g12",__FILE__,__LINE__);
  }
  //
  // geng12 correlation factor?
  //
  else if (corrfactor == "geng12" || corrfactor == "GENG12") {
    if (keyval->exists("corr_param")) {
      typedef LinearR12::GenG12CorrelationFactor::CorrelationParameters CorrParams;
      CorrParams params;
      const int num_f12 = keyval->count("corr_param");
      if (num_f12 != 0) {
        // Do I have contracted functions?
        bool contracted = (keyval->count("corr_param",0) != 0);
        if (!contracted) {
          // Primitive functions only
          for(int f=0; f<num_f12; f++) {
            double exponent = keyval->doublevalue("corr_param", f);
            LinearR12::GenG12CorrelationFactor::ContractedGeminal vtmp;
            vtmp.push_back(std::make_pair(std::make_pair(0.0,exponent),1.0));
            params.push_back(vtmp);
          }
        }
        else {
          // Contracted functions
          for(int f=0; f<num_f12; f++) {
            const int nprims = keyval->count("corr_param", f);
            if (nprims == 0)
              throw InputError("Contracted and primitive geminals cannot be mixed in the input", __FILE__, __LINE__);
            LinearR12::GenG12CorrelationFactor::ContractedGeminal vtmp;
            for(int p=0; p<nprims; p++) {
              if (keyval->count("corr_param", f, p) != 3)
                throw InputError("Invalid contracted geminal specification",__FILE__,__LINE__);
              double alpha = keyval->Va_doublevalue("corr_param", 3, f, p, 0);
              double gamma = keyval->Va_doublevalue("corr_param", 3, f, p, 1);
              double coef = keyval->Va_doublevalue("corr_param", 3, f, p, 2);
              vtmp.push_back(std::make_pair(std::make_pair(alpha,gamma),coef));
            }
            params.push_back(vtmp);
          }
        }
      }
      else {
        double exponent = keyval->doublevalue("corr_param");
        LinearR12::GenG12CorrelationFactor::ContractedGeminal vtmp;  vtmp.push_back(std::make_pair(std::make_pair(0.0,exponent),1.0));
        params.push_back(vtmp);
      }
      corrfactor_ = new LinearR12::GenG12CorrelationFactor(params);
    }
    else
      throw ProgrammingError("MBPT2_R12::MBPT2_R12() -- corr_param keyword must be given when corr_factor=g12",__FILE__,__LINE__);
  }
  //
  // stg-ng correlation factor
  //
  else if (corrfactor.find("stg") != string::npos || corrfactor.find("STG") != string::npos) {
    // how many gaussians?
    int ng12;
    {
	string::size_type pos1;
	pos1 = corrfactor.find("stg-");
	if (pos1 != 0)
	    pos1 = corrfactor.find("STG-");
	if (pos1 != 0)
	    throw InputError("Should specify Slater-type geminal correlation factor as STG-NG, where N is the number of Gaussians in the fit",__FILE__,__LINE__);
	// erage STG-
	string str1 = corrfactor.erase(0,4);
	// and trailing G also
	pos1 = corrfactor.find("G");
	if (pos1 == string::npos)
	    pos1 = corrfactor.find("g");
	if (pos1 == string::npos)
	    throw InputError("Should specify Slater-type geminal correlation factor as STG-NG, where N is the number of Gaussians in the fit",__FILE__,__LINE__);
	string ngtg_str = str1.erase(pos1,1);
	ng12 = std::atoi(ngtg_str.c_str());
	if (ng12 < 1)
	    throw InputError("Number of Gaussian geminals must be greater than 0",__FILE__,__LINE__);
    }

    if (!keyval->exists("corr_param"))
	throw InputError("keyword corr_param must be given when corrfactor=stg",__FILE__,__LINE__);

    std::vector<double> stg_exponents;
    typedef LinearR12::G12CorrelationFactor::CorrelationParameters CorrParams;
    CorrParams params;
    int num_f12 = keyval->count("corr_param");
    if (num_f12 != 0) {
        // Do I have contracted functions? Can't handle these (yet?)
        bool contracted = (keyval->count("corr_param",0) != 0);
        if (contracted)
	    throw FeatureNotImplemented("Cannot accept contracted STG correlation factors yet",__FILE__,__LINE__);
	
	// Primitive functions only
	for(int f=0; f<num_f12; f++) {
            double exponent = keyval->doublevalue("corr_param", f);
	    stg_exponents.push_back(exponent);
	}
    }
    else { // single exponent
	num_f12 = 1;
        double exponent = keyval->doublevalue("corr_param");
	stg_exponents.push_back(exponent);
    }
    // convert STGs into combinations of Gaussians
    for(int f=0; f<num_f12; f++) {
	using namespace sc::mbptr12;
#if INCLUDE_CUSP_IN_GTG_FIT
	// fit using weight exp(-0.005*r^6), which is flat to r=1, then falls slowly till r=2, then quickly decays to r=3
	PowerGaussian1D w(0.005,6,0);
#else
	// fit using Tew&Klopper's recipe: weight is r^2 exp(-2*r^2), which has maximum near r=0.75 and decays sharply near r=0 and r=1.5
	PowerGaussian1D w(2.0,2,2);
#endif
	typedef GaussianFit<Slater1D,PowerGaussian1D> GTGFit;
	GTGFit gtgfit(ng12, w, 0.0, 10.0, 1001);
	// fit r12^k exp(-gamma*r_{12})
	const int k = 0;
	const double gamma = stg_exponents[f];
	Ref<G12CorrelationFactor> cf;
	cf << stg_to_g12<G12CorrelationFactor,GTGFit>(gtgfit,gamma,k);
	params.push_back(cf->function(0));
    }

    // If stdapprox_ == C, no need for commutators
    if (stdapprox_ == LinearR12::StdApprox_C)
	corrfactor_ = new LinearR12::G12NCCorrelationFactor(params);
    else
	corrfactor_ = new LinearR12::G12CorrelationFactor(params);

    if (debug_ >= DefaultPrintThresholds::diagnostics) {
	corrfactor_->print();
    }
  }
  //
  // no explicit correlation
  //
  else if (corrfactor == "none" || corrfactor == "NONE") {
    corrfactor_ = new LinearR12::NullCorrelationFactor();
  }
  else
    throw FeatureNotImplemented("MBPT2_R12::MBPT2_R12 -- this correlation factor is not implemented",__FILE__,__LINE__);
  
  // Default is to assume GBC
  gbc_ = keyval->booleanvalue("gbc",KeyValValueboolean((int)true));
  // Default is to assume EBC
  ebc_ = keyval->booleanvalue("ebc",KeyValValueboolean((int)true));
  
  // Default is not to use Klopper-Samson approach to EBC
  const bool ks_ebcfree_default = false;
  if (!ebc_)
    ks_ebcfree_ = keyval->booleanvalue("ks_ebcfree",KeyValValueboolean((int)ks_ebcfree_default));
  else
    ks_ebcfree_ = ks_ebcfree_default;
  
  // Default is to include P in intermediate B
  omit_P_ = keyval->booleanvalue("omit_P",KeyValValueboolean((int)false));
  
  // For now the default is to use the old ABS method, of Klopper and Samson
  char* abs_method_str = keyval->pcharvalue("abs_method",KeyValValuepchar("ABS"));
  if ( !strcmp(abs_method_str,"KS") ||
       !strcmp(abs_method_str,"ks") ||
       !strcmp(abs_method_str,"ABS") ||
       !strcmp(abs_method_str,"abs") ) {
    abs_method_ = LinearR12::ABS_ABS;
  }
  else if ( !strcmp(abs_method_str,"KS+") ||
	    !strcmp(abs_method_str,"ks+") ||
            !strcmp(abs_method_str,"ABS+") ||
	    !strcmp(abs_method_str,"abs+") ) {
    abs_method_ = LinearR12::ABS_ABSPlus;
  }
  else if ( !strcmp(abs_method_str,"EV") ||
	    !strcmp(abs_method_str,"ev") ||
            !strcmp(abs_method_str,"CABS") ||
	    !strcmp(abs_method_str,"cabs") ) {
    abs_method_ = LinearR12::ABS_CABS;
  }
  else if ( !strcmp(abs_method_str,"EV+") ||
	    !strcmp(abs_method_str,"ev+") ||
            !strcmp(abs_method_str,"CABS+") ||
	    !strcmp(abs_method_str,"cabs+") ) {
    abs_method_ = LinearR12::ABS_CABSPlus;
  }
  else {
    delete[] abs_method_str;
    throw std::runtime_error("MBPT2_R12::MBPT2_R12 -- unrecognized value for abs_method");
  }
  delete[] abs_method_str;

  ansatz_ = require_dynamic_cast<LinearR12Ansatz*>(
    keyval->describedclassvalue("ansatz").pointer(),
    "MBPT2_R12::MBPT2_R12\n"
    );
  // Default constructor for LinearR12Ansatz specifies the default
  if (ansatz_.null())
    ansatz_ = new LinearR12Ansatz;
  if (ansatz()->projector() == LinearR12::Projector_1)
    throw InputError("MBPT2_R12::MBPT2_R12 -- projector 1 has not been implemented yet",__FILE__,__LINE__);
  if (ansatz()->projector() == LinearR12::Projector_3 &&
      stdapprox_ != LinearR12::StdApprox_C)
    throw InputError("MBPT2_R12::MBPT2_R12 -- projector 3 is only valid when stdapprox=C",__FILE__,__LINE__);
  
  // Default is to include all integrals
  maxnabs_ = static_cast<unsigned int>(keyval->intvalue("maxnabs",KeyValValueint(2)));
  
  spinadapted_ = false;
  if (closedshell)
    // Default is to use spin-adapted algorithm
    spinadapted_ = keyval->booleanvalue("spinadapted",KeyValValueboolean((int)true));

  // Default is to not compute MP1 energy
  include_mp1_ = keyval->booleanvalue("include_mp1",KeyValValueboolean((int)false));

  // Determine how to store MO integrals
  char *r12ints_str = keyval->pcharvalue("r12ints",KeyValValuepchar("mem-posix"));
  if (!strcmp(r12ints_str,"mem")) {
    r12ints_method_ = R12IntEvalInfo::StoreMethod::mem_only;
  }
#if HAVE_MPIIO
  else if (!strcmp(r12ints_str,"mem-mpi")) {
    r12ints_method_ = MOIntsTransformFactory::StoreMethod::mem_mpi;
  }
  else if (!strcmp(r12ints_str,"mpi")) {
    r12ints_method_ = MOIntsTransformFactory::StoreMethod::mpi;
  }
#else
  else if ( !strcmp(r12ints_str,"mem-mpi") ||
	    !strcmp(r12ints_str,"mpi") ) {
    throw std::runtime_error("MBPT2_R12::MBPT2_R12 -- the value for keyword r12ints is not valid in this environment (no MPI-I/O detected)");
  }
#endif
  else if (!strcmp(r12ints_str,"mem-posix")) {
    r12ints_method_ = R12IntEvalInfo::StoreMethod::mem_posix;
  }
  else if (!strcmp(r12ints_str,"posix")) {
    r12ints_method_ = R12IntEvalInfo::StoreMethod::posix;
  }
  else {
    delete[] r12ints_str;
    throw std::runtime_error("MBPT2_R12::MBPT2_R12 -- invalid value for keyword r12ints");
  }
  delete[] r12ints_str;
  
  // Get the prefix for the filename to store the integrals
  std::string r12ints_file_default("./");
  r12ints_file_ = keyval->stringvalue("r12ints_file",KeyValValuestring(r12ints_file_default));
  // if the last character of r12ints_file is '/' then append the default basename
  if (*(r12ints_file_.rbegin()) == '/')
    r12ints_file_ += std::string(SCFormIO::fileext_to_filename(".moints"));

  twopdm_grid_ = require_dynamic_cast<TwoBodyGrid*>(
                   keyval->describedclassvalue("twopdm_grid").pointer(),
                   "MBPT2_R12::MBPT2_R12\n"
                 );
  if (twopdm_grid_.nonnull()) {
    plot_pair_function_[0] = static_cast<unsigned int>(keyval->intvalue("plot_pair_function", 0,
									KeyValValueint(0)
									));
    plot_pair_function_[1] = static_cast<unsigned int>(keyval->intvalue("plot_pair_function", 1,
									KeyValValueint(0)
									));
  }

  safety_check_ = keyval->booleanvalue("safety_check",KeyValValueboolean((int)true));

  std::string posdef_B = keyval->stringvalue("posdef_B",KeyValValuestring("weak"));
  if (posdef_B == "no" || posdef_B == "NO" || posdef_B == "false" || posdef_B == "FALSE") {
      posdef_B_ = LinearR12::PositiveDefiniteB_no;
  }
  else if (posdef_B == "yes" || posdef_B == "YES" || posdef_B == "true" || posdef_B == "TRUE") {
      posdef_B_ = LinearR12::PositiveDefiniteB_yes;
  }
  else if (posdef_B == "weak" || posdef_B == "WEAK") {
      posdef_B_ = LinearR12::PositiveDefiniteB_weak;
  }
  else {
      throw InputError("MBPT2_R12::MBPT2_R12 -- invalid value for keyword posdef_B",__FILE__,__LINE__);
  }

  //
  //
  // Check that requested features are compatible/allowed
  //
  //

  // stdapprox must be C if corrfactor == geng12
  {
    Ref<LinearR12::GenG12CorrelationFactor> gg12ptr; gg12ptr << corrfactor_;
    if (gg12ptr.nonnull() && stdapprox_ != LinearR12::StdApprox_C) {
	throw InputError("MBPT2_R12::MBPT2_R12() -- stdapprox must be set to C when using general Geminal correlation factor",__FILE__,__LINE__);
    }
  }
  
  // Klopper and Samson's ABS method is only implemented for certain "old" methods
  // Make sure that the ABS method is available for the requested MP2-R12 energy
  const bool must_use_cabs = (!gbc_ ||
			      !ebc_ ||
			      (stdapprox_ == LinearR12::StdApprox_B && !abs_eq_obs) ||
			      (stdapprox_ == LinearR12::StdApprox_C && !abs_eq_obs) ||
			      (stdapprox_ == LinearR12::StdApprox_App && !abs_eq_obs) ||
                              !basis()->equiv(vir_basis_));
  if (must_use_cabs &&
      (abs_method_ == LinearR12::ABS_ABS || abs_method_ == LinearR12::ABS_ABSPlus))
    throw std::runtime_error("MBPT2_R12::MBPT2_R12() -- abs_method must be set to cabs or cabs+ for this MP2-R12 method");

  // Standard approximation A is not valid when gbc_ = false or ebc_ = false
  if ( (!gbc_ || !ebc_) && stdapprox_ == LinearR12::StdApprox_A )
    throw std::runtime_error("MBPT2_R12::MBPT2_R12() -- stdapprox=A is not valid when gbc_ = false or ebc_ = false");
    

  // Most calculations must use disk for now
  bool must_use_disk = true;
  // except MP2
  {
    Ref<LinearR12::NullCorrelationFactor> nullptr; nullptr << corrfactor_;
    if (nullptr.nonnull())
      must_use_disk = false;
  }

  if (must_use_disk && r12ints_method_ == R12IntEvalInfo::StoreMethod::mem_only)
    throw std::runtime_error("MBPT2_R12::MBPT2_R12 -- r12ints=mem is only possible for MP2-R12/A and MP2-R12/A' (GBC+EBC) methods");
  if (must_use_disk) {
    if (r12ints_method_ == R12IntEvalInfo::StoreMethod::mem_posix)
      r12ints_method_ = R12IntEvalInfo::StoreMethod::posix;
#if HAVE_MPIIO
    if (r12ints_method_ == R12IntEvalInfo::StoreMethod::mem_mpi)
      r12ints_method_ = R12IntEvalInfo::StoreMethod::mpi;
#endif
  }

  r12eval_ = 0;
  r12a_energy_ = 0;
  r12ap_energy_ = 0;
  r12app_energy_ = 0;
  r12b_energy_ = 0;
  r12c_energy_ = 0;
  mp2_corr_energy_ = 0.0;
  r12_corr_energy_ = 0.0;
  
}

MBPT2_R12::~MBPT2_R12()
{
}

void
MBPT2_R12::save_data_state(StateOut& s)
{
  MBPT2::save_data_state(s);
  SavableState::save_state(r12eval_.pointer(),s);
  SavableState::save_state(r12a_energy_.pointer(),s);
  SavableState::save_state(r12ap_energy_.pointer(),s);
  SavableState::save_state(r12app_energy_.pointer(),s);
  SavableState::save_state(r12b_energy_.pointer(),s);
  SavableState::save_state(r12c_energy_.pointer(),s);
  SavableState::save_state(aux_basis_.pointer(),s);
  SavableState::save_state(vir_basis_.pointer(),s);
  s.put((int)gbc_);
  s.put((int)ebc_);
  s.put((int)ks_ebcfree_);
  s.put((int)omit_P_);
  s.put((int)abs_method_);
  s.put((int)stdapprox_);
  SavableState::save_state(ansatz_.pointer(),s);
  s.put((int)spinadapted_);
  s.put((int)include_mp1_);
  s.put((int)r12ints_method_);
  s.put(r12ints_file_);
  s.put((int)safety_check_);
  s.put((int)posdef_B_);

  s.put(mp2_corr_energy_);
  s.put(r12_corr_energy_);

  SavableState::save_state(twopdm_grid_.pointer(),s);
  s.put((int)plot_pair_function_[0]);
  s.put((int)plot_pair_function_[1]);
}

void
MBPT2_R12::print(ostream&o) const
{
  o << indent << "MBPT2_R12:" << endl;
  o << incindent;

  if (!safety_check())
    o << indent << "WARNING: ---- safety checks SKIPPED ----" << endl;
  corrfactor()->print(o); o << endl;
  o << indent << "GBC assumed: " << (gbc_ ? "true" : "false") << endl;
  o << indent << "EBC assumed: " << (ebc_ ? "true" : "false") << endl;
  o << indent << "EBC-free method: " << (!ks_ebcfree_ ? "Valeev" : "Klopper and Samson") << endl;
  switch (posdef_B()) {
    case LinearR12::PositiveDefiniteB_no:     o << indent << "Do not enforce positive definiteness of B" << endl;  break;      
    case LinearR12::PositiveDefiniteB_yes:    o << indent << "Enforce positive definiteness of B" << endl;  break;      
    case LinearR12::PositiveDefiniteB_weak:   o << indent << "Enforce positive definiteness of B, but not ~B(ij)" << endl;  break;      
  }
  if (stdapprox_ == LinearR12::StdApprox_B && omit_P_) {
    o << indent << "Intermediate P is omitted" << endl;
  }
  switch(abs_method_) {
  case LinearR12::ABS_ABS :
    o << indent << "ABS method variant: ABS  (Klopper and Samson)" << endl;
    break;
  case LinearR12::ABS_ABSPlus :
    o << indent << "ABS method variant: ABS+ (Klopper and Samson using the union of OBS and ABS for RI)" << endl;
    break;
  case LinearR12::ABS_CABS :
    o << indent << "ABS method variant: CABS  (Valeev)" << endl;
    break;
  case LinearR12::ABS_CABSPlus :
    o << indent << "ABS method variant: CABS+ (Valeev using the union of OBS and ABS for RI)" << endl;
    break;
  }
  switch (stdapprox_) {
    case LinearR12::StdApprox_A :
      o << indent << "Standard Approximation: A" << endl;
    break;
    case LinearR12::StdApprox_Ap :
      o << indent << "Standard Approximation: A'" << endl;
    break;
    case LinearR12::StdApprox_App :
      o << indent << "Standard Approximation: A''" << endl;
    break;
    case LinearR12::StdApprox_B :
      o << indent << "Standard Approximation: B" << endl;
    break;
    case LinearR12::StdApprox_C :
      o << indent << "Standard Approximation: C" << endl;
    break;
  }
  ansatz()->print(o);

  o << indent << "Max # ABS indices: " << maxnabs_ << endl;

  o << indent << "Spin-adapted algorithm: " << (spinadapted_ ? "true" : "false") << endl;
  if (!vir_basis_->equiv(basis()))
    o << indent << "Compute MP1 energy: " << (include_mp1_ ? "true" : "false") << endl;

  std::string r12ints_str;
  switch (r12ints_method_) {
  case R12IntEvalInfo::StoreMethod::mem_only:
    r12ints_str = "mem"; break;
  case R12IntEvalInfo::StoreMethod::mem_posix:
    r12ints_str = "mem_posix"; break;
  case R12IntEvalInfo::StoreMethod::posix:
    r12ints_str = "posix"; break;
#if HAVE_MPIIO
  case R12IntEvalInfo::StoreMethod::mem_mpi:
    r12ints_str = "mem-mpi"; break;
  case R12IntEvalInfo::StoreMethod::mpi:
    r12ints_str = "mpi"; break;
#endif
  default:
    throw std::runtime_error("MBPT2_R12::print -- invalid value of r12ints_method_");
  }
  o << indent << "How to Store Transformed Integrals: " << r12ints_str << endl << endl;
  o << indent << "Transformed Integrals file suffix: " << r12ints_file_ << endl << endl;
  o << indent << "Auxiliary Basis Set (ABS):" << endl;
  o << incindent; aux_basis_->print(o); o << decindent << endl;
  o << indent << " Virtuals Basis Set (VBS):" << endl;
  o << incindent; vir_basis_->print(o); o << decindent << endl;
  MBPT2::print(o);
  o << decindent;
}

RefSymmSCMatrix
MBPT2_R12::density()
{
  ExEnv::out0() << "MBPT2_R12::density() is not implemented" << endl;
  abort();
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

void
MBPT2_R12::compute()
{
  if (std::string(reference_->integral()->class_name()) != integral()->class_name()) {
      FeatureNotImplemented ex(
          "cannot use a reference with a different Integral specialization",
          __FILE__, __LINE__, class_desc());
      try {
          ex.elaborate()
              << "reference uses " << reference_->integral()->class_name()
              << " but this object uses " << integral()->class_name()
              << std::endl;
        }
      catch (...) {}
      throw ex;
    }

  init_variables_();
  reference_->set_desired_value_accuracy(desired_value_accuracy()
                                         / ref_to_mp2r12_acc_);

  compute_energy_();
}

//////////////////////////////////////////////////////////////////////////////

void
MBPT2_R12::obsolete()
{
  r12eval_ = 0;
  r12a_energy_ = 0;
  r12ap_energy_ = 0;
  r12app_energy_ = 0;
  r12b_energy_ = 0;
  r12c_energy_ = 0;
  mp2_corr_energy_ = 0.0;
  r12_corr_energy_ = 0.0;
  MBPT2::obsolete();
}

//////////////////////////////////////////////////////////////////////////////

int
MBPT2_R12::gradient_implemented() const
{
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

int
MBPT2_R12::value_implemented() const
{
  return 1;
}

////////////////////////////////////////////////////////////////////////////

const Ref<LinearR12::CorrelationFactor>&
MBPT2_R12::corrfactor() const
{
  return corrfactor_;
}

////////////////////////////////////////////////////////////////////////////

void
MBPT2_R12::corrfactor(const Ref<LinearR12::CorrelationFactor>& cf)
{
  if (!corrfactor_->equiv(cf)) {
      corrfactor_ = cf;
      obsolete();
  }
}

/////////////////////////////////////////////////////////////////////////////

Ref<GaussianBasisSet>
MBPT2_R12::aux_basis() const
{
  return aux_basis_;
}

/////////////////////////////////////////////////////////////////////////////

Ref<GaussianBasisSet>
MBPT2_R12::vir_basis() const
{
  return vir_basis_;
}

/////////////////////////////////////////////////////////////////////////////

bool
MBPT2_R12::ks_ebcfree() const
{
  return ks_ebcfree_;
}

/////////////////////////////////////////////////////////////////////////////

bool
MBPT2_R12::omit_P() const
{
  return omit_P_;
}

/////////////////////////////////////////////////////////////////////////////

unsigned int
MBPT2_R12::maxnabs() const
{
  return maxnabs_;
}

/////////////////////////////////////////////////////////////////////////////

bool
MBPT2_R12::include_mp1() const
{
  return include_mp1_;
}

/////////////////////////////////////////////////////////////////////////////

bool
MBPT2_R12::gbc() const
{
  return gbc_;
}

/////////////////////////////////////////////////////////////////////////////

bool
MBPT2_R12::ebc() const
{
  return ebc_;
}

/////////////////////////////////////////////////////////////////////////////

LinearR12::ABSMethod
MBPT2_R12::abs_method() const
{
  return abs_method_;
}

/////////////////////////////////////////////////////////////////////////////

LinearR12::StandardApproximation
MBPT2_R12::stdapprox() const
{
  return stdapprox_;
}

/////////////////////////////////////////////////////////////////////////////

const Ref<LinearR12Ansatz>&
MBPT2_R12::ansatz() const
{
  return ansatz_;
}

/////////////////////////////////////////////////////////////////////////////

bool
MBPT2_R12::spinadapted() const
{
  return spinadapted_;
}

/////////////////////////////////////////////////////////////////////////////

R12IntEvalInfo::StoreMethod::type
MBPT2_R12::r12ints_method() const
{
  return r12ints_method_;
}

/////////////////////////////////////////////////////////////////////////////

const std::string&
MBPT2_R12::r12ints_file() const
{
  return r12ints_file_;
}

/////////////////////////////////////////////////////////////////////////////

bool
MBPT2_R12::safety_check() const
{
  return safety_check_;
}

/////////////////////////////////////////////////////////////////////////////

const LinearR12::PositiveDefiniteB&
MBPT2_R12::posdef_B() const
{
  return posdef_B_;
}

/////////////////////////////////////////////////////////////////////////////

double
MBPT2_R12::corr_energy()
{
  energy();
  return energy() - ref_energy();
}

/////////////////////////////////////////////////////////////////////////////

double
MBPT2_R12::r12_corr_energy()
{
  energy();
  return energy() - mp2_corr_energy_ - ref_energy();
}

/////////////////////////////////////////////////////////////////////////////

void
MBPT2_R12::init_variables_()
{
}

/////////////////////////////////////////////////////////////////////////////

void
MBPT2_R12::check_integral_factory_()
{
  // Only IntegralCints or IntegralLibint2 can be used at the moment
  bool allowed_integral_factory = false;
#if HAVE_INTEGRALCINTS
  IntegralCints* cintsintf = dynamic_cast<IntegralCints*>(integral().pointer());
  if (cintsintf) {
    allowed_integral_factory = true;
  }
#endif
#if HAVE_INTEGRALLIBINT2
  IntegralLibint2* libint2intf = dynamic_cast<IntegralLibint2*>(integral().pointer());
  if (libint2intf) {
    allowed_integral_factory = true;
  }
#endif
  if(strcmp(integral()->class_name(),"IntegralCCA") == 0) {
    allowed_integral_factory = true;
  }
  if (!allowed_integral_factory) {
    InputError ex("MBPT2_R12::check_integral_factory_(): invalid integral factory provided.",
                  __FILE__, __LINE__, 0, 0, class_desc());
    try {
      ex.elaborate() << "Try using IntegralCints, IntegralCCA, or IntegralLibint2."
                     << std::endl
                     << "The given integral factory was of type " << integral()->class_name()
                     << std::endl;
    }
    catch (...) {}
    throw ex;
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
