//
// r12technology.cc
//
// Copyright (C) 2007 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
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

#include <util/state/statein.h>
#include <util/state/stateout.h>
#include <chemistry/qc/basis/integral.h>
#  include <chemistry/qc/intv3/intv3.h>
#if HAVE_INTEGRALCINTS
#  include <chemistry/qc/cints/cints.h>
#endif
#if HAVE_INTEGRALLIBINT2
#  include <chemistry/qc/libint2/libint2.h>
#endif
#include <chemistry/qc/mbptr12/gaussianfit.h>
#include <chemistry/qc/mbptr12/gaussianfit.timpl.h>
#include <chemistry/qc/mbptr12/linearr12.timpl.h>
#include <chemistry/qc/mbptr12/print.h>
#include <chemistry/qc/mbptr12/r12technology.h>

// set to 1 to include the cusp region in the fit, otherwise use Tew&Klopper's strategy
#define INCLUDE_CUSP_IN_GTG_FIT 0

using namespace std;
using namespace sc;
using namespace sc::LinearR12;

/*--------------------------------
  MBPT2_R12
 --------------------------------*/

static ClassDesc R12Technology_cd(
  typeid(R12Technology),"R12Technology",8,"public MBPT2",
  0, 0, create<R12Technology>);

R12Technology::R12Technology(StateIn& s)
{
  int vbs_eq_obs; s.get(vbs_eq_obs); vbs_eq_obs_ = (bool)vbs_eq_obs;
  int abs_eq_obs; s.get(abs_eq_obs); abs_eq_obs_ = (bool)abs_eq_obs;

  int gbc; s.get(gbc); gbc_ = (bool)gbc;
  int ebc; s.get(ebc); ebc_ = (bool)ebc;
  int ks_ebcfree; s.get(ks_ebcfree); ks_ebcfree_ = (bool)ks_ebcfree;
  int omit_P; s.get(omit_P); omit_P_ = (bool)omit_P;
  int absmethod; s.get(absmethod); abs_method_ = (LinearR12::ABSMethod)absmethod;
  int stdapprox; s.get(stdapprox); stdapprox_ = (LinearR12::StandardApproximation) stdapprox;
  ansatz_ << SavableState::restore_state(s);
  s.get(maxnabs_);

  int safety_check; s.get(safety_check);
  safety_check_ = static_cast<bool>(safety_check);
  int posdef_B; s.get(posdef_B);
  posdef_B_ = static_cast<LinearR12::PositiveDefiniteB>(posdef_B);
}

R12Technology::R12Technology(const Ref<KeyVal>& keyval,
			     const Ref<GaussianBasisSet>& obs,
			     const Ref<GaussianBasisSet>& vbs,
			     const Ref<GaussianBasisSet>& abs)
{
  abs_eq_obs_ = abs->equiv(obs);
  vbs_eq_obs_ = vbs->equiv(obs);

  // Default is to use the R12 factor
  std::string corrfactor = keyval->stringvalue("corr_factor", KeyValValuestring("r12"));

  // Default method is MBPT2-R12/A'
  char *sa_string = keyval->pcharvalue("stdapprox",KeyValValuepchar("A'"));
  if ( !strcmp(sa_string,"A") ||
       !strcmp(sa_string,"a") ) {
    throw FeatureNotImplemented("stdapprox=A is obsolete",__FILE__,__LINE__);
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
    throw std::runtime_error("R12Technology::R12Technology() -- unrecognized value for stdapprox");
  }

  // if no explicit correlation then set stdapprox to A'
  if (sa_string == "none" || sa_string == "NONE") {
      stdapprox_ = LinearR12::StdApprox_Ap;
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
      throw ProgrammingError("R12Technology::R12Technology() -- corr_param keyword must be given when corr_factor=g12",__FILE__,__LINE__);
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
      throw ProgrammingError("R12Technology::R12Technology() -- corr_param keyword must be given when corr_factor=g12",__FILE__,__LINE__);
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
  }
  //
  // no explicit correlation
  //
  else if (corrfactor == "none" || corrfactor == "NONE") {
    corrfactor_ = new LinearR12::NullCorrelationFactor();
  }
  else
    throw FeatureNotImplemented("R12Technology::R12Technology -- this correlation factor is not implemented",__FILE__,__LINE__);
  
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
    throw std::runtime_error("R12Technology::R12Technology -- unrecognized value for abs_method");
  }
  delete[] abs_method_str;

  ansatz_ = require_dynamic_cast<LinearR12Ansatz*>(
    keyval->describedclassvalue("ansatz").pointer(),
    "R12Technology::R12Technology\n"
    );
  // Default constructor for LinearR12Ansatz specifies the default
  if (ansatz_.null())
    ansatz_ = new LinearR12Ansatz;
  if (ansatz()->projector() == LinearR12::Projector_1)
    throw InputError("R12Technology::R12Technology -- projector 1 has not been implemented yet",__FILE__,__LINE__);
  if (ansatz()->projector() == LinearR12::Projector_3 &&
      stdapprox_ != LinearR12::StdApprox_C)
    throw InputError("R12Technology::R12Technology -- projector 3 is only valid when stdapprox=C",__FILE__,__LINE__);
  
  // Default is to include all integrals, unless using A'' method
  const int default_maxnabs = (stdapprox_ == LinearR12::StdApprox_App) ? 1 : 2;
  maxnabs_ = static_cast<unsigned int>(keyval->intvalue("maxnabs",KeyValValueint(default_maxnabs)));
  
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
      throw InputError("R12Technology::R12Technology -- invalid value for keyword posdef_B",__FILE__,__LINE__);
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
	throw InputError("R12Technology::R12Technology() -- stdapprox must be set to C when using general Geminal correlation factor",__FILE__,__LINE__);
    }
  }
  
  // Klopper and Samson's ABS method is only implemented for certain "old" methods
  // Make sure that the ABS method is available for the requested MP2-R12 energy
  const bool must_use_cabs = (!gbc_ ||
			      !ebc_ ||
			      (stdapprox_ == LinearR12::StdApprox_B && !abs_eq_obs_) ||
			      (stdapprox_ == LinearR12::StdApprox_C && !abs_eq_obs_) ||
			      (stdapprox_ == LinearR12::StdApprox_App && !abs_eq_obs_) ||
                              !vbs_eq_obs_);
  if (must_use_cabs &&
      (abs_method_ == LinearR12::ABS_ABS || abs_method_ == LinearR12::ABS_ABSPlus))
    throw std::runtime_error("R12Technology::R12Technology() -- abs_method must be set to cabs or cabs+ for this MP2-R12 method");

}

R12Technology::~R12Technology()
{
}

void
R12Technology::save_data_state(StateOut& s)
{
  s.put((int)gbc_);
  s.put((int)ebc_);
  s.put((int)ks_ebcfree_);
  s.put((int)omit_P_);
  s.put((int)abs_method_);
  s.put((int)stdapprox_);
  SavableState::save_state(ansatz_.pointer(),s);
  s.put((int)safety_check_);
  s.put((int)posdef_B_);

}

void
R12Technology::print(ostream&o) const
{
  o << indent << "R12Technology:" << endl;
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

  o << decindent;
}

////////////////////////////////////////////////////////////////////////////

const Ref<LinearR12::CorrelationFactor>&
R12Technology::corrfactor() const
{
  return corrfactor_;
}

////////////////////////////////////////////////////////////////////////////

void
R12Technology::corrfactor(const Ref<LinearR12::CorrelationFactor>& cf)
{
  if (!corrfactor_->equiv(cf)) {
      corrfactor_ = cf;
  }
}

/////////////////////////////////////////////////////////////////////////////

bool
R12Technology::ks_ebcfree() const
{
  return ks_ebcfree_;
}

/////////////////////////////////////////////////////////////////////////////

bool
R12Technology::omit_P() const
{
  return omit_P_;
}

/////////////////////////////////////////////////////////////////////////////

unsigned int
R12Technology::maxnabs() const
{
  return maxnabs_;
}

/////////////////////////////////////////////////////////////////////////////

bool
R12Technology::gbc() const
{
  return gbc_;
}

/////////////////////////////////////////////////////////////////////////////

bool
R12Technology::ebc() const
{
  return ebc_;
}

/////////////////////////////////////////////////////////////////////////////

LinearR12::ABSMethod
R12Technology::abs_method() const
{
  return abs_method_;
}

/////////////////////////////////////////////////////////////////////////////

LinearR12::StandardApproximation
R12Technology::stdapprox() const
{
  return stdapprox_;
}

/////////////////////////////////////////////////////////////////////////////

const Ref<LinearR12Ansatz>&
R12Technology::ansatz() const
{
  return ansatz_;
}

/////////////////////////////////////////////////////////////////////////////

bool
R12Technology::safety_check() const
{
  return safety_check_;
}

/////////////////////////////////////////////////////////////////////////////

const LinearR12::PositiveDefiniteB&
R12Technology::posdef_B() const
{
  return posdef_B_;
}

/////////////////////////////////////////////////////////////////////////////

void
R12Technology::check_integral_factory(const Ref<Integral>& ints)
{
  // Only IntegralCints or IntegralLibint2 can be used at the moment
  bool allowed_integral_factory = false;
#if HAVE_INTEGRALCINTS
  IntegralCints* cintsintf = dynamic_cast<IntegralCints*>(ints.pointer());
  if (cintsintf) {
    allowed_integral_factory = true;
  }
#endif
#if HAVE_INTEGRALLIBINT2
  IntegralLibint2* libint2intf = dynamic_cast<IntegralLibint2*>(ints.pointer());
  if (libint2intf) {
    allowed_integral_factory = true;
  }
#endif
  if(strcmp(ints->class_name(),"IntegralCCA") == 0) {
    allowed_integral_factory = true;
  }
  if (!allowed_integral_factory) {
    InputError ex("R12Technology::check_integral_factory_(): invalid integral factory provided.",
                  __FILE__, __LINE__, 0, 0, class_desc());
    try {
      ex.elaborate() << "Try using IntegralCints, IntegralCCA, or IntegralLibint2."
                     << std::endl
                     << "The given integral factory was of type " << ints->class_name()
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
