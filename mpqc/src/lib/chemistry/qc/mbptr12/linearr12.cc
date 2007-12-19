//
// linearr12.cc
//
// Copyright (C) 2005 Edward Valeev
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

#ifdef __GNUG__
#pragma implementation
#endif

#include <strstream>
#include <sstream>
#include <util/misc/formio.h>
#include <util/class/scexception.h>
#include <util/ref/ref.h>
#include <chemistry/qc/mbptr12/linearr12.h>
#include <chemistry/qc/mbptr12/linearr12.timpl.h>
#include <chemistry/qc/mbptr12/gaussianfit.h>

using namespace sc;
using namespace LinearR12;

namespace sc {
  namespace LinearR12 {
    template <>
    bool
    CorrParamCompare<IntParamsG12>::equiv(const PrimitiveGeminal& A, const PrimitiveGeminal& B) {
      return (std::fabs(A.first-B.first) < epsilon && std::fabs(A.second-B.second) < epsilon);
    }
    
    template <>
    bool
    CorrParamCompare<IntParamsGenG12>::equiv(const PrimitiveGeminal& A, const PrimitiveGeminal& B) {
      return (std::fabs(A.first.first-B.first.first) < epsilon &&
	      std::fabs(A.first.second-B.first.second) < epsilon &&
	      std::fabs(A.second-B.second) < epsilon);
    }

    Ref<CorrelationFactor> ang_to_geng12(double alpha) {
      
      const double halfalpha = alpha/2.0;
      
      // feed to the constructor of CorrFactor
      typedef IntParamsGenG12::PrimitiveGeminal PrimitiveGeminal;
      typedef IntParamsGenG12::ContractedGeminal ContractedGeminal;
      ContractedGeminal geminal_ang;
      
      // add ang
      geminal_ang.push_back(std::make_pair(std::make_pair(halfalpha,-halfalpha),1.0));
      
      std::vector<ContractedGeminal> geminals;
      geminals.push_back(geminal_ang);
      
      Ref<CorrelationFactor> cf = new GenG12CorrelationFactor(geminals);
      return cf;
    }

    /*********************
     * GeminalDescriptor *
     *********************/
    GeminalDescriptor::GeminalDescriptor(){
      type_ = "invalid";
    }
    
    GeminalDescriptor::GeminalDescriptor(const std::string& type, const std::vector<std::string> &params){
      type_ = type;
      params_ = params;
      //compute_offsets();
    }
    
    GeminalDescriptor::GeminalDescriptor(const GeminalDescriptor& source){
      type_ = source.type_;
      params_ = source.params_;
      //compute_offsets();
    }
    
    //void GeminalDescriptor::compute_offsets() {
    //  if(type_!=std::string("invalid") && type_!=std::string("r12") && type_!=std::string("R12")){
    //    int nfunction=atoi(params_[0].c_str());
    //    offsets_=std::vector<int>(nfunction+1);
    //    int cumul_ind=nfunction+1;
    //    int nprimitve;
    //    for(int i=0; i<nfunction; i++){
    //      offsets_[i]=cumul_ind;
    //      nprimitve=atoi(params_[i+1].c_str());
    //      cumul_ind+=2*nprimitve;
    //    }
    //    offsets_[nfunction]=cumul_ind;
    //  }
    //}
    
    std::string GeminalDescriptor::type() const {
      return(type_);
    }
    
    std::vector<std::string> GeminalDescriptor::params() const {
      return(params_);
    }
    
    void GeminalDescriptor::print(std::ostream &o) {
      o << "GeminalDescriptor :" << std::endl;
      o << "type_ = " << type_ << std::endl;
      o << "params_ :" << std::endl;
      for(int i=0; i<params_.size(); i++){
        o << params_[i] << std::endl;
      }
    }
    
    bool invalid(const Ref<GeminalDescriptor>& gdesc){
      std::string type=gdesc->type();
      return((type==std::string("invalid")) ? true : false);
    }
    
    bool R12(const Ref<GeminalDescriptor>& gdesc){
      std::string type=gdesc->type();
      return(((type==std::string("R12")) || (type==std::string("r12"))) ? true : false);
    }
    
    bool STG(const Ref<GeminalDescriptor>& gdesc){
      std::string type=gdesc->type();
      return(((type==std::string("STG")) || (type==std::string("stg"))) ? true : false);
    }
    
    bool G12(const Ref<GeminalDescriptor>& gdesc){
      std::string type=gdesc->type();
      if((type==std::string("G12")) || (type==std::string("g12"))){
        return(true);
      }
      else {
        return(false);
      }
    }
    
    double single_slater_exponent(const Ref<GeminalDescriptor>& gdesc) {
      //gdesc->print();
      std::string type=gdesc->type();
      std::vector<std::string> params=gdesc->params();
      if(type!=std::string("STG")){
        throw ProgrammingError("GeminalDescriptor::single_slater_exponent(): Geminal is not not of Slater function type.",__FILE__,__LINE__);
      }
      int nfunctions=atoi(params[0].c_str());
      if(nfunctions!=1){
        throw ProgrammingError("GeminalDescriptor::single_slater_exponent(): There are more than one Slater type functions.",__FILE__,__LINE__);
      }
      int nprimitives=atoi(params[1].c_str());
      double exponent=atof(params[3].c_str());
      ExEnv::out0() << "exponent = " << exponent << std::endl;
      
      return(exponent);
    }
    
    /****************************
     * GeminalDescriptorFactory *
     ****************************/
    GeminalDescriptorFactory::GeminalDescriptorFactory() 
      : invalid_id_("invalid"),
        r12_id_("R12"),
        stg_id_("STG"),
        g12_id_("G12") {}
    
    Ref<GeminalDescriptor> GeminalDescriptorFactory::null_geminal(){
      std::vector<std::string> void_vector;
      return(Ref<GeminalDescriptor>(new GeminalDescriptor(std::string(invalid_id_),void_vector)));
    }
    
    Ref<GeminalDescriptor> GeminalDescriptorFactory::r12_geminal(){
      std::vector<std::string> void_vector;
      return(Ref<GeminalDescriptor>(new GeminalDescriptor(std::string(r12_id_),void_vector)));
    }
    
    Ref<GeminalDescriptor> GeminalDescriptorFactory::slater_geminal(double gamma){
      std::vector<std::string> param_vec(3);
      int nfunction=1;
      int nprimitive=1;
      std::strstream inout;
      inout << nfunction;
      inout >> param_vec[0];
      inout << nprimitive;
      inout >> param_vec[1];
      inout << std::setprecision(16) << gamma;
      inout >> param_vec[2];
      return(Ref<GeminalDescriptor>(new GeminalDescriptor(std::string(stg_id_),param_vec)));
    }
    
    Ref<GeminalDescriptor> GeminalDescriptorFactory::slater_geminal(const std::vector<double> &gamma) {
      int nfunction=gamma.size();
      std::vector<std::string> params(1+3*nfunction);
      std::strstream inout;
      inout << nfunction;
      inout >> params[0];
      int one=1;
      std::string one_str;
      std::ostringstream out;
      out << one;
      one_str=out.str();
      out.str("");
      double oned=1.0;
      std::string oned_str;
      out << std::setprecision(16) << oned;
      oned_str=out.str();
      out.str("");
      ExEnv::out0() << "oned_str = " << oned_str << std::endl;
      
      for(int i=0; i<nfunction; i++){
        params[i+1]=one_str;
        params[nfunction+1+2*i]=oned_str;
        out << std::setprecision(16) << gamma[i];
        params[nfunction+2+2*i]=out.str();
        out.str("");
      }
      
      return(Ref<GeminalDescriptor>(new GeminalDescriptor(std::string(stg_id_),params)));
    }
    
    Ref<GeminalDescriptor> GeminalDescriptorFactory::gaussian_geminal(double gamma){
      std::vector<std::string> param_vec(3);
      int nfunction=1;
      int nprimitive=1;
      std::strstream inout;
      inout << nfunction;
      inout >> param_vec[0];
      inout << nprimitive;
      inout >> param_vec[1];
      inout << std::setprecision(16) << gamma;
      inout >> param_vec[2];
      return(Ref<GeminalDescriptor>(new GeminalDescriptor(std::string(g12_id_),param_vec)));      
    }
    
    Ref<GeminalDescriptor> GeminalDescriptorFactory::contracted_gaussian_geminal(const std::vector<double> &coeff,
                                                                                 const std::vector<double> &gamma){
      int nfunction=1;
      unsigned int ngeminal=coeff.size();
      std::vector<std::string> param_vec(2*ngeminal+2);
      std::strstream inout;
      inout << nfunction;
      inout >> param_vec[0];
      inout << ngeminal;
      inout >> param_vec[1];
      for(int i=0; i<ngeminal; i++){
        inout << std::setprecision(16) << coeff[i];
        inout >> param_vec[2*i+2];
        inout << std::setprecision(16) << gamma[i];
        inout >> param_vec[2*i+3];
      }
      return(Ref<GeminalDescriptor>(new GeminalDescriptor(std::string(g12_id_),param_vec)));
    }
    
    Ref<GeminalDescriptor> GeminalDescriptorFactory::gaussian_geminal(const LinearR12::G12CorrelationFactor::CorrelationParameters &corrparams){
      unsigned int nfunction=corrparams.size();
      std::vector<int> offsets(nfunction+1);
      int numofparams=nfunction+1;
      int nprimitive;
      for(int i=0; i<nfunction; i++){
        offsets[i]=numofparams;
        nprimitive=corrparams[i].size();
        numofparams+=2*nprimitive;
      }
      offsets[nfunction]=numofparams;
      
      std::vector<std::string> params(numofparams);
      std::strstream inout;
      inout << nfunction;
      inout >> params[0];
      for(int i=0; i<nfunction; i++){
        nprimitive=corrparams[i].size();
        inout << nprimitive;
        inout >> params[i+1];
        for(int j=0; j<nprimitive; j++){
          double coefficient=corrparams[i][j].second;
          double exponent=corrparams[i][j].first;
          inout << std::setprecision(16) << coefficient;
          inout >> params[offsets[i]+2*j];
          inout << std::setprecision(16) << exponent;
          inout >> params[offsets[i]+2*j+1];
        }
      }
      return(Ref<GeminalDescriptor>(new GeminalDescriptor(std::string(g12_id_),params)));
    }

  };
};

CorrelationFactor::CorrelationFactor(const std::string& label, const Ref<GeminalDescriptor> &geminaldescriptor) :
  label_(label),
  geminaldescriptor_(geminaldescriptor)
{
}

CorrelationFactor::~CorrelationFactor()
{
}

CorrelationFactor::CorrelationFactor()
{
  label_=std::string("invalid");
  Ref<GeminalDescriptorFactory> gdesc_factory=new GeminalDescriptorFactory;
  geminaldescriptor_=gdesc_factory->null_geminal();
}

unsigned int
CorrelationFactor::nfunctions() const
{
  return 1;
}

unsigned int
CorrelationFactor::nprimitives(unsigned int c) const
{
  return 1;
}

const std::string&
CorrelationFactor::label() const
{
  return label_;
}

void
CorrelationFactor::print(std::ostream& os) const
{
  using std::endl;
  os << indent << "CorrelationFactor:" << endl;
  os << incindent;
  const int nfunc = nfunctions();
  for(int f=0; f<nfunc; f++) {
    os << indent << "Function " << f << ":" << endl << incindent;
    os << indent << "Functional form: " << label() << endl;
    print_params(os,f);
    //geminaldescriptor_->print(os);
    os << decindent;
  }
  os << decindent;
}

Ref<GeminalDescriptor> CorrelationFactor::geminaldescriptor() {
  return(geminaldescriptor_);
}

void
CorrelationFactor::print_params(std::ostream& os, unsigned int f) const
{
}

double
CorrelationFactor::value(unsigned int c, double r12, double r1, double r2) const
{
  return value(c,r12);
}

TwoBodyInt::tbint_type
CorrelationFactor::tbint_type_eri() const
{
  return TwoBodyInt::eri;
}

TwoBodyInt::tbint_type
CorrelationFactor::tbint_type_f12() const
{
  throw ProgrammingError("LinearR12::CorrelationFactor::tbint_type_f12() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

TwoBodyInt::tbint_type
CorrelationFactor::tbint_type_t1f12() const
{
  throw ProgrammingError("LinearR12::CorrelationFactor::tbint_type_t1f12() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

TwoBodyInt::tbint_type
CorrelationFactor::tbint_type_t2f12() const
{
  throw ProgrammingError("LinearR12::CorrelationFactor::tbint_type_t2f12() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

TwoBodyInt::tbint_type
CorrelationFactor::tbint_type_f12eri() const
{
  throw ProgrammingError("LinearR12::CorrelationFactor::tbint_type_f12eri() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

TwoBodyInt::tbint_type
CorrelationFactor::tbint_type_f12f12() const
{
  throw ProgrammingError("LinearR12::CorrelationFactor::tbint_type_f12f12() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

TwoBodyInt::tbint_type
CorrelationFactor::tbint_type_f12t1f12() const
{
  throw ProgrammingError("LinearR12::CorrelationFactor::tbint_type_f12t1f12() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

TwoBodyInt::tbint_type
CorrelationFactor::tbint_type_f12f12_anti() const
{
  throw ProgrammingError("LinearR12::CorrelationFactor::tbint_type_f12f12_anti() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

Ref<TwoBodyIntDescr>
CorrelationFactor::tbintdescr(const Ref<Integral>& IF, unsigned int f) const
{
  throw ProgrammingError("LinearR12::CorrelationFactor::tbintdescr(f) -- should not be called for this CorrelationFactor",__FILE__,__LINE__);
}

Ref<TwoBodyIntDescr>
CorrelationFactor::tbintdescr(const Ref<Integral>& IF, unsigned int fbra, unsigned int fket) const
{
  throw ProgrammingError("LinearR12::CorrelationFactor::tbintdescr(f,g) -- should not be called for this CorrelationFactor",__FILE__,__LINE__);
}

////

NullCorrelationFactor::NullCorrelationFactor() :
  CorrelationFactor()
{
}

double
NullCorrelationFactor::value(unsigned int c, double r12) const
{
  return 0.0;
}

bool
NullCorrelationFactor::equiv(const Ref<CorrelationFactor>& cf) const
{
  Ref<NullCorrelationFactor> cf_null; cf_null << cf;
  return cf_null.nonnull();
}

////

R12CorrelationFactor::R12CorrelationFactor()
{
  label_=std::string("R12");
  Ref<GeminalDescriptorFactory> gdesc_factory=new GeminalDescriptorFactory;
  geminaldescriptor_=gdesc_factory->r12_geminal();
}

TwoBodyInt::tbint_type
R12CorrelationFactor::tbint_type_f12() const
{
  return (TwoBodyInt::r12);
}

TwoBodyInt::tbint_type
R12CorrelationFactor::tbint_type_t1f12() const
{
  return (TwoBodyInt::r12t1);
}

TwoBodyInt::tbint_type
R12CorrelationFactor::tbint_type_t2f12() const
{
  return (TwoBodyInt::r12t2);
}

Ref<TwoBodyIntDescr>
R12CorrelationFactor::tbintdescr(const Ref<Integral>& IF, unsigned int f) const
{
  return new TwoBodyIntDescrR12(IF);
}

double
R12CorrelationFactor::value(unsigned int c, double r12) const
{
  return r12;
}

bool
R12CorrelationFactor::equiv(const Ref<CorrelationFactor>& cf) const
{
  Ref<R12CorrelationFactor> cf_cast; cf_cast << cf;
  return cf_cast.nonnull();
}

////

G12CorrelationFactor::G12CorrelationFactor(const CorrelationParameters& params, const Ref<GeminalDescriptor> &geminaldescriptor){
  label_=std::string("G12");
  if (geminaldescriptor.nonnull())
    geminaldescriptor_=geminaldescriptor;
  else {
    Ref<GeminalDescriptorFactory> gdesc_factory=new GeminalDescriptorFactory;
    geminaldescriptor_=gdesc_factory->gaussian_geminal(params);
  }
  params_=params;
}

unsigned int
G12CorrelationFactor::nfunctions() const
{
  return params_.size();
}

const G12CorrelationFactor::ContractedGeminal&
G12CorrelationFactor::function(unsigned int c) const
{
  return params_.at(c);
}

unsigned int
G12CorrelationFactor::nprimitives(unsigned int c) const
{
  return params_.at(c).size();
}

const G12CorrelationFactor::PrimitiveGeminal&
G12CorrelationFactor::primitive(unsigned int c, unsigned int p) const
{
  return params_.at(c).at(p);
}


TwoBodyInt::tbint_type
G12CorrelationFactor::tbint_type_f12() const
{
  return (TwoBodyInt::r12_0_g12);
}

TwoBodyInt::tbint_type
G12CorrelationFactor::tbint_type_t1f12() const
{
  return (TwoBodyInt::t1g12);
}

TwoBodyInt::tbint_type
G12CorrelationFactor::tbint_type_t2f12() const
{
  return (TwoBodyInt::t2g12);
}

TwoBodyInt::tbint_type
G12CorrelationFactor::tbint_type_f12eri() const
{
  return (TwoBodyInt::r12_m1_g12);
}

TwoBodyInt::tbint_type
G12CorrelationFactor::tbint_type_f12f12() const
{
  return (TwoBodyInt::r12_0_g12);
}

TwoBodyInt::tbint_type
G12CorrelationFactor::tbint_type_f12t1f12() const
{
  return (TwoBodyInt::g12t1g12);
}

Ref<TwoBodyIntDescr>
G12CorrelationFactor::tbintdescr(const Ref<Integral>& IF, unsigned int f) const
{
  Ref<IntParamsG12> params = new IntParamsG12(function(f));
  return new TwoBodyIntDescrG12(IF,params);
}

Ref<TwoBodyIntDescr>
G12CorrelationFactor::tbintdescr(const Ref<Integral>& IF,
                                            unsigned int fbra,
                                            unsigned int fket) const
{
  Ref<IntParamsG12> params = new IntParamsG12(function(fbra),function(fket));
  return new TwoBodyIntDescrG12(IF,params);
}

double
G12CorrelationFactor::value(unsigned int c, double r12) const
{
  double val = 0.0;
  const unsigned int nprims = nprimitives(c);
  for(unsigned int p=0; p<nprims; p++) {
    const PrimitiveGeminal& prim = primitive(c,p);
    const double expon = prim.first;
    const double coef = prim.second;
    val += coef*exp(-expon*r12*r12);
  }
  return val;
}

void
G12CorrelationFactor::print_params(std::ostream& os, unsigned int f) const
{
  using std::endl;

  os << indent << "[ Exponent Coefficient] = [ ";
  const int nprim = nprimitives(f);
  for(int p=0; p<nprim; p++) {
    const PrimitiveGeminal& prim = primitive(f,p);
    os << "[" << prim.first << " " << prim.second << "] ";
  }
  os << " ]" << endl;
}

bool
G12CorrelationFactor::equiv(const Ref<CorrelationFactor>& cf) const
{
  Ref<G12CorrelationFactor> cf_cast; cf_cast << cf;
  if (cf_cast.null()) return false;
  return CorrParamCompare<IntParamsG12>::equiv(params_,(*cf_cast).params_);
}

////

G12NCCorrelationFactor::G12NCCorrelationFactor(const CorrelationParameters& params, const Ref<GeminalDescriptor> &geminaldescriptor){
  label_=std::string("G12");
  if (geminaldescriptor.nonnull()){
  geminaldescriptor_=geminaldescriptor;
  }
  else {
    Ref<GeminalDescriptorFactory> gdesc_factory=new GeminalDescriptorFactory;
    geminaldescriptor_=gdesc_factory->gaussian_geminal(params);
  }
  params_=params;
}

unsigned int
G12NCCorrelationFactor::nfunctions() const
{
  return params_.size();
}

const G12NCCorrelationFactor::ContractedGeminal&
G12NCCorrelationFactor::function(unsigned int c) const
{
  return params_.at(c);
}

unsigned int
G12NCCorrelationFactor::nprimitives(unsigned int c) const
{
  return params_.at(c).size();
}

const G12NCCorrelationFactor::PrimitiveGeminal&
G12NCCorrelationFactor::primitive(unsigned int c, unsigned int p) const
{
  return params_.at(c).at(p);
}

TwoBodyInt::tbint_type
G12NCCorrelationFactor::tbint_type_f12() const
{
  return (TwoBodyInt::r12_0_g12);
}

TwoBodyInt::tbint_type
G12NCCorrelationFactor::tbint_type_f12eri() const
{
  return (TwoBodyInt::r12_m1_g12);
}

TwoBodyInt::tbint_type
G12NCCorrelationFactor::tbint_type_f12f12() const
{
  return (TwoBodyInt::r12_0_g12);
}

TwoBodyInt::tbint_type
G12NCCorrelationFactor::tbint_type_f12t1f12() const
{
  return (TwoBodyInt::g12t1g12);
}

TwoBodyInt::tbint_type
G12NCCorrelationFactor::tbint_type_f12f12_anti() const
{
  return (TwoBodyInt::anti_g12g12);
}

Ref<TwoBodyIntDescr>
G12NCCorrelationFactor::tbintdescr(const Ref<Integral>& IF, unsigned int f) const
{
  Ref<IntParamsG12> params = new IntParamsG12(function(f));
  return new TwoBodyIntDescrG12NC(IF,params);
}

Ref<TwoBodyIntDescr>
G12NCCorrelationFactor::tbintdescr(const Ref<Integral>& IF,
					      unsigned int fbra,
					      unsigned int fket) const
{
  Ref<IntParamsG12> params = new IntParamsG12(function(fbra),function(fket));
  return new TwoBodyIntDescrG12NC(IF,params);
}

double
G12NCCorrelationFactor::value(unsigned int c, double r12) const
{
  double val = 0.0;
  const unsigned int nprims = nprimitives(c);
  for(unsigned int p=0; p<nprims; p++) {
    const PrimitiveGeminal& prim = primitive(c,p);
    const double expon = prim.first;
    const double coef = prim.second;
    val += coef*exp(-expon*r12*r12);
  }
  return val;
}

void
G12NCCorrelationFactor::print_params(std::ostream& os, unsigned int f) const
{
  using std::endl;

  os << indent << "[ Exponent Coefficient] = [ ";
  const int nprim = nprimitives(f);
  for(int p=0; p<nprim; p++) {
    const PrimitiveGeminal& prim = primitive(f,p);
    os << "[" << prim.first << " " << prim.second << "] ";
  }
  os << " ]" << endl;
}

bool
G12NCCorrelationFactor::equiv(const Ref<CorrelationFactor>& cf) const
{
  Ref<G12NCCorrelationFactor> cf_cast; cf_cast << cf;
  if (cf_cast.null()) return false;
  return CorrParamCompare<IntParamsG12>::equiv(params_,(*cf_cast).params_);
}

G12NCCorrelationFactor::ContractedGeminal
G12NCCorrelationFactor::product(const ContractedGeminal& A,
                                       const ContractedGeminal& B)
{
  return IntParamsG12::product(A,B);
}

////

GenG12CorrelationFactor::GenG12CorrelationFactor(const CorrelationParameters& params) :
  CorrelationFactor(), params_(params)
{
  if (params_.size() != 1)
    throw ProgrammingError("GenG12CorrelationFactor::GenG12CorrelationFactor() -- only 1 general Geminal correlation factor can now be handled",__FILE__,__LINE__);
}

unsigned int
GenG12CorrelationFactor::nfunctions() const
{
  return params_.size();
}


const GenG12CorrelationFactor::ContractedGeminal&
GenG12CorrelationFactor::function(unsigned int c) const
{
  return params_.at(c);
}

unsigned int
GenG12CorrelationFactor::nprimitives(unsigned int c) const
{
  return params_.at(c).size();
}

const GenG12CorrelationFactor::PrimitiveGeminal&
GenG12CorrelationFactor::primitive(unsigned int c, unsigned int p) const
{
  return params_.at(c).at(p);
}

TwoBodyInt::tbint_type
GenG12CorrelationFactor::tbint_type_f12() const
{
  return (TwoBodyInt::r12_0_gg12);
}

TwoBodyInt::tbint_type
GenG12CorrelationFactor::tbint_type_f12eri() const
{
  return (TwoBodyInt::r12_m1_gg12);
}

TwoBodyInt::tbint_type
GenG12CorrelationFactor::tbint_type_f12f12() const
{
  return (TwoBodyInt::r12_0_gg12);
}

TwoBodyInt::tbint_type
GenG12CorrelationFactor::tbint_type_f12t1f12() const
{
  return (TwoBodyInt::gg12t1gg12);
}

Ref<TwoBodyIntDescr>
GenG12CorrelationFactor::tbintdescr(const Ref<Integral>& IF, unsigned int f) const
{
  Ref<IntParamsGenG12> params = new IntParamsGenG12(function(f));
  return new TwoBodyIntDescrGenG12(IF,params);
}

Ref<TwoBodyIntDescr>
GenG12CorrelationFactor::tbintdescr(const Ref<Integral>& IF,
					       unsigned int fbra,
					       unsigned int fket) const
{
  Ref<IntParamsGenG12> params = new IntParamsGenG12(function(fbra),function(fket));
  return new TwoBodyIntDescrGenG12(IF,params);
}

double
GenG12CorrelationFactor::value(unsigned int c, double r12) const
{
  throw ProgrammingError("GenG12CorrelationFactor::value(c,r12) -- not defined for general Geminal correlation factors",__FILE__,__LINE__);
}

double
GenG12CorrelationFactor::value(unsigned int c, double r12, double r1, double r2) const
{
  double val = 0.0;
  const unsigned int nprims = nprimitives(c);
  for(unsigned int p=0; p<nprims; p++) {
    const PrimitiveGeminal& prim = primitive(c,p);
    const double alpha = prim.first.first;
    const double gamma = prim.first.second;
    const double coef = prim.second;
    val += coef*exp( - alpha*(r1*r1 + r2*r2) - gamma*r12*r12 );
  }
  return val;
}

void
GenG12CorrelationFactor::print_params(std::ostream& os, unsigned int f) const
{
  using std::endl;

  os << indent << "[ [Alpha Gamma] Coefficient] = [ ";
  const int nprim = nprimitives(f);
  for(int p=0; p<nprim; p++) {
    const PrimitiveGeminal& prim = primitive(f,p);
    os << "[ [" << prim.first.first << " " << prim.first.second << "] " << prim.second << "] ";
  }
  os << " ]" << endl;
}

bool
GenG12CorrelationFactor::equiv(const Ref<CorrelationFactor>& cf) const
{
  Ref<GenG12CorrelationFactor> cf_cast; cf_cast << cf;
  if (cf_cast.null()) return false;
  return CorrParamCompare<IntParamsGenG12>::equiv(params_,(*cf_cast).params_);
}

