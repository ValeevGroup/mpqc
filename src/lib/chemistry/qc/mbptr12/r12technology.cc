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

#include <cassert>
#include <util/state/statein.h>
#include <util/state/stateout.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/orthog.h>
#include <chemistry/qc/basis/uncontract.h>
#include <chemistry/qc/basis/lselect.h>
#include <chemistry/qc/intv3/intv3.h>
#ifdef HAVE_LIBINT2
#  include <chemistry/qc/libint2/libint2.h>
#endif
#include <math/optimize/gaussianfit.h>
#include <math/optimize/gaussianfit.timpl.h>
#include <util/misc/print.h>
#include <chemistry/qc/mbptr12/r12technology.h>

using namespace std;
using namespace sc;

static ClassDesc R12Ansatz_cd(
  typeid(R12Technology::R12Ansatz),"R12Ansatz",3,"virtual public SavableState",
  create<R12Technology::R12Ansatz>, create<R12Technology::R12Ansatz>, create<R12Technology::R12Ansatz>);

R12Technology::R12Ansatz::R12Ansatz()
  : projector_(R12Technology::Projector_2),
    diag_(false),
    amplitudes_(R12Technology::GeminalAmplitudeAnsatz_fullopt),
    wof_(false),
    orbital_product_GG_(R12Technology::OrbProdGG_ij),
    orbital_product_gg_(R12Technology::OrbProdgg_ij) {}

R12Technology::R12Ansatz::R12Ansatz(const Ref<KeyVal>& keyval)
{
  projector_ = (R12Technology::Projector)keyval->intvalue("projector",KeyValValueint(2));
  const bool default_wof = (projector_ == R12Technology::Projector_0);
  wof_ = keyval->booleanvalue("wof",KeyValValueboolean((int)default_wof));

  diag_ = keyval->booleanvalue("diag",KeyValValueboolean((int)false));

  // amplitudes should be fixed by default if the diagonal ansatz is used
  const std::string default_amplitudes_str = diag_ ? std::string("fixed") : std::string("optimized");
  const std::string amplitudes_str = keyval->stringvalue("amplitudes",
                                                         KeyValValuestring(default_amplitudes_str));
  if (amplitudes_str == std::string("fixed"))
    amplitudes_ = R12Technology::GeminalAmplitudeAnsatz_fixed;
  else if (amplitudes_str == std::string("optimized"))
    amplitudes_ = R12Technology::GeminalAmplitudeAnsatz_fullopt;
  else if (amplitudes_str == std::string("scaledfixed"))
    amplitudes_ = R12Technology::GeminalAmplitudeAnsatz_scaledfixed;
  else
    throw InputError("Invalid value for keyword \"amplitudes\"",__FILE__,__LINE__);
  if ( (diag_==false) && (amplitudes_!=R12Technology::GeminalAmplitudeAnsatz_fullopt) ){
    throw InputError("R12Ansatz::R12Ansatz -- amplitudes can only be fixed if diag is true",__FILE__,__LINE__);
  }

  std::string op = keyval->stringvalue("orbital_product_GG",KeyValValuestring("ij"));
  if (op == "ij")
    orbital_product_GG_ = R12Technology::OrbProdGG_ij;
  else if (op == "pq")
    orbital_product_GG_ = R12Technology::OrbProdGG_pq;
  else
    throw InputError("R12Ansatz::R12Ansatz -- invalid value for orbital_product_GG",__FILE__,__LINE__);

  op = keyval->stringvalue("orbital_product_gg",KeyValValuestring("ij"));
  if (op == "ij")
    orbital_product_gg_ = R12Technology::OrbProdgg_ij;
  else if (op == "pq")
    orbital_product_gg_ = R12Technology::OrbProdgg_pq;
  else
    throw InputError("R12Ansatz::R12Ansatz -- invalid value for orbital_product_gg",__FILE__,__LINE__);
}

R12Technology::R12Ansatz::R12Ansatz(StateIn& s) :
  SavableState(s)
{
  int p; s.get(p); projector_ = (R12Technology::Projector)p;
  int d; s.get(d); diag_ = (bool)d;
  int a; s.get(a); amplitudes_ = (R12Technology::GeminalAmplitudeAnsatz)a;
  if (s.version(::class_desc<R12Ansatz>()) >= 3) {
    int w; s.get(w); wof_ = (bool)w;
  }
  if (s.version(::class_desc<R12Ansatz>()) >= 2) {
    int o; s.get(o); orbital_product_GG_ = (R12Technology::OrbitalProduct_GG)o;
    s.get(o); orbital_product_gg_ = (R12Technology::OrbitalProduct_gg)o;
  }
}

R12Technology::R12Ansatz::~R12Ansatz() {}

void
R12Technology::R12Ansatz::save_data_state(StateOut& s)
{
  s.put((int)projector_);
  s.put((int)diag_);
  s.put((int)amplitudes_);
  s.put((int)wof_);
  s.put((int)orbital_product_GG_);
  s.put((int)orbital_product_gg_);
}

void
R12Technology::R12Ansatz::print(std::ostream& o) const
{
  o << indent << "R12Ansatz:" << std::endl;
  o << incindent;

  o << indent << "Geminal orbital Product Space: ";
  switch(orbital_product_GG_) {
    case R12Technology::OrbProdGG_ij: o << "ij"; break;
    case R12Technology::OrbProdGG_pq: o << "pq"; break;
  }
  o << std::endl;

  o << indent << "Space of orbital products from which geminal substitutions are allowed: ";
  switch(orbital_product_gg_) {
    case R12Technology::OrbProdgg_ij: o << "ij"; break;
    case R12Technology::OrbProdgg_pq: o << "pq"; break;
  }
  o << std::endl;

  o << indent << "Projector: ";
  switch(projector_) {
    case R12Technology::Projector_0: o << "0  , i.e. 1"; break;
    case R12Technology::Projector_1: o << "1  , i.e. (1-P1)(1-P2)"; break;
    case R12Technology::Projector_2: o << "2  , i.e. (1-O1)(1-O2)-V1V2"; break;
    case R12Technology::Projector_3: o << "3  , i.e. 1-P1P2"; break;
  }
  o << std::endl;

  std::string amplitudes_str;
  switch (amplitudes_) {
    case R12Technology::GeminalAmplitudeAnsatz_fullopt:
      amplitudes_str = std::string("optimized"); break;
    case R12Technology::GeminalAmplitudeAnsatz_fixed:
      amplitudes_str = std::string("fixed"); break;
    case R12Technology::GeminalAmplitudeAnsatz_scaledfixed:
      amplitudes_str = std::string("scaled fixed"); break;
  }
  o << indent << "Ansatz: " << (diag_ ? "diagonal" : "orbital-invariant")
              << " with " << amplitudes_str << " amplitudes" << std::endl;
  o << indent << "WOF: " << (wof_ ? "true" : "false") << std::endl;
  o << decindent;
}

R12Technology::Projector
R12Technology::R12Ansatz::projector() const { return projector_; }

bool
R12Technology::R12Ansatz::diag() const { return diag_; }

R12Technology::GeminalAmplitudeAnsatz
R12Technology::R12Ansatz::amplitudes() const { return amplitudes_; }

bool
R12Technology::R12Ansatz::wof() const { return wof_; }

R12Technology::OrbitalProduct_GG
R12Technology::R12Ansatz::orbital_product_GG() const { return orbital_product_GG_; }

R12Technology::OrbitalProduct_gg
R12Technology::R12Ansatz::orbital_product_gg() const {
  return(orbital_product_gg_);
}

/////////////////////////////////////

namespace sc {

  template<>
  bool R12Technology::CorrParamCompare<IntParamsG12>::equiv(
                                                            const PrimitiveGeminal& A,
                                                            const PrimitiveGeminal& B) {
    return (std::fabs(A.first - B.first) < epsilon && std::fabs(A.second
        - B.second) < epsilon);
  }

}

/*********************
 * GeminalDescriptor *
 *********************/
R12Technology::GeminalDescriptor::GeminalDescriptor(){
  type_ = "invalid";
}

R12Technology::GeminalDescriptor::GeminalDescriptor(const std::string& type, const std::vector<std::string> &params){
  type_ = type;
  params_ = params;
  //compute_offsets();
}

R12Technology::GeminalDescriptor::GeminalDescriptor(const GeminalDescriptor& source){
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

std::string R12Technology::GeminalDescriptor::type() const {
  return(type_);
}

std::vector<std::string> R12Technology::GeminalDescriptor::params() const {
  return(params_);
}

void R12Technology::GeminalDescriptor::print(std::ostream &o) {
  o << indent << "GeminalDescriptor :" << std::endl;
  o << incindent << indent << "type = " << type_ << std::endl;
  o << indent << "params = [";
  for(int i=0; i<params_.size(); i++) {
    o << " " << params_[i];
  }
  o << "]" << std::endl << decindent;
}

bool R12Technology::invalid(const Ref<GeminalDescriptor>& gdesc){
  std::string type=gdesc->type();
  return((type==std::string("invalid")) ? true : false);
}

bool R12Technology::R12(const Ref<GeminalDescriptor>& gdesc){
  std::string type=gdesc->type();
  return(((type==std::string("R12")) || (type==std::string("r12"))) ? true : false);
}

bool R12Technology::STG(const Ref<GeminalDescriptor>& gdesc){
  std::string type=gdesc->type();
  return(((type==std::string("STG")) || (type==std::string("stg"))) ? true : false);
}

bool R12Technology::G12(const Ref<GeminalDescriptor>& gdesc){
  std::string type=gdesc->type();
  if((type==std::string("G12")) || (type==std::string("g12"))){
    return(true);
  }
  else {
    return(false);
  }
}

double R12Technology::single_slater_exponent(const Ref<GeminalDescriptor>& gdesc) {
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

  return(exponent);
}

/****************************
 * GeminalDescriptorFactory *
 ****************************/
R12Technology::GeminalDescriptorFactory::GeminalDescriptorFactory()
  : invalid_id_("invalid"),
    r12_id_("R12"),
    stg_id_("STG"),
    g12_id_("G12") {}

Ref<R12Technology::GeminalDescriptor>
R12Technology::GeminalDescriptorFactory::null_geminal(){
  std::vector<std::string> void_vector;
  return(Ref<GeminalDescriptor>(new GeminalDescriptor(std::string(invalid_id_),void_vector)));
}

Ref<R12Technology::GeminalDescriptor>
R12Technology::GeminalDescriptorFactory::r12_geminal(){
  std::vector<std::string> void_vector;
  return(Ref<GeminalDescriptor>(new GeminalDescriptor(std::string(r12_id_),void_vector)));
}

Ref<R12Technology::GeminalDescriptor>
R12Technology::GeminalDescriptorFactory::slater_geminal(double gamma){
  std::vector<std::string> param_vec(3);
  int nfunction=1;
  int nprimitive=1;
  std::stringstream inout;
  inout << nfunction;
  inout >> param_vec[0];
  inout << nprimitive;
  inout >> param_vec[1];
  inout << std::setprecision(16) << gamma;
  inout >> param_vec[2];
  return(Ref<GeminalDescriptor>(new GeminalDescriptor(std::string(stg_id_),param_vec)));
}

Ref<R12Technology::GeminalDescriptor>
R12Technology::GeminalDescriptorFactory::slater_geminal(const std::vector<double> &gamma) {
  const int nfunction=gamma.size();
  std::vector<std::string> params(1+3*nfunction);
  std::stringstream inout;
  inout << nfunction;
  inout >> params[0];

  for(int i=0; i<nfunction; i++){
    params[i+1] = "1";
    params[nfunction+1+2*i] = "1.0";
    std::ostringstream oss;
    oss << std::setprecision(16) << gamma[i];
    params[nfunction+2+2*i] = oss.str();
  }

  return(Ref<GeminalDescriptor>(new GeminalDescriptor(std::string(stg_id_),params)));
}

Ref<R12Technology::GeminalDescriptor>
R12Technology::GeminalDescriptorFactory::gaussian_geminal(double gamma){
  std::vector<std::string> param_vec(3);
  int nfunction=1;
  int nprimitive=1;
  std::stringstream inout;
  inout << nfunction;
  inout >> param_vec[0];
  inout << nprimitive;
  inout >> param_vec[1];
  inout << std::setprecision(16) << gamma;
  inout >> param_vec[2];
  return(Ref<GeminalDescriptor>(new GeminalDescriptor(std::string(g12_id_),param_vec)));
}

Ref<R12Technology::GeminalDescriptor>
R12Technology::GeminalDescriptorFactory::contracted_gaussian_geminal(const std::vector<double> &coeff,
                                                                     const std::vector<double> &gamma) {
  int nfunction=1;
  unsigned int ngeminal=coeff.size();
  std::vector<std::string> param_vec(2*ngeminal+2);
  std::stringstream inout;
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

Ref<R12Technology::GeminalDescriptor>
R12Technology::GeminalDescriptorFactory::gaussian_geminal(const G12CorrelationFactor::CorrelationParameters &corrparams) {
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
  std::stringstream inout;
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

R12Technology::CorrelationFactor::CorrelationFactor(const std::string& label,
                                                    const Ref<GeminalDescriptor> &geminaldescriptor) :
  label_(label),
  geminaldescriptor_(geminaldescriptor)
{
}

R12Technology::CorrelationFactor::~CorrelationFactor()
{
}

R12Technology::CorrelationFactor::CorrelationFactor()
{
  label_=std::string("invalid");
  Ref<GeminalDescriptorFactory> gdesc_factory=new GeminalDescriptorFactory;
  geminaldescriptor_=gdesc_factory->null_geminal();
}

unsigned int
R12Technology::CorrelationFactor::nfunctions() const
{
  return 1;
}

unsigned int
R12Technology::CorrelationFactor::nprimitives(unsigned int c) const
{
  return 1;
}

const std::string&
R12Technology::CorrelationFactor::label() const
{
  return label_;
}

void
R12Technology::CorrelationFactor::print(std::ostream& os) const
{
  using std::endl;
  os << indent << "CorrelationFactor:";
  if (invalid(geminaldescriptor_)) {
    os << " none" << endl;
  }
  else if (R12(geminaldescriptor_)) {
    os << " R12" << endl;
  }
  else if (STG(geminaldescriptor_)) {
    os << " STG-" << this->nprimitives(0) << "G[" << geminaldescriptor_->params()[3] << "]"<< endl << incindent;
    print_params(os,0);
    os << decindent;
  }
  else {
    os << endl << incindent;
    const int nfunc = nfunctions();
    for(int f=0; f<nfunc; f++) {
      os << indent << "Function " << f << ":" << endl << incindent;
      os << indent << "Functional form: " << label() << endl;
      print_params(os,f);
      os << decindent;
    }
    os << decindent;
  }
}

Ref<R12Technology::GeminalDescriptor> R12Technology::CorrelationFactor::geminaldescriptor() {
  return(geminaldescriptor_);
}

void
R12Technology::CorrelationFactor::print_params(std::ostream& os, unsigned int f) const
{
}

TwoBodyOper::type
R12Technology::CorrelationFactor::tbint_type_eri() const
{
  return TwoBodyOper::eri;
}

TwoBodyOper::type
R12Technology::CorrelationFactor::tbint_type_f12() const
{
  throw ProgrammingError("CorrelationFactor::tbint_type_f12() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

TwoBodyOper::type
R12Technology::CorrelationFactor::tbint_type_t1f12() const
{
  throw ProgrammingError("CorrelationFactor::tbint_type_t1f12() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

TwoBodyOper::type
R12Technology::CorrelationFactor::tbint_type_t2f12() const
{
  throw ProgrammingError("CorrelationFactor::tbint_type_t2f12() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

TwoBodyOper::type
R12Technology::CorrelationFactor::tbint_type_f12eri() const
{
  throw ProgrammingError("CorrelationFactor::tbint_type_f12eri() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

TwoBodyOper::type
R12Technology::CorrelationFactor::tbint_type_f12f12() const
{
  throw ProgrammingError("CorrelationFactor::tbint_type_f12f12() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

TwoBodyOper::type
R12Technology::CorrelationFactor::tbint_type_f12t1f12() const
{
  throw ProgrammingError("CorrelationFactor::tbint_type_f12t1f12() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

TwoBodyOper::type
R12Technology::CorrelationFactor::tbint_type_f12f12_anti() const
{
  throw ProgrammingError("CorrelationFactor::tbint_type_f12f12_anti() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

Ref<TwoBodyIntDescr>
R12Technology::CorrelationFactor::tbintdescr(const Ref<Integral>& IF, unsigned int f) const
{
  throw ProgrammingError("CorrelationFactor::tbintdescr(f) -- should not be called for this CorrelationFactor",__FILE__,__LINE__);
}

Ref<TwoBodyIntDescr>
R12Technology::CorrelationFactor::tbintdescr(const Ref<Integral>& IF, unsigned int fbra, unsigned int fket) const
{
  throw ProgrammingError("CorrelationFactor::tbintdescr(f,g) -- should not be called for this CorrelationFactor",__FILE__,__LINE__);
}

////

R12Technology::NullCorrelationFactor::NullCorrelationFactor() :
  CorrelationFactor()
{
}

double
R12Technology::NullCorrelationFactor::value(unsigned int c, double r12) const
{
  return 0.0;
}

bool
R12Technology::NullCorrelationFactor::equiv(const Ref<CorrelationFactor>& cf) const
{
  Ref<NullCorrelationFactor> cf_null; cf_null << cf;
  return cf_null;
}

Ref<TwoBodyIntDescr>
R12Technology::NullCorrelationFactor::tbintdescr(const Ref<Integral>& IF, unsigned int f) const
{
  return new TwoBodyIntDescrERI(IF);
}

Ref<TwoBodyIntDescr>
R12Technology::NullCorrelationFactor::tbintdescr(const Ref<Integral>& IF, unsigned int f, unsigned int g) const
{
  return new TwoBodyIntDescrERI(IF);
}

////

R12Technology::R12CorrelationFactor::R12CorrelationFactor()
{
  label_=std::string("R12");
  Ref<GeminalDescriptorFactory> gdesc_factory=new GeminalDescriptorFactory;
  geminaldescriptor_=gdesc_factory->r12_geminal();
}

TwoBodyOper::type
R12Technology::R12CorrelationFactor::tbint_type_f12() const
{
  return (TwoBodyOper::r12);
}

TwoBodyOper::type
R12Technology::R12CorrelationFactor::tbint_type_t1f12() const
{
  return (TwoBodyOper::r12t1);
}

TwoBodyOper::type
R12Technology::R12CorrelationFactor::tbint_type_t2f12() const
{
  return (TwoBodyOper::r12t2);
}

Ref<TwoBodyIntDescr>
R12Technology::R12CorrelationFactor::tbintdescr(const Ref<Integral>& IF, unsigned int f) const
{
  return new TwoBodyIntDescrR12(IF);
}

Ref<TwoBodyIntDescr>
R12Technology::R12CorrelationFactor::tbintdescr(const Ref<Integral>& IF, unsigned int f, unsigned int g) const
{
  return new TwoBodyIntDescrR12(IF);
}

double
R12Technology::R12CorrelationFactor::value(unsigned int c, double r12) const
{
  return r12;
}

bool
R12Technology::R12CorrelationFactor::equiv(const Ref<CorrelationFactor>& cf) const
{
  Ref<R12CorrelationFactor> cf_cast; cf_cast << cf;
  return cf_cast;
}

////

R12Technology::G12CorrelationFactor::G12CorrelationFactor(const CorrelationParameters& params, const Ref<GeminalDescriptor> &geminaldescriptor){
  label_=std::string("G12");
  if (geminaldescriptor)
    geminaldescriptor_=geminaldescriptor;
  else {
    Ref<GeminalDescriptorFactory> gdesc_factory=new GeminalDescriptorFactory;
    geminaldescriptor_=gdesc_factory->gaussian_geminal(params);
  }
  params_=params;
}

unsigned int
R12Technology::G12CorrelationFactor::nfunctions() const
{
  return params_.size();
}

const R12Technology::G12CorrelationFactor::ContractedGeminal&
R12Technology::G12CorrelationFactor::function(unsigned int c) const
{
  return params_.at(c);
}

unsigned int
R12Technology::G12CorrelationFactor::nprimitives(unsigned int c) const
{
  return params_.at(c).size();
}

const R12Technology::G12CorrelationFactor::PrimitiveGeminal&
R12Technology::G12CorrelationFactor::primitive(unsigned int c, unsigned int p) const
{
  return params_.at(c).at(p);
}


TwoBodyOper::type
R12Technology::G12CorrelationFactor::tbint_type_f12() const
{
  return (TwoBodyOper::r12_0_g12);
}

TwoBodyOper::type
R12Technology::G12CorrelationFactor::tbint_type_t1f12() const
{
  return (TwoBodyOper::t1g12);
}

TwoBodyOper::type
R12Technology::G12CorrelationFactor::tbint_type_t2f12() const
{
  return (TwoBodyOper::t2g12);
}

TwoBodyOper::type
R12Technology::G12CorrelationFactor::tbint_type_f12eri() const
{
  return (TwoBodyOper::r12_m1_g12);
}

TwoBodyOper::type
R12Technology::G12CorrelationFactor::tbint_type_f12f12() const
{
  return (TwoBodyOper::r12_0_g12);
}

TwoBodyOper::type
R12Technology::G12CorrelationFactor::tbint_type_f12t1f12() const
{
  return (TwoBodyOper::g12t1g12);
}

Ref<TwoBodyIntDescr>
R12Technology::G12CorrelationFactor::tbintdescr(const Ref<Integral>& IF, unsigned int f) const
{
  Ref<IntParamsG12> params = new IntParamsG12(function(f));
  return new TwoBodyIntDescrG12(IF,params);
}

Ref<TwoBodyIntDescr>
R12Technology::G12CorrelationFactor::tbintdescr(const Ref<Integral>& IF,
                                            unsigned int fbra,
                                            unsigned int fket) const
{
  Ref<IntParamsG12> params = new IntParamsG12(function(fbra),function(fket));
  return new TwoBodyIntDescrG12(IF,params);
}

double
R12Technology::G12CorrelationFactor::value(unsigned int c, double r12) const
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
R12Technology::G12CorrelationFactor::print_params(std::ostream& os, unsigned int f) const
{
  using std::endl;

  os << indent << "[ Exponent Coefficient] = [ ";
  const int nprim = nprimitives(f);
  for(int p=0; p<nprim; p++) {
    const PrimitiveGeminal& prim = primitive(f,p);
    os << scprintf("[ %18.12lf %18.12lf ] ", prim.first, prim.second);
  }
  os << " ]" << endl;
}

bool
R12Technology::G12CorrelationFactor::equiv(const Ref<CorrelationFactor>& cf) const
{
  Ref<G12CorrelationFactor> cf_cast; cf_cast << cf;
  if (cf_cast.null()) return false;
  return R12Technology::CorrParamCompare<IntParamsG12>::equiv(params_,(*cf_cast).params_);
}

////

R12Technology::G12NCCorrelationFactor::G12NCCorrelationFactor(const CorrelationParameters& params, const Ref<GeminalDescriptor> &geminaldescriptor){
  label_=std::string("G12");
  if (geminaldescriptor){
  geminaldescriptor_=geminaldescriptor;
  }
  else {
    Ref<GeminalDescriptorFactory> gdesc_factory=new GeminalDescriptorFactory;
    geminaldescriptor_=gdesc_factory->gaussian_geminal(params);
  }
  params_=params;
}

unsigned int
R12Technology::G12NCCorrelationFactor::nfunctions() const
{
  return params_.size();
}

const R12Technology::G12NCCorrelationFactor::ContractedGeminal&
R12Technology::G12NCCorrelationFactor::function(unsigned int c) const
{
  return params_.at(c);
}

unsigned int
R12Technology::G12NCCorrelationFactor::nprimitives(unsigned int c) const
{
  return params_.at(c).size();
}

const R12Technology::G12NCCorrelationFactor::PrimitiveGeminal&
R12Technology::G12NCCorrelationFactor::primitive(unsigned int c, unsigned int p) const
{
  return params_.at(c).at(p);
}

TwoBodyOper::type
R12Technology::G12NCCorrelationFactor::tbint_type_f12() const
{
  return (TwoBodyOper::r12_0_g12);
}

TwoBodyOper::type
R12Technology::G12NCCorrelationFactor::tbint_type_f12eri() const
{
  return (TwoBodyOper::r12_m1_g12);
}

TwoBodyOper::type
R12Technology::G12NCCorrelationFactor::tbint_type_f12f12() const
{
  return (TwoBodyOper::r12_0_g12);
}

TwoBodyOper::type
R12Technology::G12NCCorrelationFactor::tbint_type_f12t1f12() const
{
  return (TwoBodyOper::g12t1g12);
}

TwoBodyOper::type
R12Technology::G12NCCorrelationFactor::tbint_type_f12f12_anti() const
{
  return (TwoBodyOper::anti_g12g12);
}

Ref<TwoBodyIntDescr>
R12Technology::G12NCCorrelationFactor::tbintdescr(const Ref<Integral>& IF, unsigned int f) const
{
  Ref<IntParamsG12> params = new IntParamsG12(function(f));
  return new TwoBodyIntDescrG12NC(IF,params);
}

Ref<TwoBodyIntDescr>
R12Technology::G12NCCorrelationFactor::tbintdescr(const Ref<Integral>& IF,
                          unsigned int fbra,
                          unsigned int fket) const
{
  Ref<IntParamsG12> params = new IntParamsG12(function(fbra),function(fket));
  return new TwoBodyIntDescrG12NC(IF,params);
}

double
R12Technology::G12NCCorrelationFactor::value(unsigned int c, double r12) const
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
R12Technology::G12NCCorrelationFactor::print_params(std::ostream& os, unsigned int f) const
{
  using std::endl;

  os << indent << "[ Exponent Coefficient] = [ ";
  const int nprim = nprimitives(f);
  for(int p=0; p<nprim; p++) {
    const PrimitiveGeminal& prim = primitive(f,p);
    os << scprintf("[ %18.12lf %18.12lf ] ", prim.first, prim.second);
  }
  os << " ]" << endl;
}

bool
R12Technology::G12NCCorrelationFactor::equiv(const Ref<CorrelationFactor>& cf) const
{
  Ref<G12NCCorrelationFactor> cf_cast; cf_cast << cf;
  if (cf_cast.null()) return false;
  return R12Technology::CorrParamCompare<IntParamsG12>::equiv(params_,(*cf_cast).params_);
}

R12Technology::G12NCCorrelationFactor::ContractedGeminal
R12Technology::G12NCCorrelationFactor::product(const ContractedGeminal& A,
                                       const ContractedGeminal& B)
{
  return IntParamsG12::product(A,B);
}

///////////////////////////////////
static ClassDesc R12Technology_cd(
  typeid(R12Technology),"R12Technology",11,"virtual public SavableState",
  0, create<R12Technology>, create<R12Technology>);

R12Technology::R12Technology(StateIn& s)
{
  int vbs_eq_obs; s.get(vbs_eq_obs); vbs_eq_obs_ = (bool)vbs_eq_obs;
  int abs_eq_obs; s.get(abs_eq_obs); abs_eq_obs_ = (bool)abs_eq_obs;

  int gbc; s.get(gbc); gbc_ = (bool)gbc;
  int ebc; s.get(ebc); ebc_ = (bool)ebc;
  int coupling; s.get(coupling); coupling_ = (bool)coupling;
  int compute_1rdm; s.get(compute_1rdm); compute_1rdm_ = (bool)compute_1rdm;
  int coupling_1rdm_f12b; s.get(coupling_1rdm_f12b); coupling_1rdm_f12b_ = (bool)coupling_1rdm_f12b;
  int omit_P; s.get(omit_P); omit_P_ = (bool)omit_P;
  int absmethod; s.get(absmethod); abs_method_ = (ABSMethod)absmethod;
  int stdapprox; s.get(stdapprox); stdapprox_ = (StandardApproximation) stdapprox;
  ansatz_ << SavableState::restore_state(s);
  s.get(abs_nlindep_);
  s.get(abs_lindep_tol_);
  s.get(maxnabs_);

  int safety_check; s.get(safety_check);
  safety_check_ = static_cast<bool>(safety_check);
  int posdef_B; s.get(posdef_B);
  posdef_B_ = static_cast<PositiveDefiniteB>(posdef_B);

  if (s.version(::class_desc<R12Technology>()) >= 9) {
    int omit_B; s.get(omit_B); omit_B_ = (bool)omit_B;
  }
  if (s.version(::class_desc<R12Technology>()) >= 10) {
    {
      int i; s.get(i); H0_dk_approx_pauli_ = static_cast<H0_dk_approx_pauli>(i);
    }
    {
      int i; s.get(i); H0_dk_keep_ = static_cast<bool>(i);
    }
  }
}

R12Technology::R12Technology(const Ref<KeyVal>& keyval)
{
  Ref<GaussianBasisSet> obs = require_dynamic_cast<GaussianBasisSet*>(
      keyval->describedclassvalue("basis").pointer(),
      "R12Technology::R12Technology\n"
      );
  Ref<GaussianBasisSet> abs = require_dynamic_cast<GaussianBasisSet*>(
      keyval->describedclassvalue("aux_basis").pointer(),
      "R12Technology::R12Technology\n"
      );
  Ref<GaussianBasisSet> vbs = require_dynamic_cast<GaussianBasisSet*>(
      keyval->describedclassvalue("vir_basis").pointer(),
      "R12Technology::R12Technology\n"
      );

  bool abs_eq_obs = true;
  bool vbs_eq_obs = true;
  if (obs) {
    abs_eq_obs = abs ? abs->equiv(obs) : true;
    vbs_eq_obs = vbs ? vbs->equiv(obs) : true;
  }
  this->init_from_kv(keyval, abs_eq_obs, vbs_eq_obs);
}

R12Technology::R12Technology(const Ref<KeyVal>& keyval,
                             const Ref<GaussianBasisSet>& obs,
                             const Ref<GaussianBasisSet>& vbs,
                             const Ref<GaussianBasisSet>& abs)
{
  this->init_from_kv(keyval,
                     abs ? abs->equiv(obs) : true,
                     vbs ? vbs->equiv(obs) : true);
}

void
R12Technology::init_from_kv(const Ref<KeyVal>& keyval,
                            bool abs_eq_obs,
                            bool vbs_eq_obs)
{
  abs_eq_obs_ = abs_eq_obs;
  vbs_eq_obs_ = vbs_eq_obs;

  // Default is to use the R12 factor
  std::string corrfactor = keyval->stringvalue("corr_factor", KeyValValuestring("stg-6g"));

  // Default method is MBPT2-R12/A'
  std::string sa_string = keyval->stringvalue("stdapprox",KeyValValuestring("C"));
  if ( sa_string == "A" ||
       sa_string == "a" ) {
    throw FeatureNotImplemented("stdapprox=A is obsolete",__FILE__,__LINE__);
  }
  else if ( sa_string == "Ap" ||
	    sa_string == "ap" ||
	    sa_string == "A'" ||
	    sa_string == "a'" ) {
    stdapprox_ = StdApprox_Ap;
  }
  else if ( sa_string == "App" ||
	    sa_string == "app" ||
	    sa_string == "A''" ||
	    sa_string == "a''" ) {
    stdapprox_ = StdApprox_App;
  }
  else if ( sa_string == "B" ||
	    sa_string == "b" ) {
    stdapprox_ = StdApprox_B;
  }
  else if ( sa_string == "C" ||
	    sa_string == "c" ) {
    stdapprox_ = StdApprox_C;
  }
  else if ( sa_string == "Cp" ||
            sa_string == "cp" ||
            sa_string == "C'" ||
            sa_string == "c'") {
    stdapprox_ = StdApprox_Cp;
  }
  else {
    throw std::runtime_error("R12Technology::R12Technology() -- unrecognized value for stdapprox");
  }

  // if no explicit correlation then set stdapprox to A'
  if (sa_string == "none" || sa_string == "NONE") {
      stdapprox_ = StdApprox_Ap;
  }

  //
  // r12 correlation factor?
  //
  if (corrfactor == "r12" ||
      corrfactor == "R12") {
    corrfactor_ = new R12CorrelationFactor();
    throw FeatureNotImplemented("Support for corr_factor=r12 not yet implemented in IntegralLibint2 object.",
                                __FILE__, __LINE__);
  }
  //
  // g12 correlation factor?
  //
  else if (corrfactor == "g12" || corrfactor == "G12") {
    if (keyval->exists("corr_param")) {
      typedef G12CorrelationFactor::CorrelationParameters CorrParams;
      CorrParams params;
      const int num_f12 = keyval->count("corr_param");
      if (num_f12 != 0) {
        // Do I have contracted functions?
        bool contracted = (keyval->count("corr_param",0) != 0);
        if (!contracted) {
          // Primitive functions only
          for(int f=0; f<num_f12; f++) {
            double exponent = keyval->doublevalue("corr_param", f);
            G12CorrelationFactor::ContractedGeminal vtmp;
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
            G12CorrelationFactor::ContractedGeminal vtmp;
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
      // If stdapprox_ == A', A'', or B, need commutators
      if (stdapprox_ == StdApprox_Ap ||
          stdapprox_ == StdApprox_App ||
          stdapprox_ == StdApprox_B)
        corrfactor_ = new G12CorrelationFactor(params);
      else
        corrfactor_ = new G12NCCorrelationFactor(params);
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
	// erase STG-
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
    
    const double corr_param_scale = keyval->doublevalue("corr_param_scale", KeyValValuedouble(1.0));

    std::vector<double> stg_exponents;
    typedef G12CorrelationFactor::CorrelationParameters CorrParams;
    CorrParams params;
    int num_f12 = keyval->count("corr_param");
    Ref<GeminalDescriptorFactory> gdesc_factory=new GeminalDescriptorFactory;
    if (num_f12 != 0) {
        // Do I have contracted functions? Can't handle these (yet?)
        const bool contracted = (keyval->count("corr_param",0) != 0);
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
      const double exponent = keyval->doublevalue("corr_param");
      stg_exponents.push_back(exponent);
    }

    Ref<GeminalDescriptor> gdesc = gdesc_factory->slater_geminal(stg_exponents);
    // convert STGs into combinations of Gaussians
    for(int f=0; f<num_f12; f++) {
      using namespace sc::math;
      PowerExponential1D* w;
      // Default is to use TewKlopper fit
      const std::string gtg_fit_weight = keyval->stringvalue("gtg_fit_weight",KeyValValuestring(std::string("tewklopper")));
      if (gtg_fit_weight == std::string("TewKlopper") ||
	    gtg_fit_weight == std::string("TEWKLOPPER") ||
	    gtg_fit_weight == std::string("tewklopper") ) {
        // fit using Tew&Klopper's recipe: weight is r^2 exp(-2*r^2), which has maximum near r=0.75 and decays sharply near r=0 and r=1.5
        w = new PowerExponential1D(2.0,2,2);
      }
      else if (gtg_fit_weight == std::string("Cusp") ||
	    gtg_fit_weight == std::string("CUSP") ||
        gtg_fit_weight == std::string("cusp")) {
        // fit using weight exp(-0.005*r^6), which is flat to r=1, then falls slowly till r=2, then quickly decays to r=3
        w = new PowerExponential1D(0.005,6,0);
      }
      else {
        throw InputError("Invalid value for keyword gtg_fit_weight",__FILE__,__LINE__);
      }
      typedef GaussianFit<Slater1D,PowerExponential1D> GTGFit;
      GTGFit gtgfit(ng12, *w, 0.0, 10.0, 1001);
      // fit r12^k exp(-gamma*r_{12})
      const int k = 0;
      MPQC_ASSERT(k < 2 && k >= 0);
      const double gamma = stg_exponents[f];
      double scale = corr_param_scale;
      if (k == 0) // fit - e^{-\gamma r12} / \gamma
        scale *= -1.0/gamma;
      if (k == 1) // fit r12 e^{-\gamma r12}
        scale *= 1.0;
      Slater1D stg(gamma,k,scale);
      params.push_back( gtgfit(stg) );
      delete w;
    }

    // If stdapprox_ == A', A'', or B, need commutators
    if (stdapprox_ == StdApprox_Ap ||
        stdapprox_ == StdApprox_App ||
        stdapprox_ == StdApprox_B)
      corrfactor_ = new G12CorrelationFactor(params,gdesc);
    else
      corrfactor_ = new G12NCCorrelationFactor(params,gdesc);
  }
  //
  // no explicit correlation
  //
  else if (corrfactor == "none" || corrfactor == "NONE") {
    corrfactor_ = new NullCorrelationFactor();
  }
  else
    throw FeatureNotImplemented("R12Technology::R12Technology -- this correlation factor is not implemented",__FILE__,__LINE__);

  // Default is to assume GBC
  gbc_ = keyval->booleanvalue("gbc",KeyValValueboolean((int)true));
  // Default is to assume EBC
  ebc_ = keyval->booleanvalue("ebc",KeyValValueboolean((int)true));
  // Default is to not include coupling
  coupling_ = keyval->booleanvalue("coupling",KeyValValueboolean((int)false));
  // Default is to not computing MP2-R12 one-electron density
  compute_1rdm_ = keyval->booleanvalue("compute_1rdm",KeyValValueboolean((int)false));
  // Default is to not computing coupling contri. of ccsd-f12b 1e density
  coupling_1rdm_f12b_ = keyval->booleanvalue("coupling_1rdm_f12b",KeyValValueboolean((int)false));

  // Default is to include P in intermediate B
  omit_P_ = keyval->booleanvalue("omit_P",KeyValValueboolean((int)false));

  std::string abs_method_str = keyval->stringvalue("abs_method",KeyValValuestring("CABS+"));
  if ( abs_method_str == "CABS" ||
	   abs_method_str == "cabs" ) {
    abs_method_ = ABS_CABS;
  }
  else if ( abs_method_str == "CABS+" ||
      abs_method_str == "cabs+" ) {
    abs_method_ = ABS_CABSPlus;
  }
  else {
    throw std::runtime_error("R12Technology::R12Technology -- unrecognized value for abs_method");
  }

  if (!abs_eq_obs_) {  // how to get rid of linear dependencies in OBS
    abs_nlindep_ = -1;
    abs_lindep_tol_ = OverlapOrthog::default_lindep_tol();
    if (keyval->exists("abs_lindep_tol")) {
      abs_lindep_tol_ = keyval->doublevalue("abs_lindep_tol", KeyValValuedouble(OverlapOrthog::default_lindep_tol()));
    }
    else {
      abs_nlindep_ = keyval->intvalue("abs_nlindep", KeyValValueint(-1));
    }
  }

  ansatz_ = require_dynamic_cast<R12Ansatz*>(
    keyval->describedclassvalue("ansatz").pointer(),
    "R12Technology::R12Technology\n"
    );
  // Default constructor for R12Ansatz specifies the default
  if (ansatz_.null())
    ansatz_ = new R12Ansatz;
  if (ansatz()->projector() == Projector_1 && stdapprox() != StdApprox_C)
    throw InputError("R12Technology::R12Technology -- projector 1 has not been implemented yet for a standard approximation other than C",__FILE__,__LINE__);
  if (ansatz()->projector() == Projector_3)
    throw InputError("R12Technology::R12Technology -- projector 3 is obsolete",__FILE__,__LINE__);

  // Default is to include all integrals, unless using A'' method
  int default_maxnabs = (stdapprox_ == StdApprox_App) ? 1 : 2;
  // there are no ABS indices if OBS and ABS are the same
  if (abs_eq_obs_) default_maxnabs = 0;
  maxnabs_ = static_cast<unsigned int>(keyval->intvalue("maxnabs",KeyValValueint(default_maxnabs)));

  safety_check_ = keyval->booleanvalue("safety_check",KeyValValueboolean((int)true));

  std::string posdef_B = keyval->stringvalue("posdef_B",KeyValValuestring("weak"));
  if (posdef_B == "no" || posdef_B == "NO" || posdef_B == "false" || posdef_B == "FALSE") {
      posdef_B_ = PositiveDefiniteB_no;
  }
  else if (posdef_B == "yes" || posdef_B == "YES" || posdef_B == "true" || posdef_B == "TRUE") {
      posdef_B_ = PositiveDefiniteB_yes;
  }
  else if (posdef_B == "weak" || posdef_B == "WEAK") {
      posdef_B_ = PositiveDefiniteB_weak;
  }
  else {
      throw InputError("R12Technology::R12Technology -- invalid value for keyword posdef_B",__FILE__,__LINE__);
  }

  const std::string H0_dk_approx_pauli = keyval->stringvalue("H0_dk_approx_pauli", KeyValValuestring("false"));
  if (H0_dk_approx_pauli == "yes" ||
      H0_dk_approx_pauli == "YES" ||
      H0_dk_approx_pauli == "true" ||
      H0_dk_approx_pauli == "TRUE") {
    H0_dk_approx_pauli_ = H0_dk_approx_pauli_true;
  }
  else if (H0_dk_approx_pauli == "no" ||
      H0_dk_approx_pauli == "NO" ||
      H0_dk_approx_pauli == "false" ||
      H0_dk_approx_pauli == "FALSE") {
    H0_dk_approx_pauli_ = H0_dk_approx_pauli_false;
  }
  else if (H0_dk_approx_pauli == "fHf" ||
      H0_dk_approx_pauli == "fhf" ||
      H0_dk_approx_pauli == "FHF") {
    H0_dk_approx_pauli_ = H0_dk_approx_pauli_fHf;
  }
  else if (H0_dk_approx_pauli == "fHf_Q" ||
      H0_dk_approx_pauli == "fhf_q" ||
      H0_dk_approx_pauli == "FHF_Q") {
    H0_dk_approx_pauli_ = H0_dk_approx_pauli_fHf_Q;
  }
  else
    throw InputError("R12Technology::R12Technology -- invalid value for keyword H0_dk_approx_pauli",__FILE__,__LINE__);

  if (H0_dk_approx_pauli_ == H0_dk_approx_pauli_false) {
    H0_dk_keep_ = keyval->booleanvalue("H0_dk_keep",KeyValValueboolean((int)false));
  }

  //
  //
  // Check that requested features are compatible/allowed
  //
  //

  //
  // Relativistic features are only implemented for certain approximations
  //
  if ((H0_dk_approx_pauli_ == H0_dk_approx_pauli_fHf ||
       H0_dk_approx_pauli_ == H0_dk_approx_pauli_fHf_Q) &&
       stdapprox_ != StdApprox_C) {
    throw InputError("R12Technology::R12Technology -- the given value of keyword H0_dk_approx_pauli is only valid when stdapprox=C",__FILE__,__LINE__);
  }

  //
  // These are for debugging only
  //
  // Do not compute expensive parts of B?
  omit_B_ = keyval->booleanvalue("omit_B",KeyValValueboolean((int)false));
}

R12Technology::~R12Technology()
{
}

void
R12Technology::save_data_state(StateOut& s)
{
  s.put((int)gbc_);
  s.put((int)ebc_);
  s.put((int)coupling_);
  s.put((int)compute_1rdm_);
  s.put((int)omit_P_);
  s.put((int)abs_method_);
  s.put((int)stdapprox_);
  SavableState::save_state(ansatz_.pointer(),s);
  s.put(abs_nlindep_);
  s.put(abs_lindep_tol_);
  s.put(maxnabs_);
  s.put((int)safety_check_);
  s.put((int)posdef_B_);
  s.put((int)omit_B_);
  s.put((int)H0_dk_approx_pauli_);
  s.put((int)H0_dk_keep_);
  s.put((int)omit_B_);
}

void
R12Technology::print(ostream&o) const
{
  o << indent << "R12Technology:" << endl;
  o << incindent;

  if (!safety_check())
    o << indent << "WARNING: ---- safety checks SKIPPED ----" << endl;
  corrfactor()->print(o); o << endl;
  o << indent << "Coupling included: " << (coupling_ ? "true" : "false") << endl;
  o << indent << "GBC assumed: " << (gbc_ ? "true" : "false") << endl;
  o << indent << "EBC assumed: " << (ebc_ ? "true" : "false") << endl;
  o << indent << "EBC-free method: Valeev" << endl;
  switch (posdef_B()) {
    case PositiveDefiniteB_no:     o << indent << "Do not enforce positive definiteness of B" << endl;  break;
    case PositiveDefiniteB_yes:    o << indent << "Enforce positive definiteness of B" << endl;  break;
    case PositiveDefiniteB_weak:   o << indent << "Enforce positive definiteness of B, but not ~B(ij)" << endl;  break;
  }
  if (stdapprox_ == StdApprox_B && omit_P_) {
    o << indent << "Intermediate P is omitted" << endl;
  }
  switch(abs_method_) {
  case ABS_CABS :
    o << indent << "RI method for many-electron integrals: CABS" << endl;
    break;
  case ABS_CABSPlus :
    o << indent << "RI method for many-electron integrals: CABS+ (CABS using the union of OBS and ABS for RI)" << endl;
    break;
  }
  if (!this->abs_eq_obs_) {
    if (this->abs_nlindep() != -1)
      o << indent << "# of linearly depenendent vectors to be removed from ABS = " << this->abs_nlindep() << endl;
    else
      o << indent << "ABS linear depenendence tolerance = " << this->abs_lindep_tol() << endl;
  }
  switch (stdapprox_) {
    case StdApprox_Ap :
      o << indent << "Standard Approximation: A'" << endl;
    break;
    case StdApprox_App :
      o << indent << "Standard Approximation: A''" << endl;
    break;
    case StdApprox_B :
      o << indent << "Standard Approximation: B" << endl;
    break;
    case StdApprox_C :
      o << indent << "Standard Approximation: C" << endl;
    break;
    case StdApprox_Cp :
      o << indent << "Standard Approximation: C'" << endl;
    break;
  }
  ansatz()->print(o);

  o << indent << "H0_dk_approx_pauli: ";
  switch(H0_dk_approx()) {
    case H0_dk_approx_pauli_true:    o << "true" << endl;  break;
    case H0_dk_approx_pauli_false:   o << "false" << endl; break;
    case H0_dk_approx_pauli_fHf:     o << "fHf" << endl;   break;
    case H0_dk_approx_pauli_fHf_Q:   o << "fHf_Q" << endl; break;
  }
  if (H0_dk_approx() == H0_dk_approx_pauli_false &&
      stdapprox() == StdApprox_App) {
    o << indent << "H0_dk_keep: " << (this->H0_dk_keep() ? "true" : "false") << endl;
  }

  o << indent << "Max # ABS indices: " << maxnabs_ << endl;

  o << decindent;
}

////////////////////////////////////////////////////////////////////////////

const Ref<R12Technology::CorrelationFactor>&
R12Technology::corrfactor() const
{
  return corrfactor_;
}

////////////////////////////////////////////////////////////////////////////

void
R12Technology::corrfactor(const Ref<CorrelationFactor>& cf)
{
  if (!corrfactor_->equiv(cf)) {
      corrfactor_ = cf;
  }
}

/////////////////////////////////////////////////////////////////////////////

bool
R12Technology::omit_P() const
{
  return omit_P_;
}

/////////////////////////////////////////////////////////////////////////////

R12Technology::H0_dk_approx_pauli
R12Technology::H0_dk_approx() const
{
  return H0_dk_approx_pauli_;
}

/////////////////////////////////////////////////////////////////////////////

bool
R12Technology::H0_dk_keep() const
{
  return H0_dk_keep_;
}

/////////////////////////////////////////////////////////////////////////////

unsigned int
R12Technology::maxnabs() const
{
  return maxnabs_;
}

/////////////////////////////////////////////////////////////////////////////

int
R12Technology::abs_nlindep() const
{
  return abs_nlindep_;
}

/////////////////////////////////////////////////////////////////////////////

double
R12Technology::abs_lindep_tol() const
{
  return abs_lindep_tol_;
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

bool
R12Technology::coupling() const
{
  return coupling_;
}

/////////////////////////////////////////////////////////////////////////////

bool
R12Technology::compute_1rdm() const
{
  return compute_1rdm_;
}

bool
R12Technology::coupling_1rdm_f12b() const
{
  return coupling_1rdm_f12b_;
}

/////////////////////////////////////////////////////////////////////////////

R12Technology::ABSMethod
R12Technology::abs_method() const
{
  return abs_method_;
}

/////////////////////////////////////////////////////////////////////////////

R12Technology::StandardApproximation
R12Technology::stdapprox() const
{
  return stdapprox_;
}

/////////////////////////////////////////////////////////////////////////////

const Ref<R12Technology::R12Ansatz>&
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

R12Technology::PositiveDefiniteB
R12Technology::posdef_B() const
{
  return posdef_B_;
}

/////////////////////////////////////////////////////////////////////////////

bool
R12Technology::omit_B() const
{
  return omit_B_;
}

/////////////////////////////////////////////////////////////////////////////

void
R12Technology::check_integral_factory(const Ref<Integral>& ints)
{
  // any factory can support pure MP2 calculations! Thus only test for nontrivial corr factors
  Ref<NullCorrelationFactor> nullcf; nullcf << corrfactor();
  if (nullcf.null()) {
    // Only IntegralLibint2 can be used at the moment
    bool allowed_integral_factory = false;
#ifdef HAVE_LIBINT2
    IntegralLibint2* libint2intf = dynamic_cast<IntegralLibint2*>(ints.pointer());
    if (libint2intf) {
      allowed_integral_factory = true;
    }
#endif
    if (!allowed_integral_factory) {
      InputError ex("R12Technology::check_integral_factory_(): invalid integral factory provided.",
                  __FILE__, __LINE__, 0, 0, class_desc());
      try {
        ex.elaborate() << "Try using IntegralLibint2."
                       << std::endl
                       << "The given integral factory was of type " << ints->class_name()
                       << std::endl;
      }
      catch (...) {}
      throw ex;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////

std::string
R12Technology::default_cabs_name(const std::string& obs_name) {
  std::string cabs_name;
  if (obs_name == "cc-pVDZ-F12")
    cabs_name = "cc-pVDZ-F12-CABS";
  else if (obs_name == "cc-pVTZ-F12")
    cabs_name = "cc-pVTZ-F12-CABS";
  else if (obs_name == "cc-pVQZ-F12")
    cabs_name = "cc-pVQZ-F12-CABS";
  else if (obs_name == "aug-cc-pVDZ")
    cabs_name = "aug-cc-pVDZ-CABS";
  else if (obs_name == "aug-cc-pVTZ")
    cabs_name = "aug-cc-pVTZ-CABS";
  else if (obs_name == "aug-cc-pVQZ")
    cabs_name = "aug-cc-pVQZ-CABS";
  else if (obs_name == "aug-cc-pV5Z")
    cabs_name = "aug-cc-pV5Z-CABS";
  return cabs_name;
}

Ref<GaussianBasisSet>
R12Technology::make_auto_cabs(const Ref<GaussianBasisSet>& bs) {
  Ref<GaussianBasisSet> cabs;

  // if bs has a known name, use a hardwired CABS name
  const std::string cabs_name = default_cabs_name(std::string(bs->label()));

  if (cabs_name.empty()) { // no CABS in the library -- make a relatively conservative CABS
    // default CABS = Uncontracted(aug-cc-pV5Z) limited up to 3 L_occ + 1
    Ref<GaussianBasisSet> cabs_core;
    {
    Ref<AssignedKeyVal> tmpkv = new AssignedKeyVal;
    tmpkv->assign("name", "aug-cc-pV5Z");
    tmpkv->assign("puream", "true");
    tmpkv->assign("molecule", bs->molecule().pointer());
    Ref<KeyVal> kv = tmpkv;
    cabs_core = new GaussianBasisSet(kv);
    }
    Ref<GaussianBasisSet> cabs_core_uncontr;
    {
    Ref<AssignedKeyVal> tmpkv = new AssignedKeyVal;
    tmpkv->assign("basis", cabs_core.pointer());
    Ref<KeyVal> kv = tmpkv;
    cabs_core_uncontr = new UncontractedBasisSet(kv);
    }
    {
    Ref<AssignedKeyVal> tmpkv = new AssignedKeyVal;
    tmpkv->assign("basis", cabs_core_uncontr.pointer());
    int lmax;
    if (bs->molecule()->max_z() <= 20) // H - Ca
    lmax = 3;
    else if (bs->molecule()->max_z() <= 56) // Sc - Ba
    lmax = 6;
    else
    // La - higher
    lmax = 9;
    tmpkv->assign("lmax", lmax);
    Ref<KeyVal> kv = tmpkv;
    cabs = new LSelectBasisSet(kv);
    }
  }
  else { // have an appropriate CABS in the library
    Ref<AssignedKeyVal> tmpkv = new AssignedKeyVal;
    tmpkv->assign("name", cabs_name.c_str());
    tmpkv->assign("molecule", bs->molecule().pointer());
    Ref<KeyVal> kv = tmpkv;
    cabs = new GaussianBasisSet(kv);
  }
  return cabs;
}

double
R12Technology::default_stg_exponent(const std::string & obs_name)
{
  double result = 0.0;
  /// recommended values from K.A. Peterson, T.B. Adler, H.J. Werner, J. Chem. Phys. 128 (2008) 084102.
  if (obs_name == "cc-pVDZ-F12")
    result = 0.9;
  else if (obs_name == "cc-pVTZ-F12")
    result = 1.0;
  else if (obs_name == "cc-pVQZ-F12")
    result = 1.1;
  else if (obs_name == "aug-cc-pVDZ")
    result = 1.1;
  else if (obs_name == "aug-cc-pVTZ")
    result = 1.2;
  else if (obs_name == "aug-cc-pVQZ")
    result = 1.4;
  else if (obs_name == "aug-cc-pV5Z")
    result = 1.4;
  return result;
}



/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
