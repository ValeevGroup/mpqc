//
// intcca.cc
//
// Copyright (C) 2004 Sandia National Laboratories
//
// Author: Joe Kenny <jpkenny@sandia.gov>
// Maintainer: Joe Kenny
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

#include <util/state/stateio.h>
#include <util/misc/ccaenv.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/intcca/intcca.h>
#include <chemistry/qc/intcca/obintcca.h>
#include <chemistry/qc/intcca/tbintcca.h>
#include <util/class/scexception.h>
#ifdef INTV3_ORDER
  #include <chemistry/qc/intv3/cartitv3.h>
  #include <chemistry/qc/intv3/tformv3.h>
#else
  #include <chemistry/qc/intcca/cartit.h>
  #include <chemistry/qc/intcca/tform.h>
#endif

using namespace std;
using namespace sc;
using namespace Chemistry::QC::GaussianBasis;

static ClassDesc IntegralCCA_cd(
  typeid(IntegralCCA),"IntegralCCA",1,"public Integral",
  0, create<IntegralCCA>, create<IntegralCCA>);

extern Ref<Integral> default_integral;

/* may need to add optional "eval_factory" argument to this method in
integral class to get this capability
Integral*
Integral::get_default_integral()
{
  if (default_integral.null())
    default_integral = new IntegralCCA();

  return default_integral;
}
*/

IntegralCCA::IntegralCCA(IntegralEvaluatorFactory eval_factory, bool use_opaque,
                         const Ref<GaussianBasisSet> &b1,
                         const Ref<GaussianBasisSet> &b2,
                         const Ref<GaussianBasisSet> &b3,
                         const Ref<GaussianBasisSet> &b4):
  Integral(b1,b2,b3,b4), eval_factory_(eval_factory)
{
  use_opaque_ = use_opaque;
  initialize_transforms();
}

IntegralCCA::IntegralCCA(const Ref<KeyVal> &keyval):
  Integral(keyval)
{

  initialize_transforms();

  string buffer = keyval->stringvalue("integral_buffer");
  if ( keyval->error() != KeyVal::OK ) buffer = "opaque";
  if ( buffer == "opaque" ) use_opaque_ = 1;
  else if ( buffer == "array" ) use_opaque_ = 0;
  else {
    InputError ex( "integral_buffer must be either opaque or array",
                   __FILE__, __LINE__,"integral_buffer",
                   buffer.c_str(),class_desc() );
    throw ex;
  }

  factory_type_ = keyval->stringvalue("evaluator_factory");
  if ( keyval->error() != KeyVal::OK ) {
    factory_type_ = string("MPQC.IntegralEvaluatorFactory");
  }
#ifdef INTV3_ORDER
/* this has gotten more complicated
  if(package_ == "cints") {
    InputError ex("using intv3 ordering, can't use cints",__FILE__, __LINE__);
    try { ex.elaborate() << "INTV3_ORDER=yes in LocalMakefile,"
                         << " this option is for development use only";
    }
    catch (...) {}
    throw ex;
  }
*/
#endif

  sc_molecule_ << keyval->describedclassvalue("molecule");
  if (sc_molecule_.null())
    throw InputError("molecule is required",__FILE__,__LINE__);

  gov::cca::Services &services = *CCAEnv::get_services();
  gov::cca::ports::BuilderService &bs = *CCAEnv::get_builder_service();
  gov::cca::TypeMap &type_map = *CCAEnv::get_type_map();
  gov::cca::ComponentID &my_id = *CCAEnv::get_component_id();

  /////////////////////////////////
  // get and configure eval factory
  /////////////////////////////////

  //get the factory
  fac_id_ = bs.createInstance("evaluator_factory",factory_type_,type_map);
  services.registerUsesPort("IntegralEvaluatorFactory",
                            "Chemistry.QC.GaussianBasis.IntegralEvaluatorFactory",
                            type_map);
  fac_con_ = bs.connect(my_id,"IntegralEvaluatorFactory",
                        fac_id_,"IntegralEvaluatorFactory");
  eval_factory_ = services.getPort("IntegralEvaluatorFactory");

  // create config objects
  ob_config_ = 
    Chemistry::Chemistry_QC_GaussianBasis_ObIntEvalConfig::_create();
  tb_config_ = 
    Chemistry::Chemistry_QC_GaussianBasis_TbIntEvalConfig::_create();

  // set default package on config objects
  package_ = keyval->stringvalue("default_package");
  if ( keyval->error() == KeyVal::OK ) {
    ExEnv::out0() << indent << "Default integral package: " << package_ 
		  << std::endl;
    if( package_ == "intv3" ) {
      ob_config_.set_default_pkg(Package_INTV3);
      tb_config_.set_default_pkg(Package_INTV3);
    }
    else if( package_ == "cints" ) {
      ob_config_.set_default_pkg(Package_CINTS);
      tb_config_.set_default_pkg(Package_CINTS);
    }
    else {
      InputError ex("Unrecognized package",__FILE__, __LINE__,
                  "default_package",package_.c_str(),class_desc());
      throw ex;
    }
  }

  // fill in config objects
  int ntype = keyval->count("type");
  int npkg = keyval->count("package");
  string tp, pkg;
  if( ntype != npkg ) throw InputError("ntype != npackage",__FILE__,__LINE__);
  for( int i=0; i<ntype; ++i) {
    tp = keyval->stringvalue("type",i);
    pkg = keyval->stringvalue("package",i);
    set_config(tp,pkg);
    ExEnv::out0() << indent << "Integral type " << tp << ": " 
		  << pkg << std::endl;
  }

  // pass configs to factory
  eval_factory_.set_obint_config( ob_config_ );
  eval_factory_.set_tbint_config( tb_config_ );

  // set molecule on factory
  molecule_ = Chemistry::Chemistry_Molecule::_create();
  molecule_.initialize(sc_molecule_->natom(),"bohr");
  for( int i=0; i<sc_molecule_->natom(); ++i ) {
    molecule_.set_atomic_number( i, sc_molecule_->Z(i) );
    for( int j=0; j<3; ++j )
      molecule_.set_cart_coor( i, j, sc_molecule_->r(i,j) );
  }
  eval_factory_.set_molecule(molecule_);

}

IntegralCCA::IntegralCCA(StateIn& s) :
  Integral(s)
{
  initialize_transforms();
}

void
IntegralCCA::save_data_state(StateOut& s)
{
  Integral::save_data_state(s);
}

IntegralCCA::~IntegralCCA()
{
  free_transforms();
}

Integral*
IntegralCCA::clone()
{
  // ???
  return new IntegralCCA(eval_factory_,use_opaque_);
  // this wouldn't take much work
  //throw FeatureNotImplemented("clone not implemented",
  //                            __FILE__,__LINE__);
}

CartesianIter *
IntegralCCA::new_cartesian_iter(int l)
{
#ifdef INTV3_ORDER
  return new CartesianIterV3(l);
#else
  return new CartesianIterCCA(l);
  //return new CartesianIterV3(l);
#endif
}

RedundantCartesianIter *
IntegralCCA::new_redundant_cartesian_iter(int l)
{
#ifdef INTV3_ORDER
  return new RedundantCartesianIterV3(l);
#else
  return new RedundantCartesianIterCCA(l);
  //return new RedundantCartesianIterV3(l);
#endif
}

RedundantCartesianSubIter *
IntegralCCA::new_redundant_cartesian_sub_iter(int l)
{
#ifdef INTV3_ORDER
  return new RedundantCartesianSubIterV3(l);
#else
  return new RedundantCartesianSubIterCCA(l);
  //return new RedundantCartesianSubIterV3(l);
#endif
}

SphericalTransformIter *
IntegralCCA::new_spherical_transform_iter(int l, int inv, int subl)
{
#ifdef INTV3_ORDER

  if (l>maxl_ || l<0)
      throw ProgrammingError("new_spherical_transform_iter: bad l",
                             __FILE__,__LINE__);
  if (subl == -1) subl = l;
  if (subl < 0 || subl > l || (l-subl)%2 != 0)
      throw ProgrammingError("new_spherical_transform_iter: bad subl",
                             __FILE__,__LINE__);
  if (inv)
      return new SphericalTransformIter(ist_[l][(l-subl)/2]);
  return new SphericalTransformIter(st_[l][(l-subl)/2]);

#else

  // CINTS version
  if (l>maxl_ || l<0)
      throw ProgrammingError("new_spherical_transform_iter: bad l",
                             __FILE__,__LINE__);
  if (subl == -1) subl = l;
  if (subl < 0 || subl > l || (l-subl)%2 != 0)
      throw ProgrammingError("new_spherical_transform_iter: bad subl",
                             __FILE__,__LINE__);
  if (inv)
      return new SphericalTransformIter(ist_[l][(l-subl)/2]);
  return new SphericalTransformIter(st_[l][(l-subl)/2]);

#endif

}

const SphericalTransform *
IntegralCCA::spherical_transform(int l, int inv, int subl)
{
#ifdef INTV3_ORDER

  // INTV3 version
  if (l>maxl_ || l<0)
      throw ProgrammingError("spherical_transform_iter: bad l",
                             __FILE__,__LINE__);
  if (subl == -1) subl = l;
  if (subl < 0 || subl > l || (l-subl)%2 != 0)
      throw ProgrammingError("spherical_transform_iter: bad subl",
                             __FILE__,__LINE__);
  if (inv)
      return ist_[l][(l-subl)/2];
  return st_[l][(l-subl)/2];

#else

  // CINTS version
  if (l>maxl_ || l<0)
      throw ProgrammingError("spherical_transform_iter: bad l",
                             __FILE__,__LINE__);
  if (subl == -1) subl = l;
  if (subl < 0 || subl > l || (l-subl)%2 != 0)
      throw ProgrammingError("spherical_transform_iter: bad subl",
                             __FILE__,__LINE__);
  if (inv)
      return ist_[l][(l-subl)/2];
  return st_[l][(l-subl)/2];

#endif

}

Ref<OneBodyInt>
IntegralCCA::overlap()
{
  return new OneBodyIntCCA(this, bs1_, bs2_, eval_factory_, 
			   &Int1eCCA::overlap, use_opaque_ );
}

Ref<OneBodyInt>
IntegralCCA::kinetic()
{
  return new OneBodyIntCCA(this, bs1_, bs2_, eval_factory_, 
			   &Int1eCCA::kinetic, use_opaque_ );
}

Ref<OneBodyInt>
IntegralCCA::nuclear()
{
  return new OneBodyIntCCA(this, bs1_, bs2_, eval_factory_, 
			   &Int1eCCA::nuclear, use_opaque_ );
}

Ref<OneBodyInt>
IntegralCCA::hcore()
{
  return new OneBodyIntCCA(this, bs1_, bs2_, eval_factory_, &Int1eCCA::hcore, 
                           use_opaque_ );
}

Ref<OneBodyInt>
IntegralCCA::point_charge(const Ref<PointChargeData>& dat)
{
  throw FeatureNotImplemented("point_charge not implemented",
                              __FILE__,__LINE__);
}

Ref<OneBodyInt>
IntegralCCA::efield_dot_vector(const Ref<EfieldDotVectorData>&dat)
{
  throw FeatureNotImplemented("efield_dot_vector not iplemented",
                              __FILE__,__LINE__);
}

Ref<OneBodyInt>
IntegralCCA::dipole(const Ref<DipoleData>& dat)
{
  throw FeatureNotImplemented("dipole not implemented",
                              __FILE__,__LINE__);
}

Ref<OneBodyInt>
IntegralCCA::quadrupole(const Ref<DipoleData>& dat)
{
  throw FeatureNotImplemented("quadrupole not implemented",
                              __FILE__,__LINE__);
}

Ref<OneBodyDerivInt>
IntegralCCA::overlap_deriv()
{
   return new OneBodyDerivIntCCA(this, bs1_, bs2_, eval_factory_, use_opaque_, 
                                 ObIntEvalType_OVERLAP);
}

Ref<OneBodyDerivInt>
IntegralCCA::kinetic_deriv()
{
   return new OneBodyDerivIntCCA(this, bs1_, bs2_, eval_factory_, use_opaque_, 
                                 ObIntEvalType_KINETIC);
}

Ref<OneBodyDerivInt>
IntegralCCA::nuclear_deriv()
{
   return new OneBodyDerivIntCCA(this, bs1_, bs2_, eval_factory_, use_opaque_, 
                                 ObIntEvalType_NUCLEAR);
}

Ref<OneBodyDerivInt>
IntegralCCA::hcore_deriv()
{
   return new OneBodyDerivIntCCA(this, bs1_, bs2_, eval_factory_, use_opaque_,
                                 ObIntEvalType_OEHAM);
}

Ref<TwoBodyInt>
IntegralCCA::electron_repulsion()
{
   return new TwoBodyIntCCA(this, bs1_, bs2_, bs3_, bs4_, 
                            storage_, eval_factory_, use_opaque_, 
			    TbIntEvalType_ERI4 );
}

Ref<TwoBodyDerivInt>
IntegralCCA::electron_repulsion_deriv()
{
   return new TwoBodyDerivIntCCA(this, bs1_, bs2_, bs3_, bs4_, 
                                 storage_, eval_factory_, 
				 use_opaque_, TbIntEvalType_ERI4 );
}

void
IntegralCCA::set_basis(const Ref<GaussianBasisSet> &b1,
                       const Ref<GaussianBasisSet> &b2,
                       const Ref<GaussianBasisSet> &b3,
                       const Ref<GaussianBasisSet> &b4)
{
  free_transforms();
  Integral::set_basis(b1,b2,b3,b4);
  initialize_transforms();
}

void
IntegralCCA::free_transforms()
{
#ifdef INTV3_ORDER

  // INTV3 version
  int i,j;
  for (i=0; i<=maxl_; i++) {
      for (j=0; j<=i/2; j++) {
          delete st_[i][j];
          delete ist_[i][j];
        }
      delete[] st_[i];
      delete[] ist_[i];
    }
  delete[] st_;
  delete[] ist_;
  st_ = 0;
  ist_ = 0;

#else

  // CINTS version
  int i,j;
  for (i=0; i<=maxl_; i++) {
      for (j=0; j<=i/2; j++) {
          delete st_[i][j];
          delete ist_[i][j];
        }
      delete[] st_[i];
      delete[] ist_[i];
    }
  if (maxl_ >= 0) {
    delete[] st_;
    delete[] ist_;
  }
  st_ = NULL;
  ist_ = NULL;

#endif

}

 void
 IntegralCCA::initialize_transforms()
 {
   maxl_ = -1;
   int maxam;
   maxam = bs1_.nonnull()?bs1_->max_angular_momentum():-1;
   if (maxl_ < maxam) maxl_ = maxam;
   maxam = bs2_.nonnull()?bs2_->max_angular_momentum():-1;
   if (maxl_ < maxam) maxl_ = maxam;
   maxam = bs3_.nonnull()?bs3_->max_angular_momentum():-1;
   if (maxl_ < maxam) maxl_ = maxam;
   maxam = bs4_.nonnull()?bs4_->max_angular_momentum():-1;
   if (maxl_ < maxam) maxl_ = maxam;
 
#ifdef INTV3_ORDER

   // INTV3 version
   st_ = new SphericalTransform**[maxl_+1];
   ist_ = new ISphericalTransform**[maxl_+1];
   int i,j;
   for (i=0; i<=maxl_; i++) {
       st_[i] = new SphericalTransform*[i/2+1];
       ist_[i] = new ISphericalTransform*[i/2+1];
       for (j=0; j<=i/2; j++) {
           st_[i][j] = new SphericalTransformV3(i,i-2*j);
           ist_[i][j] = new ISphericalTransformV3(i,i-2*j);
         }
     }
#else

  // CINTS version
  if (maxl_ >= 0) {
    st_ = new SphericalTransform**[maxl_+1];
    ist_ = new ISphericalTransform**[maxl_+1];
    int i,j;
    for (i=0; i<=maxl_; i++) {
      st_[i] = new SphericalTransform*[i/2+1];
      ist_[i] = new ISphericalTransform*[i/2+1];
      for (j=0; j<=i/2; j++) {
        st_[i][j] = new SphericalTransformCCA(i,i-2*j);
        ist_[i][j] = new ISphericalTransformCCA(i,i-2*j);
      }
    }
  }
  else {
    st_ = NULL;
    ist_ = NULL;
  }

#endif

}

void 
IntegralCCA::set_config( string type, string package) 
{

  // set package
  Package pkg;
  if( package == "intv3" )
    pkg = Package_INTV3;
  else if( package == "cints" )
    pkg = Package_CINTS;
  else 
    throw InputError("Unrecognized package",__FILE__, __LINE__,
		     "package",package_.c_str(),class_desc());
    
  // one-body types
  if( type == "overlap" )
    ob_config_.set_pkg_config(ObIntEvalType_OVERLAP,pkg);
  else if( type == "kinetic" )
    ob_config_.set_pkg_config(ObIntEvalType_KINETIC,pkg);
  else if( type == "nuclear" )
    ob_config_.set_pkg_config(ObIntEvalType_NUCLEAR,pkg);
  else if( type == "hcore" )
    ob_config_.set_pkg_config(ObIntEvalType_OEHAM,pkg);
  else if( type == "pointcharge1" )
    ob_config_.set_pkg_config(ObIntEvalType_POINTCHARGE1,pkg);
  else if( type == "pointcharge2" )
    ob_config_.set_pkg_config(ObIntEvalType_POINTCHARGE2,pkg);
  else if( type == "efield_dot_vector" )
    ob_config_.set_pkg_config(ObIntEvalType_EFIELD_DOT_VECTOR,pkg);
  else if( type == "dipole" )
    ob_config_.set_pkg_config(ObIntEvalType_DIPOLE,pkg);
  else if( type == "quadrupole" )
    ob_config_.set_pkg_config(ObIntEvalType_QUADRUPOLE,pkg);
  else if( type == "dk" )
    ob_config_.set_pkg_config(ObIntEvalType_DK,pkg);
  else if( type == "ecp" )
    ob_config_.set_pkg_config(ObIntEvalType_ECP,pkg);
  else if( type == "projmpole" )
    ob_config_.set_pkg_config(ObIntEvalType_PROJMPOLE,pkg);
  
  // two-body types
  else if( type == "eri2" )
    tb_config_.set_pkg_config(TbIntEvalType_ERI2,pkg);
  else if( type == "eri3" )
    tb_config_.set_pkg_config(TbIntEvalType_ERI3,pkg);
  else if( type == "eri4" )
    tb_config_.set_pkg_config(TbIntEvalType_ERI4,pkg);
  else if( type == "grt" )
    tb_config_.set_pkg_config(TbIntEvalType_GRT,pkg);

  // punt
  else 
    throw InputError("Unrecognized integral type",__FILE__, __LINE__,
		     "type",type.c_str(),class_desc());  

}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
