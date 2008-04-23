//
// intcca.h
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

#ifndef _chemistry_cca_int_intcca_h
#define _chemistry_cca_int_intcca_h

#include <gov_cca.hxx>
#include <chemistry/qc/basis/integral.h>
#include <Chemistry_QC_GaussianBasis_IntegralEvaluatorFactoryInterface.hxx>
#include <ChemistryCXX_Molecule.hxx>
#include <Chemistry_QC_GaussianBasis_DerivCentersInterface.hxx>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/transform.h>
#include <chemistry/qc/basis/cartiter.h>
#include "obintcca.h"
#include "tbintcca.h"

namespace sc {

/** IntegralCCA provides an SC client for CCA IntegralEvaluator components. */

  class IntegralCCA : public Integral {

  private:
    
    //---------------------------------------------------------------------
    // function object code for evaluator generation 
    //---------------------------------------------------------------------

    class onebody_generator {

    private:

      Integral* integral_;
      Chemistry::QC::GaussianBasis::IntegralEvaluatorFactoryInterface factory_;
      bool reorder_;
      Ref<GaussianBasisSet> bs1_, bs2_;

    public:
      
      onebody_generator( ) { }
      
      onebody_generator(
          Integral* integral,
          Chemistry::QC::GaussianBasis::IntegralEvaluatorFactoryInterface fac, 
          bool reorder ):
	integral_(integral), factory_(fac), reorder_(reorder) 
      { }

      void set_basis( Ref<GaussianBasisSet> bs1,
		      Ref<GaussianBasisSet> bs2 ) 
      { bs1_ = bs1, bs2_ = bs2; }
      
      Ref<OneBodyIntCCA> generate(
          Chemistry::QC::GaussianBasis::CompositeIntegralDescrInterface cdesc )
      {
	Ref<OneBodyIntCCA> eval;
	eval = new OneBodyIntCCA( integral_, bs1_, bs2_, 
				  factory_, cdesc, reorder_ );
	return eval;
      }

    };

    //---------------------------------------------------------------------

    class onebody_deriv_generator {

    private:

      Integral* integral_;
      Chemistry::QC::GaussianBasis::IntegralEvaluatorFactoryInterface factory_;
      bool reorder_;
      Ref<GaussianBasisSet> bs1_, bs2_;

    public:
      
      onebody_deriv_generator( ) { }

      onebody_deriv_generator(
          Integral* integral,
          Chemistry::QC::GaussianBasis::IntegralEvaluatorFactoryInterface fac, 
          bool reorder ):
	integral_(integral), factory_(fac), reorder_(reorder)
      { }

      void set_basis( Ref<GaussianBasisSet> bs1,
		      Ref<GaussianBasisSet> bs2 ) 
      { bs1_ = bs1, bs2_ = bs2; }

      Ref<OneBodyDerivIntCCA> generate(
          Chemistry::QC::GaussianBasis::CompositeIntegralDescrInterface cdesc )
      {
	Ref<OneBodyDerivIntCCA> eval;
	eval = new OneBodyDerivIntCCA( integral_, bs1_, bs2_, 
				       factory_, cdesc, reorder_ );
	return eval;
      }

    };

    //------------------------------------------------------------------------

    class twobody_generator {

    private:

      Ref<GaussianBasisSet> bs1_, bs2_, bs3_, bs4_;
      Integral* integral_;
      Chemistry::QC::GaussianBasis::IntegralEvaluatorFactoryInterface factory_;

    public:
      
      twobody_generator( ) { }

      twobody_generator(
          Integral* integral,
          Chemistry::QC::GaussianBasis::IntegralEvaluatorFactoryInterface fac ):
	integral_(integral), factory_(fac)
      { }

      void set_basis( Ref<GaussianBasisSet> bs1,
		      Ref<GaussianBasisSet> bs2,
		      Ref<GaussianBasisSet> bs3,
		      Ref<GaussianBasisSet> bs4 ) 
      { bs1_ = bs1; bs2_ = bs2; bs3_ = bs3; bs4_ = bs4; }

      Ref<TwoBodyIntCCA> generate(
          Chemistry::QC::GaussianBasis::CompositeIntegralDescrInterface cdesc )
      {
	Ref<TwoBodyIntCCA> eval;
	eval = new TwoBodyIntCCA( integral_, bs1_, bs2_, bs3_, bs4_,
				  factory_, cdesc );
	return eval;
      }

    };

    //------------------------------------------------------------------------

    class twobody_deriv_generator {

    private:

      Ref<GaussianBasisSet> bs1_, bs2_, bs3_, bs4_;
      Integral* integral_;
      Chemistry::QC::GaussianBasis::IntegralEvaluatorFactoryInterface factory_;

    public:
      
      twobody_deriv_generator( ) { }

      twobody_deriv_generator(
          Integral* integral,
          Chemistry::QC::GaussianBasis::IntegralEvaluatorFactoryInterface fac ):
	integral_(integral), factory_(fac)
      { }

      void set_basis( Ref<GaussianBasisSet> bs1,
		      Ref<GaussianBasisSet> bs2,
		      Ref<GaussianBasisSet> bs3,
		      Ref<GaussianBasisSet> bs4 ) 
      { bs1_ = bs1; bs2_ = bs2; bs3_ = bs3; bs4_ = bs4; }

      Ref<TwoBodyDerivIntCCA> generate(
          Chemistry::QC::GaussianBasis::CompositeIntegralDescrInterface cdesc )
      {
	Ref<TwoBodyDerivIntCCA> eval;
	eval = new TwoBodyDerivIntCCA( integral_, bs1_, bs2_, bs3_, bs4_,
				       factory_, cdesc );
	return eval;
      }

    };

    //----------------------------------------------------------------------
	
    // the function object

    template< typename eval_type, typename generator_type >
    class sc_eval_factory {

    private:
      
      generator_type generator_;

    public:

      sc_eval_factory() { }
      
      sc_eval_factory( generator_type generator ):
      generator_(generator)
      { }

      Ref<eval_type> operator() (
          Chemistry::QC::GaussianBasis::CompositeIntegralDescrInterface cdesc )
      {
	return generator_.generate( cdesc );
      }

    };

    //----------------------------------------------------------------------
    
    int maxl_;
    bool intv3_order_;
    bool use_superfac_;
    gov::cca::ComponentID fac_id_;
    gov::cca::ConnectionID fac_con_;
    Ref<Molecule> sc_molecule_;
    ChemistryCXX::Molecule molecule_;
    std::vector<Chemistry::QC::GaussianBasis::DerivCentersInterface> cca_dcs_;
    std::vector<Chemistry::QC::GaussianBasis::IntegralDescrInterface> descs_;
    std::string buffer_;
    std::string default_subfactory_;
    Chemistry::QC::GaussianBasis::IntegralEvaluatorFactoryInterface eval_factory_;
    Chemistry::QC::GaussianBasis::CompositeIntegralDescrInterface eval_req_;

    sidl::array<std::string> types_;
    sidl::array<std::string> derivs_;
    sidl::array<std::string> sfacs_;
    std::vector<std::string> std_types_;
    std::vector<std::string> std_derivs_;
    std::vector<std::string> std_sfacs_;

    onebody_generator obgen_;
    onebody_deriv_generator obdgen_;
    twobody_generator tbgen_;
    twobody_deriv_generator tbdgen_;
    sc_eval_factory<OneBodyInt,onebody_generator> get_onebody;
    sc_eval_factory<OneBodyDerivInt,onebody_deriv_generator> get_onebody_deriv;
    sc_eval_factory<TwoBodyInt,twobody_generator> get_twobody;
    sc_eval_factory<TwoBodyDerivInt,twobody_deriv_generator> get_twobody_deriv;
    
    SphericalTransform ***st_;
    ISphericalTransform ***ist_;
    
    void init_factory();
    void init_generators();
    void free_transforms();
    void initialize_transforms();
    void free_transformsV3();
    void initialize_transformsV3();
    CartesianIter* new_cartesian_iterV3(int l);
    RedundantCartesianIter* new_redundant_cartesian_iterV3(int l);
    RedundantCartesianSubIter* new_redundant_cartesian_sub_iterV3(int l);
    SphericalTransformIter* new_spherical_transform_iterV3(int l, int inv, int subl);
    const SphericalTransform* spherical_transformV3(int l, int inv, int subl);
    
  public:

    /** The KeyVal constructor.
        This constructor is used when the framework is embedded.
        The following keywords are read:

        <dl>
        <dt><tt>evaluator_factory</tt><dd> This gives the symbol name of a 
        CCA IntegralEvaluatorFactory component.  This symbol name should
        also appear in the cca-load argument.  The default is
        <tt>MPQC.IntegralEvaluatorFactory</tt>.

        <dt><tt>integral_package</tt><dd> If the default 
        <tt>MPQC.IntegralEvaluatorFactory</tt> is used, then this 
        option may be used to specify the integrals package to use
        (<tt>intv3</tt>, <tt>cints</tt>, or <tt>libint2</tt>).  The default is <tt>intv3</tt>.

        <dt><tt>molecule</tt><dd> This gives a molecule object, it is required.
        </dl>
    */

    IntegralCCA(const Ref<KeyVal>&);

    IntegralCCA( const Ref<GaussianBasisSet> &b1,
                 const Ref<GaussianBasisSet> &b2,
                 const Ref<GaussianBasisSet> &b3,
                 const Ref<GaussianBasisSet> &b4,
                 std::string default_sf,
                 bool use_superfac,
                 sidl::array<std::string> types,
                 sidl::array<std::string> derivs,
                 sidl::array<std::string> sfacs,
                 bool intv3_order
                );

    IntegralCCA(StateIn&);

    ~IntegralCCA();

    void save_data_state(StateOut&);

    Integral* clone();

    void set_storage(size_t i) 
    { 
      storage_=i; 
      eval_factory_.set_storage( storage_ );
    };

    CartesianIter * new_cartesian_iter(int);
    RedundantCartesianIter * new_redundant_cartesian_iter(int);
    RedundantCartesianSubIter * new_redundant_cartesian_sub_iter(int);
    SphericalTransformIter * new_spherical_transform_iter(int l,
                                                          int inv=0,
                                                          int subl=-1);
    const SphericalTransform * spherical_transform(int l,
                                                   int inv=0, int subl=-1);

    Ref<OneBodyInt> overlap();

    Ref<OneBodyInt> p_dot_nuclear_p();

    Ref<OneBodyInt> p_cross_nuclear_p();

    Ref<OneBodyInt> kinetic();

    Ref<OneBodyInt> point_charge(const Ref<PointChargeData>& =0);

    Ref<OneBodyInt> nuclear();

    Ref<OneBodyInt> hcore();

    Ref<OneBodyInt> efield_dot_vector(const Ref<EfieldDotVectorData>& =0);

    Ref<OneBodyInt> dipole(const Ref<DipoleData>& =0);

    Ref<OneBodyInt> quadrupole(const Ref<DipoleData>& =0);

    Ref<OneBodyDerivInt> overlap_deriv();
                                     
    Ref<OneBodyDerivInt> kinetic_deriv();
                                     
    Ref<OneBodyDerivInt> nuclear_deriv();
                                     
    Ref<OneBodyDerivInt> hcore_deriv();
                                     
    Ref<TwoBodyInt> electron_repulsion();

    Ref<TwoBodyDerivInt> electron_repulsion_deriv();

    Ref<TwoBodyInt> grt();

    void set_basis(const Ref<GaussianBasisSet> &b1,
                   const Ref<GaussianBasisSet> &b2 = 0,
                   const Ref<GaussianBasisSet> &b3 = 0,
                   const Ref<GaussianBasisSet> &b4 = 0);

Integral::CartesianOrdering cartesian_ordering() const;
		    
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
