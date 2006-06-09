#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/tbint.h>
#include <Chemistry_QC_GaussianBasis_IntegralDescr.hh>
#include <Chemistry_QC_GaussianBasis_CompositeIntegralDescr.hh>
#include <Chemistry_QC_GaussianBasis_DerivCenters.hh>
#include <Chemistry_QC_GaussianBasis_Molecular.hh>
#include <Chemistry_CompositeIntegralDescr.hh>
#include <limits.h>
#include <vector>
#include <utility>
#include <sidl_SIDLException.hh>

namespace MpqcCca {

  typedef Chemistry::QC::GaussianBasis::IntegralDescr QC_IntegralDescr;
  typedef Chemistry::QC::GaussianBasis::DerivCenters QC_DerivCenters;
  typedef Chemistry::QC::GaussianBasis::CompositeIntegralDescr
    QC_CompIntegralDescr;

  class onebody_onecenter_computer {

  private:
    int sh1_;
    Chemistry::QC::GaussianBasis::Molecular bs1_;

  public:
    onebody_onecenter_computer() { }

    void set_shells( int sh1 )
    { sh1_ = sh1; }

    void compute(
        sc::OneBodyOneCenterInt* eval,
        QC_DerivCenters* dc )
    { eval->compute_shell( sh1_ ); }

    double compute_bounds(
        sc::OneBodyOneCenterInt* eval )
    {  
      sidl::SIDLException ex = sidl::SIDLException::_create();
      try {
        ex.setNote("MPQC doesn't support one body in bounds");
        ex.add(__FILE__, __LINE__,"");
      }
      catch(...) { }
      throw ex;
    }

  };


  class onebody_onecenter_deriv_computer {

  private:
    int sh1_;
    Chemistry::QC::GaussianBasis::Molecular bs1_;

  public:
    onebody_onecenter_deriv_computer() { }

    void set_shells( int sh1 )
    { sh1_ = sh1; }

    void compute(
        sc::OneBodyOneCenterDerivInt* eval,
        QC_DerivCenters* dc )
    {
    }

    double compute_bounds(
        sc::OneBodyOneCenterDerivInt* eval )
    {
      sidl::SIDLException ex = sidl::SIDLException::_create();
      try {
        ex.setNote("MPQC doesn't support one body int bounds");
        ex.add(__FILE__, __LINE__,"");
      }
      catch(...) { }
      throw ex;
    }

  };


  class onebody_computer {

  private:
    int sh1_, sh2_;

  public:
    onebody_computer() { }

    void set_shells( int sh1, int sh2 )
    { sh1_=sh1; sh2_=sh2; }

    void compute( sc::OneBodyInt* eval, QC_DerivCenters* dc )
    { 
      eval->compute_shell( sh1_, sh2_ ); 
    }

    double compute_bounds( sc::OneBodyInt* eval )
    {
      sidl::SIDLException ex = sidl::SIDLException::_create();
      try {
        ex.setNote("MPQC doesn't support one body in bounds");
        ex.add(__FILE__, __LINE__,"");
      }
      catch(...) { }
      throw ex;
    }
  };


  class onebody_deriv_computer {

  private:
    int sh1_, sh2_;

  public:
    onebody_deriv_computer() { }

    void set_shells( int sh1, int sh2 )
    { sh1_=sh1; sh2_=sh2; }

    void compute( sc::OneBodyDerivInt* eval, QC_DerivCenters* dc )
    {
      eval->compute_shell( sh1_, sh2_, dc->get_deriv_atom() );
    }

    double compute_bounds( sc::OneBodyDerivInt* eval )
    {
      sidl::SIDLException ex = sidl::SIDLException::_create();
      try {
        ex.setNote("MPQC doesn't support one body int bounds");
        ex.add(__FILE__, __LINE__,"");
      }
      catch(...) { }
      throw ex;
    }
  };


  class twobody_threecenter_computer {

  private:
    int sh1_, sh2_, sh3_;

  public:
    twobody_threecenter_computer() { }

    void set_shells( int sh1, int sh2, int sh3 )
    { sh1_=sh1; sh2_=sh2; sh3_=sh3; }

    void compute( sc::TwoBodyThreeCenterInt* eval, QC_DerivCenters* dc )
    { eval->compute_shell( sh1_, sh2_, sh3_ ); }

    double compute_bounds( sc::TwoBodyThreeCenterInt* eval )
    { return eval->shell_bound( sh1_, sh2_, sh3_); }
  };


  class twobody_threecenter_deriv_computer {

  private:
    int sh1_, sh2_, sh3_;

  public:
    twobody_threecenter_deriv_computer() { }

    void set_shells( int sh1, int sh2, int sh3 )
    { sh1_=sh1; sh2_=sh2; sh3_=sh3; }

    void compute( sc::TwoBodyThreeCenterDerivInt* eval, QC_DerivCenters* dc )
    {
    }

    double compute_bounds( sc::TwoBodyThreeCenterDerivInt* eval )
    { return eval->shell_bound( sh1_, sh2_, sh3_); }
  };


  class twobody_computer {

  private:
    int sh1_, sh2_, sh3_, sh4_;

  public:
    twobody_computer() { }

    void set_shells( int sh1, int sh2, int sh3, int sh4 )
    { sh1_=sh1; sh2_=sh2; sh3_=sh3; sh4_=sh4; };

    void compute( sc::TwoBodyInt* eval, QC_DerivCenters* dc )
    { eval->compute_shell( sh1_, sh2_, sh3_, sh4_ ); }

    double compute_bounds( sc::TwoBodyInt* eval )
    {
      return eval->shell_bound( sh1_, sh2_, sh3_, sh4_ );
    }
  };


  class twobody_deriv_computer {

  private:
    int sh1_, sh2_, sh3_, sh4_;
 
  public:
    sc::DerivCenters sc_dc_;

  public:
    twobody_deriv_computer() { }

    void set_shells( int sh1, int sh2, int sh3, int sh4 )
    { sh1_=sh1; sh2_=sh2; sh3_=sh3; sh4_=sh4; };

    void compute( sc::TwoBodyDerivInt* eval, QC_DerivCenters* dc )
    {
      sc_dc_.clear();
      eval->compute_shell(sh1_,sh2_,sh3_,sh4_,sc_dc_);

      dc->clear();
      if( sc_dc_.has_omitted_center() )
        dc->add_omitted(sc_dc_.omitted_center(),sc_dc_.omitted_atom());
      for( int i=0; i<sc_dc_.n(); ++i)
        dc->add_center(sc_dc_.center(i),sc_dc_.atom(i));

    }

    double compute_bounds( sc::TwoBodyDerivInt* eval )
    {
      return eval->shell_bound( sh1_, sh2_, sh3_, sh4_ );
    }
  };


  template< typename eval_type, typename computer_type >
  class IntegralEvaluator {

  public:
    
    IntegralEvaluator( ) { }
    
    ~IntegralEvaluator( ) { }
    
  private:
    
    std::vector< std::pair<eval_type*,QC_IntegralDescr> > evals_;
    std::vector< QC_DerivCenters > dcs_;
    
  public:
    
    void add_evaluator ( void* eval, QC_IntegralDescr desc ) 
    {
      eval_type* eval_ptr;   
      eval_ptr = static_cast<eval_type*>(eval);
      std::pair<eval_type*,QC_IntegralDescr> p(eval_ptr,desc);
      evals_.push_back( p );
      dcs_.push_back( p.second.get_deriv_centers() );
    }
    
    double* get_buffer ( QC_IntegralDescr desc ) 
    {
      for( int i=0; i<evals_.size(); ++i)
	if( desc.get_type() == evals_[i].second.get_type() &&
	    desc.get_deriv_lvl() == evals_[i].second.get_deriv_lvl() ) {
	  return const_cast<double*>( evals_[i].first->buffer() );
	}
      return NULL;
    }
    
    Chemistry::QC::GaussianBasis::CompositeIntegralDescr get_descriptor ()
    {
      Chemistry::QC::GaussianBasis::CompositeIntegralDescr cdesc = 
        Chemistry::CompositeIntegralDescr::_create();
      for( int i=0; i<evals_.size(); ++i)
        cdesc.add_descr( evals_[i].second );
      return cdesc;
    }

    void compute ( computer_type* computer ) 
    {
      for( int i=0; i<evals_.size(); ++i) {
        if( evals_[i].second.get_deriv_lvl() == 0 )
          computer->compute( evals_[i].first, NULL );
        else
          computer->compute( evals_[i].first, &(dcs_[i]) );
      }
    }

    double compute_bounds( computer_type* computer )
    {
      // this is obviously not going to work for multiple evals
      // that will require interface work
      if( evals_.size() ) 
        for( int i=0; i<evals_.size(); ++i)
          return computer->compute_bounds( evals_[i].first );

      return 0.0;
    }
      
      
    sidl::array<double> compute_array ( computer_type* computer ) 
    {
      /*
      compute( shellnum1, shellnum2, deriv_level, deriv_atom );
  
      int lower[1] = {0};
      int upper[1]; upper[0] = max_nshell2_-1;
      int stride[1] = {1};
      sidl_buffer_.borrow( const_cast<double*>(sc_buffer_), 1, 
                           lower, upper, stride);
      return sidl_buffer_;
      */
    }
      
  };

  template< typename eval_type, typename computer_type >
  class CompositeIntegralEvaluator {

  public:

    CompositeIntegralEvaluator( ) { }

    ~CompositeIntegralEvaluator( ) { }

  private:

    std::vector< std::pair<eval_type*,QC_CompIntegralDescr> > evals_;

  public:

    void add_evaluator ( void* eval, QC_CompIntegralDescr cdesc )
    {
      eval_type* eval_ptr;
      eval_ptr = static_cast<eval_type*>(eval);
      std::pair<eval_type*,QC_CompIntegralDescr> p(eval_ptr,cdesc);
      evals_.push_back( p );
    }

    double* get_buffer ( QC_IntegralDescr desc, 
                         sc::TwoBodyInt::tbint_type tbt )
    {
      for( int i=0; i<evals_.size(); ++i)
        for( int j=0; i<evals_[i].second.get_n_descr(); ++ j)
          if( desc.get_type() == evals_[i].second.get_descr(j).get_type() &&
              desc.get_deriv_lvl() == 
                evals_[i].second.get_descr(j).get_deriv_lvl() )
            return const_cast<double*>( evals_[i].first->buffer( tbt ) );
      return NULL;
    }

    Chemistry::QC::GaussianBasis::DerivCenters get_deriv_centers ()
    {
      // later
    }

    Chemistry::QC::GaussianBasis::CompositeIntegralDescr get_descriptor ()
    {
      Chemistry::QC::GaussianBasis::CompositeIntegralDescr cdesc =
        Chemistry::CompositeIntegralDescr::_create();
      for( int i=0; i<evals_.size(); ++i)
        for( int j=0; j<evals_[i].second.get_n_descr(); ++j )
          cdesc.add_descr( evals_[i].second.get_descr(j) );
      return cdesc;
    }

    void compute ( computer_type* computer )
    {
      for( int i=0; i<evals_.size(); ++i)
        computer->compute( evals_[i].first, 0 );
    }

    double compute_bounds( computer_type* computer )
    {
      // this is obviously not going to work for multiple evals
      // that will require interface work
      if( evals_.size() )
        for( int i=0; i<evals_.size(); ++i)
          return computer->compute_bounds( evals_[i].first );

      return 0.0;
    }


    sidl::array<double> compute_array ( computer_type* computer )
    {
      /*
      compute( shellnum1, shellnum2, deriv_level, deriv_atom );

      int lower[1] = {0};
      int upper[1]; upper[0] = max_nshell2_-1;
      int stride[1] = {1};
      sidl_buffer_.borrow( const_cast<double*>(sc_buffer_), 1,
                           lower, upper, stride);
      return sidl_buffer_;
      */
    }

  };
}
