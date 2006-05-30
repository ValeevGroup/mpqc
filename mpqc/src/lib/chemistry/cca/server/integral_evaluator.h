#include <chemistry/qc/basis/integral.h>
#include <Chemistry_QC_GaussianBasis_IntegralDescr.hh>
#include <Chemistry_QC_GaussianBasis_CompositeIntegralDescr.hh>
#include <Chemistry_QC_GaussianBasis_DerivCenters.hh>
#include <Chemistry_QC_GaussianBasis_Molecular.hh>
#include <Chemistry_CompositeIntegralDescr.hh>
#include <limits.h>
#include <vector>
#include <utility>
#include <sidl_SIDLException.hh>
using namespace std;
using namespace sc;
using namespace Chemistry::QC::GaussianBasis;

namespace MpqcCca {

  typedef Chemistry::QC::GaussianBasis::IntegralDescr QC_IntegralDescr;
  typedef Chemistry::QC::GaussianBasis::DerivCenters QC_DerivCenters;

  class onebody_onecenter_computer {

  private:
    int sh1_;
    Molecular bs1_;

  public:
    onebody_onecenter_computer() { }

    void set_shells( int sh1 )
    { sh1_ = sh1; }

    void compute( OneBodyOneCenterInt* eval )
    { eval->compute_shell( sh1_ ); }

    double compute_bounds( OneBodyOneCenterInt* eval )
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

  class onebody_computer {

  private:
    int sh1_, sh2_;

  public:
    onebody_computer() { }

    void set_shells( int sh1, int sh2 )
    { sh1_=sh1; sh2_=sh2; }

    void compute( OneBodyInt* eval )
    { eval->compute_shell( sh1_, sh2_ ); }

    double compute_bounds( OneBodyInt* eval )
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

  class twobody_threecenter_computer {

  private:
    int sh1_, sh2_, sh3_;

  public:
    twobody_threecenter_computer() { }

    void set_shells( int sh1, int sh2, int sh3 )
    { sh1_=sh1; sh2_=sh2; sh3_=sh3; }

    void compute( TwoBodyThreeCenterInt* eval )
    { eval->compute_shell( sh1_, sh2_, sh3_ ); }

    double compute_bounds( TwoBodyThreeCenterInt* eval )
    { return eval->shell_bound( sh1_, sh2_, sh3_); }
  };

  class twobody_computer {

  private:
    int sh1_, sh2_, sh3_, sh4_;

  public:
    twobody_computer() { }

    void set_shells( int sh1, int sh2, int sh3, int sh4 )
    { sh1_=sh1; sh2_=sh2; sh3_=sh3; sh4_=sh4; };

    void compute( TwoBodyInt* eval )
    { eval->compute_shell( sh1_, sh2_, sh3_, sh4_ ); }

    double compute_bounds( TwoBodyInt* eval )
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
    
    vector< pair<eval_type*,QC_IntegralDescr> > evals_;
    
  public:
    
    void add_evaluator ( void* eval, IntegralDescr desc ) 
    {
      eval_type* eval_ptr;   
      eval_ptr = static_cast<eval_type*>(eval);
      pair<eval_type*,IntegralDescr> p(eval_ptr,desc);
      evals_.push_back( p );
    }
    
    double* get_buffer ( IntegralDescr desc ) 
    {
      for( int i=0; i<evals_.size(); ++i)
	if( desc.get_type() == evals_[i].second.get_type() &&
	    desc.get_deriv_lvl() == evals_[i].second.get_deriv_lvl() ) {
	  return const_cast<double*>( evals_[i].first->buffer() );
	}
      return NULL;
    }
    
    Chemistry::QC::GaussianBasis::DerivCenters get_deriv_centers ()
    {
      // later
    }
    
    Chemistry::QC::GaussianBasis::CompositeIntegralDescr get_descriptor ()
    {
      CompositeIntegralDescr cdesc = 
        Chemistry::CompositeIntegralDescr::_create();
      for( int i=0; i<evals_.size(); ++i)
        cdesc.add_descr( evals_[i].second );
      return cdesc;
    }

    void compute ( computer_type* computer ) 
    {
      for( int i=0; i<evals_.size(); ++i) 
        computer->compute( evals_[i].first );
    }

    double compute_bounds( computer_type* computer )
    {
      // this is obviously not going to work for multiple evals
      // that will require interface work
      for( int i=0; i<evals_.size(); ++i)
        return computer->compute_bounds( evals_[i].first );
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
