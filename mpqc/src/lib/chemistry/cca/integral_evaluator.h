#include <chemistry/qc/basis/integral.h>
#include <Chemistry_QC_GaussianBasis_IntegralDescr.hh>
#include <Chemistry_QC_GaussianBasis_CompositeIntegralDescr.hh>
#include <Chemistry_QC_GaussianBasis_DerivCenters.hh>
#include <Chemistry_QC_GaussianBasis_Molecular.hh>
#include <vector>
#include <utility>
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
  };


  template< typename eval_type, typename computer_type >
  class IntegralEvaluator {

  public:

    IntegralEvaluator( ) { };
    
    ~IntegralEvaluator( )
    {
      /*
      if( reorder_ ) {
	delete [] temp_buffer_;
	for( int i=0; i<=maxam_; ++i)
	  delete [] reorder_[i];
	delete [] reorder_;
      }
      */
    }
    
    void add_evaluator ( void* eval, IntegralDescr desc );
    void set_basis( vector<Molecular> );
    void set_reorder (bool reorder) { reorder_ = reorder; } 
    void* get_buffer ( IntegralDescr desc ); 
    QC_DerivCenters get_deriv_centers ();
    CompositeIntegralDescr get_descriptor ();
    void compute ( computer_type* computer ); 
    sidl::array<double> compute_array ( computer_type* computer ); 

  private:

    bool reorder_;

    vector< pair<eval_type*,QC_IntegralDescr> > evals_;

    void initialize_reorder(); 
    void reorder( int shellnum1, int shellnum2 );
    void reorder_doublet( sc::GaussianShell* s1, 
			  sc::GaussianShell* s2,
			  int nc1, int nc2, int is_deriv );
  
     };

};
