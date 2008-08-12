#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/tbint.h>
#include <Chemistry_QC_GaussianBasis_DescrInterface.hxx>
#include <Chemistry_QC_GaussianBasis_CompositeDescrInterface.hxx>
#include <Chemistry_QC_GaussianBasis_DerivCentersInterface.hxx>
#include <Chemistry_QC_GaussianBasis_MolecularInterface.hxx>
#include <ChemistryDescrCXX_CompositeDescr.hxx>
#include <limits.h>
#include <vector>
#include <utility>
#include <sidl_SIDLException.hxx>

namespace MpqcCca {

  typedef Chemistry::QC::GaussianBasis::DescrInterface 
    QC_IntegralDescr;
  typedef Chemistry::QC::GaussianBasis::DerivCentersInterface
    QC_DerivCenters;
  typedef Chemistry::QC::GaussianBasis::CompositeDescrInterface
    QC_CompIntegralDescr;

  class onebody_onecenter_computer {

  private:
    int sh1_;

  public:
    onebody_onecenter_computer():
      sh1_(-1) { }

    void set_shells( int sh1 )
    {
      sh1_=sh1;
    }


    void compute(
        sc::OneBodyOneCenterInt* eval,
        QC_DerivCenters* dc )
    {
      eval->compute_shell( sh1_ );
    }

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

  public:
    onebody_onecenter_deriv_computer():
      sh1_(-1) { }

    void set_shells( int sh1 )
    {
      sh1_=sh1;
    }


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
    onebody_computer():
      sh1_(-1), sh2_(-1) { }

    void set_shells( int sh1, int sh2 )
    {
      sh1_=sh1; sh2_=sh2;
    }

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
    onebody_deriv_computer():
      sh1_(-1), sh2_(-1) { }

    void set_shells( int sh1, int sh2 )
    {
      sh1_=sh1; sh2_=sh2;
    }


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
    twobody_threecenter_computer():
      sh1_(-1), sh2_(-1), sh3_(-1) { }

    void set_shells( int sh1, int sh2, int sh3 )
    {
      sh1_=sh1; sh2_=sh2; sh3_=sh3;
    }


    void compute( sc::TwoBodyThreeCenterInt* eval, QC_DerivCenters* dc )
    {
      eval->compute_shell( sh1_, sh2_, sh3_ );
    }

    double compute_bounds( sc::TwoBodyThreeCenterInt* eval )
    { return eval->shell_bound( sh1_, sh2_, sh3_); }
  };


  class twobody_threecenter_deriv_computer {

  private:
    int sh1_, sh2_, sh3_;

  public:
    twobody_threecenter_deriv_computer():
      sh1_(-1), sh2_(-1), sh3_(-1) { }

    void set_shells( int sh1, int sh2, int sh3 )
    {
      sh1_=sh1; sh2_=sh2; sh3_=sh3;
    }

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
    twobody_computer():
      sh1_(-1), sh2_(-1), sh3_(-1), sh4_(-1) { }

    void set_shells( int sh1, int sh2, int sh3, int sh4 )
    {
      sh1_=sh1; sh2_=sh2; sh3_=sh3; sh4_=sh4;
    }


    void compute( sc::TwoBodyInt* eval, QC_DerivCenters* dc )
    { 
      eval->compute_shell( sh1_, sh2_, sh3_, sh4_ );
    }

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
    twobody_deriv_computer():
    sh1_(-1), sh2_(-1), sh3_(-1), sh4_(-1) { }

    void set_shells( int sh1, int sh2, int sh3, int sh4 )
    { 
      sh1_=sh1; sh2_=sh2; sh3_=sh3; sh4_=sh4; 
    }

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

  
  //---------------------------------------------------------------------------


  template< typename eval_type, typename computer_type >
  class IntegralEvaluator {

  public:
    
    IntegralEvaluator( ) { }
    
    ~IntegralEvaluator( ) {
    }
    
  private:
    
    std::vector< std::pair<eval_type*,QC_IntegralDescr> > evals_;
    std::vector< QC_DerivCenters > dcs_;
    std::vector< int > deriv_lvls_;
    std::vector< std::string > types_;
    std::vector< double* > buffers_;
    sidl::array<double> sidl_buffer_;
    
  public:
    
    void add_evaluator ( void* eval, QC_IntegralDescr desc ) 
    {
      eval_type* eval_ptr;   
      eval_ptr = static_cast<eval_type*>(eval);
      std::pair<eval_type*,QC_IntegralDescr> p(eval_ptr,desc);
      evals_.push_back( p );
      dcs_.push_back( p.second.get_deriv_centers() );
      deriv_lvls_.push_back( p.second.get_deriv_lvl() );
      types_.push_back( p.second.get_type() );
      buffers_.push_back( const_cast<double*>( p.first->buffer()) );
    }
    
    double* get_buffer ( QC_IntegralDescr desc ) 
    {
      for( int i=0; i<evals_.size(); ++i)
	if( desc.get_type() == types_[i] &&
	    desc.get_deriv_lvl() == deriv_lvls_[i] ) {
	  return const_cast<double*>( buffers_[i] );
	}
      return NULL;
    }

    sidl::array<double> get_array ( QC_IntegralDescr desc,
                                    int buffer_size )
    {
      for( int i=0; i<evals_.size(); ++i)
        if( desc.get_type() == types_[i] &&
            desc.get_deriv_lvl() == deriv_lvls_[i] ) {

          int lower[1] = {0};
          int upper[1];
          upper[0] = buffer_size - 1;
          int stride[1] = {1};

          sidl_buffer_.borrow( buffers_[i], 1, lower, upper, stride);
          return sidl_buffer_;
        }
    }
    
    Chemistry::QC::GaussianBasis::CompositeDescrInterface
    get_descriptor ()
    {
      Chemistry::QC::GaussianBasis::CompositeDescrInterface cdesc = 
        ChemistryDescrCXX::CompositeDescr::_create();
      for( int i=0; i<evals_.size(); ++i)
        cdesc.add_descr( evals_[i].second );
      return cdesc;
    }

    void compute ( computer_type* computer ) 
    {
      for( int i=0; i<evals_.size(); ++i) {
        if( deriv_lvls_[i] == 0 )
          computer->compute( evals_[i].first, NULL );
        else
          computer->compute( evals_[i].first, &(dcs_[i]) );
      }
    }

    double compute_bounds( computer_type* computer )
    {
      double bounds=0;
      if( evals_.size() ) 
        for( int i=0; i<evals_.size(); ++i)
          bounds = std::max( computer->compute_bounds(evals_[i].first), 
                             bounds );

      return bounds;
    }
      
      
    sidl::array<double> compute_array ( computer_type* computer,
                                        std::string type,
                                        int deriv_lvl,
                                        int buffer_size ) 
    {
      int lower[1] = {0};
      int upper[1];
      upper[0] = buffer_size - 1;
      int stride[1] = {1};

      for( int i=0; i<evals_.size(); ++i)
        if( types_[i] == type && deriv_lvls_[i] == deriv_lvl ) {

          if( deriv_lvls_[i] == 0 )
            computer->compute( evals_[i].first, NULL );
          else
            computer->compute( evals_[i].first, &(dcs_[i]) );

          sidl_buffer_.borrow( const_cast<double*>(buffers_[i]),
                               1, lower, upper, stride);
          return sidl_buffer_;
        }

      return NULL;
    }
      
  };


  //-------------------------------------------------------------------------
  template< typename eval_type, typename computer_type >
  class CompositeIntegralEvaluator {

  public:

    CompositeIntegralEvaluator( ): 
      have_new_shells_(true), sh1_(-1), sh2_(-1), sh3_(-1), sh4_(-1) { }

    ~CompositeIntegralEvaluator( ) {
    }

  private:

    std::vector< std::pair<eval_type*,QC_CompIntegralDescr> > evals_;
    std::vector< sidl::array<double> > sidl_buffers_;
    int sh1_, sh2_, sh3_, sh4_;
    bool have_new_shells_;
    std::vector< double* > buffers_;
    std::vector< std::string > types_;
    std::vector< int > deriv_lvls_;

  public:

    void set_shells( int sh1, int sh2, int sh3, int sh4 )
    {
      have_new_shells_ = false;
      if( sh1_ != sh1 ) {
        sh1_ = sh1;
        have_new_shells_ = true;
      }
      if( sh2 != -1 && sh2_ != sh2 ) {
        sh2_ = sh2;
        have_new_shells_ = true;
      }
      if( sh3 != -1 && sh3_ != sh3 ) {
        sh3_ = sh3;
        have_new_shells_ = true;
      }
      if( sh4 != -1 && sh4_ != sh4 ) {
        sh4_ = sh4;
        have_new_shells_ = true;
      }
    }

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
      std::string type = desc.get_type();
      int deriv_lvl = desc.get_deriv_lvl();
      for( int i=0; i<evals_.size(); ++i)
        for( int j=0; i<evals_[i].second.get_n_descr(); ++ j)
          if( type == evals_[i].second.get_descr(j).get_type() &&
              deriv_lvl == evals_[i].second.get_descr(j).get_deriv_lvl() ) {
            const double *b = evals_[i].first->buffer( tbt );
            types_.push_back( type );
            deriv_lvls_.push_back( deriv_lvl );
            buffers_.push_back( const_cast<double*>( b ) );
            return buffers_.back();
          }
      return NULL;
    }

    sidl::array<double> get_array ( QC_IntegralDescr desc,
                                    sc::TwoBodyInt::tbint_type tbt,
                                    int buffer_size )
    {
      std::string type = desc.get_type();
      int deriv_lvl = desc.get_deriv_lvl();
      for( int i=0; i<evals_.size(); ++i) {
        for( int j=0; j<evals_[i].second.get_n_descr(); ++ j) {
          if( type == evals_[i].second.get_descr(j).get_type() &&
              deriv_lvl == evals_[i].second.get_descr(j).get_deriv_lvl() ) {
            const double *b = evals_[i].first->buffer( tbt );
            types_.push_back( type );
            deriv_lvls_.push_back( deriv_lvl );
            buffers_.push_back( const_cast<double*>( b ) );
          }
        }
      }

      int lower[1] = {0};
      int upper[1];
      upper[0] = buffer_size - 1;
      int stride[1] = {1};

      for( int i=0; i<buffers_.size(); ++i)
        if( types_[i] == type && deriv_lvls_[i] == deriv_lvl ) {
          sidl::array<double> sa;
          sa.borrow( buffers_[i], 1, lower, upper, stride);
          sidl_buffers_.push_back(sa);
          return sidl_buffers_.back();
        }

    }

    Chemistry::QC::GaussianBasis::CompositeDescrInterface
    get_descriptor ()
    {
      Chemistry::QC::GaussianBasis::CompositeDescrInterface cdesc =
        ChemistryDescrCXX::CompositeDescr::_create();
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
      double bounds=0;
      if( evals_.size() )
        for( int i=0; i<evals_.size(); ++i)
          bounds = std::max( computer->compute_bounds(evals_[i].first),
                             bounds );

      return bounds;
    }


    sidl::array<double> compute_array ( computer_type* computer,
                                        std::string type,
                                        int deriv_lvl,
                                        int buffer_size )
    {
      if( have_new_shells_ )
        compute( computer );

      for( int i=0; i<buffers_.size(); ++i)
        if( types_[i] == type && deriv_lvls_[i] == deriv_lvl ) {
          return sidl_buffers_[i];
        }

      return NULL;
    } 

  };
}
