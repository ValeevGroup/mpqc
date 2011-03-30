#include <vector>
#include <Chemistry_QC_GaussianBasis_DescrInterface.hxx>
#include <chemistry/qc/basis/gaussbas.h>
#include <chemistry/qc/intv3/cartitv3.h>
#include <chemistry/qc/intv3/intv3.h>
#include <sidl_SIDLException.hxx>
#include <ccaiter.h>
#include <buffer_size.h>

namespace MpqcCca {

  class ReorderEngine {

  public:
    
    ReorderEngine( ) { }
    
    ~ReorderEngine( ) { }
    
  private:
    int n_center_, maxam_, size_, max_deriv_lvl_;
    int **reorder_;
    std::vector<int> deriv_lvls_;
    std::vector<int> segments_;
    std::vector<double*> buffers_;
    double* temp_buffer_;
    BufferSize buffer_size_;

    sc::Ref<sc::GaussianBasisSet> basis_sets_[4];
    int shell_ids_[4];
    sc::GaussianShell* shells_[4];
    
  public:
    
    void init( int n,
               sc::Ref<sc::GaussianBasisSet>,
               sc::Ref<sc::GaussianBasisSet>,
               sc::Ref<sc::GaussianBasisSet>,
               sc::Ref<sc::GaussianBasisSet> ); 
    
    void check_temp_buffer( int deriv_lvl, int n_segment );

    void
    add_buffer ( 
      double* buffer,
      Chemistry::QC::GaussianBasis::DescrInterface desc 
    );
  
    void do_it( int s1, int s2, int s3, int s4 );
 
    void reorder_1c( double* buf, double* tbuf, int offset );
    void reorder_2c( double* buf, double* tbuf, int offset, bool is_deriv );
    void reorder_3c( double* buf, double* tbuf, int offset );
    void reorder_4c( double* buf, double* tbuf, int offset );
    
  };
}
