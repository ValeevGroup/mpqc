#include <vector>
#include <Chemistry_QC_GaussianBasis_IntegralDescr.hh>
#include <Chemistry_QC_GaussianBasis_Molecular.hh>
#include <chemistry/qc/basis/gaussbas.h>
#include <chemistry/qc/intv3/cartitv3.h>
#include <chemistry/qc/intv3/intv3.h>
#include <sidl_SIDLException.hh>
#pragma implementation "ccaiter.h"
#include <ccaiter.h>

sc::Ref<sc::GaussianBasisSet> basis_cca_to_sc(
    Chemistry::QC::GaussianBasis::Molecular &cca_basis );

namespace MpqcCca {

  class ReorderEngine {

  public:
    
    ReorderEngine( ) { }
    
    ~ReorderEngine( ) { }
    
  private:
    int n_center_, maxam_, max_deriv_lvl_, max_segment_size_;
    int **reorder_;
    std::vector<int> deriv_lvls_;
    std::vector<double*> buffers_;
    double* temp_buffer_;
    sc::Ref<sc::GaussianBasisSet> basis_sets_[4];
    int shell_ids_[4];
 
    //Ref<GaussianShell> shells_[4];
    sc::GaussianShell* shells_[4];
    
  public:
    
    void init( int n,
               sc::Ref<sc::GaussianBasisSet>,
               sc::Ref<sc::GaussianBasisSet>,
               sc::Ref<sc::GaussianBasisSet>,
               sc::Ref<sc::GaussianBasisSet> ); 

    void add_buffer ( double* buffer,
                      Chemistry::QC::GaussianBasis::IntegralDescr desc );
  
    void do_it( int s1, int s2, int s3, int s4 );
 
    void reorder_1c( double* buf, double* tbuf, int offset );
    void reorder_2c( double* buf, double* tbuf, int is_deriv );
    void reorder_3c( double* buf, double* tbuf, int offset );
    void reorder_4c( double* buf, double* tbuf, int offset );
    
  };
}
