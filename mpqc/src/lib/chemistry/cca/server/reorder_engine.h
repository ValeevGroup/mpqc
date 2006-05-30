#include <vector>
#include <Chemistry_QC_GaussianBasis_IntegralDescr.hh>
#include <Chemistry_QC_GaussianBasis_Molecular.hh>
#include <chemistry/qc/basis/gaussbas.h>
#include <chemistry/qc/intv3/cartitv3.h>
#include <chemistry/qc/intv3/intv3.h>
#include <sidl_SIDLException.hh>
#pragma implementation "ccaiter.h"
#include <ccaiter.h>

using namespace std;
using namespace sc;
using namespace Chemistry::QC::GaussianBasis;

Ref<GaussianBasisSet> basis_cca_to_sc( Molecular &cca_basis );

namespace MpqcCca {

  class ReorderEngine {

  public:
    
    ReorderEngine( ) { }
    
    ~ReorderEngine( ) { }
    
  private:
    n_center_;
    maxam_;
    map<IntegralDescr,double*> deriv_lvls_;
    map<IntegralDescr,double*> buffers_;
    map<IntegralDescr,double*> temp_buffers_;
    vector< Ref<GaussianBasisSet> > basis_sets_;
    vector< int > shell_ids_;
    vector< Ref<GaussianShell> > shells_;
    
  public:
    
    void init( int n, Ref<GaussianBasisSet>, Ref<GaussianBasisSet>,
               Ref<GaussianBasisSet>, Ref<GaussianBasisSet> ); 

    void add_buffer ( double* buffer, IntegralDescr desc );
  
    void do_it( int s1, int s2, int s3, int s4 );
    
  };
}
