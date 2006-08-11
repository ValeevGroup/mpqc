#include <vector>
#include <chemistry/qc/basis/gaussbas.h>
#include <sidl_SIDLException.hh>

/**
     This utility class keeps track of what size temporary buffers should be.
*/

namespace MpqcCca {

  class BufferSize {

  public:
    
    BufferSize( ) { }
    
    ~BufferSize( ) { }
    
  private:
    int n_center_, max_deriv_lvl_, max_segment_size_, max_n_segment_, size_;
    sc::Ref<sc::GaussianBasisSet> basis_sets_[4];
    
  public:
    
    void init( int n,
               sc::Ref<sc::GaussianBasisSet>,
               sc::Ref<sc::GaussianBasisSet>,
               sc::Ref<sc::GaussianBasisSet>,
               sc::Ref<sc::GaussianBasisSet> ); 
    
    void update( int deriv_lvl, int n_segment );

    int size() { return size_; }
    
  };
}
