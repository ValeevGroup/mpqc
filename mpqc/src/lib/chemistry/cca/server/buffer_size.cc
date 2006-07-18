#include <buffer_size.h>

using namespace std;
using namespace sc;

using namespace MpqcCca;

void
BufferSize::init( int n, 
                  Ref<GaussianBasisSet> bs1,
                  Ref<GaussianBasisSet> bs2,
                  Ref<GaussianBasisSet> bs3,
                  Ref<GaussianBasisSet> bs4  ) 
{
  if( 1 <= n <= 4 )
    n_center_ = n;
  else {
    sidl::SIDLException ex = sidl::SIDLException::_create();
    try {
      ex.setNote("invalid n centers");
      ex.add(__FILE__, __LINE__,"");
    }
    catch(...) { }
    throw ex;
  }

  basis_sets_[0] = bs1;
  basis_sets_[1] = bs2;
  basis_sets_[2] = bs3;
  basis_sets_[3] = bs4;

  // check that we have needed basis sets
  for( int i=0; i<n_center_; ++i )
    if( basis_sets_[i].null() ) {
      sidl::SIDLException ex = sidl::SIDLException::_create();
      try {
        ex.setNote("Problem with basis sets in reorder engine");
        ex.add(__FILE__, __LINE__,"");
      }
      catch(...) { }
      throw ex;
    }

  // terminology:
  // Segment: regarding basis function ids, the smallest repeating 
  // portion of a buffer.
  // A simple 0th order buffer has one segment.
  // A 1st order nuclear derivative buffer has three segements (dx,dy,dx), etc.

  // compute max segment size
  max_segment_size_ = basis_sets_[0]->max_ncartesian_in_shell();
  for( int i=1; i<n_center_; ++i)
    max_segment_size_ *= basis_sets_[i]->max_ncartesian_in_shell();

  // until we discover otherwise
  max_deriv_lvl_ = 0;
  max_n_segment_ = 1;
  size_ = max_segment_size_;

}

void
BufferSize::update( int deriv_lvl, int n_segment )
{
  if( deriv_lvl > max_deriv_lvl_ )
    max_deriv_lvl_ = deriv_lvl;
  if( n_segment > max_n_segment_ )
    max_n_segment_ = n_segment;

  int nderiv;
  if( max_deriv_lvl_ == 1 ) {
    nderiv = 3;
    if( n_center_ > 2 )
      nderiv *= n_center_ - 1;
  }
  else 
    nderiv = 1;

  size_ = max_segment_size_ * max_n_segment_ * nderiv;
}

