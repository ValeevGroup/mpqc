#include <reorder_engine.h>
using namespace MpqcCca

void
ReorderEngine::init( int n, 
              Ref<GuassianBasisSet> bs1,
              Ref<GaussianBasisSet> bs2,
              Ref<GaussianBasisSet> bs3,
              Ref<GaussianBasisSet> bs4  ) 
{
  n_center_ = n;

  basis_sets_.push_back(bs1);
  basis_sets_.push_back(bs2);
  basis_sets_.push_back(bs3);
  basis_sets_.push_back(bs4);

  // check that we have needed basis sets
  for( int i=0; i<n_center_; ++i )
    if( basis_set_[i].null() ) {
      sidl::SIDLException ex = sidl::SIDLException::_create();
      try {
        ex.setNote("Problem with basis sets in reorder engine");
        ex.add(__FILE__, __LINE__,"");
      }
      catch(...) { }
      throw ex;
    }

  // terminology:
  // segment = regarding basis function ids, the smallest repeating 
  //           portion of a buffer
  // A simple 0th order buffer has one segments.
  // A 1st order nuclear derivative buffer has three segements (dx,dy,dx), etc.

  // compute max segment size and max am
  max_segment_size_ = basis_sets_[0]->max_ncartesian_in_shell();
  maxam_ = basis_sets_[0]->max_angular_momentum();
  for( int i=1; i<n_center_; ++i) {
    max_segment_size_ *= basis_sets_[i]->max_ncartesian_in_shell();
    maxam_ = max( maxam_, basis_sets_[i]->max_angular_momentum() );
  }

  // construct reorder arrays
  for( int i=1; i<=maxam_; ++i) {

    Integral integral = new IntegralV3;
    CartesianIter *v3iter = integral->new_cartesian_iter(i);
    MPQC::CartesianIterCCA iter(i);
    MPQC::CartesianIterCCA *ccaiter = &iter;
    ccaiter->start();
    int ncf = ccaiter->n();

    reorder_[i] = new int[ncf];
    v3iter->start();
    for( int j=0; j<ncf; ++j) {
      ccaiter->start();
      for( int k=0; k<ncf; ++k) {
        if( v3iter->a() == ccaiter->a() &&
            v3iter->b() == ccaiter->b() &&
            v3iter->c() == ccaiter->c() ) {
          reorder_[i][j] = k;
          k=ncf; //break k loop
        }
        else ccaiter->next();
      }
      v3iter->next();
    }
  }

}

void
ReorderEngine::add_buffer ( double* buffer, IntegralDescr desc )
{
  buffers_[desc] = buffer;
  int deriv_lvl = desc.get_deriv_lvl();
  deriv_lvls_[desc] = deriv_lvl;
  if( deriv_lvl == 0 )
    temp_buffers_[desc] = new double[max_segment_size_];
  else if( deriv_lvl == 1 )
    temp_buffers_[desc] = new double[3*max_segment_size_];
  else {
    sidl::SIDLException ex = sidl::SIDLException::_create();
    try {
      ex.setNote("invalid deriv_lvl in reorder engine");
      ex.add(__FILE__, __LINE__,"");
    }
    catch(...) { }
    throw ex;
  }
}
  
void
ReorderEngine::do_it( int s1, int s2, int s3, int s4 )
{

  shells_ids_[0] = s1;
  shells_ids_[1] = s2;
  shells_ids_[2] = s3;
  shells_ids_[3] = s4;

  for( int i=0; i<n_center_; ++i )
    shells_[i] = basis_sets_[i]->shell( shell_ids[i] );

  int segment_size=1;
  for( int i=0; i<n_center_; ++i )
    segment_size *= shells_[i]->nfunction();

  map<IntegralDescr,double*>::iterator buff_it = buffers_.begin();
  for( ; buff_it != buffers_.end(); ++buff_it ) {

    int real_size;
    if( deriv_lvls_[ (*buff_it).first ] == 0 )
      real_size = segment_size;
    else if( deriv_lvls_[ (*buff_it).first ] == 1 )
      real_size = segment_size * 3;

    double* buf = (*buff_it).second;
    double* tbuf = temp_buffers_[ (*buff_it).first ];
    for( int i=0; i<real_size; ++i)
      tbuf[i] = buf[i];

    if( n_center_ == 1 )
      reorder_1c();
    else if( n_center_ == 2 )
      reorder_2c();
    else if( n_center_ == 3 )
      reorder_3c();
    else if( n_center_ == 4 )
      reorder_4c();
  }
}

void
ReorderEngine::reorder_1c()
{
}

void
ReorderEngine::reorder_2c()
{
}

void
ReorderEngine::reorder_3c()
{
}

void
ReorderEngine::reorder_4c()
{
}




}
