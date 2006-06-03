#include <reorder_engine.h>
using namespace MpqcCca;

void
ReorderEngine::init( int n, 
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
      ex.setNote("invalid n centers in reorder engine");
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

  // compute max segment size and max am
  int max_segment_size_ = basis_sets_[0]->max_ncartesian_in_shell();
  maxam_ = basis_sets_[0]->max_angular_momentum();
  for( int i=1; i<n_center_; ++i) {
    max_segment_size_ *= basis_sets_[i]->max_ncartesian_in_shell();
    maxam_ = max( maxam_, basis_sets_[i]->max_angular_momentum() );
  }

  max_deriv_lvl_ = 0; // until we discover otherwise
  temp_buffer_ = new double[max_segment_size_*3];

  // construct reorder arrays

  reorder_ = new int*[maxam_+1];
  reorder_[0] = new int[1];
  reorder_[0][0] = 0;

  for( int i=1; i<=maxam_; ++i) {

    IntegralV3 integral;
    CartesianIter *v3iter = integral.new_cartesian_iter(i);
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
  buffers_.push_back( buffer );
  int deriv_lvl = desc.get_deriv_lvl();
  deriv_lvls_.push_back( deriv_lvl );
  if( max_deriv_lvl_ < deriv_lvl ) {
    max_deriv_lvl_ = deriv_lvl;
    //delete temp_buffer_; //mistmatched ???
    //temp_buffer_ = new double[max_segment_size_*3];
  } 
  if( max_deriv_lvl_ > 1 ) {
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
  shell_ids_[0] = s1;
  shell_ids_[1] = s2;
  shell_ids_[2] = s3;
  shell_ids_[3] = s4;

  for( int i=0; i<n_center_; ++i )
    shells_[i] = &( basis_sets_[i]->shell(shell_ids_[i]) );

  int segment_size=1;
  for( int i=0; i<n_center_; ++i )
    segment_size *= shells_[i]->nfunction();

  vector<double*>::iterator buff_it = buffers_.begin();
  vector<int>::iterator lvl_it = deriv_lvls_.begin();
  for( ; buff_it != buffers_.end(); ++buff_it ) {

    int real_size;
    int deriv_lvl = *lvl_it;
    if( deriv_lvl == 0 )
      real_size = segment_size;
    else if( deriv_lvl == 1 )
      real_size = segment_size * 3;

    double* buf = *buff_it;
    for( int i=0; i<real_size; ++i)
      temp_buffer_[i] = buf[i];

    int deriv_offset;
    if( n_center_ == 1 ) {
      if( deriv_lvl == 0 )
        reorder_1c( buf, temp_buffer_, 0 );
    }
    else if( n_center_ == 2 ) {
      if( deriv_lvl == 0 )
        reorder_2c( buf, temp_buffer_, 0 );
      else if( deriv_lvl == 1 ) {
        // a 2-center intv3 derivative buffer is composed of 
        // nfunc triplets (dx,dy,dz) which must be repacked
        // into 3 (dx,dy,dz) shell doublets of nfunc entries
        reorder_2c( buf, temp_buffer_, 1 );
      }
    }
    else if( n_center_ == 3 ) {
      if( deriv_lvl == 0 )
        reorder_3c( buf, temp_buffer_, 0 );
    }
    else if( n_center_ == 4 ) {
      if( deriv_lvl == 0 )
        reorder_4c( buf, temp_buffer_, 0 );
      else if( deriv_lvl == 1 )
        for(int i=0; i<3; ++i) {
          deriv_offset = i*segment_size;
          reorder_4c( buf, temp_buffer_, deriv_offset );
        }
    }
  }

}

void
ReorderEngine::reorder_1c( double* buf, double* tbuf, int offset )
{
}

void
ReorderEngine::reorder_2c( double* buf, double* tbuf, int is_deriv )
{
  GaussianShell* s1 = shells_[0];
  GaussianShell* s2 = shells_[1];
  int nc1 = s1->ncontraction();
  int nc2 = s2->ncontraction();

  int index=0, con2_offset=0, local2_offset, 
      c1_base, c2_base;

  int temp;
  con2_offset = s2->nfunction();

  int s1_is_cart, s2_is_cart, nfunc, s1_nfunc, s2_nfunc; 
  nfunc = s1->nfunction() * s2->nfunction();

  c1_base = 0;
  for( int c1=0; c1<nc1; ++c1 ) {

    //c1_base = index;
    if(c1>0) c1_base += s1->nfunction(c1-1) * con2_offset;

    s1_is_cart = s1->is_cartesian(c1);
    s1_nfunc = s1->nfunction(c1);

    for( int fc1=0; fc1<s1_nfunc; ++fc1 ) {

      if( s1_is_cart )
        c2_base = c1_base + reorder_[s1->am(c1)][fc1] * con2_offset;
      else
        c2_base = c1_base + fc1 * con2_offset;

      local2_offset = 0;
      for( int c2=0; c2<nc2; ++c2 ) {
        if( c2>0 ) local2_offset += s2->nfunction(c2-1);
        s2_is_cart = s2->is_cartesian(c2);
        s2_nfunc = s2->nfunction(c2);

        if(!is_deriv) {
          if( s2_is_cart )
            for( int fc2=0; fc2<s2_nfunc; ++fc2 ) {
              buf[ c2_base + local2_offset + reorder_[s2->am(c2)][fc2] ]
                = tbuf[index];
              ++index;
            }
          else
            for( int fc2=0; fc2<s2_nfunc; ++fc2 ) {
              buf[ c2_base + local2_offset + fc2 ] = tbuf[index];
              ++index;
            }
        }
        else {
          if( s2_is_cart )
            for( int fc2=0; fc2<s2_nfunc; ++fc2 ) {
              for(int di=0; di<3; ++di) {
                buf[ c2_base + local2_offset + reorder_[s2->am(c2)][fc2]
                        + di*nfunc ]
                  = tbuf[index];
                ++index;
              }
            }
          else
            for( int fc2=0; fc2<s2_nfunc; ++fc2 ) {
              for(int di=0; di<3; ++di) {
                buf[ c2_base + local2_offset + fc2 + di*nfunc ] 
                  = tbuf[index];
                ++index;
              }
            }
        }

      }
    }
  }

}

void
ReorderEngine::reorder_3c( double* buf, double* tbuf, int offset )
{
}

void
ReorderEngine::reorder_4c( double* buf, double* tbuf, int offset )
{
  GaussianShell* s1 = shells_[0];
  GaussianShell* s2 = shells_[1];
  GaussianShell* s3 = shells_[2];
  GaussianShell* s4 = shells_[3];
  int nc1 = s1->ncontraction();
  int nc2 = s2->ncontraction();
  int nc3 = s3->ncontraction();
  int nc4 = s4->ncontraction();

  int index=offset, con2_offset=0, con3_offset=0, con4_offset=0, con_offset,
    local2_offset, local3_offset, local4_offset,
    c1_base, c2_base, c3_base, c4_base;                                       

  int temp;                                                                       for( int c4=0; c4<nc4; ++c4 )
    con4_offset += s4->nfunction(c4);                                           
  temp = 0;                                                                       con3_offset = con4_offset;
  for( int c3=0; c3<nc3; ++c3 )                                                     temp += s3->nfunction(c3);
  con3_offset *= temp;                                                          
  temp = 0;                                                                       con2_offset = con3_offset;
  for( int c2=0; c2<nc2; ++c2 )                                                     temp += s2->nfunction(c2);
  con2_offset *= temp;                                                          

  int s1_is_cart, s2_is_cart, s3_is_cart, s4_is_cart,                                 s1_nfunc, s2_nfunc, s3_nfunc, s4_nfunc;
                                                                                  for( int c1=0; c1<nc1; ++c1 ) {
                                                                                    c1_base = index;
    s1_is_cart = s1->is_cartesian(c1);                                              s1_nfunc = s1->nfunction(c1);
                                                                                    for( int fc1=0; fc1<s1_nfunc; ++fc1 ) {
                                                                                      if( s1_is_cart )
        c2_base = c1_base + reorder_[s1->am(c1)][fc1] * con2_offset;                  else
        c2_base = c1_base + fc1 * con2_offset;                                  
      local2_offset = 0;                                                              for( int c2=0; c2<nc2; ++c2 ) {

        if( c2>0 ) local2_offset += s2->nfunction(c2-1);
        s2_is_cart = s2->is_cartesian(c2);                                              s2_nfunc = s2->nfunction(c2);
                                                                                        for( int fc2=0; fc2<s2_nfunc; ++fc2 ) {

          if( s2_is_cart )
            c3_base = c2_base + (local2_offset + reorder_[s2->am(c2)][fc2])
                        * con3_offset;
          else
            c3_base = c2_base + (local2_offset + fc2) * con3_offset;

          local3_offset = 0;
          for( int c3=0; c3<nc3; ++c3 ) {

            if( c3>0 ) local3_offset += s3->nfunction(c3-1);
            s3_is_cart = s3->is_cartesian(c3);
            s3_nfunc = s3->nfunction(c3);

            for( int fc3=0; fc3<s3_nfunc; ++fc3 ) {

              if( s3_is_cart )
                c4_base = c3_base + (local3_offset + reorder_[s3->am(c3)][fc3])
                            * con4_offset;
              else
                c4_base = c3_base + (local3_offset + fc3) * con4_offset;

              local4_offset = 0;
              for( int c4=0; c4<nc4; ++c4 ) {

                if( c4>0 ) local4_offset += s4->nfunction(c4-1);
                s4_is_cart = s4->is_cartesian(c4);
                s4_nfunc = s4->nfunction(c4);

                if( s4_is_cart )
                  for( int fc4=0; fc4<s4_nfunc; ++fc4 ) {
                    buf[ c4_base + local4_offset + reorder_[s4->am(c4)][fc4] ]
                      = temp_buffer_[index];
                    ++index;
                  }
                else
                  for( int fc4=0; fc4<s4_nfunc; ++fc4 ) {
                    buf[ c4_base + local4_offset + fc4 ] = temp_buffer_[index];
                    ++index;
                  }
              }
            }
          }
        }
      }
    }
  }

}

